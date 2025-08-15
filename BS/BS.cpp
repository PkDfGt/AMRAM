#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cmath>
#include <omp.h>
#include "../AESN/AMAES.h"
#include "BS.h"
#include <openssl/evp.h>

using ByteVec = std::vector<unsigned char>;

int total_comps = 0;
int curr_comps = 0;
int ctr = 0;
constexpr int INT_BYTE_LEN = 4;

void printProgressBar(float progress) {
    int barWidth = 40;
    std::cout << "\r[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << static_cast<int>(progress * 100.0) << " %";
    std::cout.flush();
}

ByteVec encodeInt(int value, int size) {
    if (size < INT_BYTE_LEN) throw std::invalid_argument("size too small for encoding int");
    ByteVec buffer(size, 0x00);
    buffer[size - 4] = (value >> 24) & 0xFF;
    buffer[size - 3] = (value >> 16) & 0xFF;
    buffer[size - 2] = (value >> 8) & 0xFF;
    buffer[size - 1] = value & 0xFF;
    return buffer;
}

int decodeInt(const ByteVec& data) {
    if (data.size() < INT_BYTE_LEN) throw std::invalid_argument("data too small to decode int");
    int offset = data.size() - 4;
    return (static_cast<uint32_t>(data[offset]) << 24) |
           (static_cast<uint32_t>(data[offset + 1]) << 16) |
           (static_cast<uint32_t>(data[offset + 2]) << 8)  |
           (static_cast<uint32_t>(data[offset + 3]));
}

// Bar
int calcTotalComparisons(int n) {
    int levels = static_cast<int>(std::log2(n));
    int total = 0;
    
    for (int i = 1; i <= levels; ++i) {
        total += (n >> 1) * i;
    }
    return total;
}

void normalCompAndSwap(std::vector<ByteVec>& A, int i, int j, bool ascending,
                 const ByteVec& normal_key) {
    int val_i = decodeInt(NormalDec(A[i], normal_key));
    int val_j = decodeInt(NormalDec(A[j], normal_key));
    bool need_swap = (ascending && val_i > val_j) || (!ascending && val_i < val_j);
    int mask = -static_cast<int>(need_swap); // 0xFFFFFFFF if true, 0x00000000 if false
    int tmp = val_i ^ val_j;
    tmp &= mask;
    val_i ^= tmp;
    val_j ^= tmp;
    A[i] = NormalEnc(encodeInt(val_i, BSize), normal_key, generateRandomBytes(16));
    A[j] = NormalEnc(encodeInt(val_j, BSize), normal_key, generateRandomBytes(16));

    // Update bar
    /*
    #pragma omp atomic
    ++curr_comps;
    int step = total_comps / 100;
    if (step == 0) step = 1;
    if (curr_comps % step == 0 || curr_comps == total_comps) {
        printProgressBar(static_cast<float>(curr_comps) / total_comps);
    }  
    */
}

// Iterative Bitonic Sort (non-parallel)
/*
void normalBitonicSortFromSequence(std::vector<ByteVec>& A, int startIndex, int lastIndex,
                                   bool ascending, const ByteVec& normal_key) {
    int noOfElements = lastIndex - startIndex + 1;

    for (int j = noOfElements / 2; j > 0; j /= 2) {
        int counter = 0;
        for (int i = startIndex; i + j <= lastIndex; ++i) {
            if (counter < j) {
                normalCompAndSwap(A, i, i + j, ascending, normal_key);
                ++counter;
            } else {
                counter = 0;
                i += j - 1;
            }
        }
    }
}

void normalBitonicSort(std::vector<ByteVec>& A, bool ascending, const ByteVec& normal_key) {
    int N = A.size();
    assert((N & (N - 1)) == 0);
    for (int j = 2; j <= N; j <<= 1) {
        #pragma omp parallel for
        for (int start = 0; start < N; start += j) {
            bool dir = (((start / j) % 2) == 0) ? ascending : !ascending;
            normalBitonicSortFromSequence(A, start, start + j - 1, dir, normal_key);
        }
    }
}
*/

void normalBitonicSort(std::vector<ByteVec>& A, bool ascending, const ByteVec& normal_key) {
    int N = A.size();
    assert((N & (N - 1)) == 0); // Power2
    for (int k = 2; k <= N; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < N; ++i) {
                int ixj = i ^ j;
                if (ixj > i) {
                    bool dir = ((i & k) == 0) ? ascending : !ascending;
                    normalCompAndSwap(A, i, ixj, dir, normal_key);
                }
            }
        }
    }
}

void anamorphiCompAndSwap(std::vector<ByteVec>& A, int i, int j,
                 bool normal_asc, bool covert_asc,
                 const ByteVec& normal_key, const ByteVec& covert_key,
                 std::vector<int>& write_count) {
    ByteVec iv_i = deriveIV(normal_key, i, write_count[i]);
    ByteVec iv_j = deriveIV(normal_key, j, write_count[j]);
    ++write_count[i];
    ++write_count[j];

    auto [ni, ci] = AnamorphicDec(A[i], normal_key, covert_key, iv_i);
    auto [nj, cj] = AnamorphicDec(A[j], normal_key, covert_key, iv_j);

    int val_i = decodeInt(ni); 
    int val_j = decodeInt(nj);
    int cov_i = decodeInt(ci);
    int cov_j = decodeInt(cj);

    // Swap normal values (branchless)
    bool normal_need_swap = (normal_asc && val_i > val_j) || (!normal_asc && val_i < val_j);
    int normal_mask = -static_cast<int>(normal_need_swap);
    int tmp = val_i ^ val_j;
    tmp &= normal_mask;
    val_i ^= tmp;
    val_j ^= tmp;

    // Swap covert values (branchless)
    bool covert_need_swap = (covert_asc && cov_i > cov_j) || (!covert_asc && cov_i < cov_j);
    int covert_mask = -static_cast<int>(covert_need_swap);
    tmp = cov_i ^ cov_j;
    tmp &= covert_mask;
    cov_i ^= tmp;
    cov_j ^= tmp;

    ByteVec new_iv_i = deriveIV(normal_key, i, write_count[i]);
    ByteVec new_iv_j = deriveIV(normal_key, j, write_count[j]);

    A[i] = AnamorphicEnc(encodeInt(val_i, BSize), normal_key, encodeInt(cov_i, 16), covert_key, new_iv_i);
    A[j] = AnamorphicEnc(encodeInt(val_j, BSize), normal_key, encodeInt(cov_j, 16), covert_key, new_iv_j);    
    
    // Update bar
    /*
    #pragma omp atomic
    ++curr_comps;
    int step = total_comps / 100;
    if (step == 0) step = 1;
    if (curr_comps % step == 0 || curr_comps == total_comps) {
        printProgressBar(static_cast<float>(curr_comps) / total_comps);
    }
    */ 
}

/*
void anamorphicBitonicSortFromSequence(std::vector<ByteVec>& A, int startIndex, int lastIndex,
                                       bool normal_asc, bool covert_asc,
                                       const ByteVec& normal_key, const ByteVec& covert_key,
                                       std::vector<int>& write_count) {
    int noOfElements = lastIndex - startIndex + 1;

    for (int j = noOfElements / 2; j > 0; j /= 2) {
        int counter = 0;
        for (int i = startIndex; i + j <= lastIndex; ++i) {
            if (counter < j) {
                anamorphiCompAndSwap(A, i, i + j, normal_asc, covert_asc, normal_key, covert_key, write_count);
                ++counter;
            } else {
                counter = 0;
                i += j - 1;
            }
        }
    }
}

void anamorphicBitonicSort(std::vector<ByteVec>& A, 
                                        bool normal_asc, bool covert_asc,
                                        const ByteVec& normal_key, const ByteVec& covert_key,
                                        std::vector<int>& write_count) {
    int N = A.size();

    for (int j = 2; j <= N; j <<= 1) {
        #pragma omp parallel for schedule(static)
        for (int start = 0; start < N; start += j) {
            bool dir_normal = (((start / j) % 2) == 0) ? normal_asc : !normal_asc;
            bool dir_covert = (((start / j) % 2) == 0) ? covert_asc : !covert_asc;

            anamorphicBitonicSortFromSequence(A, start, start + j - 1, dir_normal, dir_covert, normal_key, covert_key, write_count);
        }
    }
}

*/

void anamorphicBitonicSort(std::vector<ByteVec>& A, bool normal_asc, bool covert_asc, 
                        const ByteVec& normal_key, const ByteVec& covert_key,
                        std::vector<int>& write_count) {
    int N = A.size();
    assert((N & (N - 1)) == 0);  // N must be a power of 2
    for (int k = 2; k <= N; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < N; ++i) {
                int ixj = i ^ j;
                if (ixj > i) {
                    bool dir_normal = ((i & k) == 0) ? normal_asc : !normal_asc;
                    bool dir_covert = ((i & k) == 0) ? covert_asc : !covert_asc;
                    anamorphiCompAndSwap(A, i, ixj, dir_normal, dir_covert, normal_key, covert_key, write_count);
                }
            }
        }
    }
}

EVP_CIPHER_CTX* getThreadLocalPRFContext() {
    static thread_local EVP_CIPHER_CTX* ctx = [] {
        EVP_CIPHER_CTX* c = EVP_CIPHER_CTX_new();
        if (!c) throw std::runtime_error("Failed to alloc thread-local PRF context");
        if (EVP_EncryptInit_ex(c, EVP_aes_256_ecb(), nullptr, nullptr, nullptr) != 1)
            throw std::runtime_error("Failed to init PRF context");
        EVP_CIPHER_CTX_set_padding(c, 0);
        return c;
    }();
    return ctx;
}

ByteVec deriveIV(const ByteVec& key, int addr, int count) {
    ByteVec input(16, 0x00);
    input[0] = (addr >> 24) & 0xFF;
    input[1] = (addr >> 16) & 0xFF;
    input[2] = (addr >> 8) & 0xFF;
    input[3] = addr & 0xFF;
    input[4] = (count >> 24) & 0xFF;
    input[5] = (count >> 16) & 0xFF;
    input[6] = (count >> 8) & 0xFF;
    input[7] = count & 0xFF;

    ByteVec iv(16);  // Output = AES-ECB(key, input)
    int out_len = 0;

    EVP_CIPHER_CTX* ctx = getThreadLocalPRFContext();

    // Reset + reinit key
    if (EVP_EncryptInit_ex(ctx, nullptr, nullptr, key.data(), nullptr) != 1)
        throw std::runtime_error("PRF Reset failed");

    if (EVP_EncryptUpdate(ctx, iv.data(), &out_len, input.data(), input.size()) != 1)
        throw std::runtime_error("PRF Encrypt failed");

    return iv;
}

/*
ByteVec deriveIV(const ByteVec& key, int addr, int count) {
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) throw std::runtime_error("EVP_CIPHER_CTX_new failed");

    ByteVec input(16, 0x00);
    input[0] = (addr >> 24) & 0xFF;
    input[1] = (addr >> 16) & 0xFF;
    input[2] = (addr >> 8) & 0xFF;
    input[3] = addr & 0xFF;
    input[4] = (count >> 24) & 0xFF;
    input[5] = (count >> 16) & 0xFF;
    input[6] = (count >> 8) & 0xFF;
    input[7] = count & 0xFF;

    ByteVec iv(16); // 128-bit IV

    int out_len = 0;
    if (EVP_EncryptInit_ex(ctx, EVP_aes_256_ecb(), nullptr, key.data(), nullptr) != 1)
        throw std::runtime_error("PRF EncryptInit failed");

    EVP_CIPHER_CTX_set_padding(ctx, 0);

    if (EVP_EncryptUpdate(ctx, iv.data(), &out_len, input.data(), input.size()) != 1)
        throw std::runtime_error("PRF EncryptUpdate failed");

    EVP_CIPHER_CTX_free(ctx);

    return iv;
}
*/