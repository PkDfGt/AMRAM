#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cmath>
#include <omp.h>
#include "../AESN/AMAES.h"
#include "TC.h"
#include <openssl/evp.h>

using ByteVec = std::vector<unsigned char>;

int total_comps = 0;
int curr_comps = 0;
int ctr = 0;

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

constexpr int INT_BYTE_LEN = 4;

ByteVec encodeInt(int value, int mark, int size) {
    if (size < 5) throw std::invalid_argument("Size must be at least 5 bytes");

    ByteVec buffer(size, 0x00);
    // value
    buffer[size - 5] = (value >> 24) & 0xFF;
    buffer[size - 4] = (value >> 16) & 0xFF;
    buffer[size - 3] = (value >> 8) & 0xFF;
    buffer[size - 2] = value & 0xFF;

    // mark
    buffer[size - 1] = mark & 0xFF;

    return buffer;
}

std::pair<int, int> decodeInt(const ByteVec& data) {
    if (data.size() < 5)
        throw std::invalid_argument("Data too small to decode int+mark");

    int offset = data.size() - 5;
    int value = (static_cast<uint32_t>(data[offset]) << 24) |
                (static_cast<uint32_t>(data[offset + 1]) << 16) |
                (static_cast<uint32_t>(data[offset + 2]) << 8) |
                (static_cast<uint32_t>(data[offset + 3]));

    int mark = static_cast<int>(data[data.size() - 1]);

    return {value, mark};
}

int calcTotalComparisons(int n) {
    int levels = static_cast<int>(std::log2(n));
    int total = 0;
    for (int i = 1; i <= levels; ++i) {
        total += (n >> 1) * i;
    }
    return total;
}

int normalCompact(std::vector<ByteVec>& A, const ByteVec& normal_key,
                   int z, int offset, int len) {
    assert((len & (len - 1)) == 0);
    if (len == 2) {
        auto [xi, mi] = decodeInt(NormalDec(A[offset], normal_key));
        auto [xj, mj] = decodeInt(NormalDec(A[offset + 1], normal_key));

        int b = ((1 - mi) * mj) ^ z;
        int mask = -b;  // if b==1 → mask = -1 (all bits 1), if b==0 → mask = 0
        // Branchless swap
        int tmp_x = (xi ^ xj) & mask;
        xi ^= tmp_x;
        xj ^= tmp_x;
        int tmp_m = (mi ^ mj) & mask;
        mi ^= tmp_m;
        mj ^= tmp_m;

        A[offset] = NormalEnc(encodeInt(xi, mi, BSize), normal_key, generateRandomBytes(16));
        A[offset + 1] = NormalEnc(encodeInt(xj, mj, BSize), normal_key, generateRandomBytes(16));
        return mi + mj;
    }

    int half = len / 2;
    int mleft = normalCompact(A, normal_key, z % half, offset, half);
    int mright = normalCompact(A, normal_key, (z + mleft) % half, offset + half, half);
    int s = (((z % half) + mleft >= half) ? 1 : 0) ^ ((z >= half) ? 1 : 0);
    #pragma omp parallel for
    for (int i = 0; i < half; ++i) {
        int i1 = offset + i;
        int i2 = offset + i + half;

        auto [xi, mi] = decodeInt(NormalDec(A[i1], normal_key));
        auto [xj, mj] = decodeInt(NormalDec(A[i2], normal_key));

        int cond = s ^ (i >= ((z + mleft) % half));
        int mask = -cond;  // cond == 1 → mask = -1 (0xFFFFFFFF), cond == 0 → mask = 0

        // Branchless swap
        int tmp_x = (xi ^ xj) & mask;
        xi ^= tmp_x;
        xj ^= tmp_x;
        int tmp_m = (mi ^ mj) & mask;
        mi ^= tmp_m;
        mj ^= tmp_m;

        A[i1] = NormalEnc(encodeInt(xi, mi, BSize), normal_key, generateRandomBytes(16));
        A[i2] = NormalEnc(encodeInt(xj, mj, BSize), normal_key, generateRandomBytes(16));
    }

    return mleft + mright;
}

std::pair<int, int> anamorphicCompact(std::vector<ByteVec>& A,
                                      const ByteVec& normal_key,
                                      const ByteVec& covert_key,
                                      int nz,  // normal z
                                      int cz,  // covert z
                                      int offset,
                                      int len,
                                      std::vector<int>& write_count) {
    assert((len & (len - 1)) == 0);

    if (len == 2) {
        int i = offset;
        int j = offset + 1;

        ByteVec iv_i = deriveIV(normal_key, i, write_count[i]);
        ByteVec iv_j = deriveIV(normal_key, j, write_count[j]);

        auto [ni, ci] = AnamorphicDec(A[i], normal_key, covert_key, iv_i);
        auto [nj, cj] = AnamorphicDec(A[j], normal_key, covert_key, iv_j);
        ++write_count[i];
        ++write_count[j];

        auto [nv_i, nm_i] = decodeInt(ni);
        auto [nv_j, nm_j] = decodeInt(nj);
        auto [cv_i, cm_i] = decodeInt(ci);
        auto [cv_j, cm_j] = decodeInt(cj);

        int nb = ((1 - nm_i) * nm_j) ^ nz;
        int nb_mask = -nb;

        int tmp_nv = (nv_i ^ nv_j) & nb_mask;
        nv_i ^= tmp_nv;
        nv_j ^= tmp_nv;

        int tmp_nm = (nm_i ^ nm_j) & nb_mask;
        nm_i ^= tmp_nm;
        nm_j ^= tmp_nm;

        int cb = ((1 - cm_i) * cm_j) ^ cz;
        int cb_mask = -cb;

        int tmp_cv = (cv_i ^ cv_j) & cb_mask;
        cv_i ^= tmp_cv;
        cv_j ^= tmp_cv;

        int tmp_cm = (cm_i ^ cm_j) & cb_mask;
        cm_i ^= tmp_cm;
        cm_j ^= tmp_cm;


        ByteVec new_iv_i = deriveIV(normal_key, i, write_count[i]);
        ByteVec new_iv_j = deriveIV(normal_key, j, write_count[j]);

        A[i] = AnamorphicEnc(encodeInt(nv_i, nm_i, BSize), normal_key,
                             encodeInt(cv_i, cm_i, 16), covert_key, new_iv_i);
        A[j] = AnamorphicEnc(encodeInt(nv_j, nm_j, BSize), normal_key,
                             encodeInt(cv_j, cm_j, 16), covert_key, new_iv_j);

        return {nm_i + nm_j, cm_i + cm_j};
    }

    int half = len / 2;
    auto [ml1, cl1] = anamorphicCompact(A, normal_key, covert_key, nz % half, cz % half, offset, half, write_count);
    auto [ml2, cl2] = anamorphicCompact(A, normal_key, covert_key, (nz + ml1) % half, (cz + cl1) % half, offset + half, half, write_count);

    int normal_total = ml1 + ml2;
    int covert_total = cl1 + cl2;

    int sn = (((nz % half) + ml1 >= half) ? 1 : 0) ^ ((nz >= half) ? 1 : 0);
    int sc = (((cz % half) + cl1 >= half) ? 1 : 0) ^ ((cz >= half) ? 1 : 0);

    #pragma omp parallel for
    for (int i = 0; i < half; ++i) {
        int i1 = offset + i;
        int i2 = offset + i + half;

        ByteVec iv1 = deriveIV(normal_key, i1, write_count[i1]);
        ByteVec iv2 = deriveIV(normal_key, i2, write_count[i2]);

        auto [n1, c1] = AnamorphicDec(A[i1], normal_key, covert_key, iv1);
        auto [n2, c2] = AnamorphicDec(A[i2], normal_key, covert_key, iv2);
        ++write_count[i1];
        ++write_count[i2];

        auto [nv1, nm1] = decodeInt(n1);
        auto [nv2, nm2] = decodeInt(n2);
        auto [cv1, cm1] = decodeInt(c1);
        auto [cv2, cm2] = decodeInt(c2);

        int nb = sn ^ (i >= ((nz + ml1) % half));
        int nb_mask = -nb;

        int tmp_nv = (nv1 ^ nv2) & nb_mask;
        nv1 ^= tmp_nv;
        nv2 ^= tmp_nv;

        int tmp_nm = (nm1 ^ nm2) & nb_mask;
        nm1 ^= tmp_nm;
        nm2 ^= tmp_nm;

        int cb = sc ^ (i >= ((cz + cl1) % half));
        int cb_mask = -cb;

        int tmp_cv = (cv1 ^ cv2) & cb_mask;
        cv1 ^= tmp_cv;
        cv2 ^= tmp_cv;

        int tmp_cm = (cm1 ^ cm2) & cb_mask;
        cm1 ^= tmp_cm;
        cm2 ^= tmp_cm;

        ByteVec new_iv1 = deriveIV(normal_key, i1, write_count[i1]);
        ByteVec new_iv2 = deriveIV(normal_key, i2, write_count[i2]);

        A[i1] = AnamorphicEnc(encodeInt(nv1, nm1, BSize), normal_key,
                              encodeInt(cv1, cm1, 16), covert_key, new_iv1);
        A[i2] = AnamorphicEnc(encodeInt(nv2, nm2, BSize), normal_key,
                              encodeInt(cv2, cm2, 16), covert_key, new_iv2);
    }

    return {normal_total, covert_total};
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