#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <chrono>
#include <omp.h>
#include <fstream>
#include <filesystem>
#include "../AESN/AMAES.h"
#include "TC.h"

using ByteVec = std::vector<unsigned char>;
using namespace std;
using namespace std::chrono;

int BSize = 16;

bool verifyNormalCompact(const std::vector<ByteVec>& A, const ByteVec& key) {
    bool seen_zero = false;
    for (size_t i = 0; i < A.size(); ++i) {
        auto [val, mark] = decodeInt(NormalDec(A[i], key));
        if (seen_zero && mark == 1) {
            std::cerr << "NormalCompact: mark=1 after mark=0 at index " << i << std::endl;
            return false;
        }
        if (mark == 0) seen_zero = true;
    }
    return true;
}

bool verifyAnamorphicCompact(const std::vector<ByteVec>& A,
                             const ByteVec& normal_key,
                             const ByteVec& covert_key,
                             const std::vector<int>& write_count) {
    bool seen_normal_zero = false;
    bool seen_covert_zero = false;

    for (size_t i = 0; i < A.size(); ++i) {
        ByteVec iv = deriveIV(normal_key, i, write_count[i]);
        auto [n, c] = AnamorphicDec(A[i], normal_key, covert_key, iv);
        auto [nv, nm] = decodeInt(n);
        auto [cv, cm] = decodeInt(c);

        if (seen_normal_zero && nm == 1) {
            std::cerr << "AnamorphicCompact: normal mark=1 after normal mark=0 at index " << i << std::endl;
            return false;
        }
        if (seen_covert_zero && cm == 1) {
            std::cerr << "AnamorphicCompact: covert mark=1 after covert mark=0 at index " << i << std::endl;
            return false;
        }

        if (nm == 0) seen_normal_zero = true;
        if (cm == 0) seen_covert_zero = true;
    }

    return true;
}

void runTestCompact_VaryN_FixedB(const std::vector<int>& ns, int fixed_bsize,
                                  const std::vector<int>& thread_counts,
                                  const ByteVec& normal_key, const ByteVec& covert_key) {
    std::ofstream normal_file("Result/compact_normal_N_vary.csv");
    std::ofstream anam_file("Result/compact_anamorphic_N_vary.csv");
    normal_file << "N,threads,time_ms\n";
    anam_file << "N,threads,time_ms\n";

    for (int N_exp : ns) {
        int length = 1 << N_exp;
        BSize = fixed_bsize;
        for (int threads : thread_counts) {
            std::cout << "[Compact] N=2^" << N_exp << " (" << length
                      << "), Threads=" << threads << ", BSize=" << BSize << "\n";
            omp_set_num_threads(threads);

            // Normal data
            std::vector<ByteVec> normal_data(length);
            for (int i = 0; i < length; ++i) {
                int val = rand() % (2 * length + 1) - length;
                int mark = rand() % 2;
                normal_data[i] = NormalEnc(encodeInt(val, mark, BSize), normal_key, generateRandomBytes(16));
            }
            auto start = high_resolution_clock::now();
            normalCompact(normal_data, normal_key, 0, 0, length);
            auto end = high_resolution_clock::now();
            double time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
            normal_file << length << "," << threads << "," << fixed << setprecision(3) << time_ms << "\n";

            // Anamorphic data
            std::vector<ByteVec> anam_data(length);
            std::vector<int> write_count(length, 0);
            for (int i = 0; i < length; ++i) {
                int nv = rand() % (2 * length + 1) - length;
                int cv = rand() % (2 * length + 1) - length;
                int nm = rand() % 2;
                int cm = rand() % 2;
                ByteVec iv = deriveIV(normal_key, i, write_count[i]);
                anam_data[i] = AnamorphicEnc(encodeInt(nv, nm, BSize), normal_key,
                                             encodeInt(cv, cm, 16), covert_key, iv);
            }
            start = high_resolution_clock::now();
            anamorphicCompact(anam_data, normal_key, covert_key, 0, 0, 0, length, write_count);
            end = high_resolution_clock::now();
            time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
            anam_file << length << "," << threads << "," << fixed << setprecision(3) << time_ms << "\n";
            assert(verifyNormalCompact(normal_data, normal_key) && "NormalCompact failed");
            assert(verifyAnamorphicCompact(anam_data, normal_key, covert_key, write_count) && "AnamorphicCompact failed");
    
        }
    }

    normal_file.close();
    anam_file.close();
}


void runTestCompact_VaryB_FixedN(const std::vector<int>& bsizes, int fixed_n,
                                  const std::vector<int>& thread_counts,
                                  const ByteVec& normal_key, const ByteVec& covert_key) {
    std::ofstream normal_file("Result/compact_normal_B_N20_vary.csv");
    std::ofstream anam_file("Result/compact_anamorphic_B_N20_vary.csv");
    normal_file << "BSize,threads,time_ms\n";
    anam_file << "BSize,threads,time_ms\n";

    for (int bsize : bsizes) {
        BSize = bsize;
        for (int threads : thread_counts) {
            std::cout << "[Compact] BSize=" << bsize << ", Threads=" << threads
                      << ", N=" << fixed_n << "\n";
            omp_set_num_threads(threads);

            // Normal data
            std::vector<ByteVec> normal_data(fixed_n);
            for (int i = 0; i < fixed_n; ++i) {
                int val = rand() % (2 * fixed_n + 1) - fixed_n;
                int mark = rand() % 2;
                normal_data[i] = NormalEnc(encodeInt(val, mark, BSize), normal_key, generateRandomBytes(16));
            }
            auto start = high_resolution_clock::now();
            normalCompact(normal_data, normal_key, 0, 0, fixed_n);
            auto end = high_resolution_clock::now();
            double time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
            normal_file << bsize << "," << threads << "," << fixed << setprecision(3) << time_ms << "\n";

            // Anamorphic data
            std::vector<ByteVec> anam_data(fixed_n);
            std::vector<int> write_count(fixed_n, 0);
            for (int i = 0; i < fixed_n; ++i) {
                int nv = rand() % (2 * fixed_n + 1) - fixed_n;
                int cv = rand() % (2 * fixed_n + 1) - fixed_n;
                int nm = rand() % 2;
                int cm = rand() % 2;
                ByteVec iv = deriveIV(normal_key, i, write_count[i]);
                anam_data[i] = AnamorphicEnc(encodeInt(nv, nm, BSize), normal_key,
                                             encodeInt(cv, cm, 16), covert_key, iv);
            }
            start = high_resolution_clock::now();
            anamorphicCompact(anam_data, normal_key, covert_key, 0, 0, 0, fixed_n, write_count);
            end = high_resolution_clock::now();
            time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
            anam_file << bsize << "," << threads << "," << fixed << setprecision(3) << time_ms << "\n";
            assert(verifyNormalCompact(normal_data, normal_key) && "NormalCompact failed");
            assert(verifyAnamorphicCompact(anam_data, normal_key, covert_key, write_count) && "AnamorphicCompact failed");
    
        }
    }

    normal_file.close();
    anam_file.close();
}

int main() {
    srand(static_cast<unsigned>(time(nullptr)));
    std::filesystem::create_directory("Result");

    ByteVec normal_key(32, 0x0f);
    ByteVec covert_key(32, 0x01);

    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 32};

    // Test 1: Fix BSize = 16, Vary N = 2^15~2^20
    std::vector<int> ns = {15, 16, 17, 18, 19, 20};
    int fixed_bsize = 16;
    runTestCompact_VaryN_FixedB(ns, fixed_bsize, thread_counts, normal_key, covert_key);

    // Test 2: Fix N = 2^20, Vary BSize
    int fixed_n = 1 << 20;
    std::vector<int> bsizes = {1024 * 4, 1024, 256, 64, 16};
    runTestCompact_VaryB_FixedN(bsizes, fixed_n, thread_counts, normal_key, covert_key);

    return 0;
}
