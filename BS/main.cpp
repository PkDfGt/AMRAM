#include <fstream>
#include <iomanip>  // for std::setprecision
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <iostream>
#include <chrono>
#include <omp.h>
#include "../AESN/AMAES.h"
#include "BS.h"

using ByteVec = std::vector<unsigned char>;
using namespace std;
using namespace std::chrono;
int BSize =  16;

bool verifySortedNormal(const std::vector<ByteVec>& A, const ByteVec& key, bool ascending) {
    for (size_t i = 1; i < A.size(); ++i) {
        int prev = decodeInt(NormalDec(A[i - 1], key));
        int curr = decodeInt(NormalDec(A[i], key));
        if ((ascending && prev > curr) || (!ascending && prev < curr)) return false;
    }
    return true;
}

bool verifySortedAnamorphic(const std::vector<ByteVec>& A,
                            const ByteVec& normal_key,
                            const ByteVec& covert_key,
                            const std::vector<int>& write_count,
                            bool normal_asc,
                            bool covert_asc) {
    for (size_t i = 1; i < A.size(); ++i) {
        ByteVec iv1 = deriveIV(normal_key, i - 1, write_count[i - 1]);
        ByteVec iv2 = deriveIV(normal_key, i, write_count[i]);

        auto [n1, c1] = AnamorphicDec(A[i - 1], normal_key, covert_key, iv1);
        auto [n2, c2] = AnamorphicDec(A[i], normal_key, covert_key, iv2);

        int norm1 = decodeInt(n1), norm2 = decodeInt(n2);
        int cov1 = decodeInt(c1), cov2 = decodeInt(c2);
        // std::cout << norm1 << " " << cov1 << "\n";

        if ((normal_asc && norm1 > norm2) || (!normal_asc && norm1 < norm2)) return false;
        if ((covert_asc && cov1 > cov2) || (!covert_asc && cov1 < cov2)) return false;
    }
    return true;
}

void runTestVaryN_FixedB_SizeThread(const std::vector<int>& ns, int fixed_bsize,
                                   const std::vector<int>& thread_counts,
                                   const ByteVec& normal_key, const ByteVec& covert_key) {
    std::ofstream normal_file("Result/normal_N_vary.csv");
    std::ofstream anam_file("Result/anamorphic_N_vary.csv");
    normal_file << "N,threads,time_ms\n";
    anam_file << "N,threads,time_ms\n";

    for (int N_exp : ns) {
        int length = 1 << N_exp;
        BSize = fixed_bsize;
        for (int threads : thread_counts) {
            std::cout << "[BS] N=2^" << N_exp << " (" << length
                      << "), Threads=" << threads << ", BSize=" << BSize << "\n";
            omp_set_num_threads(threads);

            // Anamorphic
            std::vector<ByteVec> anam_data(length);
            std::vector<int> write_count(length, 0);
            for (int i = 0; i < length; ++i) {
                int norm_val = rand() % (2 * length + 1) - length;
                int cov_val = rand() % (2 * length + 1) - length;
                anam_data[i] = AnamorphicEnc(encodeInt(norm_val, BSize), normal_key,
                                            encodeInt(cov_val, 16), covert_key,
                                            deriveIV(normal_key, i, 0));
            }
            total_comps = calcTotalComparisons(length);
            curr_comps = 0;
            auto start2 = std::chrono::high_resolution_clock::now();
            anamorphicBitonicSort(anam_data, true, false, normal_key, covert_key, write_count);
            auto end2 = std::chrono::high_resolution_clock::now();
            double anam_time = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000.0;
            // std::cout << "Anamorphic Sort time: " << anam_time << " ms\n";
            anam_file << length << "," << threads << "," << std::fixed << std::setprecision(3) << anam_time << "\n";

            // Normal
            std::vector<ByteVec> normal_data(length);
            srand(static_cast<unsigned>(time(nullptr)));
            for (int i = 0; i < length; ++i) {
                int val = rand() % (2 * length + 1) - length;
                normal_data[i] = NormalEnc(encodeInt(val, BSize), normal_key, generateRandomBytes(16));
            }
            total_comps = calcTotalComparisons(length);
            curr_comps = 0;
            auto start = std::chrono::high_resolution_clock::now();
            normalBitonicSort(normal_data, true, normal_key);
            auto end = std::chrono::high_resolution_clock::now();
            double normal_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
            // std::cout << "Normal Sort time: " << normal_time << " ms\n";
            normal_file << length << "," << threads << "," << std::fixed << std::setprecision(3) << normal_time << "\n";
            // assert(verifySortedNormal(normal_data, normal_key, true) && "Normal sort failed");
            // assert(verifySortedAnamorphic(anam_data, normal_key, covert_key, write_count, true, false) && "Anamorphic sort failed");
        }
    }
    normal_file.close();
    anam_file.close();
}

void runTestVaryB_FixedN_Thread(const std::vector<int>& bsizes, int fixed_n,
                               const std::vector<int>& thread_counts,
                               const ByteVec& normal_key, const ByteVec& covert_key) {
    std::ofstream normal_file("Result/normal_B_N20_vary.csv");
    std::ofstream anam_file("Result/anamorphic_B_N20_vary.csv");
    normal_file << "BSize,threads,time_ms\n";
    anam_file << "BSize,threads,time_ms\n";

    for (int bsize : bsizes) {
        for (int threads : thread_counts) {
            std::cout << "[BS] BSize=" << bsize << ", Threads=" << threads
                      << ", N=" << fixed_n << "\n";
            BSize = bsize;
            omp_set_num_threads(threads);

            // Anamorphic
            std::vector<ByteVec> anam_data(fixed_n);
            std::vector<int> write_count(fixed_n, 0);
            for (int i = 0; i < fixed_n; ++i) {
                int norm_val = rand() % (2 * fixed_n + 1) - fixed_n;
                int cov_val = rand() % (2 * fixed_n + 1) - fixed_n;
                anam_data[i] = AnamorphicEnc(encodeInt(norm_val, BSize), normal_key,
                                            encodeInt(cov_val, 16), covert_key,
                                            deriveIV(normal_key, i, 0));
            }
            total_comps = calcTotalComparisons(fixed_n);
            curr_comps = 0;
            auto start2 = std::chrono::high_resolution_clock::now();
            anamorphicBitonicSort(anam_data, true, false, normal_key, covert_key, write_count);
            auto end2 = std::chrono::high_resolution_clock::now();
            double anam_time = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2).count() / 1000.0;
            std::cout << "Anamorphic Sort time: " << anam_time << " ms\n";
            anam_file << bsize << "," << threads << "," << std::fixed << std::setprecision(3) << anam_time << "\n";

            // Normal
            std::vector<ByteVec> normal_data(fixed_n);
            srand(static_cast<unsigned>(time(nullptr)));
            for (int i = 0; i < fixed_n; ++i) {
                int val = rand() % (2 * fixed_n + 1) - fixed_n;
                normal_data[i] = NormalEnc(encodeInt(val, BSize), normal_key, generateRandomBytes(16));
            }
            total_comps = calcTotalComparisons(fixed_n);
            curr_comps = 0;
            auto start = std::chrono::high_resolution_clock::now();
            normalBitonicSort(normal_data, true, normal_key);
            auto end = std::chrono::high_resolution_clock::now();
            double normal_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
            std::cout << "Normal Sort time: " << normal_time << " ms\n";
            normal_file << bsize << "," << threads << "," << std::fixed << std::setprecision(3) << normal_time << "\n";
            // assert(verifySortedNormal(normal_data, normal_key, true) && "Normal sort failed");
            // assert(verifySortedAnamorphic(anam_data, normal_key, covert_key, write_count, true, false) && "Anamorphic sort failed");
        }
    }
    normal_file.close();
    anam_file.close();
}

int main() {
    ByteVec normal_key(32, 0x0f);
    ByteVec covert_key(32, 0x01);

    // Vary ThreadCount
    std::vector<int> thread_counts = {1, 2, 4, 8, 16, 32};
    // Vary BlockSize
    std::vector<int> bsizes = {1024 * 4, 1024, 256, 64, 16};
    int fixed_n = 1 << 20;
    runTestVaryB_FixedN_Thread(bsizes, fixed_n, thread_counts, normal_key, covert_key);
    // Vary DBSize
    std::vector<int> ns = {15, 16, 17, 18, 19, 20};
    int fixed_bsize = 16;
    runTestVaryN_FixedB_SizeThread(ns, fixed_bsize, thread_counts, normal_key, covert_key);
    return 0;
    
}

/*
void runTestDBSizes(const std::vector<int>& sizes, const ByteVec& normal_key, const ByteVec& covert_key) {
    std::ofstream normal_file("normal.csv");
    std::ofstream anam_file("anamorphic.csv");
    normal_file << "size,time_ms\n";
    anam_file << "size,time_ms\n";

    for (int power2 : sizes) {
        int length = 1 << power2;
        std::cout << "Testing size = 2^" << power2 << " = " << length << std::endl;
        omp_set_num_threads(8); 

        // Anamorphic
        std::vector<ByteVec> anam_data(length);
        std::vector<int> write_count(length, 0);
        
        for (int i = 0; i < length; ++i) {
            int norm_val = rand() % (2 * length + 1) - length;
            int cov_val = rand() % (2 * length + 1) - length;
            anam_data[i] = AnamorphicEnc(encodeInt(norm_val,BSize), normal_key, encodeInt(cov_val,BSize), covert_key, deriveIV(normal_key, i, 0));
        }

        total_comps = calcTotalComparisons(length);
        curr_comps = 0;
        auto start2 = high_resolution_clock::now();
        anamorphicBitonicSort(anam_data, true, false, normal_key, covert_key, write_count);
        auto end2 = high_resolution_clock::now();
        double anam_time = duration_cast<microseconds>(end2 - start2).count() / 1000.0;
        std::cout << "Anamorphic Sort time: " << anam_time << " ms\n";
        anam_file << length << "," << std::fixed << std::setprecision(3) << anam_time << "\n";

        // Normal
        std::vector<ByteVec> normal_data(length);
        srand(static_cast<unsigned>(time(nullptr)));
        for (int i = 0; i < length; ++i) {
            int val = rand() % (2 * length + 1) - length;
            normal_data[i] = NormalEnc(encodeInt(val,BSize), normal_key, generateRandomBytes(16));
        }

        total_comps = calcTotalComparisons(length);
        curr_comps = 0;
        auto start = high_resolution_clock::now();
        normalBitonicSort(normal_data, true, normal_key);
        auto end = high_resolution_clock::now();
        double normal_time = duration_cast<microseconds>(end - start).count() / 1000.0;
        std::cout << "Normal Sort time: " << normal_time << " ms\n";
        normal_file << length << "," << std::fixed << std::setprecision(3) << normal_time << "\n";
    }
    normal_file.close();
    anam_file.close();
}


void runTestBlockSizes(int length, const std::vector<int>& bsizes, const ByteVec& normal_key, const ByteVec& covert_key) {
    std::ofstream normal_file("Result/normal_BSize.csv");
    std::ofstream anam_file("Result/anamorphic_BSize.csv");
    normal_file << "BSize,time_ms\n";
    anam_file << "BSize,time_ms\n";

    for (int bsize : bsizes) {
        std::cout << "Testing BSize = " << bsize << std::endl;
        BSize = bsize;
        omp_set_num_threads(8);

        // anamorphic
        std::vector<ByteVec> anam_data(length);
        std::vector<int> write_count(length, 0);
        for (int i = 0; i < length; ++i) {
            int norm_val = rand() % (2 * length + 1) - length;
            int cov_val = rand() % (2 * length + 1) - length;
            anam_data[i] = AnamorphicEnc(encodeInt(norm_val, BSize), normal_key, encodeInt(cov_val, BSize), covert_key, deriveIV(normal_key, i, 0));
        }
        total_comps = calcTotalComparisons(length);
        curr_comps = 0;
        auto start2 = high_resolution_clock::now();
        anamorphicBitonicSort(anam_data, true, false, normal_key, covert_key, write_count);
        auto end2 = high_resolution_clock::now();
        double anam_time = duration_cast<microseconds>(end2 - start2).count() / 1000.0;
        std::cout << "Anamorphic Sort time: " << anam_time << " ms\n";
        anam_file << bsize << "," << std::fixed << std::setprecision(3) << anam_time << "\n";

        // normal
        std::vector<ByteVec> normal_data(length);
        srand(static_cast<unsigned>(time(nullptr)));
        for (int i = 0; i < length; ++i) {
            int val = rand() % (2 * length + 1) - length;
            normal_data[i] = NormalEnc(encodeInt(val, BSize), normal_key, generateRandomBytes(16));
        }
        total_comps = calcTotalComparisons(length);
        curr_comps = 0;
        auto start = high_resolution_clock::now();
        normalBitonicSort(normal_data, true, normal_key);
        auto end = high_resolution_clock::now();
        double normal_time = duration_cast<microseconds>(end - start).count() / 1000.0;
        std::cout << "Normal Sort time: " << normal_time << " ms\n";
        normal_file << bsize << "," << std::fixed << std::setprecision(3) << normal_time << "\n";

        // assert(verifySortedNormal(normal_data, normal_key, true) && "Normal sort failed");
        // assert(verifySortedAnamorphic(anam_data, normal_key, covert_key, write_count, true, false) && "Anamorphic sort failed");
    }

    normal_file.close();
    anam_file.close();
}

void runTestThreadCounts(int length,
                         const ByteVec& normal_key, const ByteVec& covert_key,
                         const std::vector<int>& thread_counts) {
    std::ofstream normal_file("Result/normal_threads.csv");
    std::ofstream anam_file("Result/anamorphic_threads.csv");
    normal_file << "threads,time_ms\n";
    anam_file << "threads,time_ms\n";

    for (int threads : thread_counts) {
        std::cout << "Testing threads = " << threads << std::endl;
        omp_set_num_threads(threads);


        // anamorphic
        std::vector<ByteVec> anam_data(length);
        std::vector<int> write_count(length, 0);
        for (int i = 0; i < length; ++i) {
            int norm_val = rand() % (2 * length + 1) - length;
            int cov_val = rand() % (2 * length + 1) - length;
            // std::cout << norm_val << " " << cov_val << "\n";
            anam_data[i] = AnamorphicEnc(encodeInt(norm_val, BSize), normal_key,
                                         encodeInt(cov_val, BSize), covert_key,
                                        deriveIV(normal_key, i, 0));
        }

        total_comps = calcTotalComparisons(length);
        curr_comps = 0;
        auto start2 = high_resolution_clock::now();
        anamorphicBitonicSort(anam_data, true, false, normal_key, covert_key, write_count);
        auto end2 = high_resolution_clock::now();
        double anam_time = duration_cast<microseconds>(end2 - start2).count() / 1000.0;
        std::cout << "Anamorphic Sort time: " << anam_time << " ms\n";
        anam_file << threads << "," << std::fixed << std::setprecision(3) << anam_time << "\n"; 

        // normal
        std::vector<ByteVec> normal_data(length);
        srand(static_cast<unsigned>(time(nullptr)));
        for (int i = 0; i < length; ++i) {
            int val = rand() % (2 * length + 1) - length;
            normal_data[i] = NormalEnc(encodeInt(val, BSize), normal_key, generateRandomBytes(16));
        }
        total_comps = calcTotalComparisons(length);
        curr_comps = 0;
        auto start = high_resolution_clock::now();
        normalBitonicSort(normal_data, true, normal_key);
        auto end = high_resolution_clock::now();
        double normal_time = duration_cast<microseconds>(end - start).count() / 1000.0;
        std::cout << "Normal Sort time: " << normal_time << " ms\n";
        normal_file << threads << "," << std::fixed << std::setprecision(3) << normal_time << "\n";
        
        // assert(verifySortedNormal(normal_data, normal_key, true) && "Normal sort failed");
        // assert(verifySortedAnamorphic(anam_data, normal_key, covert_key, write_count, true, false) && "Anamorphic sort failed");
    }

    normal_file.close();
    anam_file.close();
}

int main() {
    ByteVec normal_key(32, 0x0f);
    ByteVec covert_key(32, 0x01);
    // Test DBSize, with #threads = 8, BlockSize = 16 bytes
    // std::vector<int> sizes = {5};//, 6};//, 7, 8, 9, 10};
    // std::vector<int> DBlog2sizes = {15, 16, 17, 18, 19, 20};
    // runTestDBSizes(DBlog2sizes, normal_key, covert_key);

    // Test BlockSize，with #threads = 8, DBSize = 2^15
    // int length = 1 << 5;
    // std::vector<int> bsizes = {256};//, 64, 16};
    // std::vector<int> bsizes = {1024*256, 1024*64, 1024*16, 1024*4, 1024, 256, 64, 16};
    // runTestBlockSizes(length, bsizes, normal_key, covert_key);

    // Test #threads, with DBSize = 2^15, BlockSize = 16 bytes
    int thread_test_len = 1 << 15;
    std::vector<int> thread_counts = {1, 2, 4, 8, 16}; //, 32};
    runTestThreadCounts(thread_test_len, normal_key, covert_key, thread_counts);
    return 0;
}
*/