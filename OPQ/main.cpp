#include "OPQ.h"
#include "PQ.h"
#include "../AESN/AMAES.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <numeric>
#include <cassert>
#include <fstream>
#include <filesystem>
#include <omp.h>

using namespace std::chrono;

void print_progress_bar(int current, int total, const std::string& phase) {
    const int bar_width = 50;
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(bar_width * progress);

    std::cout << "\r[" << phase << "] [";
    for (int i = 0; i < bar_width; ++i)
        std::cout << (i < pos ? '=' : (i == pos ? '>' : ' '));
    std::cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")" << std::flush;

    if (current == total) std::cout << "\n";
}

void benchmark_pq_power2(const std::string& label) {
    const int fixed_bsize = 16;
    ByteVec key(32, 0x01);
    std::filesystem::create_directory("Result");
    std::ofstream fout("Result/" + label + "_power2.csv");
    fout << "power2,N,avg_insert_ms,avg_extract_ms\n";

    for (int power2 = 15; power2 <= 20; ++power2) {
        int N = 1 << power2;
        BinaryHeapPQ pq(N, fixed_bsize, key);

        std::vector<int> prios(N);
        std::iota(prios.begin(), prios.end(), 2);
        std::shuffle(prios.begin(), prios.end(), std::mt19937{std::random_device{}()});

        auto start_insert = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Insert");
            pq.insert(prios[i], 100 + i);
        }
        auto end_insert = std::chrono::high_resolution_clock::now();

        std::sort(prios.begin(), prios.end());

        auto start_extract = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Extract");
            PlainEle e = pq.extract_min();
            assert(e.p == prios[i]);
        }
        auto end_extract = std::chrono::high_resolution_clock::now();

        double insert_ms = std::chrono::duration<double, std::milli>(end_insert - start_insert).count() / N;
        double extract_ms = std::chrono::duration<double, std::milli>(end_extract - start_extract).count() / N;

        fout << power2 << "," << N << "," << insert_ms << "," << extract_ms << "\n";
        std::cout << "power2 = " << power2 << " (N = " << N << ")"
                  << ", avg insert = " << insert_ms << " ms"
                  << ", avg extract = " << extract_ms << " ms\n";
    }

    fout.close();
}

void benchmark_pq_bsize(const std::string& label) {
    const int power2 = 20;
    const int N = 1 << power2;
    std::vector<int> bsizes = {1024 * 4, 1024, 256, 64, 16};
    ByteVec key(32, 0x01);

    std::filesystem::create_directory("Result");
    std::ofstream fout("Result/N20_" + label + "_bsize.csv");
    fout << "bsize,N,avg_insert_ms,avg_extract_ms\n";

    for (int bsize : bsizes) {
        BinaryHeapPQ pq(N, bsize, key);

        std::vector<int> prios(N);
        std::iota(prios.begin(), prios.end(), 2);
        std::shuffle(prios.begin(), prios.end(), std::mt19937{std::random_device{}()});

        auto start_insert = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Insert");
            pq.insert(prios[i], 100 + i);
        }
        auto end_insert = std::chrono::high_resolution_clock::now();

        std::sort(prios.begin(), prios.end());

        auto start_extract = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Extract");
            PlainEle e = pq.extract_min();
            assert(e.p == prios[i]);
        }
        auto end_extract = std::chrono::high_resolution_clock::now();

        double insert_ms = std::chrono::duration<double, std::milli>(end_insert - start_insert).count() / N;
        double extract_ms = std::chrono::duration<double, std::milli>(end_extract - start_extract).count() / N;

        fout << bsize << "," << N << "," << insert_ms << "," << extract_ms << "\n";
        std::cout << "bsize = " << bsize
                  << ", avg insert = " << insert_ms << " ms"
                  << ", avg extract = " << extract_ms << " ms\n";
    }

    fout.close();
}

void benchmark_opq(int mode, const std::string& label) {
    const int fixed_bsize = 16;
    std::cout << "\n===== [" << label << " Mode: Varying power2, fixed bsize=16] =====\n";
    ByteVec key1(32, 0x01);
    ByteVec key2(32, 0x02);
    ByteVec ivkey(32, 0x02);

    std::filesystem::create_directory("Result");
    std::ofstream fout("Result/" + label + "_power2.csv");
    fout << "power2,N,avg_insert_ms,avg_extract_ms\n";
    omp_set_num_threads(8);

    for (int power2 = 15; power2 <= 20; ++power2) {
        int N = 1 << power2;
        int Z = 4, Zr = power2;
        OPQ opq(N, Z, Zr, fixed_bsize, mode, key1, key2, ivkey);

        std::vector<int> p1(N), p2(N);
        std::iota(p1.begin(), p1.end(), 2);
        std::iota(p2.begin(), p2.end(), 2);
        std::shuffle(p1.begin(), p1.end(), std::mt19937{std::random_device{}()});
        std::shuffle(p2.begin(), p2.end(), std::mt19937{std::random_device{}()});

        auto start_insert = high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Insert");
            if (mode == 0) {
                opq.insert(p1[i], 10 + i, p1[i] % N);
            } else {
                opq.ainsert(p1[i], 10 + i, p1[i] % N, p2[i], 20 + i, p2[i] % N);
            }
        }
        auto end_insert = high_resolution_clock::now();

        std::sort(p1.begin(), p1.end());

        auto start_extract = high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Extract");
            if (mode == 0) {
                auto ele = opq.extract_min();
                assert(ele.p == p1[i]);
            } else {
                auto [ne, ae] = opq.aextract_min();
                assert(ne.p == p1[i] || ae.p == p1[i]);
            }
        }
        auto end_extract = high_resolution_clock::now();

        double insert_ms = duration<double, std::milli>(end_insert - start_insert).count() / N;
        double extract_ms = duration<double, std::milli>(end_extract - start_extract).count() / N;

        fout << power2 << "," << N << "," << insert_ms << "," << extract_ms << "\n";

        std::cout << "power2 = " << power2 << " (N = " << N << ")"
                  << ", avg insert = " << insert_ms << " ms"
                  << ", avg extract = " << extract_ms << " ms\n";
    }

    fout.close();
}

void benchmark_bsize_opq(int mode, const std::string& label) {
    std::cout << "\n===== [" << label << " Mode: Varying bsize, fixed power2=15] =====\n";
    std::vector<int> bsizes = {
        1024 * 4, 1024, 256, 64, 16
    };
    int power2 = 20;
    int N = 1 << power2;
    int Z = 4, Zr = power2;
    ByteVec key1(32, 0x01);
    ByteVec key2(32, 0x02);
    ByteVec ivkey(32, 0x02);

    std::filesystem::create_directory("Result");
    std::ofstream fout("Result/N20_" + label + "_bsize.csv");
    fout << "bsize,N,avg_insert_ms,avg_extract_ms\n";

    for (int bsize : bsizes) {
        if (mode == 1 && bsize < 16) continue;

        OPQ opq(N, Z, Zr, bsize, mode, key1, key2, ivkey);

        std::vector<int> p1(N), p2(N);
        std::iota(p1.begin(), p1.end(), 2);
        std::iota(p2.begin(), p2.end(), 2);
        std::shuffle(p1.begin(), p1.end(), std::mt19937{std::random_device{}()});
        std::shuffle(p2.begin(), p2.end(), std::mt19937{std::random_device{}()});

        auto start_insert = high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Insert");
            if (mode == 0) {
                opq.insert(p1[i], 10 + i, p1[i] % N);
            } else {
                opq.ainsert(p1[i], 10 + i, p1[i] % N, p2[i], 20 + i, p2[i] % N);
            }
        }
        auto end_insert = high_resolution_clock::now();

        std::sort(p1.begin(), p1.end());

        auto start_extract = high_resolution_clock::now();
        for (int i = 0; i < N; ++i) {
            print_progress_bar(i + 1, N, "Extract");
            if (mode == 0) {
                auto ele = opq.extract_min();
                assert(ele.p == p1[i]);
            } else {
                auto [ne, ae] = opq.aextract_min();
                assert(ne.p == p1[i] || ae.p == p1[i]);
            }
        }
        auto end_extract = high_resolution_clock::now();

        double insert_ms = duration<double, std::milli>(end_insert - start_insert).count() / N;
        double extract_ms = duration<double, std::milli>(end_extract - start_extract).count() / N;

        fout << bsize << "," << N << "," << insert_ms << "," << extract_ms << "\n";

        std::cout << "bsize = " << bsize
                  << ", avg insert = " << insert_ms << " ms"
                  << ", avg extract = " << extract_ms << " ms\n";
    }

    fout.close();
}

int main() {
    benchmark_opq(0, "Normal");
    benchmark_opq(0, "Anamorphic");
    benchmark_pq_power2("HeapPQ");
    
    benchmark_bsize_opq(0, "Normal");
    benchmark_bsize_opq(1, "Anamorphic");
    benchmark_pq_bsize("HeapPQ");
    return 0;
}


/*

#include "OPQ.h"
#include "../AESN/AMAES.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <tuple>
#include <cassert>

void test_normal_opq() {
    const int power2 = 10;
    const int N = 1 << power2;
    const int Z = 4;
    const int Zr = power2;
    const int bsize = 64;
    const int mode = 0;
    ByteVec key1(32, 0x01);
    ByteVec key2(32, 0x02);
    ByteVec ivkey(32, 0x02);
    OPQ opq(N, Z, Zr, bsize, mode, key1, key2, ivkey);

    const int m = N;
    std::vector<int> priorities(m);
    for (int i = 0; i < m; ++i) priorities[i] = i + 2;
    std::shuffle(priorities.begin(), priorities.end(), std::mt19937{std::random_device{}()});

    std::cout << "Normal Mode: Shuffled priorities:\n";
    // for (int i = 0; i < m; ++i) std::cout << priorities[i] << " ";
    // std::cout << "\n";

    for (int i = 0; i < m; ++i) {
        opq.insert(priorities[i], 10 + i, priorities[i] % N);
    }
    // opq.print_tree();

    std::sort(priorities.begin(), priorities.end()); // expect ascending extract

    for (int i = 0; i < m; ++i) {
        auto ele = opq.extract_min();
        if (ele.p != priorities[i]) {
            std::cerr << "❌ Normal Error: expected p=" << priorities[i] << ", got p=" << ele.p << "\n";
            assert(false);
        }
    }

    std::cout << "✅ Normal test passed: all priorities extracted in correct order.\n\n";
}

void test_anamorphic_opq() {
    const int power2 = 10;
    const int N = 1 << power2;
    const int Z = 4;
    const int Zr = power2;
    const int bsize = 64;
    const int mode = 1;
    ByteVec key1(32, 0x01);
    ByteVec key2(32, 0x02);
    ByteVec ivkey(32, 0x02);
    OPQ opq(N, Z, Zr, bsize, mode, key1, key2, ivkey);

    const int m = N;
    std::vector<int> priorities1(m);
    std::vector<int> priorities2(m);
    for (int i = 0; i < m; ++i) priorities1[i] = i + 2;
    for (int i = 0; i < m; ++i) priorities2[i] = i + 2;
    std::shuffle(priorities1.begin(), priorities1.end(), std::mt19937{std::random_device{}()});
    std::shuffle(priorities2.begin(), priorities2.end(), std::mt19937{std::random_device{}()});

    std::cout << "Anamorphic Mode: Shuffled priorities:\n";
    // for (int i = 0; i < m; ++i) std::cout << priorities[i] << " ";
    // std::cout << "\n";

    for (int i = 0; i < m; ++i) {
        int np = priorities1[i];
        int ap = priorities2[i];
        int nv = 10 + i;
        int av = 20 + i;
        int npos = np % N;
        int apos = ap % N;
        opq.ainsert(np, nv, npos, ap, av, apos);
    }
    // opq.print_anamorphic_tree();

    std::sort(priorities1.begin(), priorities1.end()); // expect ascending extract

    for (int i = 0; i < m; ++i) {
        auto [ne, ae] = opq.aextract_min();
        if (ne.p != priorities1[i] && ae.p != priorities1[i]) {
            std::cerr << "❌ Anamorphic Error: expected p=" << priorities1[i]
                      << ", got normal=" << ne.p << ", anamorphic=" << ae.p << "\n";
            assert(false);
        }       
    }
    std::cout << "✅ Anamorphic test passed: all priorities extracted in correct order.\n\n";
}

int main() {
    test_normal_opq();
    test_anamorphic_opq();
    return 0;
}

*/