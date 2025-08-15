#include "OORAM.h"
#include "Utils.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cassert>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <thread>


namespace fs = std::filesystem;
using ByteVec = std::vector<unsigned char>;

// benchmark_access: dbsize
void benchmark_access(int dbsize, int bsize) {
    fs::create_directories("Result");
    std::string filename = "Result/ooram_BS_dbsize_" + std::to_string(dbsize) + "_bsize_" + std::to_string(bsize) + ".csv";

    std::ofstream fout(filename);
    fout << "dbsize,acsize,init_ms,access_ms,total_ms\n";
    
    ByteVec normal_key(32, 0xAA), covert_key(32, 0xBB),
            iv_key(32, 0xCC), pos_key(32, 0xDD);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> distDATA(0, 1 << 10);

    std::vector<int> DB(dbsize);
    std::vector<int> AS;
    for (int i = 0; i < dbsize; ++i) {
        DB[i] = distDATA(gen);
        // std::cout << "DB[i]: " << DB[i] << "\n";
    }
    generateBSAddr(dbsize, AS);

    int Z = 4;
    int acsize = AS.size();
    int Zr = 2 * static_cast<int>(std::ceil(std::log2(dbsize + acsize)));
    int  mode = 0;
    // int power2 = static_cast<int>(std::ceil(std::log2(dbsize)));

    OORAM ooram(dbsize, acsize, Z, Zr, bsize, mode,
                normal_key, covert_key, iv_key, pos_key);

    auto t0 = std::chrono::high_resolution_clock::now();
    ooram.init(DB, AS);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<int> cur_vals = DB;
    auto t2 = std::chrono::high_resolution_clock::now();
    size_t real_actimes = 100;
    for (size_t i = 0; i < real_actimes; ++i) {
        print_progress(i + 1, real_actimes, "Access");
        ooram.access();
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << std::endl;
    // ooram.verify_BS_result();

    auto init_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto access_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    fout << dbsize << "," << real_actimes << "," << init_ms << "," << access_ms << "," << (init_ms + access_ms) << "\n";

    fout.close();
    std::cout << "Finished benchmark_access, saved to " << filename << std::endl;
}

// benchmark_aaccess: dbsize
void benchmark_aaccess(int dbsize, int bsize) {
    fs::create_directories("Result");
    std::string filename = "Result/ooram_BSTC_dbsize_" + std::to_string(dbsize) + "_bsize_" + std::to_string(bsize) + ".csv";

    std::ofstream fout(filename);
    fout << "dbsize,acsize,init_ms,access_ms,total_ms\n";

    ByteVec normal_key(32, 0xAA), covert_key(32, 0xBB),
            iv_key(32, 0xCC), pos_key(32, 0xDD);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> distDATA(0, 1 << 10);
    std::uniform_int_distribution<> distM(0, 1);

    // Init input
    std::vector<int> DB(dbsize), aDB(dbsize);
    std::vector<int> AS;
    std::vector<int> aAS;
    for (int i = 0; i < dbsize; ++i) {
        DB[i] = distDATA(gen);
        aDB[i] = distDATA(gen);
    }
    generateBSAddr(dbsize, AS);
    generateTCAddr(0, dbsize, aAS);
    int first_dim = AS.size() / aAS.size();
    int total_len = first_dim * aAS.size();
    std::vector<bool> aM_flat(total_len);
    for (int i = 0; i < total_len; ++i) {
        aM_flat[i] = static_cast<bool>(distM(gen));
        // std::cout << aM_flat[i] << "\n";
    }

    int Z = 4;
    int acsize = AS.size();
    int Zr = 2 * static_cast<int>(std::ceil(std::log2(dbsize + acsize)));
    int mode = 1;
    // int power2 = static_cast<int>(std::ceil(std::log2(dbsize)));


    OORAM ooram(dbsize, acsize, Z, Zr, bsize, mode,
                normal_key, covert_key, iv_key, pos_key);

    auto t0 = std::chrono::high_resolution_clock::now();
    ooram.ainit(DB, AS, aDB, aAS, aM_flat);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<int> cur_vals = DB, acur_vals = aDB;
    size_t real_actimes = 100; // (AS.size() >> 1)
    auto t2 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < real_actimes; ++i) {
        print_progress(i + 1, real_actimes, "aAccess");
        ooram.aaccess();
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << std::endl;
    // ooram.verify_BSTC_result();
    

    auto init_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto access_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    fout << dbsize << "," << real_actimes << "," << init_ms << "," << access_ms << "," << (init_ms + access_ms) << "\n";

    fout.close();
    std::cout << "Finished benchmark_aaccess, saved to " << filename << std::endl;
}

void run_test_task(int dbsize, int bsize) {
    benchmark_aaccess(dbsize, bsize);
    benchmark_access(dbsize, bsize);
}

int main() {
    const int max_threads = 1;

    // === Vary dbsize ===
    // std::vector<int> dbsizes = {1<<19, 1<<18, 1<<17};// 1<<20, 
    // int bsize_fixed = 16;

    /*
    std::cout << "=== Running benchmarks: varying dbsize, fixed bsize = " << bsize_fixed << " ===\n";
    for (size_t i = 0; i < dbsizes.size(); i += max_threads) {
        std::vector<std::thread> batch_threads;
        for (size_t j = i; j < i + max_threads && j < dbsizes.size(); ++j) {
            int dbsize = dbsizes[j];
            batch_threads.emplace_back(run_test_task, dbsize, bsize_fixed);
        }
        for (auto& t : batch_threads) {
            t.join();
        }
    }
    */

    // === Vary bsize ===
    int dbsize_fixed = 1 << 20;
    std::vector<int> bsizes = {4096, 1024, 256, 64, 16};
    std::cout << "\n=== Running benchmarks: fixed dbsize = " << dbsize_fixed << ", varying bsize ===\n";
    for (size_t i = 0; i < bsizes.size(); i += max_threads) {
        std::vector<std::thread> batch_threads;
        for (size_t j = i; j < i + max_threads && j < bsizes.size(); ++j) {
            int bsize = bsizes[j];
            batch_threads.emplace_back(run_test_task, dbsize_fixed, bsize);
        }
        for (auto& t : batch_threads) {
            t.join();
        }
    }

    std::cout << "\nAll benchmarks finished.\n";
    return 0;
}
