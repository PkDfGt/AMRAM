#include "OORAM.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cassert>
#include <chrono>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;
using ByteVec = std::vector<unsigned char>;

// benchmark_access：Fix dbsize，vary acsize
void benchmark_access(int dbsize, int acsize, int bsize) {
    fs::create_directories("Result");
    std::string filename = "Result/ooram_dbsize_" + std::to_string(dbsize)
                         + "_acsize_" + std::to_string(acsize) + "_bsize_" + std::to_string(bsize) + ".csv";

    std::ofstream fout(filename);

    int Z = 4;
    int Zr = 2 * static_cast<int>(std::ceil(std::log2(dbsize + acsize)));
    int mode = 0;

    ByteVec normal_key(32, 0xAA), covert_key(32, 0xBB),
            iv_key(32, 0xCC), pos_key(32, 0xDD);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> distAS(0, dbsize - 1);
    std::uniform_int_distribution<> distOP(0, 1);
    std::uniform_int_distribution<> distDATA(1, 1 << 10);

    std::vector<int> DB(dbsize), AS(acsize);
    for (int i = 0; i < dbsize; ++i) DB[i] = distDATA(gen);
    for (int& a : AS) a = distAS(gen);

    OORAM ooram(dbsize, acsize, Z, Zr, bsize, mode,
                normal_key, covert_key, iv_key, pos_key);

    auto t0 = std::chrono::high_resolution_clock::now();
    ooram.init(DB, AS);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<int> cur_vals = DB;
    auto t2 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < AS.size(); ++i) {
        print_progress(i + 1, AS.size(), "Access");
        int addr = AS[i];
        int op = distOP(gen);
        int new_val = distDATA(gen);
        int val = ooram.access(op, addr, new_val);
        if (op) cur_vals[addr] = new_val;
        assert(val == cur_vals[addr]);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << std::endl;

    auto init_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    auto access_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    auto access_enc_dec_time = ooram.get_enc_dec_time();
    std::cout << "normal_access_enc_dec_time: " << access_enc_dec_time << "\n";

    double access_logic_time = access_ms - access_enc_dec_time;

    fout << "dbsize,acsize,init_ms,access_ms,enc_dec_ms,logic_ms,total_ms\n";

    fout << dbsize << "," << acsize << "," 
        << init_ms << "," 
        << access_ms << ","
        << access_enc_dec_time << ","
        << access_logic_time << ","
        << (init_ms + access_ms) << "\n";

    fout.close();
    std::cout << "Finished benchmark_access, saved to " << filename << std::endl;
}

// benchmark_aaccess：Fix dbsize，Vary acsize
void benchmark_aaccess(int dbsize, int acsize, int bsize) {
    fs::create_directories("Result");
    std::string filename = "Result/ooram_aaccess_dbsize_" + std::to_string(dbsize)
                         + "_acsize_" + std::to_string(acsize) + "_bsize_" + std::to_string(bsize) + ".csv";

    std::ofstream fout(filename);

    int Z = 4;
    int Zr = 2 * static_cast<int>(std::ceil(std::log2(dbsize + acsize)));
    int mode = 1;

    ByteVec normal_key(32, 0xAA), covert_key(32, 0xBB),
            iv_key(32, 0xCC), pos_key(32, 0xDD);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> distAS(0, dbsize - 1);
    std::uniform_int_distribution<> distOP(0, 1);
    std::uniform_int_distribution<> distDATA(1, 1 << 10);

    std::vector<int> DB(dbsize), aDB(dbsize);
    for (int i = 0; i < dbsize; ++i) {
        DB[i] = distDATA(gen);
        aDB[i] = distDATA(gen);
    }

    std::vector<int> AS(acsize), aAS(acsize);
    for (int i = 0; i < acsize; ++i) {
        AS[i] = distAS(gen);
        aAS[i] = distAS(gen);
    }

    OORAM ooram(dbsize, acsize, Z, Zr, bsize, mode,
                normal_key, covert_key, iv_key, pos_key);

    auto t0 = std::chrono::high_resolution_clock::now();
    ooram.ainit(DB, AS, aDB, aAS);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<int> cur_vals = DB, acur_vals = aDB;
    auto t2 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < AS.size(); ++i) {
        print_progress(i + 1, AS.size(), "aAccess");
        int addr = AS[i], aaddr = aAS[i];
        int op = distOP(gen), aop = distOP(gen);
        int new_val = distDATA(gen), anew_val = distDATA(gen);
        auto [val, aval] = ooram.aaccess(op, addr, new_val, aop, aaddr, anew_val);
        if (op) cur_vals[addr] = new_val;
        if (aop) acur_vals[aaddr] = anew_val;
        assert(val == cur_vals[addr]);
        assert(aval == acur_vals[aaddr]);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << std::endl;

    double init_ms   = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    double access_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count();
    double access_enc_dec_time = ooram.get_enc_dec_time();
    std::cout << "anamorphic_access_enc_dec_time: " << access_enc_dec_time << "\n";

    double access_logic_time = access_ms - access_enc_dec_time;

    fout << "dbsize,acsize,init_ms,access_ms,enc_dec_ms,logic_ms,total_ms\n";

    fout << dbsize << "," << acsize << "," 
        << init_ms << "," 
        << access_ms << ","
        << access_enc_dec_time << ","
        << access_logic_time << ","
        << (init_ms + access_ms) << "\n";

    fout.close();
    std::cout << "Finished benchmark_aaccess, saved to " << filename << std::endl;
}

int main() {
    
    // 第一组：固定 block size=16 bytes
    int block_size_fixed = 16;

    // dbsize=2^15
    int dbsize1 = 1 << 15;
    std::vector<int> acsizes1;
    for (int i = 15; i <= 20; ++i) {
        acsizes1.push_back(1 << i);
    }

    // dbsize=2^20
    int dbsize2 = 1 << 20;
    std::vector<int> acsizes2;
    for (int i = 20; i <= 25; ++i) {
        acsizes2.push_back(1 << i);
    }

    std::cout << "=== First task: fix block size=16 bytes ===" << std::endl;

    for (auto acsize : acsizes1) {
        benchmark_access(dbsize1, acsize, block_size_fixed);
        benchmark_aaccess(dbsize1, acsize, block_size_fixed);
    }
    for (auto acsize : acsizes2) {
        benchmark_access(dbsize2, acsize, block_size_fixed);
        benchmark_aaccess(dbsize2, acsize, block_size_fixed);
    }
    


    // Second task: fix dbsize=2^20, acsize=2^20；vary blocksize 
    int dbsize_fixed = 1 << 20; // 1 << 15;
    int acsize_fixed = dbsize_fixed;
    std::vector<int> block_sizes = {4096, 1024, 256, 64, 16};
    std::vector<int> block_sizes2 = {4096, 1024, 256, 64, 16}; //, 64, 256, 1024, 4096};

    std::cout << "=== Second task: fix dbsize=2^20, acsize=2^20 ===" << std::endl;

    for (auto bsize : block_sizes) {
        benchmark_aaccess(dbsize_fixed, acsize_fixed, bsize);
    }

    for (auto bsize : block_sizes2) {
        benchmark_access(dbsize_fixed, acsize_fixed, bsize);
    }

    return 0;
}
