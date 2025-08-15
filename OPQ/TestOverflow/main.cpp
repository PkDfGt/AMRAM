#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cassert>
#include <string>
#include <cmath>
#include <functional>
#include <thread>
#include "nOPQ.h"
#include "aOPQ.h"

int get_pos_from_p(int p, int dbsize) {
    std::hash<int> hasher;
    size_t h = hasher(p);
    return h % dbsize;
}

void run_opq_overflow_test(int dbsize, int total_ops, int Z, int Zr, const std::string& filename) {
    nOPQ nopq(dbsize, Z, Zr);

    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, dbsize - 1);

    int p_counter = 1;

    // maxBroot range from 2^0 to 2^10
    std::vector<int> maxBroot;
    for (int e = 0; e <= 10; ++e) {
        maxBroot.push_back(1 << e);
    }
    std::vector<int> overflow_counts(maxBroot.size(), 0);

    // Initialize OPQ
    for (int i = 0; i < dbsize; ++i) {
        int pos = dis(gen);
        nopq.insert(p_counter++, pos);
    }

    for (int i = 0; i < total_ops; ++i) {
        auto min_ele = nopq.extract_min();
        assert(min_ele.p == p_counter - dbsize);

        int new_pos = dis(gen);
        nopq.insert(p_counter++, new_pos);

        int root_used = nopq.get_rootused();
        for (size_t j = 0; j < maxBroot.size(); ++j) {
            if (root_used > maxBroot[j]) {
                overflow_counts[j]++;
            }
        }

        if (i % 10000000 == 0) {
            std::cout << "After " << i / 1000000 << " million operations:\n";
            for (size_t j = 0; j < maxBroot.size(); ++j) {
                std::cout << "  maxBroot " << maxBroot[j] << ": "
                          << overflow_counts[j] << " overflows, "
                          << static_cast<double>(overflow_counts[j]) / (i + 1) << " probability\n";
            }
        }
    }

    // Write
    std::ofstream ofs(filename);
    ofs << "maxBroot,overflow_count,overflow_probability\n";
    for (size_t j = 0; j < maxBroot.size(); ++j) {
        double prob = static_cast<double>(overflow_counts[j]) / total_ops;
        ofs << maxBroot[j] << "," << overflow_counts[j] << "," << prob << "\n";
    }
    ofs.close();

    std::cout << "Test completed. Overflow stats saved to " << filename << std::endl;
}

void run_aopq_overflow_test(int dbsize, int total_ops, int Z, int Zr, const std::string& filename) {
    aOPQ aopq(dbsize, Z, Zr);

    std::mt19937 gen(std::random_device{}());

    // Init
    std::vector<int> priorities(dbsize);
    for (int i = 0; i < dbsize; ++i) priorities[i] = i + 1;

    auto normal_prios = priorities;
    auto anam_prios = priorities;

    std::shuffle(normal_prios.begin(), normal_prios.end(), gen);
    std::shuffle(anam_prios.begin(), anam_prios.end(), gen);
    int np, ap, npos, apos;
    // Insert
    for (int i = 0; i < dbsize; ++i) {
        np = normal_prios[i];
        ap = anam_prios[i];
        npos = get_pos_from_p(np, dbsize);
        apos = get_pos_from_p(ap, dbsize);
        aopq.ainsert(np, npos, ap, apos);
    }

    // Root overflow
    std::vector<int> maxBroot;
    for (int e = 0; e <= 10; ++e) maxBroot.push_back(1 << e);
    std::vector<int> overflow_counts(maxBroot.size(), 0);

    int inserted_count = dbsize;
    int normal_idx = 0;
    int anam_idx = 0;
    int offset = dbsize;

    for (int i = 0; i < total_ops; ++i) {
        auto [min_normal, min_anam] = aopq.aextract_min();

        if (min_normal.p != min_anam.p) {
            std::cout << "optimes: " << i << " ";
            std::cout << "min_normal.p: " << min_normal.p << " min_anam.p " << min_anam.p << "\n";
        }
        assert(min_normal.p == min_anam.p);

        np = normal_prios[normal_idx++] + offset;
        ap = anam_prios[anam_idx++] + offset;
        npos = get_pos_from_p(np, dbsize);
        apos = get_pos_from_p(ap, dbsize);
        aopq.ainsert(np, npos, ap, apos);

        inserted_count++;
        if (inserted_count % dbsize == 0) {
            // Shuffle normal/anam order
            std::shuffle(normal_prios.begin(), normal_prios.end(), gen);
            std::shuffle(anam_prios.begin(), anam_prios.end(), gen);
            normal_idx = 0;
            anam_idx = 0;
            offset += dbsize;
        }

        // Overflow
        int root_used = aopq.get_rootused();
        for (size_t j = 0; j < maxBroot.size(); ++j)
            if (root_used > maxBroot[j]) overflow_counts[j]++;

        if (i % 10000000 == 0 && i > 0) {
            std::cout << "After " << (i / 1000000) << " million operations:\n";
            for (size_t j = 0; j < maxBroot.size(); ++j) {
                std::cout << "  maxBroot " << maxBroot[j] << ": "
                          << overflow_counts[j] << " overflows, "
                          << static_cast<double>(overflow_counts[j]) / (i + 1) << " probability\n";
            }
        }
    }

    // Write
    std::ofstream ofs(filename);
    ofs << "maxBroot,overflow_count,overflow_probability\n";
    for (size_t j = 0; j < maxBroot.size(); ++j) {
        double prob = static_cast<double>(overflow_counts[j]) / total_ops;
        ofs << maxBroot[j] << "," << overflow_counts[j] << "," << prob << "\n";
    }
    ofs.close();

    std::cout << "AOPQ test completed. Overflow stats saved to " << filename << std::endl;
}


void run_aopq_worst_overflow_test(int dbsize, int total_ops, int Z, int Zr, const std::string& filename) {
    aOPQ aopq(dbsize, Z, Zr);

    std::mt19937 gen(std::random_device{}());

    std::vector<int> priorities(dbsize);
    for (int i = 0; i < dbsize; ++i) priorities[i] = i + 1;

    auto normal_prios = priorities;
    auto anam_prios = priorities;

    int np, ap, npos, apos;
    for (int i = 0; i < dbsize; ++i) {
        np = normal_prios[i];
        ap = anam_prios[i];
        npos = get_pos_from_p(np, dbsize);
        apos = get_pos_from_p(ap, dbsize);
        aopq.ainsert(np, npos, ap, apos);
    }

    std::vector<int> maxBroot;
    for (int e = 0; e <= 10; ++e) maxBroot.push_back(1 << e);
    std::vector<int> overflow_counts(maxBroot.size(), 0);

    int inserted_count = dbsize;
    int normal_idx = 0;
    int anam_idx = 0;
    int offset = dbsize;

    for (int i = 0; i < total_ops; ++i) {
        auto [min_normal, min_anam] = aopq.aextract_min();

        if (min_normal.p != min_anam.p) {
            std::cout << "optimes: " << i << " ";
            std::cout << "min_normal.p: " << min_normal.p << " min_anam.p " << min_anam.p << "\n";
        }
        assert(min_normal.p == min_anam.p);

        np = normal_prios[normal_idx++] + offset;
        ap = anam_prios[anam_idx++] + offset;
        npos = get_pos_from_p(np, dbsize);
        apos = get_pos_from_p(ap, dbsize);
        aopq.ainsert(np, npos, ap, apos);

        inserted_count++;
        if (inserted_count % dbsize == 0) {
            normal_idx = 0;
            anam_idx = 0;
            offset += dbsize;
        }

        int root_used = aopq.get_rootused();
        for (size_t j = 0; j < maxBroot.size(); ++j)
            if (root_used > maxBroot[j]) overflow_counts[j]++;

        if (i % 10000000 == 0 && i > 0) {
            std::cout << "After " << (i / 1000000) << " million operations:\n";
            for (size_t j = 0; j < maxBroot.size(); ++j) {
                std::cout << "  maxBroot " << maxBroot[j] << ": "
                          << overflow_counts[j] << " overflows, "
                          << static_cast<double>(overflow_counts[j]) / (i + 1) << " probability\n";
            }
        }
    }

    std::ofstream ofs(filename);
    ofs << "maxBroot,overflow_count,overflow_probability\n";
    for (size_t j = 0; j < maxBroot.size(); ++j) {
        double prob = static_cast<double>(overflow_counts[j]) / total_ops;
        ofs << maxBroot[j] << "," << overflow_counts[j] << "," << prob << "\n";
    }
    ofs.close();

    std::cout << "AOPQ test completed. Overflow stats saved to " << filename << std::endl;
}

void run_test_task(int dbsize, int total_ops, int Z, int Zr, const std::string& prefix) {
    std::string n_filename = prefix + "_nOPQ_dbsize_" + std::to_string(dbsize) +
                             "_Z_" + std::to_string(Z) + "_Zr_" + std::to_string(Zr) + ".csv";
    run_opq_overflow_test(dbsize, total_ops, Z, Zr, n_filename);

    std::string a_filename = prefix + "_aOPQ_dbsize_" + std::to_string(dbsize) +
                             "_Z_" + std::to_string(Z) + "_Zr_" + std::to_string(Zr) + ".csv";
    run_aopq_overflow_test(dbsize, total_ops, Z, Zr, a_filename);

    // std::string a_w_filename = prefix + "_aOPQ_worst_dbsize_" + std::to_string(dbsize) +
    //                          "_Z_" + std::to_string(Z) + "_Zr_" + std::to_string(Zr) + ".csv";
    // run_aopq_worst_overflow_test(dbsize, total_ops, Z, Zr, a_w_filename);
}

int main() {
    const int total_ops = 1000000000;
    const int fixed_dbsize = 1 << 15;
    const int fixed_Zr = 15;

    std::vector<std::tuple<int, int, int, int, std::string>> tasks;

    // First task
    
    for (int Z = 2; Z <= 2; ++Z) {
        tasks.emplace_back(fixed_dbsize, total_ops, Z, fixed_Zr, "group1");
    }
    
    // Second task
    
    for (int Z : {4,3}) {
        for (int pow = 20; pow >= 15; --pow) {
            int dbsize = 1 << pow;
            int Zr = pow;
            tasks.emplace_back(dbsize, total_ops, Z, Zr, "group2");
        }
    } 

    const int max_threads = 4;  //4 thread

    for (size_t i = 0; i < tasks.size(); i += max_threads) {
        std::vector<std::thread> batch_threads;
        for (size_t j = i; j < i + max_threads && j < tasks.size(); ++j) {
            auto& [dbsize, total_ops, Z, Zr, prefix] = tasks[j];
            batch_threads.emplace_back(run_test_task, dbsize, total_ops, Z, Zr, prefix);
        }
        for (auto& t : batch_threads) {
            t.join();
        }
    }

    return 0;
}


