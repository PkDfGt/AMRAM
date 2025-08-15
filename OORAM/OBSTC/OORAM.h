// OORAM.h
#ifndef OORAM_H
#define OORAM_H

#include <vector>
#include <cassert>
#include <stdexcept>
#include <cstdint>
#include <climits>
#include <memory>
#include <iostream>
#include <stack>
#include "Utils.h"
#include "../../OPQ/OPQ.h"

using ByteVec = std::vector<unsigned char>;

inline void print_progress(size_t current, size_t total, const std::string& phase) {
    static const int bar_width = 50;
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(bar_width * progress);
    std::cout << "\r[" << phase << "] [";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "%\r";
    std::cout.flush();
}

class OORAM {
public:
    OORAM(int dbsize_, int acsize_, int Z_, int Zr_, int bsize_, int mode_,
          const ByteVec& nkey_, const ByteVec& ckey_, const ByteVec& ikey_, const ByteVec& pkey_)
        : dbsize(dbsize_), acsize(acsize_), Z(Z_), Zr(Zr_), bsize(bsize_),
          mode(mode_), ctrleaf(round_to_pow2(dbsize + acsize)), ascending(true),
          BS_i(0), BS_j(1), BS_k(2),
          normal_key(nkey_), covert_key(ckey_), iv_key(ikey_), pos_key(pkey_),
          op_ctr(0)
    {
        stk.push({0, dbsize, 0, 0,0,0,0});
    }

    void init(std::vector<int> DB, std::vector<int> AS);
    void access();
    void ainit(std::vector<int> DB, std::vector<int> AS, 
        std::vector<int> aDB, std::vector<int> aAS, std::vector<bool> aM);
    void aaccess();

    void verify_BS_result();
    void verify_BSTC_result();
    void print_tree();
    void print_anamorphic_tree();

private:
    int dbsize, acsize, Z, Zr, bsize, mode, ctrleaf;
    bool ascending;
    int BS_i, BS_j, BS_k;
    ByteVec normal_key, covert_key, iv_key, pos_key;
    int op_ctr;
    std::stack<StackFrame> stk;
    std::unique_ptr<OPQ> opq;
    int tcop_size, toa_size, tc_result;
    int dummy_ctr;

    int round_to_pow2(int x);
    int encode_val_ts(int timestamp);
    int decode_val_ts(int val);
    int encode_val_data(const std::vector<bool>& M, int data);
    std::pair<bool, int> decode_val_data(int m_id, int val);
    bool is_ts(int val);
    int generate_random_leaf();
    int generate_pseudorandom_leaf(int p);
};

#endif