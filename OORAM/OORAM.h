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
#include "../OPQ/OPQ.h"

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
          mode(mode_), ctrleaf(round_to_pow2(dbsize + acsize)), 
          normal_key(nkey_), covert_key(ckey_), iv_key(ikey_), pos_key(pkey_)
    {}

    void init(std::vector<int> DB, std::vector<int> AS);
    int access(bool op, int addr, int data);
    void ainit(std::vector<int> DB, std::vector<int> AS, std::vector<int> aDB, std::vector<int> aAS);
    std::pair<int, int> aaccess(bool op, int addr, int data, bool aop, int aaddr, int adata);
    double get_enc_dec_time();

private:
    int dbsize, acsize, Z, Zr, bsize, mode, ctrleaf;
    ByteVec normal_key, covert_key, iv_key, pos_key;
    std::unique_ptr<OPQ> opq;

    int round_to_pow2(int x);
    int encode_val_ts(int timestamp);
    int encode_val_data(int addr, int data);
    bool is_ts(int val);
    int decode_val_ts(int val);
    std::pair<int, int> decode_val_data(int val);
    int generate_random_leaf();
    int generate_pseudorandom_leaf(int p);
};

#endif