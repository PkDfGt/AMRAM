// OPQ.h
#ifndef OPQ_H
#define OPQ_H

#include <vector>
#include <cassert>
#include <stdexcept>
#include <cstdint>
#include <climits>
#include <chrono>
#include <iostream>

using ByteVec = std::vector<unsigned char>;

class OPQ {
public:
    double total_enc_dec_time = 0.0;
    struct PlainEle { int p, v, pos, tau; }; // Element

    OPQ(int N_, int Z_, int Zr_, int bsize_, int mode_, const ByteVec& nkey_, const ByteVec& ckey_, const ByteVec& ikey_);
    void init(std::vector<int> Prios, std::vector<int> Vals);
    void insert(int p, int v, int pos);
    PlainEle extract_min();
    void ainit(std::vector<int> nPrios, std::vector<int> nVals, std::vector<int> aPrios, std::vector<int> aVals, ByteVec pos_key);
    void ainsert(int np, int nv, int npos, int ap, int av, int apos);
    std::pair<PlainEle, PlainEle> aextract_min();

    void print_tree();
    void print_anamorphic_tree();
    void print_elements(PlainEle e);
    std::vector<std::vector<std::pair<PlainEle, PlainEle>>> export_anamorphic_tree_snapshot();

    template<typename Func>
    auto measure_time(Func&& func) {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = func();
        auto end = std::chrono::high_resolution_clock::now();
        double us = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        total_enc_dec_time += us;
        return result;
    }

private:
    int N, L, Z, Zr, bsize, tree_size, tau, mode;
    ByteVec normal_key, covert_key, iv_key, aIV;
    std::vector<int> write_count;
    std::vector<std::vector<ByteVec>> tree;
    std::vector<ByteVec> cmins;

    std::vector<int> get_path_indices(int leaf) const;
    std::vector<int> get_two_path_indices(int leaf1, int leaf2) const;
    bool check_common_path(int leaf, int idx) const;

    ByteVec enc_element(int p, int v, int pos, int tau);
    ByteVec enc_element(PlainEle e);
    PlainEle dec_element(const ByteVec& enc_ele);
    ByteVec aenc_element(int p1, int v1, int pos1, int tau1,
                          int p2, int v2, int pos2, int tau2, ByteVec iv); 
    ByteVec aenc_element(PlainEle ne, PlainEle ae, ByteVec iv);
    std::pair<PlainEle, PlainEle> adec_element(const ByteVec& enc_ele, ByteVec iv);
    ByteVec gen_dummy_enc();
    ByteVec gen_dummy_aenc(int bucketID, int slotID);
    int generate_random_leaf(int b, int e);
    int generate_pseudorandom_leaf(int p, ByteVec pos_key);

    static void write_int(ByteVec& buf, int offset, int val);
    static int read_int(const ByteVec& buf, int offset);

    std::vector<PlainEle> read_elements_from_path(const std::vector<int>& path, bool skip_match = false, const PlainEle& match = {});
    void evict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& elements);
    std::pair<std::vector<PlainEle>, std::vector<PlainEle>> aread_elements_from_path(const std::vector<int>& path, bool skip_match = false, const PlainEle& match = {}); 
    void aevict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& normal_elements, std::vector<PlainEle>& anamorphic_elements);

    ByteVec deriveIV(int addr, int count) const;

    void reset_time();
};

#endif