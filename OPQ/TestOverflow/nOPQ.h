#ifndef AOPQ_H
#define AOPQ_H

#include <vector>
#include <iostream>
#include <climits>
#include <cassert>
#include <cmath>

class nOPQ {
public:
    struct PlainEle { int p, pos; };
    nOPQ(int N_, int Z_, int Zr_);
    void insert(int p, int pos);
    PlainEle extract_min();
    void ainsert(int p, int pos, int ap, int apos);
    std::pair<PlainEle, PlainEle> aextract_min();
    void print_tree() const;
    int get_rootused();

private:

    int N, L, tree_size, Z, Zr;
    std::vector<std::vector<PlainEle>> tree;
    std::vector<PlainEle> cmins;
    
    int generate_random_leaf(int b, int e);
    std::vector<int> get_path_indices(int leaf) const;
    std::vector<int> get_two_path_indices(int leaf1, int leaf2) const;

    bool check_common_path(int pos, int idx) const;
    std::vector<PlainEle> read_elements_from_path(const std::vector<int>& path, bool skip_match=false, const PlainEle& match={});
    void evict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& elements);
};

#endif
