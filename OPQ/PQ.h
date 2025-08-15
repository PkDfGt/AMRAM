#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include "../AESN/AMAES.h"

using ByteVec = std::vector<unsigned char>;

struct PlainEle {
    int p, v;
};

class BinaryHeapPQ {
public:
    BinaryHeapPQ(int N_,  int bsize_, const ByteVec& key);

    void insert(int p, int v);
    PlainEle extract_min();
    bool empty() const;
    void print_decrypted();

private:
    int N, bsize;
    ByteVec key;
    std::vector<ByteVec> heap;

    void heapify_up(int idx);
    void heapify_down(int idx);

    void write_int(ByteVec& buf, int offset, int val);
    int read_int(const ByteVec& buf, int offset);
    ByteVec enc_element(int p, int v);
    ByteVec enc_element(PlainEle e);
    PlainEle dec_element(const ByteVec& enc_ele);

    int parent(int i) { return (i - 1) / 2; }
    int left(int i) { return 2 * i + 1; }
    int right(int i) { return 2 * i + 2; }
};
