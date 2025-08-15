#include "PQ.h"
#include "../AESN/AMAES.h"
#include <cassert>

BinaryHeapPQ::BinaryHeapPQ(int N_,  int bsize_, const ByteVec& key): N(N_), bsize(bsize_), key(key) {}

bool BinaryHeapPQ::empty() const {
    return heap.empty();
}

void BinaryHeapPQ::write_int(ByteVec& buf, int offset, int val) {
    if (offset + 4 > static_cast<int>(buf.size()))
        throw std::runtime_error("buffer overflow in write_int");

    buf[offset] = static_cast<unsigned char>((val >> 24) & 0xFF);
    buf[offset + 1] = static_cast<unsigned char>((val >> 16) & 0xFF);
    buf[offset + 2] = static_cast<unsigned char>((val >> 8) & 0xFF);
    buf[offset + 3] = static_cast<unsigned char>(val & 0xFF);
}

int BinaryHeapPQ::read_int(const ByteVec& buf, int offset) {
    return (static_cast<uint32_t>(buf[offset]) << 24) |
           (static_cast<uint32_t>(buf[offset + 1]) << 16) |
           (static_cast<uint32_t>(buf[offset + 2]) << 8) |
           (static_cast<uint32_t>(buf[offset + 3]));
}

ByteVec BinaryHeapPQ::enc_element(int p, int v) {
    assert(bsize >= 16);
    ByteVec plaintext(bsize, 0);
    write_int(plaintext, 0, p);
    write_int(plaintext, 4, v);
    ByteVec iv = generateRandomBytes(16);
    return NormalEnc(plaintext, key, iv);
}

ByteVec BinaryHeapPQ::enc_element(PlainEle e) {
    return enc_element(e.p, e.v);
}

PlainEle BinaryHeapPQ::dec_element(const ByteVec& enc_ele) {
    ByteVec dec = NormalDec(enc_ele, key);
    int p_ = read_int(dec, 0);
    int v_ = read_int(dec, 4);
    return {p_, v_};
}

void BinaryHeapPQ::insert(int p, int v) {
    PlainEle e{p, v};
    heap.push_back(enc_element(e));
    heapify_up(heap.size() - 1);
}

PlainEle BinaryHeapPQ::extract_min() {
    if (heap.empty())
        throw std::runtime_error("PQ is empty");

    PlainEle res = dec_element(heap[0]);
    heap[0] = heap.back();
    heap.pop_back();
    if (!heap.empty()) heapify_down(0);
    return res;
}

void BinaryHeapPQ::heapify_up(int idx) {
    PlainEle e_idx = dec_element(heap[idx]);
    while (idx > 0) {
        int par = parent(idx);
        PlainEle e_par = dec_element(heap[par]);
        if (e_idx.p < e_par.p) {
            std::swap(heap[idx], heap[par]);
            idx = par;
        } else {
            break;
        }
    }
}

void BinaryHeapPQ::heapify_down(int idx) {
    int size = heap.size();
    PlainEle e_idx = dec_element(heap[idx]);

    while (true) {
        int l = left(idx), r = right(idx);
        int smallest = idx;
        PlainEle e_smallest = e_idx;

        if (l < size) {
            PlainEle e_l = dec_element(heap[l]);
            if (e_l.p < e_smallest.p) {
                smallest = l;
                e_smallest = e_l;
            }
        }

        if (r < size) {
            PlainEle e_r = dec_element(heap[r]);
            if (e_r.p < e_smallest.p) {
                smallest = r;
                e_smallest = e_r;
            }
        }

        if (smallest != idx) {
            std::swap(heap[idx], heap[smallest]);
            idx = smallest;
        } else {
            break;
        }
    }
}



void BinaryHeapPQ::print_decrypted() {
    std::cout << "Encrypted Heap (Decrypted View):\n";
    for (int i = 0; i < (int)heap.size(); ++i) {
        auto e = dec_element(heap[i]);
        std::cout << "[" << i << "] p=" << e.p << ", v=" << e.v << "\n";
    }
}

/*
int main() {
    ByteVec key(32, 0x00);  // 256-bit key (all zeros for test)

    BinaryHeapPQ pq(4, 64, key);
    pq.insert(5, 50);
    pq.insert(3, 30);
    pq.insert(7, 70);
    pq.insert(1, 10);
    pq.insert(75, 50);
    pq.insert(334, 30);
    pq.insert(2, 70);
    pq.insert(6, 10);

    pq.print_decrypted();

    while (!pq.empty()) {
        PlainEle e = pq.extract_min();
        std::cout << "Extracted: p=" << e.p << ", v=" << e.v << "\n";
    }

    return 0;
}
*/