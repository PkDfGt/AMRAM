// OPQ.cpp
#include "../AESN/AMAES.h"
#include "OPQ.h"
#include <climits>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <openssl/evp.h>
#include <omp.h>

using ByteVec = std::vector<unsigned char>;

OPQ::OPQ(int N_, int Z_, int Zr_, int bsize_, int mode_, const ByteVec& nkey_, const ByteVec& ckey_, const ByteVec& ikey_)
    : N(N_), Z(Z_), Zr(Zr_), bsize(bsize_), tau(1), mode(mode_), normal_key(nkey_), covert_key(ckey_), iv_key(ikey_) {
    aIV = generateRandomBytes(16);
    L = static_cast<int>(std::ceil(std::log2(N_))) + 1;
    tree_size = (1 << L) - 1;

    write_count.resize(tree_size, 0);
    tree.resize(tree_size);
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        tree[i].resize(cap);
    }
    cmins.resize(tree_size);

    // Non-root
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            if (mode == 0) tree[i][j] = gen_dummy_enc();
            else tree[i][j] = gen_dummy_aenc(i, j);
        }
        if (mode == 0) cmins[i] = gen_dummy_enc();
        else cmins[i] = gen_dummy_aenc(i, cap);
    }
}

void OPQ::init(std::vector<int> Prios, std::vector<int> Vals) {
    std::vector<std::vector<PlainEle>> plain_tree;
    std::vector<PlainEle> plain_cmins;

    assert(Prios.size() == Vals.size());
    int total = (int)Prios.size();

    // 1. Initialization
    plain_tree.resize(tree_size);
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        plain_tree[i].reserve(cap);
    }
    plain_cmins.resize(tree_size, {INT_MAX, 0, 0, 0});

    // 2. Randomly insert
    for (int i = 0; i < total; ++i) {
        int p = Prios[i];
        int v = Vals[i];
        int pos = generate_random_leaf(0, N - 1);
        PlainEle ele = {p, v, pos, tau++};

        auto path = get_path_indices(pos);
        for (int idx : path) {
            int cap = (idx == 0) ? Zr : Z;
            if ((int)plain_tree[idx].size() < cap) {
                plain_tree[idx].push_back(ele);
                break;
            }
        }
    }

    // Update plain_cmins
    for (int i = tree_size - 1; i >= 0; --i) {
        PlainEle min_ele = {INT_MAX, 0, 0, 0};
        for (const PlainEle& e : plain_tree[i]) {
            if (e.p != 0 && e.p < min_ele.p) {
                min_ele = e;
            }
        }
        if (i * 2 + 1 < tree_size) {
            if (plain_cmins[i*2+1].p != 0 && plain_cmins[i*2+1].p < min_ele.p) {
                min_ele = plain_cmins[i*2+1];
            }
        }
        if (i * 2 + 2 < tree_size) {
            if (plain_cmins[i*2+2].p != 0 && plain_cmins[i*2+2].p < min_ele.p) {
                min_ele = plain_cmins[i*2+2];
            }
        }
        plain_cmins[i] = min_ele;
    }

    // 3. Insert all
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;

        int inserted = 0;
        for (const PlainEle& e : plain_tree[i]) {
            if (inserted < cap) {
                tree[i][inserted++] = enc_element(e);
            }
        }
        // Pad dummy
        for (; inserted < cap; ++inserted) {
            tree[i][inserted] = gen_dummy_enc();
        }

        cmins[i] = enc_element(plain_cmins[i]);
    }

    plain_tree.clear();
    plain_cmins.clear();

    reset_time();
}

void OPQ::ainit(std::vector<int> nPrios, std::vector<int> nVals,
                std::vector<int> aPrios, std::vector<int> aVals, ByteVec poskey) {
    assert(nPrios.size() == nVals.size());
    assert(aPrios.size() == aVals.size());

    int total_normal = (int)nPrios.size();
    int total_anam = (int)aPrios.size();

    std::vector<std::vector<std::pair<PlainEle, PlainEle>>> plain_tree(tree_size);
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        plain_tree[i].resize(cap, { {0,0,0,0}, {0,0,0,0} });  // 全 dummy
    }

    std::vector<std::pair<PlainEle, PlainEle>> plain_cmins(tree_size, { {INT_MAX,0,0,0}, {INT_MAX,0,0,0} });

    // Insert normal elements
    for (int i = 0; i < total_normal; ++i) {
        int p = nPrios[i];
        int v = nVals[i];
        int pos = generate_pseudorandom_leaf(p, poskey);
        PlainEle ele = {p, v, pos, tau++};

        auto path = get_path_indices(pos);
        bool inserted = false;
        for (int idx : path) {
            int cap = (idx == 0) ? Zr : Z;
            for (int j = 0; j < cap; ++j) {
                if (plain_tree[idx][j].first.p == 0) {  // normal
                    plain_tree[idx][j].first = ele;
                    inserted = true;
                    break;
                }
            }
            if (inserted) break;
        }
    }

    // Insert anamorphic elements
    for (int i = 0; i < total_anam; ++i) {
        int p = aPrios[i];
        int v = aVals[i];
        int pos = generate_pseudorandom_leaf(p, poskey);
        PlainEle ele = {p, v, pos, tau++};

        auto path = get_path_indices(pos);
        bool inserted = false;
        for (int idx : path) {
            int cap = (idx == 0) ? Zr : Z;
            for (int j = 0; j < cap; ++j) {
                if (plain_tree[idx][j].first.p == 0 && plain_tree[idx][j].second.p == 0) {
                    plain_tree[idx][j].second = ele;
                    inserted = true;
                    break;
                }
            }
            if (inserted) break;
        }
    }

    // ✨ Step: Update plain_cmins
    for (int i = tree_size - 1; i >= 0; --i) {
        PlainEle min_normal = {INT_MAX, 0, 0, 0};
        PlainEle min_anam   = {INT_MAX, 0, 0, 0};

        int cap = (i == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            if (plain_tree[i][j].first.p != 0 && plain_tree[i][j].first.p < min_normal.p) {
                min_normal = plain_tree[i][j].first;
            }
            if (plain_tree[i][j].second.p != 0 && plain_tree[i][j].second.p < min_anam.p) {
                min_anam = plain_tree[i][j].second;
            }
        }
        // Children cmins
        int left = i * 2 + 1, right = i * 2 + 2;
        if (left < tree_size) {
            if (plain_cmins[left].first.p != 0 && plain_cmins[left].first.p < min_normal.p)
                min_normal = plain_cmins[left].first;
            if (plain_cmins[left].second.p != 0 && plain_cmins[left].second.p < min_anam.p)
                min_anam = plain_cmins[left].second;
        }
        if (right < tree_size) {
            if (plain_cmins[right].first.p != 0 && plain_cmins[right].first.p < min_normal.p)
                min_normal = plain_cmins[right].first;
            if (plain_cmins[right].second.p != 0 && plain_cmins[right].second.p < min_anam.p)
                min_anam = plain_cmins[right].second;
        }
        plain_cmins[i].first = min_normal;
        plain_cmins[i].second = min_anam;
    }

    // Write to tree
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            auto& [ne, ae] = plain_tree[i][j];
            ByteVec iv = deriveIV(((i + 1) << 8) | j, write_count[i]);
            tree[i][j] = aenc_element(ne, ae, iv);
        }

        // cmins encryption
        auto& [ne, ae] = plain_cmins[i];
        ByteVec iv = deriveIV(((i + 1) << 8) | cap, write_count[i]);
        cmins[i] = aenc_element(ne, ae, iv);
    }
    
    reset_time();
}

void OPQ::insert(int p, int v, int pos) {
    // Randomly read and select two disjoint paths except at root
    int leaf1 = generate_random_leaf(0, (N >> 1) - 1);
    int leaf2 = generate_random_leaf(N >> 1, N - 1);
    std::vector<int> combined_nodes = get_two_path_indices(leaf1, leaf2);
    std::vector<PlainEle> elements = read_elements_from_path(combined_nodes);
    elements.push_back({p, v, pos, tau++});

    //Evict and Updatemin
    evict_and_update_min(combined_nodes, elements);
}

OPQ::PlainEle OPQ::extract_min() {
    // ReadNRm
    PlainEle min_ele = dec_element(cmins[0]);
    std::vector<int> path = get_path_indices(min_ele.pos);
    std::vector<PlainEle> elements = read_elements_from_path(path, true, min_ele);
    // Evict and Updatemin
    evict_and_update_min(path, elements);  
    return min_ele;
}

void OPQ::ainsert(int np, int nv, int npos, int ap, int av, int apos) {
    // Randomly read and select two disjoint paths except at root
    int leaf1 = generate_random_leaf(0, N / 2 - 1);
    int leaf2 = generate_random_leaf(N / 2, N - 1);
    std::vector<int> combined_nodes = get_two_path_indices(leaf1, leaf2);
    auto [normal_elements, anamorphic_elements] = aread_elements_from_path(combined_nodes); 
    normal_elements.push_back({np, nv, npos, tau});
    anamorphic_elements.push_back({ap, av, apos, tau++});

    //Evict and Updatemin
    aevict_and_update_min(combined_nodes, normal_elements, anamorphic_elements);
}

std::pair<OPQ::PlainEle, OPQ::PlainEle> OPQ::aextract_min() {
    // ReadNRm
    auto [min_normal, min_anamorphic] = adec_element(cmins[0], deriveIV(((0 + 1) << 8) | Zr, write_count[0]));
    std::vector<int> path = get_path_indices(min_normal.pos);
    auto [normal_elements, anamorphic_elements] = aread_elements_from_path(path, true, min_normal);

    // Evict and Updatemin
    aevict_and_update_min(path, normal_elements, anamorphic_elements); 
    return {min_normal, min_anamorphic};
}

bool OPQ::check_common_path(int leaf, int idx) const {
    int index = leaf + N - 1;
    while (index >= 0) {
        if (index == idx) return true;
        else if (index < idx) return false;
        index = (index - 1) >> 1;
    }
    return false;
}

std::vector<int> OPQ::get_path_indices(int leaf) const {
    int index = leaf + N - 1;
    int i = 0;
    std::vector<int> path(L);
    while (index > 0) {
        path[i++] = index;
        index = (index - 1) >> 1;
    }
    path[i++] = index;
    assert(i == L);
    return path;
}

std::vector<int> OPQ::get_two_path_indices(int leaf1, int leaf2) const {
    int index1 = leaf1 + N - 1;
    int index2 = leaf2 + N - 1;
    int i = 0;
    std::vector<int> path((L << 1) - 1);
    while (index1 > 0) {
        path[i++] = index1;
        path[i++] = index2;
        index1 = (index1 - 1) >> 1;
        index2 = (index2 - 1) >> 1;
    }
    path[i++] = index1;
    assert(i == (L << 1) - 1);
    return path;
}

ByteVec OPQ::enc_element(int p, int v, int pos, int tau) {
    assert(bsize >= 16);
    ByteVec plaintext(bsize, 0);

    write_int(plaintext, 0, p);
    write_int(plaintext, 4, v);
    write_int(plaintext, 8, pos);
    write_int(plaintext, 12, tau);

    ByteVec iv = generateRandomBytes(16);

    return measure_time([&]() {
        return NormalEnc(plaintext, normal_key, iv);
    });
    // return NormalEnc(plaintext, normal_key, iv);
}

ByteVec OPQ::enc_element(PlainEle e) {
    return enc_element(e.p, e.v, e.pos, e.tau);
}

OPQ::PlainEle OPQ::dec_element(const ByteVec& enc_ele) {
    ByteVec dec = measure_time([&]() {
        return NormalDec(enc_ele, normal_key);
    });
    // ByteVec dec = NormalDec(enc_ele, normal_key);

    int p_ = read_int(dec, 0);
    int v_ = read_int(dec, 4);
    int pos_ = read_int(dec, 8);
    int tau_ = read_int(dec, 12);
    return {p_, v_, pos_, tau_};
}

ByteVec OPQ::aenc_element(int p1, int v1, int pos1, int tau1,
                          int p2, int v2, int pos2, int tau2, ByteVec iv) {
    assert(bsize >= 16);
    ByteVec normal_plain(bsize, 0);
    ByteVec covert_plain(16, 0);

    write_int(normal_plain, 0, p1);
    write_int(normal_plain, 4, v1);
    write_int(normal_plain, 8, pos1);
    write_int(normal_plain, 12, tau1);

    write_int(covert_plain, 0, p2);
    write_int(covert_plain, 4, v2);
    write_int(covert_plain, 8, pos2);
    write_int(covert_plain, 12, tau2);

    // ByteVec iv = generateRandomBytes(16);  // or use aIV
    
    return measure_time([&]() {
        return AnamorphicEnc(normal_plain, normal_key, covert_plain, covert_key, iv);
    });
    // return AnamorphicEnc(normal_plain, normal_key, covert_plain, covert_key, iv);
}

ByteVec OPQ::aenc_element(PlainEle ne, PlainEle ae, ByteVec iv) {
    return aenc_element(ne.p, ne.v, ne.pos, ne.tau,
                        ae.p, ae.v, ae.pos, ae.tau, iv);
}

std::pair<OPQ::PlainEle, OPQ::PlainEle> OPQ::adec_element(const ByteVec& enc_ele, ByteVec iv) {
    auto [normal_plain, covert_plain] = measure_time([&]() {
        return AnamorphicDec(enc_ele, normal_key, covert_key, iv);
    });
    // auto [normal_plain, covert_plain] = AnamorphicDec(enc_ele, normal_key, covert_key, iv);

    int p1 = read_int(normal_plain, 0);
    int v1 = read_int(normal_plain, 4);
    int pos1 = read_int(normal_plain, 8);
    int tau1 = read_int(normal_plain, 12);

    int p2 = read_int(covert_plain, 0);
    int v2 = read_int(covert_plain, 4);
    int pos2 = read_int(covert_plain, 8);
    int tau2 = read_int(covert_plain, 12);

    return {{p1, v1, pos1, tau1}, {p2, v2, pos2, tau2}};
}

void OPQ::write_int(ByteVec& buf, int offset, int val) {
    if (offset + 4 > static_cast<int>(buf.size()))
        throw std::runtime_error("buffer overflow in write_int");

    buf[offset] = static_cast<unsigned char>((val >> 24) & 0xFF);
    buf[offset + 1] = static_cast<unsigned char>((val >> 16) & 0xFF);
    buf[offset + 2] = static_cast<unsigned char>((val >> 8) & 0xFF);
    buf[offset + 3] = static_cast<unsigned char>(val & 0xFF);
}

int OPQ::read_int(const ByteVec& buf, int offset) {
    return (static_cast<uint32_t>(buf[offset]) << 24) |
           (static_cast<uint32_t>(buf[offset + 1]) << 16) |
           (static_cast<uint32_t>(buf[offset + 2]) << 8) |
           (static_cast<uint32_t>(buf[offset + 3]));
}

int OPQ::generate_random_leaf(int b, int e) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(b, e);
    return dis(gen);
}

int OPQ::generate_pseudorandom_leaf(int p, ByteVec pos_key) {
    unsigned char input[16] = {0};
    input[0] = (p >> 24) & 0xFF;
    input[1] = (p >> 16) & 0xFF;
    input[2] = (p >> 8) & 0xFF;
    input[3] = p & 0xFF;

    unsigned char output[16];
    int outlen = 0;

    static thread_local EVP_CIPHER_CTX* ctx = [] {
        EVP_CIPHER_CTX* c = EVP_CIPHER_CTX_new();
        if (!c) throw std::runtime_error("Failed to alloc thread-local PRF context");
        if (EVP_EncryptInit_ex(c, EVP_aes_256_ecb(), nullptr, nullptr, nullptr) != 1)
            throw std::runtime_error("Failed to init PRF context");
        EVP_CIPHER_CTX_set_padding(c, 0);
        return c;
    }();
    EVP_EncryptInit_ex(ctx, EVP_aes_256_ecb(), NULL, pos_key.data(), NULL);
    EVP_CIPHER_CTX_set_padding(ctx, 0);
    EVP_EncryptUpdate(ctx, output, &outlen, input, 16);

    uint32_t val = (output[0]<<24) | (output[1]<<16) | (output[2]<<8) | output[3];

    return val % N;
}


ByteVec OPQ::gen_dummy_enc() {
    return measure_time([&]() {
        return NormalEnc(ByteVec(bsize, 0), normal_key, generateRandomBytes(16));
    });
    // return NormalEnc(ByteVec(bsize, 0), normal_key, generateRandomBytes(16));
}

ByteVec OPQ::gen_dummy_aenc(int bucketID, int slotID) {
    return measure_time([&]() {
        return AnamorphicEnc(ByteVec(bsize, 0), normal_key, 
                    ByteVec(16, 0), covert_key, 
                    deriveIV(((bucketID + 1) << 8) | slotID, write_count[bucketID]));
    });
    // return AnamorphicEnc(ByteVec(bsize, 0), normal_key, 
    //                 ByteVec(16, 0), covert_key, 
    //                 deriveIV(((bucketID + 1) << 8) | slotID, write_count[bucketID]));
}

std::vector<OPQ::PlainEle> OPQ::read_elements_from_path(const std::vector<int>& path, bool skip_match, const PlainEle& match) {
    std::vector<PlainEle> elements;
    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            const ByteVec& slot = tree[idx][j];
            PlainEle e = dec_element(slot);
            if (e.p != 0) {
                if (!skip_match || e.p != match.p) {
                    elements.push_back(e);
                }
            }
        }
    }
    return elements;
    
    /*
    std::vector<PlainEle> elements;
    #pragma omp parallel
    {
        std::vector<PlainEle> local_elements;
        #pragma omp for nowait
        for (int i = 0; i < (int)path.size(); ++i) {
            int idx = path[i];
            int cap = (idx == 0) ? Zr : Z;
            for (int j = 0; j < cap; ++j) {
                const ByteVec& slot = tree[idx][j];
                PlainEle e = dec_element(slot);
                if (e.p != 0) {
                    if (!skip_match || e.p != match.p) {
                        local_elements.push_back(e);
                    }
                }
            }
        }
        #pragma omp critical
        elements.insert(elements.end(), local_elements.begin(), local_elements.end());
    }
    return elements;
    */
}

void OPQ::evict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& elements) {
    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        int count = 0;
        std::vector<PlainEle> bucket_elements;
        // 1. Evict
        int write_pos = 0, read_pos = 0;
        for (; read_pos < (int)elements.size() && count < cap; ++read_pos) {
            const PlainEle& e = elements[read_pos];
            if (check_common_path(e.pos, idx)) {
                bucket_elements.push_back(e);
                count++;
            } else {
                elements[write_pos++] = e; // Non-matching
            }
        }
        // Remaining
        for (int k = read_pos; k < (int)elements.size(); ++k) {
            elements[write_pos++] = elements[k];
        }
        elements.resize(write_pos);

        // 2. Updatemin
        PlainEle minEle = {INT_MAX, 0, 0, 0};
        int left = (idx << 1) + 1, right = (idx << 1) + 2;

        if (left < tree_size) {
            PlainEle tmp = dec_element(cmins[left]);
            if (tmp.p != 0 && tmp.p < minEle.p) minEle = tmp;
        }
        if (right < tree_size) {
            PlainEle tmp = dec_element(cmins[right]);
            if (tmp.p != 0 && tmp.p < minEle.p) minEle = tmp;
        }

        int inserted = 0;
        for (const auto& e : bucket_elements) {
            if (e.p != 0 && e.p < minEle.p) minEle = e;
            tree[idx][inserted++] = enc_element(e);
        }

        for (; inserted < cap; ++inserted) {
            tree[idx][inserted] = gen_dummy_enc();
        }

        cmins[idx] = enc_element(minEle);
    }
    assert(elements.size() == 0);
}

std::pair<std::vector<OPQ::PlainEle>, std::vector<OPQ::PlainEle>> OPQ::aread_elements_from_path(const std::vector<int>& path, bool skip_match, const PlainEle& match) {
    std::vector<PlainEle> normal_elements;
    std::vector<PlainEle> anamorphic_elements;

    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            ByteVec& enc_slot = tree[idx][j];
            auto [ne, ae] = adec_element(enc_slot, deriveIV(((idx + 1) << 8) | j, write_count[idx]));
            if (ne.p != 0) {
                if (!skip_match || ne.p != match.p) {
                    normal_elements.push_back(ne);
                }
            }
            if (ae.p != 0) {
                if (!skip_match || ae.p != match.p) {
                    anamorphic_elements.push_back(ae);
                }
            }
        }
        write_count[idx]++;
    }
    return {normal_elements, anamorphic_elements};
}

void OPQ::aevict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& normal_elements, std::vector<PlainEle>& anamorphic_elements) {
    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        int count = 0;
        std::vector<PlainEle> normal_bucket, anamorphic_bucket;

        // normal
        int write_pos = 0, read_pos = 0;
        for (; read_pos < (int)normal_elements.size() && count < cap; ++read_pos) {
            const PlainEle& e = normal_elements[read_pos];
            if (check_common_path(e.pos, idx)) {
                normal_bucket.push_back(e);
                ++count;
            } else {
                normal_elements[write_pos++] = e;
            }
        }
        for (int k = read_pos; k < (int)normal_elements.size(); ++k) {
            normal_elements[write_pos++] = normal_elements[k];
        }
        normal_elements.resize(write_pos);

        // anamorphic
        write_pos = read_pos = 0;
        for (; read_pos < (int)anamorphic_elements.size() && count < cap; ++read_pos) {
            const PlainEle& e = anamorphic_elements[read_pos];
            if (check_common_path(e.pos, idx)) {
                anamorphic_bucket.push_back(e);
                ++count;
            } else {
                anamorphic_elements[write_pos++] = e;
            }
        }
        for (int k = read_pos; k < (int)anamorphic_elements.size(); ++k) {
            anamorphic_elements[write_pos++] = anamorphic_elements[k];
        }
        anamorphic_elements.resize(write_pos);

        // Find left/right min
        PlainEle min_n = {INT_MAX, 0, 0, 0};
        PlainEle min_a = {INT_MAX, 0, 0, 0};
        int left = 2 * idx + 1, right = 2 * idx + 2;
        if (left < tree_size) {
            auto [ln, la] = adec_element(cmins[left], deriveIV(((left + 1) << 8) | Z, write_count[left]));
            if (ln.p != 0 && ln.p < min_n.p) min_n = ln;
            if (la.p != 0 && la.p < min_a.p) min_a = la;
        }
        if (right < tree_size) {
            auto [rn, ra] = adec_element(cmins[right], deriveIV(((right + 1) << 8) | Z, write_count[right]));
            if (rn.p != 0 && rn.p < min_n.p) min_n = rn;
            if (ra.p != 0 && ra.p < min_a.p) min_a = ra;
        }

        // Writeback
        int inserted = 0;
        for (const auto& e : normal_bucket) {
            if (e.p != 0 && e.p < min_n.p) min_n = e;
            tree[idx][inserted] = aenc_element(e, {0, 0, 0, 0}, deriveIV(((idx + 1) << 8) | inserted, write_count[idx]));
            inserted++;
        }
        for (const auto& e : anamorphic_bucket) {
            if (e.p != 0 && e.p < min_a.p) min_a = e;
            tree[idx][inserted] = aenc_element({0, 0, 0, 0}, e, deriveIV(((idx + 1) << 8) | inserted, write_count[idx]));
            inserted++;
        }
        for (; inserted < cap; ++inserted) {
            tree[idx][inserted] = gen_dummy_aenc(idx, inserted);
        }

        cmins[idx] = aenc_element(min_n, min_a, deriveIV(((idx + 1) << 8) | cap, write_count[idx]));
    }

    assert(normal_elements.size() == 0);
    assert(anamorphic_elements.size() == 0);
}

EVP_CIPHER_CTX* getThreadLocalPRFContext() {
    static thread_local EVP_CIPHER_CTX* ctx = [] {
        EVP_CIPHER_CTX* c = EVP_CIPHER_CTX_new();
        if (!c) throw std::runtime_error("Failed to alloc thread-local PRF context");
        if (EVP_EncryptInit_ex(c, EVP_aes_256_ecb(), nullptr, nullptr, nullptr) != 1)
            throw std::runtime_error("Failed to init PRF context");
        EVP_CIPHER_CTX_set_padding(c, 0);
        return c;
    }();
    return ctx;
}

ByteVec OPQ::deriveIV(int addr, int count) const {
    ByteVec input(16, 0x00);
    input[0] = (addr >> 24) & 0xFF;
    input[1] = (addr >> 16) & 0xFF;
    input[2] = (addr >> 8) & 0xFF;
    input[3] = addr & 0xFF;
    input[4] = (count >> 24) & 0xFF;
    input[5] = (count >> 16) & 0xFF;
    input[6] = (count >> 8) & 0xFF;
    input[7] = count & 0xFF;

    ByteVec iv(16);  // Output = AES-ECB(key, input)
    int out_len = 0;

    EVP_CIPHER_CTX* ctx = getThreadLocalPRFContext();

    // Reset + reinit key
    if (EVP_EncryptInit_ex(ctx, nullptr, nullptr, iv_key.data(), nullptr) != 1)
        throw std::runtime_error("PRF Reset failed");

    if (EVP_EncryptUpdate(ctx, iv.data(), &out_len, input.data(), input.size()) != 1)
        throw std::runtime_error("PRF Encrypt failed");

    return iv;
}

void OPQ::print_tree() {
    reset_time();
    std::cout << "\n====== OPQ Tree Structure ======\n";
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        std::cout << "Node " << i << " [";
        for (int j = 0; j < cap; ++j) {
            PlainEle e = dec_element(tree[i][j]);
            if (e.p == 0) std::cout << "(d)";
            else std::cout << "(p=" << e.p << ", v=" << e.v << ", pos=" << e.pos << ", tau=" << e.tau << ")";
            if (j != cap - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    std::cout << "\n====== cmins (per node min) ======\n";
    for (int i = 0; i < tree_size; ++i) {
        PlainEle e = dec_element(cmins[i]);
        std::cout << "cmins[" << i << "]: ";
        if (e.p == 0) std::cout << "(d)";
        else std::cout << "(p=" << e.p << ", v=" << e.v << ", pos=" << e.pos << ", tau=" << e.tau << ")";
        std::cout << "\n";
    }
    std::cout << "==================================\n";
    reset_time();
}

std::vector<std::vector<std::pair<OPQ::PlainEle, OPQ::PlainEle>>> OPQ::export_anamorphic_tree_snapshot() {
    reset_time();
    std::vector<std::vector<std::pair<PlainEle, PlainEle>>> snapshot(tree_size);
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            auto [ne, ae] = adec_element(tree[i][j], deriveIV(((i + 1) << 8) | j, write_count[i]));
            snapshot[i].emplace_back(ne, ae);
        }
    }
    reset_time();
    return snapshot;
}

void OPQ::print_anamorphic_tree() {
    reset_time();
    std::cout << "\n====== OPQ Dual Tree Structure (Normal | Anamorphic) ======\n";
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        std::cout << "Node " << i << " [";

        for (int j = 0; j < cap; ++j) {
            auto [ne, ae] = adec_element(tree[i][j], deriveIV(((i + 1) << 8) | j, write_count[i]));
            std::cout << "(";
            if (ne.p == 0) ;//std::cout << "d";
            else std::cout << "p=" << ne.p << ",v=" << ne.v << ",pos=" << ne.pos << ",tau=" << ne.tau;
            std::cout << " | ";
            if (ae.p == 0) ;//std::cout << "d";
            else std::cout << "p=" << ae.p << ",v=" << ae.v << ",pos=" << ae.pos << ",tau=" << ae.tau;
            std::cout << ")";
            if (j != cap - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    std::cout << "\n====== cmins (Normal | Anamorphic) ======\n";
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        auto [ne, ae] = adec_element(cmins[i], deriveIV(((i + 1) << 8) | cap, write_count[i]));
        std::cout << "cmins[" << i << "]: (";
        std::cout << "p=" << ne.p << ",v=" << ne.v << ",pos=" << ne.pos << ",tau=" << ne.tau;
        std::cout << " | ";
        std::cout << "p=" << ae.p << ",v=" << ae.v << ",pos=" << ae.pos << ",tau=" << ae.tau;
        std::cout << ")\n";
    }
    std::cout << "===========================================================\n";
    reset_time();
}

void OPQ::print_elements(PlainEle e) {
    std::cout << e.p  << " " << e.v << " " << e.pos << " " << e.tau << "\n";
}

void OPQ::reset_time() {
    total_enc_dec_time = 0.0;
}
