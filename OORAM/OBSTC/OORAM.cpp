#include "OORAM.h"
#include <unordered_map>
#include <memory>
#include <utility>
#include <random>
#include <iostream>
#include <openssl/evp.h>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <stack>
#include "Utils.h"

using ByteVec = std::vector<unsigned char>;

void OORAM::init(std::vector<int> DB, std::vector<int> AS) {
    /// 1. Initialize OPQ with size dbsize + acsize
    opq = std::make_unique<OPQ>(
        ctrleaf, Z, Zr, bsize, mode,
        normal_key, covert_key, iv_key
    );

    int n = (int)AS.size();
    std::vector<int> prios, vals;
    prios.reserve(n + dbsize);
    vals.reserve(n + dbsize);

    int addr, cur_ts, prio, val;
    std::vector<bool> zero_M(10, 0);  // M 全 0

    // 2. Collect AS
    for (int i = 0; i < n; ++i) {
        addr = AS[i];
        cur_ts = i + 1;
        if ((DB[addr] >> 31) == 0) {
            prio = (cur_ts << 1) | 1;  // odd
            val = encode_val_data(zero_M, DB[addr]);
        } else {
            prio = (DB[addr] & 0x7FFFFFFF) << 1;  // even
            val = encode_val_ts(cur_ts);
        }
        DB[addr] = cur_ts | (1 << 31);

        prios.push_back(prio);
        vals.push_back(val);
    }

    // 3. Collect dbsize
    for (int addr = 0; addr < dbsize; ++addr) {
        if ((DB[addr] >> 31) == 0) {
            prio = (acsize + addr) << 2;
            val = encode_val_data(zero_M, DB[addr]);
        } else {
            prio = (DB[addr] & 0x7FFFFFFF) << 1;
            val = encode_val_ts(dbsize + acsize + addr);
        }
        prios.push_back(prio);
        vals.push_back(val);
    }

    // 4. Insert
    opq->init(prios, vals);
    /*
    for (size_t i = 0; i < prios.size(); ++i) {
        print_progress(i + 1, prios.size(), "Init");
        int pos = generate_random_leaf();
        opq->insert(prios[i], vals[i], pos);
    }   
    */
    // opq->print_tree();
}

// CompreNSwap
void OORAM::access() {
    // ==========================
    // Extract for first address
    // ==========================
    OPQ::PlainEle e_ts1 = opq->extract_min();
    int ts1 = decode_val_ts(e_ts1.v);
    OPQ::PlainEle e_data1 = opq->extract_min();
    auto [_, val1] = decode_val_data(0, e_data1.v);

    // ==========================
    // Repeat for second address
    // ==========================
    OPQ::PlainEle e_ts2 = opq->extract_min();
    int ts2 = decode_val_ts(e_ts2.v);
    OPQ::PlainEle e_data2 = opq->extract_min();
    auto [__, val2] = decode_val_data(0, e_data2.v);

    int swap_val1 = val1;
    int swap_val2 = val2;
    bool need_swap = BSIte(dbsize, swap_val1, swap_val2, BS_i, BS_j, BS_k, ascending);

    // Step 4: re-insert
    int prio1 = (ts1 << 1) | 1;  // odd priority
    int pos1 = generate_random_leaf();
    int prio2 = (ts2 << 1) | 1;
    int pos2 = generate_random_leaf();
    if (need_swap) {
        opq->insert(prio1, e_data2.v, pos1);
        opq->insert(prio2, e_data1.v, pos2);
    } else {
        opq->insert(prio1, e_data1.v, pos1);
        opq->insert(prio2, e_data2.v, pos2);
    }
}

void OORAM::verify_BS_result() {
    std::vector<int> res;
    for (int i = 0; i < dbsize; i++) {
        OPQ::PlainEle e_data = opq->extract_min();
        auto [_, val] = decode_val_data(0, e_data.v);
        res.push_back(val);
    }
    verifySorted(res, ascending);
    /*
    for (int i = 0; i < dbsize; i++) {
        std::cout << res[i] << "\n";
    }
    */
}

void OORAM::ainit(std::vector<int> DB, std::vector<int> AS,
                  std::vector<int> aDB, std::vector<int> aAS, std::vector<bool> aM) {
    assert(!stk.empty());
    /// 1. Initialize OPQ with size dbsize + acsize
    opq = std::make_unique<OPQ>(ctrleaf, Z, Zr, bsize, mode,
                                 normal_key, covert_key, iv_key);
    tcop_size = aAS.size();
    tc_result = 0;
    // 2. Collect AS
    int cur_ts;
    int naddr, nprio, nval;
    int aaddr, aprio, aval;

    int n = (int)AS.size();
    std::vector<int> nprios, nvals, aprios, avals;

    nprios.reserve(n + dbsize);
    nvals.reserve(n + dbsize);
    aprios.reserve(n + dbsize);
    avals.reserve(n + dbsize);
    std::vector<bool> zero_M(10, 0); 
    int M_index = 0;
    int first_dim = AS.size() / aAS.size();
    toa_size = first_dim * aAS.size();
    dummy_ctr = 0x0FFFFFFF;

    std::vector<bool> curM;

    for (int i = 0; i < n; ++i) {
        cur_ts = i + 1;

        naddr = AS[i];
        if ((DB[naddr] >> 31) == 0) {
            nprio = (cur_ts << 1) | 1;
            nval = encode_val_data(zero_M, DB[naddr]);
        } else {
            nprio = (DB[naddr] & 0x7FFFFFFF) << 1;
            nval = encode_val_ts(cur_ts);
        }
        DB[naddr] = cur_ts | (1 << 31);

        aaddr = aAS[i % aAS.size()];
        if ((aDB[aaddr] >> 31) == 0) {
            aprio = (cur_ts << 1) | 1;
            curM.assign(aM.begin() + M_index, aM.begin() + M_index + first_dim);
            aval = encode_val_data(curM, aDB[aaddr]);
        } else {
            aprio = (aDB[aaddr] & 0x7FFFFFFF) << 1;
            aval = encode_val_ts(cur_ts);
        }
        aDB[aaddr] = cur_ts | (1 << 31);

        nprios.push_back(nprio);
        nvals.push_back(nval);
        aprios.push_back(aprio);
        avals.push_back(aval);
        
        M_index += first_dim;
        // std::cout << M_index << "\n;";
    }

    // 3. Collect dbsize
    for (int addr = 0; addr < dbsize; ++addr) {
        if ((DB[addr] >> 31) == 0) {
            nprio = (acsize + addr) << 2;
            nval = encode_val_data(zero_M, DB[addr]);
        } else {
            nprio = (DB[addr] & 0x7FFFFFFF) << 1;
            nval = encode_val_ts(dbsize + acsize + addr);
        }

        if ((aDB[addr] >> 31) == 0) {
            aprio = (acsize + addr) << 2;
            aval = encode_val_data(zero_M, aDB[addr]);
        } else {
            aprio = (aDB[addr] & 0x7FFFFFFF) << 1;
            aval = encode_val_ts(dbsize + acsize + addr);
        }

        nprios.push_back(nprio);
        nvals.push_back(nval);
        aprios.push_back(aprio);
        avals.push_back(aval);
    }

    int total = (int)nprios.size();
    std::random_device rd;
    std::mt19937 g(rd());
    std::vector<int> idx(total);
    for (int i = 0; i < total; ++i) idx[i] = i;
    std::shuffle(idx.begin(), idx.end(), g);

    // 4. Insert
    opq->ainit(nprios, nvals, aprios, avals, pos_key);
    /*
    for (int i = 0; i < total; ++i) {
        print_progress(i + 1, total, "aInit");
        int pos0 = generate_pseudorandom_leaf(nprios[i]);
        int pos1 = generate_pseudorandom_leaf(aprios[idx[i]]);
        opq->ainsert(nprios[i], nvals[i], pos0, aprios[idx[i]], avals[idx[i]], pos1);
    }
    */
    // print_anamorphic_tree();
}

void OORAM::aaccess() {
    
    // First
    auto [e_ts_n1, e_ts_a1] = opq->aextract_min();
    int ts_n1 = decode_val_ts(e_ts_n1.v);
    int ts_a1 = decode_val_ts(e_ts_a1.v);

    int m_id = op_ctr / tcop_size;
    op_ctr += 2;

    auto [e_data_n1, e_data_a1] = opq->aextract_min();
    auto [_, val_n1] = decode_val_data(0, e_data_n1.v);         // normal 用 m_id=0
    auto [mark1, val_a1] = decode_val_data(m_id, e_data_a1.v); // anamorphic 用 m_id1

    // Second
    auto [e_ts_n2, e_ts_a2] = opq->aextract_min();
    int ts_n2 = decode_val_ts(e_ts_n2.v);
    int ts_a2 = decode_val_ts(e_ts_a2.v);

    auto [e_data_n2, e_data_a2] = opq->aextract_min();
    auto [__, val_n2] = decode_val_data(0, e_data_n2.v);         // normal
    auto [mark2, val_a2] = decode_val_data(m_id, e_data_a2.v);  // anamorphic

    int mark1_var = mark1;  // 这里mark1是之前解码得到的int或bool变量
    int mark2_var = mark2;


    // Step 2: generate common pos from prio
    int prio_n1 = (ts_n1 << 1) | 1;
    int prio_a1 = (ts_a1 << 1) | 1;
    int pos_n1  = generate_pseudorandom_leaf(prio_n1);
    int pos_a1  = generate_pseudorandom_leaf(prio_a1);

    int prio_n2 = (ts_n2 << 1) | 1;
    int prio_a2 = (ts_a2 << 1) | 1;
    int pos_n2  = generate_pseudorandom_leaf(prio_n2);
    int pos_a2  = generate_pseudorandom_leaf(prio_a2);

    bool bs_need_swap = BSIte(dbsize, val_n1, val_n2, BS_i, BS_j, BS_k, ascending);
    assert(!stk.empty());
    bool tc_need_swap = compactIte(val_a1, mark1_var, val_a2, mark2_var, stk, tc_result);

    // 编码成 2 bit：bs_need_swap 在高位，tc_need_swap 在低位
    int swap_case = (bs_need_swap << 1) | (tc_need_swap ? 1 : 0);

    switch (swap_case) {
        case 0b11: // bs=1, tc=1, 都 swap
            opq->ainsert(prio_n1, e_data_n2.v, pos_n1,
                        prio_a1, e_data_a2.v, pos_a1);
            opq->ainsert(prio_n2, e_data_n1.v, pos_n2,
                        prio_a2, e_data_a1.v, pos_a2);
            break;
        case 0b10: // bs=1, tc=0, 只 swap normal
            opq->ainsert(prio_n1, e_data_n2.v, pos_n1,
                        prio_a1, e_data_a1.v, pos_a1);
            opq->ainsert(prio_n2, e_data_n1.v, pos_n2,
                        prio_a2, e_data_a2.v, pos_a2);
            break;
        case 0b01: // bs=0, tc=1, 只 swap anamorphic
            opq->ainsert(prio_n1, e_data_n1.v, pos_n1,
                        prio_a1, e_data_a2.v, pos_a1);
            opq->ainsert(prio_n2, e_data_n2.v, pos_n2,
                        prio_a2, e_data_a1.v, pos_a2);
            break;
        case 0b00: // bs=0, tc=0, 都不 swap
            opq->ainsert(prio_n1, e_data_n1.v, pos_n1,
                        prio_a1, e_data_a1.v, pos_a1);
            opq->ainsert(prio_n2, e_data_n2.v, pos_n2,
                        prio_a2, e_data_a2.v, pos_a2);
            break;
    }
    if (op_ctr % tcop_size == 0) {
        while (!stk.empty()) {
            stk.pop();
        }
        stk.push({0, dbsize, 0, 0,0,0,0});
    }
    // print_anamorphic_tree();
}

void OORAM::verify_BSTC_result() {
    std::vector<int> res_n;  // normal (BS) tree
    std::vector<int> rm_a;  // anamorphic (TC) tree
    int m_id = op_ctr / tcop_size - 1;

    for (int i = 0; i < dbsize; i++) {
        auto [e_data_n, e_data_a] = opq->aextract_min();
        
        auto [_, val_n] = decode_val_data(0, e_data_n.v);  // normal tree取第0 bit
        auto [mark_a, __] = decode_val_data(m_id, e_data_a.v); // anamorphic tree同样取
        
        res_n.push_back(val_n);
        rm_a.push_back(mark_a);
    }

    // 验证 normal tree 是否已排序
    verifySorted(res_n, ascending);
    int power2 = static_cast<int>(std::ceil(std::log2(dbsize)));
    if (power2 >> 1 << 1 != power2) {    
        verifyTC(rm_a);
    }
}


void OORAM::print_tree() {
    opq->print_tree();
}

void OORAM::print_anamorphic_tree() {
    auto snapshot = opq->export_anamorphic_tree_snapshot();

    std::cout << "\n====== OPQ Dual Tree Structure (Normal | Anamorphic) ======\n";
    for (size_t i = 0; i < snapshot.size(); ++i) {
        std::cout << "Node " << i << " [";
        for (size_t j = 0; j < snapshot[i].size(); ++j) {
            auto& [ne, ae] = snapshot[i][j];
            std::cout << "(";

            // normal tree
            if (ne.p == 0) ;// std::cout << "d";
            else {
                std::cout << "p=" << ne.p << ",pos=" << ne.pos << ",tau=" << ne.tau;
                if ((ne.v >> 31) == 0) {
                    int ts = decode_val_ts(ne.v);
                    std::cout << ",ts=" << ts;
                } else {
                    auto [_, data] = decode_val_data(0, ne.v);
                    std::cout << ",data=" << data << ",M=[";
                    for (int m_id = 0; m_id < 10; ++m_id) {
                        auto [m_bit, __] = decode_val_data(m_id, ne.v);
                        std::cout << m_bit;
                        if (m_id != 9) std::cout << ",";
                    }
                    std::cout << "]";
                }
            }

            std::cout << " | ";

            // anamorphic tree
            if (ae.p == 0) ;// std::cout << "d";
            else {
                std::cout << "p=" << ae.p << ",pos=" << ae.pos << ",tau=" << ae.tau;
                if ((ae.v >> 31) == 0) {
                    int ts = decode_val_ts(ae.v);
                    std::cout << ",ts=" << ts;
                } else {
                    auto [_, data] = decode_val_data(0, ae.v);
                    std::cout << ",data=" << data << ",M=[";
                    for (int m_id = 0; m_id < 10; ++m_id) {
                        auto [m_bit, __] = decode_val_data(m_id, ae.v);
                        std::cout << m_bit;
                        if (m_id != 9) std::cout << ",";
                    }
                    std::cout << "]";
                }
            }

            std::cout << ")";
            if (j != snapshot[i].size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "===========================================================\n";
}


int OORAM::round_to_pow2(int x) {
    if (x <= 1) return 1;
    int power = 1;
    while (power < x) power <<= 1;
    return power;
}

int OORAM::encode_val_ts(int timestamp) {
    return timestamp & 0x7FFFFFFF;  // MSB = 0
}

int OORAM::decode_val_ts(int val) {
    return val & 0x7FFFFFFF;
}

int OORAM::encode_val_data(const std::vector<bool>& M, int data) {
    assert(M.size() <= 10);
    assert((data & 0x1FFFFF) == data);  // data 占 21 bit

    int m_bits = 0;
    for (size_t i = 0; i < M.size(); ++i) {
        m_bits |= (static_cast<int>(M[i]) & 0x1) << i;
    }

    // 编码：MSB=1 | M<<21 | data
    return (1 << 31) | (m_bits << 21) | data;
}

std::pair<bool, int> OORAM::decode_val_data(int m_id, int val) {
    assert(m_id >= 0 && m_id < 10);

    int data = val & 0x1FFFFF;                  // bit 0~20
    int m_bits = (val >> 21) & 0x3FF;            // bit 21~30，最多10 bit
    bool m_bit = (m_bits >> m_id) & 0x1;          // 取第 m_id 位

    return {m_bit, data};
}

bool OORAM::is_ts(int val) {
    return (val & (1 << 31)) == 0;
}

int OORAM::generate_random_leaf() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, ctrleaf - 1);
    return dis(gen);
}

int OORAM::generate_pseudorandom_leaf(int p) {
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

    return val % ctrleaf;
}



