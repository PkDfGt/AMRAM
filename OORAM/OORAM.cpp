#include "OORAM.h"
#include "OPQ.h"
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

    // 2. Collect AS
    for (int i = 0; i < n; ++i) {
        print_progress(i + 1, n, "Init");
        addr = AS[i];
        cur_ts = i + 1;
        if ((DB[addr] >> 31) == 0) {
            prio = (cur_ts << 1) | 1;  // odd
            val = encode_val_data(addr, DB[addr]);
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
        print_progress(addr + 1, dbsize, "InitAddr");
        if ((DB[addr] >> 31) == 0) {
            prio = (acsize + addr) << 2;
            val = encode_val_data(addr, DB[addr]);
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

int OORAM::access(bool op, int addr, int data) {
    // Step 1: extract timestamp
    OPQ::PlainEle e_ts = opq->extract_min();
    // std::cout << "opq::::" << opq->total_enc_dec_time << "\n";
    int ts = decode_val_ts(e_ts.v);

    // Step 2: extract data entry
    OPQ::PlainEle e_data = opq->extract_min();
    auto [a, val] = decode_val_data(e_data.v);
    assert(a == addr);  // Ensure correctness

    // Step 3: if write operation, update value
    if (op) {
        val = data;
    }

    // Step 4: re-insert updated value
    int prio = (ts << 1) | 1;  // odd priority
    int pos = generate_random_leaf();
    int new_val = encode_val_data(addr, val);
    opq->insert(prio, new_val, pos);

    return val;
}

void OORAM::ainit(std::vector<int> DB, std::vector<int> AS,
                  std::vector<int> aDB, std::vector<int> aAS) {
    /// 1. Initialize OPQ with size dbsize + acsize
    opq = std::make_unique<OPQ>(ctrleaf, Z, Zr, bsize, mode,
                                 normal_key, covert_key, iv_key);
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

    for (int i = 0; i < n; ++i) {
        print_progress(i + 1, n, "Init");
        cur_ts = i + 1;

        naddr = AS[i];
        if ((DB[naddr] >> 31) == 0) {
            nprio = (cur_ts << 1) | 1;
            nval = encode_val_data(naddr, DB[naddr]);
        } else {
            nprio = (DB[naddr] & 0x7FFFFFFF) << 1;
            nval = encode_val_ts(cur_ts);
        }
        DB[naddr] = cur_ts | (1 << 31);

        aaddr = aAS[i];
        if ((aDB[aaddr] >> 31) == 0) {
            aprio = (cur_ts << 1) | 1;
            aval = encode_val_data(aaddr, aDB[aaddr]);
        } else {
            aprio = (aDB[aaddr] & 0x7FFFFFFF) << 1;
            aval = encode_val_ts(cur_ts);
        }
        aDB[aaddr] = cur_ts | (1 << 31);

        nprios.push_back(nprio);
        nvals.push_back(nval);
        aprios.push_back(aprio);
        avals.push_back(aval);
    }

    // 3. Collect dbsize
    for (int addr = 0; addr < dbsize; ++addr) {
        print_progress(addr + 1, dbsize, "InitAddr");
        if ((DB[addr] >> 31) == 0) {
            nprio = (acsize + addr) << 2;
            nval = encode_val_data(addr, DB[addr]);
        } else {
            nprio = (DB[addr] & 0x7FFFFFFF) << 1;
            nval = encode_val_ts(dbsize + acsize + addr);
        }

        if ((aDB[addr] >> 31) == 0) {
            aprio = (acsize + addr) << 2;
            aval = encode_val_data(addr, aDB[addr]);
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
    //std::shuffle(aprios.begin(), aprios.end(), g);
    //std::shuffle(avals.begin(), avals.end(), g);

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

    // opq->print_anamorphic_tree();
}

std::pair<int, int> OORAM::aaccess(bool op, int addr, int data,
                                   bool aop, int aaddr, int adata) {
    // Step 1: extract two timestamp entries
    auto [e_ts_n, e_ts_a] = opq->aextract_min();
    // std::cout << "opq::::" << opq->total_enc_dec_time << "\n";

    int ts_n = decode_val_ts(e_ts_n.v);
    int ts_a = decode_val_ts(e_ts_a.v);

    // Step 2: extract two data entries
    auto [e_data_n, e_data_a] = opq->aextract_min();

    auto [addr_n, val_n] = decode_val_data(e_data_n.v);
    auto [addr_a, val_a] = decode_val_data(e_data_a.v);
 
    assert(addr_n == addr);
    assert(addr_a == aaddr);

    // Step 3: write if needed
    if (op)   val_n = data;
    if (aop)  val_a = adata;

    // Step 4: generate common pos from prio
    int prio_n = (ts_n << 1) | 1;
    int prio_a = (ts_a << 1) | 1;
    int pos_n  = generate_pseudorandom_leaf(prio_n);
    int pos_a  = generate_pseudorandom_leaf(prio_a);

    // Step 5: encode and insert
    int val_encoded_n = encode_val_data(addr_n, val_n);
    int val_encoded_a = encode_val_data(addr_a, val_a);
    opq->ainsert(prio_n, val_encoded_n, pos_n,
            prio_a, val_encoded_a, pos_a);

    return {val_n, val_a};  // return both values
}

double OORAM::get_enc_dec_time() {
    return opq->total_enc_dec_time / 1000 / 1000;
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

int OORAM::encode_val_data(int addr, int data) {
    assert((addr & 0xFFFFF) == addr);  // 20 bits
    assert((data & 0x7FF) == data);    // 11 bits
    return (1 << 31) | (data << 20) | addr;  // MSB = 1
}

bool OORAM::is_ts(int val) {
    return (val & (1 << 31)) == 0;
}

int OORAM::decode_val_ts(int val) {
    return val & 0x7FFFFFFF;
}

std::pair<int, int> OORAM::decode_val_data(int val) {
    int addr = val & 0xFFFFF;
    int data = (val >> 20) & 0x7FF;
    return {addr, data};
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



