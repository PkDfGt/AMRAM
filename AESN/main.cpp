#include <openssl/aes.h>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <chrono>
#include <cassert>
#include <iomanip>
#include "AMAES.h"

using namespace std;
using namespace std::chrono;

using ByteVec = std::vector<unsigned char>;

int main() {
    const int iterations = 1 << 10;
    ByteVec key_normal(32, 0x11);
    ByteVec key_covert(32, 0x22);
    ByteVec IV1(16, 0x77);

    cout << "===== Total timing for Normal vs Anamorphic encryption =====" << endl;
    cout << "Size(Bytes)\tBandwidthRate\tNormal(us)\tAnam(us)" << endl;

    for (int exp = 4; exp <= 20; exp += 2) {
        size_t normal_size = 1ULL << exp;
        ByteVec plain_normal(normal_size, 0xAA);
        ByteVec plain_covert(16, 0xBB);
        double bandwidth_rate = static_cast<double>(plain_covert.size()) / normal_size;

        // === Normal timing ===
        auto start_normal = high_resolution_clock::now();
        for (int i = 0; i < iterations; ++i) {
            ByteVec iv = IV1;
            ByteVec enc = NormalEnc(plain_normal, key_normal, iv);
            ByteVec dec = NormalDec(enc, key_normal);
            assert(dec == plain_normal);
        }
        auto end_normal = high_resolution_clock::now();
        // double total_normal_us = duration_cast<microseconds>(end_normal - start_normal).count() / iterations;
        double total_normal_us  = duration_cast<nanoseconds>(end_normal - start_normal).count() / 1000.0 / iterations;

        // === Anamorphic timing ===
        auto start_ana = high_resolution_clock::now();
        for (int i = 0; i < iterations; ++i) {
            ByteVec enc = AnamorphicEnc(plain_normal, key_normal, plain_covert, key_covert, IV1);
            auto [dec_normal, dec_covert] = AnamorphicDec(enc, key_normal, key_covert, IV1);
            assert(dec_normal == plain_normal && dec_covert == plain_covert);
        }
        auto end_ana = high_resolution_clock::now();
        // double total_ana_us = duration_cast<microseconds>(end_normal - start_normal).count() / iterations;
        double total_ana_us = duration_cast<nanoseconds>(end_ana - start_ana).count() / 1000.0 / iterations;

        // Print result
        cout << normal_size << "\t\t"
             << fixed << setprecision(6) << bandwidth_rate << "\t"
             << fixed << setprecision(3) << total_normal_us  << "\t\t"
             << fixed << setprecision(3) << total_ana_us << endl;
    }

    return 0;
}
