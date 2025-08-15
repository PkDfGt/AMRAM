#include <iostream>
#include <vector>
#include <chrono>
#include <cassert>
#include <iomanip>
#include <openssl/evp.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/bn.h>
#include <openssl/sha.h>
#include "../AESN/AMAES.h"

using namespace std;
using namespace std::chrono;
using ByteVec = vector<unsigned char>;

// MGF1 for SHA256
ByteVec MGF1(const ByteVec& seed, size_t length) {
    ByteVec output;
    for (uint32_t counter = 0; output.size() < length; ++counter) {
        ByteVec data(seed);
        for (int i = 3; i >= 0; --i)
            data.push_back((counter >> (8 * i)) & 0xFF);
        unsigned char hash[SHA256_DIGEST_LENGTH];
        SHA256(data.data(), data.size(), hash);
        output.insert(output.end(), hash, hash + SHA256_DIGEST_LENGTH);
    }
    output.resize(length);
    return output;
}

// OAEP Encode
ByteVec OAEPEncode(const ByteVec& msg, const ByteVec& r, size_t k) {
    size_t hLen = SHA256_DIGEST_LENGTH;
    assert(r.size() == hLen);
    assert(msg.size() <= k - 2 * hLen - 2);

    unsigned char lHash[hLen];
    SHA256(nullptr, 0, lHash);

    ByteVec PS(k - msg.size() - 2 * hLen - 2, 0x00);
    ByteVec DB(lHash, lHash + hLen);
    DB.insert(DB.end(), PS.begin(), PS.end());
    DB.push_back(0x01);
    DB.insert(DB.end(), msg.begin(), msg.end());

    ByteVec dbMask = MGF1(r, k - hLen - 1);
    ByteVec maskedDB(DB.size());
    for (size_t i = 0; i < DB.size(); ++i)
        maskedDB[i] = DB[i] ^ dbMask[i];

    ByteVec seedMask = MGF1(maskedDB, hLen);
    ByteVec maskedSeed(hLen);
    for (size_t i = 0; i < hLen; ++i)
        maskedSeed[i] = r[i] ^ seedMask[i];

    ByteVec EM = {0x00};
    EM.insert(EM.end(), maskedSeed.begin(), maskedSeed.end());
    EM.insert(EM.end(), maskedDB.begin(), maskedDB.end());
    return EM;
}

// OAEP Decode
pair<ByteVec, ByteVec> OAEPDecode(const ByteVec& EM, size_t k) {
    size_t hLen = SHA256_DIGEST_LENGTH;
    assert(EM.size() == k);
    assert(EM[0] == 0x00);

    ByteVec maskedSeed(EM.begin() + 1, EM.begin() + 1 + hLen);
    ByteVec maskedDB(EM.begin() + 1 + hLen, EM.end());

    ByteVec seedMask = MGF1(maskedDB, hLen);
    ByteVec seed(hLen);
    for (size_t i = 0; i < hLen; ++i)
        seed[i] = maskedSeed[i] ^ seedMask[i];

    ByteVec dbMask = MGF1(seed, k - hLen - 1);
    ByteVec DB(maskedDB.size());
    for (size_t i = 0; i < DB.size(); ++i)
        DB[i] = maskedDB[i] ^ dbMask[i];

    unsigned char lHash[hLen];
    SHA256(nullptr, 0, lHash);
    if (!equal(DB.begin(), DB.begin() + hLen, lHash)) {
        throw runtime_error("lHash mismatch");
    }

    size_t i = hLen;
    while (i < DB.size() && DB[i] == 0x00)
        ++i;
    if (i == DB.size() || DB[i] != 0x01)
        throw runtime_error("0x01 separator not found");
    ++i;

    ByteVec msg(DB.begin() + i, DB.end());
    return {msg, seed};
}

// EVP for RSA_NOPADDING Decryption
ByteVec EVP_RSA_NoPadding_Encrypt(EVP_PKEY* pkey, const ByteVec& plaintext) {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new(pkey, nullptr);
    if (!ctx) throw runtime_error("EVP_PKEY_CTX_new failed");

    if (EVP_PKEY_encrypt_init(ctx) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_encrypt_init failed");
    }

    if (EVP_PKEY_CTX_set_rsa_padding(ctx, RSA_NO_PADDING) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_CTX_set_rsa_padding failed");
    }

    size_t outlen = 0;
    if (EVP_PKEY_encrypt(ctx, nullptr, &outlen, plaintext.data(), plaintext.size()) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_encrypt length query failed");
    }

    ByteVec outbuf(outlen);
    if (EVP_PKEY_encrypt(ctx, outbuf.data(), &outlen, plaintext.data(), plaintext.size()) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_encrypt failed");
    }
    outbuf.resize(outlen);

    EVP_PKEY_CTX_free(ctx);
    return outbuf;
}

// EVP for RSA_NOPADDING Encryption
ByteVec EVP_RSA_NoPadding_Decrypt(EVP_PKEY* pkey, const ByteVec& ciphertext) {
    EVP_PKEY_CTX* ctx = EVP_PKEY_CTX_new(pkey, nullptr);
    if (!ctx) throw runtime_error("EVP_PKEY_CTX_new failed");

    if (EVP_PKEY_decrypt_init(ctx) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_decrypt_init failed");
    }

    if (EVP_PKEY_CTX_set_rsa_padding(ctx, RSA_NO_PADDING) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_CTX_set_rsa_padding failed");
    }

    size_t outlen = 0;
    if (EVP_PKEY_decrypt(ctx, nullptr, &outlen, ciphertext.data(), ciphertext.size()) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_decrypt length query failed");
    }

    ByteVec outbuf(outlen);
    if (EVP_PKEY_decrypt(ctx, outbuf.data(), &outlen, ciphertext.data(), ciphertext.size()) <= 0) {
        EVP_PKEY_CTX_free(ctx);
        throw runtime_error("EVP_PKEY_decrypt failed");
    }
    outbuf.resize(outlen);

    EVP_PKEY_CTX_free(ctx);
    return outbuf;
}

void test(size_t msg_len, bool use_ana_r) {
    const int N = 1 << 10;

    double total_time_base = 0.0;
    double total_time_total = 0.0;

    ByteVec msg(msg_len, 0xAB);
    ByteVec covert_msg(16, 0xCC);
    ByteVec covert_key(32, 0x22);
    ByteVec iv(16, 0x77);

    // Key generation
    RSA* rsa = RSA_new();
    BIGNUM* bn = BN_new();
    BN_set_word(bn, RSA_F4);
    RSA_generate_key_ex(rsa, 2048, bn, nullptr);
    BN_free(bn);

    EVP_PKEY* pkey = EVP_PKEY_new();
    EVP_PKEY_assign_RSA(pkey, rsa);

    size_t k = EVP_PKEY_size(pkey);

    for (int i = 0; i < N; ++i) {
        ByteVec r;
        high_resolution_clock::time_point start_r_enc, end_r_enc, start_r_dec, end_r_dec;
        double normal_enc_time = 0, normal_dec_time = 0;


        auto start_enc = high_resolution_clock::now();
        if (use_ana_r) {
            start_r_enc = high_resolution_clock::now();
            ByteVec enc_result = NormalEnc(covert_msg, covert_key, iv);
            r.clear();
            r.insert(r.end(), enc_result.begin(), enc_result.begin() + 16);
            r.insert(r.end(), enc_result.begin() + 16, enc_result.begin() + 32);
            assert(r.size() == SHA256_DIGEST_LENGTH);
            end_r_enc = high_resolution_clock::now();
        } else {
            r = ByteVec(SHA256_DIGEST_LENGTH, 0x42);
        }
        ByteVec encoded = OAEPEncode(msg, r, k);
        ByteVec cipher = EVP_RSA_NoPadding_Encrypt(pkey, encoded);
        auto end_enc = high_resolution_clock::now();


        auto start_dec = high_resolution_clock::now();
        ByteVec decrypted = EVP_RSA_NoPadding_Decrypt(pkey, cipher);
        auto [decoded_msg, extracted_r] = OAEPDecode(decrypted, k);
        if (use_ana_r) {
            ByteVec iv(extracted_r.begin(), extracted_r.begin() + 16);
            ByteVec ct(extracted_r.begin() + 16, extracted_r.end());

            ByteVec recovered(ct.size());
            EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
            EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), nullptr, covert_key.data(), iv.data());
            EVP_CIPHER_CTX_set_padding(ctx, 0);

            int len = 0, total = 0;
            EVP_DecryptUpdate(ctx, recovered.data(), &len, ct.data(), ct.size()); total += len;
            EVP_DecryptFinal_ex(ctx, recovered.data() + total, &len); total += len;
            EVP_CIPHER_CTX_free(ctx);
            recovered.resize(total);
            assert(recovered == covert_msg);
        }
        auto end_dec = high_resolution_clock::now();


        assert(decoded_msg == msg);
        assert(extracted_r == r);

        double time_total = duration_cast<nanoseconds>(end_enc - start_enc + end_dec - start_dec).count();
        total_time_base += time_total;
    }

    EVP_PKEY_free(pkey); // Release

    string r_type = use_ana_r ? "anam" : "normal";
    cout << r_type << "\t" << msg_len << "\t"
         << fixed << setprecision(3) << total_time_base / N / 1000
         << " us" << endl;
}

int main() {
    test(64, true);
    test(64, false);
    test(16, true);
    test(16, false);
    return 0;
}
