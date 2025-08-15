#include <openssl/aes.h>
#include <openssl/evp.h>
#include <cstring>
#include <stdexcept>
#include <chrono>
#include <iostream>
#include <random>
#include "AMAES.h"

using namespace std;
using namespace std::chrono;
using ByteVec = std::vector<unsigned char>;

ByteVec generateRandomBytes(int length) {
    thread_local std::mt19937 gen(std::random_device{}());
    thread_local std::uniform_int_distribution<> dis(0, 255);

    ByteVec bytes(length);
    for (int i = 0; i < length; ++i) {
        bytes[i] = static_cast<unsigned char>(dis(gen));
    }
    return bytes;
}

// Helper: Initialize thread_local AES_KEY
void init_encrypt_key(AES_KEY& aes_key, bool& initialized, const ByteVec& key) {
    if (!initialized) {
        if (AES_set_encrypt_key(key.data(), 256, &aes_key) < 0)
            throw std::runtime_error("AES_set_encrypt_key failed");
        initialized = true;
    }
}

void init_decrypt_key(AES_KEY& aes_key, bool& initialized, const ByteVec& key) {
    if (!initialized) {
        if (AES_set_decrypt_key(key.data(), 256, &aes_key) < 0)
            throw std::runtime_error("AES_set_decrypt_key failed");
        initialized = true;
    }
}

// NormalEnc: output = IV || ciphertext
ByteVec NormalEnc(const ByteVec& plaintext, const ByteVec& key, const ByteVec& iv) {
    if (plaintext.size() % AES_BLOCK_SIZE != 0)
        throw std::runtime_error("NormalEnc: plaintext must be AES_BLOCK_SIZE aligned");

    thread_local AES_KEY aes_key;
    thread_local bool aes_key_initialized = false;

    init_encrypt_key(aes_key, aes_key_initialized, key);

    ByteVec ciphertext(plaintext.size());
    ByteVec iv_copy = iv;

    AES_cbc_encrypt(plaintext.data(), ciphertext.data(), plaintext.size(), &aes_key, iv_copy.data(), AES_ENCRYPT);

    ByteVec output = iv;
    output.insert(output.end(), ciphertext.begin(), ciphertext.end());
    return output;
}

// NormalDec: input = IV || ciphertext
ByteVec NormalDec(const ByteVec& input, const ByteVec& key) {
    if (input.size() < AES_BLOCK_SIZE || (input.size() - AES_BLOCK_SIZE) % AES_BLOCK_SIZE != 0)
        throw std::runtime_error("NormalDec: invalid input size");

    thread_local AES_KEY aes_key;
    thread_local bool aes_key_initialized = false;

    init_decrypt_key(aes_key, aes_key_initialized, key);

    ByteVec iv(input.begin(), input.begin() + AES_BLOCK_SIZE);
    ByteVec ciphertext(input.begin() + AES_BLOCK_SIZE, input.end());
    ByteVec plaintext(ciphertext.size());

    AES_cbc_encrypt(ciphertext.data(), plaintext.data(), ciphertext.size(), &aes_key, iv.data(), AES_DECRYPT);
    return plaintext;
}

// AnamorphicEnc: output = IV2 || ciphertext
ByteVec AnamorphicEnc(const ByteVec& normal_plain, const ByteVec& normal_key,
                      const ByteVec& covert_plain, const ByteVec& covert_key,
                      const ByteVec& IV1) {
    if (normal_plain.size() % AES_BLOCK_SIZE != 0 || covert_plain.size() % AES_BLOCK_SIZE != 0)
        throw std::runtime_error("AnamorphicEnc: inputs must be AES block aligned");

    thread_local AES_KEY covert_aes_key;
    thread_local bool covert_aes_key_initialized = false;
    init_encrypt_key(covert_aes_key, covert_aes_key_initialized, covert_key);

    thread_local AES_KEY normal_aes_key;
    thread_local bool normal_aes_key_initialized = false;
    init_encrypt_key(normal_aes_key, normal_aes_key_initialized, normal_key);

    // Step 1: covert_plain → IV2 = Enc_CBC(covert_key, IV1)
    ByteVec IV2(covert_plain.size());
    ByteVec iv1_copy = IV1;
    AES_cbc_encrypt(covert_plain.data(), IV2.data(), covert_plain.size(), &covert_aes_key, iv1_copy.data(), AES_ENCRYPT);

    // Step 2: normal_plain → ciphertext = Enc_CBC(normal_key, IV2)
    ByteVec ciphertext(normal_plain.size());
    ByteVec iv2_copy = IV2;
    AES_cbc_encrypt(normal_plain.data(), ciphertext.data(), normal_plain.size(), &normal_aes_key, iv2_copy.data(), AES_ENCRYPT);

    // Output = IV2 || ciphertext
    ByteVec output = IV2;
    output.insert(output.end(), ciphertext.begin(), ciphertext.end());
    return output;
}

// AnamorphicDec: input = IV2 || ciphertext
std::pair<ByteVec, ByteVec> AnamorphicDec(const ByteVec& input,
                                          const ByteVec& normal_key,
                                          const ByteVec& covert_key,
                                          const ByteVec& IV1) {
    if (input.size() < AES_BLOCK_SIZE || input.size() % AES_BLOCK_SIZE != 0)
        throw std::runtime_error("AnamorphicDec: invalid input size");

    thread_local AES_KEY normal_aes_key;
    thread_local bool normal_aes_key_initialized = false;
    init_decrypt_key(normal_aes_key, normal_aes_key_initialized, normal_key);

    thread_local AES_KEY covert_aes_key;
    thread_local bool covert_aes_key_initialized = false;
    init_decrypt_key(covert_aes_key, covert_aes_key_initialized, covert_key);

    // Step 1: extract IV2 and ciphertext
    ByteVec IV2(input.begin(), input.begin() + AES_BLOCK_SIZE);
    ByteVec ciphertext(input.begin() + AES_BLOCK_SIZE, input.end());

    // Step 2: ciphertext → normal_plain using AES-CBC(normal_key, IV2)
    ByteVec normal_plain(ciphertext.size());
    ByteVec iv2_copy = IV2;
    AES_cbc_encrypt(ciphertext.data(), normal_plain.data(), ciphertext.size(), &normal_aes_key, iv2_copy.data(), AES_DECRYPT);

    // Step 3: IV2 → covert_plain using AES-CBC(covert_key, IV1)
    ByteVec covert_plain(IV2.size());
    ByteVec iv1_copy = IV1;
    AES_cbc_encrypt(IV2.data(), covert_plain.data(), IV2.size(), &covert_aes_key, iv1_copy.data(), AES_DECRYPT);

    return {normal_plain, covert_plain};
}


/*
ByteVec NormalEnc(const ByteVec& plaintext, const ByteVec& key, const ByteVec& iv) {
    thread_local EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (!ctx) throw std::runtime_error("EVP_CIPHER_CTX_new failed");
    ByteVec ciphertext(plaintext.size() + AES_BLOCK_SIZE);
    int len = 0, total_len = 0;
    EVP_CIPHER_CTX_reset(ctx);
    if (EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), nullptr, key.data(), iv.data()) != 1)
        throw std::runtime_error("EVP_EncryptInit_ex failed");

    if (EVP_EncryptUpdate(ctx, ciphertext.data(), &len, plaintext.data(), plaintext.size()) != 1)
        throw std::runtime_error("EVP_EncryptUpdate failed");
    total_len = len;

    if (EVP_EncryptFinal_ex(ctx, ciphertext.data() + len, &len) != 1)
        throw std::runtime_error("EVP_EncryptFinal_ex failed");
    total_len += len;

    ciphertext.resize(total_len);
    ByteVec output = iv;
    output.insert(output.end(), ciphertext.begin(), ciphertext.end());
    return output;
}

ByteVec NormalDec(const ByteVec& input, const ByteVec& key) {
    thread_local EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    if (input.size() < AES_BLOCK_SIZE)
        throw std::runtime_error("Input too short");
    ByteVec iv(input.begin(), input.begin() + AES_BLOCK_SIZE);
    ByteVec ciphertext(input.begin() + AES_BLOCK_SIZE, input.end());
    if (!ctx) throw std::runtime_error("EVP_CIPHER_CTX_new failed");
    ByteVec plaintext(ciphertext.size()); // slightly overallocated
    int len = 0, total_len = 0;
    EVP_CIPHER_CTX_reset(ctx);
    if (EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), nullptr, key.data(), iv.data()) != 1)
        throw std::runtime_error("EVP_DecryptInit_ex failed");

    if (EVP_DecryptUpdate(ctx, plaintext.data(), &len, ciphertext.data(), ciphertext.size()) != 1)
        throw std::runtime_error("EVP_DecryptUpdate failed");
    total_len = len;

    if (EVP_DecryptFinal_ex(ctx, plaintext.data() + len, &len) != 1)
        throw std::runtime_error("EVP_DecryptFinal_ex failed");
    total_len += len;

    plaintext.resize(total_len);
    return plaintext;
}

ByteVec AnamorphicEnc(const ByteVec& normal_plain, const ByteVec& normal_key,
                      const ByteVec& covert_plain, const ByteVec& covert_key,
                      const ByteVec& IV1) {
    thread_local EVP_CIPHER_CTX* ctx_covert = EVP_CIPHER_CTX_new();
    thread_local EVP_CIPHER_CTX* ctx_normal = EVP_CIPHER_CTX_new();

    // Step 1: Encrypt covert_plain with AES-CTR using covert_key and IV1 -> IV2
    ByteVec IV2(covert_plain.size(), 0);
    int len = 0;
    EVP_CIPHER_CTX_reset(ctx_covert);
    EVP_EncryptInit_ex(ctx_covert, EVP_aes_256_ctr(), nullptr, covert_key.data(), IV1.data());
    EVP_CIPHER_CTX_set_padding(ctx_covert, 0);
    EVP_EncryptUpdate(ctx_covert, IV2.data(), &len, covert_plain.data(), covert_plain.size());

    // Step 2: Encrypt normal_plain with AES-CBC using normal_key and IV2
    ByteVec ciphertext(normal_plain.size(), 0);
    EVP_CIPHER_CTX_reset(ctx_normal);
    EVP_EncryptInit_ex(ctx_normal, EVP_aes_256_cbc(), nullptr, normal_key.data(), IV2.data());
    EVP_CIPHER_CTX_set_padding(ctx_normal, 0);
    EVP_EncryptUpdate(ctx_normal, ciphertext.data(), &len, normal_plain.data(), normal_plain.size());

    // Combine
    ByteVec output = IV2;
    output.insert(output.end(), ciphertext.begin(), ciphertext.end());
    return output;
}

std::pair<ByteVec, ByteVec> AnamorphicDec(const ByteVec& input,
                                          const ByteVec& normal_key,
                                          const ByteVec& covert_key,
                                          const ByteVec& IV1) {
    thread_local EVP_CIPHER_CTX* ctx_normal = EVP_CIPHER_CTX_new();
    thread_local EVP_CIPHER_CTX* ctx_covert = EVP_CIPHER_CTX_new();

    ByteVec IV2(input.begin(), input.begin() + AES_BLOCK_SIZE);
    ByteVec ciphertext(input.begin() + AES_BLOCK_SIZE, input.end());

    // Step 1: Decrypt ciphertext with AES-CBC using normal_key and IV2
    ByteVec normal_plain(ciphertext.size(), 0);
    int len = 0;
    EVP_CIPHER_CTX_reset(ctx_normal);
    EVP_DecryptInit_ex(ctx_normal, EVP_aes_256_cbc(), nullptr, normal_key.data(), IV2.data());
    EVP_CIPHER_CTX_set_padding(ctx_normal, 0);
    EVP_DecryptUpdate(ctx_normal, normal_plain.data(), &len, ciphertext.data(), ciphertext.size());

    // Step 2: Decrypt IV2 with AES-CTR using covert_key and IV1
    ByteVec covert_plain(AES_BLOCK_SIZE, 0);
    EVP_CIPHER_CTX_reset(ctx_covert);
    EVP_DecryptInit_ex(ctx_covert, EVP_aes_256_ctr(), nullptr, covert_key.data(), IV1.data());
    EVP_CIPHER_CTX_set_padding(ctx_covert, 0);
    EVP_DecryptUpdate(ctx_covert, covert_plain.data(), &len, IV2.data(), AES_BLOCK_SIZE);

    return {normal_plain, covert_plain};
}
*/

/*
ByteVec AnamorphicEnc(const ByteVec& normal_plain, const ByteVec& normal_key,
                      const ByteVec& covert_plain, const ByteVec& covert_key,
                      const ByteVec& IV1) {
    // Step 1: Encrypt covert_plain with AES-CBC using covert_key and IV1 -> IV2
    ByteVec IV2(covert_plain.size());
    EVP_CIPHER_CTX* ctx1 = EVP_CIPHER_CTX_new();
    EVP_EncryptInit_ex(ctx1, EVP_aes_256_ctr(), nullptr, covert_key.data(), IV1.data());
    EVP_CIPHER_CTX_set_padding(ctx1, 0); // Disable padding

    int len1 = 0;
    EVP_EncryptUpdate(ctx1, IV2.data(), &len1, covert_plain.data(), covert_plain.size());
    EVP_EncryptFinal_ex(ctx1, IV2.data() + len1, &len1);
    EVP_CIPHER_CTX_free(ctx1);

    // Step 2: Encrypt normal_plain with AES-CBC using normal_key and IV2
    ByteVec ciphertext(normal_plain.size());
    EVP_CIPHER_CTX* ctx2 = EVP_CIPHER_CTX_new();
    EVP_EncryptInit_ex(ctx2, EVP_aes_256_cbc(), nullptr, normal_key.data(), IV2.data());
    EVP_CIPHER_CTX_set_padding(ctx2, 0);

    int len2 = 0;
    EVP_EncryptUpdate(ctx2, ciphertext.data(), &len2, normal_plain.data(), normal_plain.size());
    EVP_EncryptFinal_ex(ctx2, ciphertext.data() + len2, &len2);
    EVP_CIPHER_CTX_free(ctx2);

    ByteVec output = IV2;
    output.insert(output.end(), ciphertext.begin(), ciphertext.end());
    return output;
}

std::pair<ByteVec, ByteVec> AnamorphicDec(const ByteVec& input,
                                          const ByteVec& normal_key,
                                          const ByteVec& covert_key,
                                          const ByteVec& IV1) {
    ByteVec IV2(input.begin(), input.begin() + AES_BLOCK_SIZE);
    ByteVec ciphertext(input.begin() + AES_BLOCK_SIZE, input.end());

    // Step 1: Decrypt ciphertext with AES-CBC using normal_key and IV2
    ByteVec normal_plain(ciphertext.size());
    EVP_CIPHER_CTX* ctx1 = EVP_CIPHER_CTX_new();
    EVP_DecryptInit_ex(ctx1, EVP_aes_256_cbc(), nullptr, normal_key.data(), IV2.data());
    EVP_CIPHER_CTX_set_padding(ctx1, 0);

    int len1 = 0;
    EVP_DecryptUpdate(ctx1, normal_plain.data(), &len1, ciphertext.data(), ciphertext.size());
    EVP_DecryptFinal_ex(ctx1, normal_plain.data() + len1, &len1);
    EVP_CIPHER_CTX_free(ctx1);

    // Step 2: Decrypt IV2 with AES-CBC using covert_key and IV1 to get covert_plain
    ByteVec covert_plain(AES_BLOCK_SIZE);
    EVP_CIPHER_CTX* ctx2 = EVP_CIPHER_CTX_new();
    EVP_DecryptInit_ex(ctx2, EVP_aes_256_ctr(), nullptr, covert_key.data(), IV1.data());
    EVP_CIPHER_CTX_set_padding(ctx2, 0);

    int len2 = 0;
    EVP_DecryptUpdate(ctx2, covert_plain.data(), &len2, IV2.data(), AES_BLOCK_SIZE);
    EVP_DecryptFinal_ex(ctx2, covert_plain.data() + len2, &len2);
    EVP_CIPHER_CTX_free(ctx2);

    return {normal_plain, covert_plain};
}
*/