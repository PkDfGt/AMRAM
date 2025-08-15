#ifndef AMAES_H
#define AMAES_H

#include <vector>
#include <utility>

using ByteVec = std::vector<unsigned char>;

ByteVec generateRandomBytes(int length);
ByteVec NormalEnc(const ByteVec& plaintext, const ByteVec& key, const ByteVec& iv);
ByteVec NormalDec(const ByteVec& input, const ByteVec& key);

ByteVec AnamorphicEnc(const ByteVec& normal_plain, const ByteVec& normal_key,
                         const ByteVec& covert_plain, const ByteVec& covert_key,
                         const ByteVec& IV1);

std::pair<ByteVec, ByteVec> AnamorphicDec(const ByteVec& input,
                                             const ByteVec& normal_key,
                                             const ByteVec& covert_key,
                                             const ByteVec& IV1);

#endif // AMAES_H
