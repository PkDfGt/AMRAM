#ifndef NORMAL_BS_H
#define NORMAL_BS_H

#include <vector>
#include <cstdint>

using ByteVec = std::vector<unsigned char>;

extern int total_comps;
extern int curr_comps;
extern int BSize;

void printProgressBar(float progress);
std::pair<int, int> decodeInt(const ByteVec& data);
ByteVec encodeInt(int value, int mark, int size);
int calcTotalComparisons(int n);
int normalCompact(std::vector<ByteVec>& A, const ByteVec& normal_key,
                    int z, int offset, int len);
std::pair<int, int> anamorphicCompact(std::vector<ByteVec>& A,
                                      const ByteVec& normal_key,
                                      const ByteVec& covert_key,
                                      int nz,  // normal z
                                      int cz,  // covert z
                                      int offset,
                                      int len,
                                      std::vector<int>& write_count);
ByteVec deriveIV(const ByteVec& key, int addr, int count);
#endif  // NORMAL_BS_H
