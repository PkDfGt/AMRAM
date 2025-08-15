#ifndef NORMAL_BS_H
#define NORMAL_BS_H

#include <vector>
#include <cstdint>

using ByteVec = std::vector<unsigned char>;

extern int total_comps;
extern int curr_comps;
extern int BSize;

void printProgressBar(float progress);
int decodeInt(const ByteVec& data);
ByteVec encodeInt(int value, int size);
int calcTotalComparisons(int n);
void normalCompAndSwap(std::vector<ByteVec>& A, int i, int j, bool ascending,
                 const ByteVec& key);
void normalBitonicSort(std::vector<ByteVec>& A, bool ascending, const ByteVec& key);
void anamorphiCompAndSwap(std::vector<ByteVec>& A, int i, int j,
                 bool normal_asc, bool covert_asc,
                 const ByteVec& normal_key, const ByteVec& covert_key,
                 std::vector<int>& write_count);
void anamorphicBitonicSort(std::vector<ByteVec>& A, bool normal_asc,
                 bool covert_asc, const ByteVec& normal_key, 
                 const ByteVec& covert_key, std::vector<int>& write_count);
ByteVec deriveIV(const ByteVec& key, int addr, int count);
#endif  // NORMAL_BS_H
