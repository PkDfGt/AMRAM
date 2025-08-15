// Utils.cpp
#include "Utils.h"
#include <iostream>
#include <cassert>
#include <algorithm>
#include <stack>

// =====================
// Tree-based Compact (TC)
// =====================

void generateTCAddr(int offset, int len, std::vector<int>& output) {
    if (len == 2) {
        output.push_back(offset);
        output.push_back(offset+1);
        return;
    }

    int half = len / 2;
    generateTCAddr(offset, half, output);
    generateTCAddr(offset+half, half, output);
    for (int i=0; i<half; ++i) {
        output.push_back(offset+i);
        output.push_back(offset+half+i);
    }
}

bool verifyTC(const std::vector<int>& compactedM) {
    bool seen_zero = false;
    for (size_t i = 0; i < compactedM.size(); ++i) {
        if (compactedM[i] == 0) {
            seen_zero = true;
        } else if (seen_zero && compactedM[i] == 1) {
            std::cout << "❌ Verify failed at position " << i << "\n";
            return false;
        }
    }
    std::cout << "✅ Verify passed: M has all 1s in front and 0s at the end.\n";
    return true;
}

bool compactIte(int& ele1, int& mark1, int& ele2, int& mark2, std::stack<StackFrame>& stack, int& result) {
    bool cond = false;
    if (stack.empty()) {
        return cond;
    }
    while (!stack.empty()) {
        auto& frame = stack.top();
        
        if (frame.len == 2) {
            // Base condition
            int xi = ele1, mi = mark1;
            int xj = ele2, mj = mark2;

            cond = ((1 - mi) * mj) ^ frame.z;
            int mask = -cond;

            // Branchless swap
            int tmp_x = (xi ^ xj) & mask;
            xi ^= tmp_x;
            xj ^= tmp_x;

            int tmp_m = (mi ^ mj) & mask;
            mi ^= tmp_m;
            mj ^= tmp_m;

            ele1 = xi;  ele2 = xj;
            mark1 = mi;  mark2 = mj;

            result = mi + mj;
            stack.pop();
            return cond;
        }
        
        int half = frame.len / 2;
        
        if (frame.stage == 0) {
            // Left half
            frame.stage = 1;
            stack.push({frame.offset, half, frame.z % half, 0, 0, 0, 0});
            continue;
        } else if (frame.stage == 1) {
            // Right half
            frame.mleft = result;
            frame.stage = 2;
            stack.push({frame.offset + half, half, (frame.z + frame.mleft) % half, 0, 0, 0, 0});
            continue;
        } else if (frame.stage == 2) {
            // Merge
            if (frame.merge_i == 0) {
                frame.mright = result; // Set for first
            }
            
            if (frame.merge_i < half) {
                int i = frame.merge_i++;

                int xi = ele1, mi = mark1;
                int xj = ele2, mj = mark2;

                int s = (((frame.z % half) + frame.mleft >= half) ? 1 : 0) ^ ((frame.z >= half) ? 1 : 0);
                cond = s ^ (i >= ((frame.z + frame.mleft) % half));
                int mask = -cond;

                // Branchless swap
                int tmp_x = (xi ^ xj) & mask;
                xi ^= tmp_x;
                xj ^= tmp_x;

                int tmp_m = (mi ^ mj) & mask;
                mi ^= tmp_m;
                mj ^= tmp_m;

                ele1 = xi;  ele2 = xj;
                mark1 = mi;  mark2 = mj; 

                if (frame.merge_i == half) {
                    result = frame.mleft + frame.mright;
                    stack.pop();
                }
                return cond;
            }
        }
        
        stack.pop();
    }
    
    return cond;
}

// =====================
// Bitonic Sort
// =====================

void generateBSAddr(int N, std::vector<int>& output) {
    assert((N & (N - 1)) == 0);
    for (int k = 2; k <= N; k <<= 1) {
        for (int j = k >> 1; j > 0; j >>= 1) {
            for (int i = 0; i < N; ++i) {
                int ixj = i ^ j;
                if (ixj > i) {
                    output.push_back(i);
                    output.push_back(ixj);
                }
            }
        }
    }
}

bool BSIte(int N, int& ele1, int& ele2, int& i, int& j, int& k, bool ascending) {
    assert((N & (N - 1)) == 0 && N > 1);

    // Compare
    bool dir = ((i & k) == 0) ? ascending : !ascending;
    bool need_swap = (dir && ele1 > ele2) || (!dir && ele1 < ele2);
    int mask = -static_cast<int>(need_swap);
    int tmp = ele1 ^ ele2;
    tmp &= mask;
    ele1 ^= tmp;
    ele2 ^= tmp;
    
    do {
        i++;
        if (i == N) {
            i = 0;
            j >>= 1;
            if (j == 0) {
                k <<= 1;
                j = k >> 1;
                if (k > N) {
                    k = 2;
                    j = 1;
                    i = 0;
                }
            }
        }
    } while ((i ^ j) <= i && k <= N);

    return need_swap;
}

bool verifySorted(const std::vector<int>& A, bool ascending) {
    for (size_t i = 1; i < A.size(); ++i) {
        if ((A[i] < A[i - 1] && ascending) || (A[i] > A[i - 1] && !ascending)) {
            std::cout << "❌ Verify failed at position " << i
                      << ": previous=" << A[i - 1] << ", current=" << A[i] << "\n";
            return false;
        }
    }
    std::cout << "✅ Verify passed: sequence is sorted in "
              << (ascending ? "ascending" : "descending") << " order.\n";
    return true;
}
