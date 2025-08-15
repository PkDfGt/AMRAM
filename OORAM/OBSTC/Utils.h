// Utils.h
#pragma once
#include <vector>
#include <stack>

struct StackFrame {
    int offset;
    int len;
    int z;
    int stage; // 0:Init, 1:Left, 2:Right
    int mleft;
    int mright;
    int merge_i; // Trace merge stage
};

// Tree-based Compact (TC)
void generateTCAddr(int offset, int len, std::vector<int>& output);
bool verifyTC(const std::vector<int>& compactedM);
bool compactIte(int& ele1, int& mark1, int& ele2, int& mark2, std::stack<StackFrame>& stack, int& result);

// Bitonic Sort
void generateBSAddr(int N, std::vector<int>& output);
bool BSIte(int N, int& ele1, int& ele2, int& i, int& j, int& k, bool ascending);
bool verifySorted(const std::vector<int>& A, bool ascending);
