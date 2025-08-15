#include "nOPQ.h"
#include <random>

nOPQ::nOPQ(int N_, int Z_, int Zr_)
    : N(N_), Z(Z_), Zr(Zr_) 
{
    L = static_cast<int>(std::ceil(std::log2(N_))) + 1;
    tree_size = (1 << L) - 1;

    tree.resize(tree_size);
    cmins.resize(tree_size, {INT_MAX, 0});

    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        tree[i].resize(cap, {0, 0});
    }
}

void nOPQ::insert(int p, int pos) {
    int leaf1 = nOPQ::generate_random_leaf(0, (N >> 1) - 1);
    int leaf2 = nOPQ::generate_random_leaf(N >> 1, N - 1);

    std::vector<int> combined_nodes = get_two_path_indices(leaf1, leaf2);
    std::vector<PlainEle> elements = read_elements_from_path(combined_nodes);
    elements.push_back({p, pos});

    //Evict and Updatemin
    evict_and_update_min(combined_nodes, elements);
}

nOPQ::PlainEle nOPQ::extract_min() {
    PlainEle min_ele = cmins[0];
    std::vector<int> path = get_path_indices(min_ele.pos);

    std::vector<PlainEle> elements = read_elements_from_path(path, true, min_ele);

    evict_and_update_min(path, elements);
    return min_ele;
}

int nOPQ::generate_random_leaf(int b, int e) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(b, e);
    return dis(gen);
}

std::vector<int> nOPQ::get_path_indices(int leaf) const {
    int index = leaf + N - 1;
    int i = 0;
    std::vector<int> path(L);
    while (index > 0) {
        path[i++] = index;
        index = (index - 1) >> 1;
    }
    path[i++] = index;
    assert(i == L);
    return path;
}

std::vector<int> nOPQ::get_two_path_indices(int leaf1, int leaf2) const {
    int index1 = leaf1 + N - 1;
    int index2 = leaf2 + N - 1;
    int i = 0;
    std::vector<int> path((L << 1) - 1);
    while (index1 > 0) {
        path[i++] = index1;
        path[i++] = index2;
        index1 = (index1 - 1) >> 1;
        index2 = (index2 - 1) >> 1;
    }
    path[i++] = index1;
    assert(i == (L << 1) - 1);
    return path;
}

bool nOPQ::check_common_path(int leaf, int idx) const {
    int index = leaf + N - 1;
    while (index >= 0) {
        if (index == idx) return true;
        else if (index < idx) return false;
        index = (index - 1) >> 1;
    }
    return false;
}

std::vector<nOPQ::PlainEle> nOPQ::read_elements_from_path(const std::vector<int>& path, bool skip_match, const PlainEle& match) {
    std::vector<PlainEle> elements;
    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            PlainEle e = tree[idx][j];
            if (e.p != 0 && (!skip_match || e.p != match.p)) {
                elements.push_back(e);
            }
        }
    }
    return elements;
}

void nOPQ::evict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& elements) {
    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;

        std::vector<PlainEle> bucket;
        int write_pos = 0;
        for (PlainEle& e : elements) {
            if ((int)bucket.size() < cap && check_common_path(e.pos, idx)) {
                bucket.push_back(e);
            } else {
                elements[write_pos++] = e;
            }
        }
        elements.resize(write_pos);

        PlainEle min_e = {INT_MAX, 0};
        int left = 2 * idx + 1, right = 2 * idx + 2;
        if (left < tree_size && cmins[left].p < min_e.p && cmins[left].p != 0)
            min_e = cmins[left];
        if (right < tree_size && cmins[right].p < min_e.p && cmins[right].p != 0)
            min_e = cmins[right];

        int inserted = 0;
        for (const auto& e : bucket) {
            if (e.p != 0 && e.p < min_e.p) min_e = e;
            tree[idx][inserted++] = e; 
        }
        for (; inserted < cap; ++inserted) {
            tree[idx][inserted] = {0, 0};  // dummy
        }

        cmins[idx] = min_e;
    }

    if (!elements.empty()) {
        int old_Zr = Zr;
        Zr += elements.size();
        tree[0].resize(Zr);
        int k = old_Zr;
        for (const auto& e : elements) {
            tree[0][k++] = e;
        }
        elements.clear();

        PlainEle min_e = cmins[0];
        for (int i = old_Zr; i < Zr; ++i) {
            if (tree[0][i].p < min_e.p && tree[0][i].p != 0)
                min_e = tree[0][i];
        }
        cmins[0] = min_e;
    }
}

int nOPQ::get_rootused() { 
    int root_used = 0;
    for (const auto& e : tree[0]) {
        if (e.p != 0) ++root_used;
    }
    return root_used;
}

void nOPQ::print_tree() const {
    std::cout << "====== nOPQ Tree ======\n";
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        std::cout << "Node " << i << ": [";
        for (int j = 0; j < cap; ++j) {
            const PlainEle& e = tree[i][j];
            if (e.p == 0) std::cout << "d";
            else std::cout << "(p=" << e.p << ",pos=" << e.pos << ")";
            if (j != cap - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "=== cmins ===\n";
    for (int i = 0; i < tree_size; ++i) {
        const PlainEle& e = cmins[i];
        if (e.p == 0 || e.p == INT_MAX) std::cout << "cmins[" << i << "]: d\n";
        else std::cout << "cmins[" << i << "]: (p=" << e.p << ",pos=" << e.pos << ")\n";
    }
}
