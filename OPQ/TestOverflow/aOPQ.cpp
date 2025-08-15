#include "aOPQ.h"
#include <random>

aOPQ::aOPQ(int N_, int Z_, int Zr_)
    : N(N_), Z(Z_), Zr(Zr_) 
{
    L = static_cast<int>(std::ceil(std::log2(N_))) + 1;
    tree_size = (1 << L) - 1;

    tree.resize(tree_size);
    cmins.resize(tree_size, {PlainEle{0, 0}, PlainEle{0, 0}});

    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z; 
        tree[i].resize(cap, {PlainEle{0, 0}, PlainEle{0, 0}});  // {normal, anamorphic}
    }
}

void aOPQ::ainsert(int np, int npos, int ap, int apos) {
    int leaf1 = generate_random_leaf(0, N / 2 - 1);
    int leaf2 = generate_random_leaf(N / 2, N - 1);
    std::vector<int> combined_nodes = get_two_path_indices(leaf1, leaf2);

    auto [normal_elements, anamorphic_elements] = aread_elements_from_path(combined_nodes); 
    normal_elements.push_back({np, npos});
    anamorphic_elements.push_back({ap, apos});

    aevict_and_update_min(combined_nodes, normal_elements, anamorphic_elements);
}

std::pair<aOPQ::PlainEle, aOPQ::PlainEle> aOPQ::aextract_min() {
    auto [min_normal, min_anamorphic] = cmins[0];
    std::vector<int> path = get_path_indices(min_normal.pos);
    auto [normal_elements, anamorphic_elements] = aread_elements_from_path(path, true, min_normal);

    aevict_and_update_min(path, normal_elements, anamorphic_elements);

    return {min_normal, min_anamorphic};
}

std::pair<std::vector<aOPQ::PlainEle>, std::vector<aOPQ::PlainEle>> aOPQ::aread_elements_from_path(
    const std::vector<int>& path, bool skip_match, const PlainEle& match) 
{
    std::vector<PlainEle> normal_elements;
    std::vector<PlainEle> anamorphic_elements;

    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        for (int j = 0; j < cap; ++j) {
            const PlainEle& normal_elem = tree[idx][j].first;
            if (normal_elem.p != 0) {
                if (!skip_match || normal_elem.p != match.p) {
                    normal_elements.push_back(normal_elem);
                }
            }
            const PlainEle& anamorphic_elem = tree[idx][j].second;
            if (anamorphic_elem.p != 0) {
                if (!skip_match || anamorphic_elem.p != match.p) {
                    anamorphic_elements.push_back(anamorphic_elem);
                }
            }
        }
    }

    return {normal_elements, anamorphic_elements};
}

void aOPQ::aevict_and_update_min(const std::vector<int>& path, std::vector<PlainEle>& normal_elements, std::vector<PlainEle>& anamorphic_elements) {
    for (int idx : path) {
        int cap = (idx == 0) ? Zr : Z;
        int count = 0;
        std::vector<PlainEle> normal_bucket, anamorphic_bucket;

        int write_pos = 0, read_pos = 0;
        for (; read_pos < (int)normal_elements.size() && count < cap; ++read_pos) {
            const PlainEle& e = normal_elements[read_pos];
            if (check_common_path(e.pos, idx)) {
                normal_bucket.push_back(e);
                ++count;
            } else {
                normal_elements[write_pos++] = e;
            }
        }
        for (int k = read_pos; k < (int)normal_elements.size(); ++k) {
            normal_elements[write_pos++] = normal_elements[k];
        }
        normal_elements.resize(write_pos);

        write_pos = read_pos = 0;
        for (; read_pos < (int)anamorphic_elements.size() && count < cap; ++read_pos) {
            const PlainEle& e = anamorphic_elements[read_pos];
            if (check_common_path(e.pos, idx)) {
                anamorphic_bucket.push_back(e);
                ++count;
            } else {
                anamorphic_elements[write_pos++] = e;
            }
        }
        for (int k = read_pos; k < (int)anamorphic_elements.size(); ++k) {
            anamorphic_elements[write_pos++] = anamorphic_elements[k];
        }
        anamorphic_elements.resize(write_pos);

        PlainEle min_n = {INT_MAX, 0};
        PlainEle min_a = {INT_MAX, 0};
        int left = 2 * idx + 1, right = 2 * idx + 2;

        if (left < tree_size) {
            const PlainEle& ln = cmins[left].first;
            if (ln.p != 0 && ln.p < min_n.p) min_n = ln;

            const PlainEle& la = cmins[left].second;
            if (la.p != 0 && la.p < min_a.p) min_a = la;
        }

        if (right < tree_size) {
            const PlainEle& rn = cmins[right].first;
            if (rn.p != 0 && rn.p < min_n.p) min_n = rn;

            const PlainEle& ra = cmins[right].second;
            if (ra.p != 0 && ra.p < min_a.p) min_a = ra;
        }

        int inserted = 0;
        for (const auto& e : normal_bucket) {
            if (e.p != 0 && e.p < min_n.p) min_n = e;
            tree[idx][inserted++] = {e, PlainEle{0, 0}};
        }

        for (const auto& e : anamorphic_bucket) {
            if (e.p != 0 && e.p < min_a.p) min_a = e;
            tree[idx][inserted++] = {PlainEle{0, 0}, e};
        }

        for (; inserted < cap; ++inserted) {
            tree[idx][inserted] = {PlainEle{0, 0}, PlainEle{0, 0}};  // dummy value
        }

        cmins[idx] = {min_n, min_a};
    }

    if (!normal_elements.empty() || !anamorphic_elements.empty()) {
        int old_Zr = Zr;
        Zr += normal_elements.size() + anamorphic_elements.size();
        tree[0].resize(Zr);
        int k = old_Zr;
        PlainEle min_n = cmins[0].first;
        PlainEle min_a = cmins[0].second;
        for (const auto& e : normal_elements) {
            if (e.p != 0 && e.p < min_n.p) min_n = e;
            tree[0][k++] = {e, PlainEle{0, 0}};
        }
        for (const auto& e : anamorphic_elements) {
            if (e.p != 0 && e.p < min_a.p) min_a = e;
            tree[0][k++] = {PlainEle{0, 0}, e};
        }
        normal_elements.clear();
        anamorphic_elements.clear();

        cmins[0] = {min_n, min_a};
    }
}

int aOPQ::generate_random_leaf(int b, int e) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(b, e);
    return dis(gen);
}

std::vector<int> aOPQ::get_path_indices(int leaf) const {
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

std::vector<int> aOPQ::get_two_path_indices(int leaf1, int leaf2) const {
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

bool aOPQ::check_common_path(int leaf, int idx) const {
    int index = leaf + N - 1;
    while (index >= 0) {
        if (index == idx) return true;
        else if (index < idx) return false;
        index = (index - 1) >> 1;
    }
    return false;
}

int aOPQ::get_rootused() {
    int root_used = 0;
    for (const auto& elem_pair : tree[0]) {
        if (elem_pair.first.p != 0 || elem_pair.second.p != 0) ++root_used;     // anamorphic
    }
    return root_used;
}

void aOPQ::print_tree() const {
    std::cout << "====== aOPQ Tree ======\n";
    for (int i = 0; i < tree_size; ++i) {
        int cap = (i == 0) ? Zr : Z;
        std::cout << "Node " << i << ": [";
        for (int j = 0; j < cap; ++j) {
            const auto& elem_pair = tree[i][j];
            const PlainEle& normal = elem_pair.first;
            const PlainEle& anam = elem_pair.second;

            std::cout << "(";
            if (normal.p == 0) std::cout << "d";
            else std::cout << "n(p=" << normal.p << ",pos=" << normal.pos << ")";
            std::cout << "|";
            if (anam.p == 0) std::cout << "d";
            else std::cout << "a(p=" << anam.p << ",pos=" << anam.pos << ")";
            std::cout << ")";

            if (j != cap - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    std::cout << "=== cmins ===\n";
    for (int i = 0; i < tree_size; ++i) {
        const auto& min_pair = cmins[i];
        const PlainEle& normal_min = min_pair.first;
        const PlainEle& anam_min = min_pair.second;

        std::cout << "cmins[" << i << "]: ";
        if (normal_min.p == 0 || normal_min.p == INT_MAX)
            std::cout << "n: d";
        else
            std::cout << "n: (p=" << normal_min.p << ",pos=" << normal_min.pos << ")";

        std::cout << " | ";

        if (anam_min.p == 0 || anam_min.p == INT_MAX)
            std::cout << "a: d";
        else
            std::cout << "a: (p=" << anam_min.p << ",pos=" << anam_min.pos << ")";

        std::cout << "\n";
    }
}
