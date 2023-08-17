#include <vector>

template <typename T>
class DST {
public:
    DST(std::vector<T> & input, T neutral) : input(std::move(input)),
        tree(std::__lg(this->input.size()), std::vector<T>(this->input.size(), neutral)),
        h(std::__lg(this->input.size()) - 1)
    {
        buildDST(0, this->input.size() - 1);
    }

    T queryDST(int l, int r) {
        int level = h - std::__lg(l ^ r);
        return std::min(tree[level][l], tree[level][r]);
    }

private:
    std::vector<T> input;
    std::vector<std::vector<T>> tree;
    int h;

    void buildDST(int l, int r, int level = 0) {
        int mid = (l + r) / 2;
        for (int i = mid; i > l - 1; i--) {
            tree[level][i] = std::min(tree[level][i + 1], input[i]);
        }
        tree[level][mid + 1] = input[mid + 1];
        for (int i = mid + 2; i < r + 1; i++) {
            tree[level][i] = std::min(tree[level][i - 1], input[i]);
        }
        if (l < mid) {
            buildDST(l, mid, level + 1);
            buildDST(mid + 1, r, level + 1);
        }
    }
};
