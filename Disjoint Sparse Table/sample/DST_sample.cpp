#include <iostream>
#include <vector>

void buildDST(std::vector<std::vector<int>> & tree, std::vector<int> & input, int l, int r, int level = 0) {
    int mid = (l + r) / 2;
    for (int i = mid; i > l - 1; i--) {
        tree[level][i] = std::min(tree[level][i + 1], input[i]);
    }
    tree[level][mid + 1] = input[mid + 1];
    for (int i = mid + 2; i < r + 1; i++) {
        tree[level][i] = std::min(tree[level][i - 1], input[i]);
    }
    if (l < mid) {
        buildDST(tree, input, l, mid, level + 1);
        buildDST(tree, input, mid + 1, r, level + 1);
    }
}

int get(std::vector<std::vector<int>> & tree, int l, int r) {
    int level = 2 - std::__lg(l ^ r);
    return std::min(tree[level][l], tree[level][r]);
}

int main() {
    std::vector<int> a = {3, 8, 7, -2, 9, 5, 18, 11};
    std::vector<std::vector<int>> tree(3, std::vector<int>(8, INT32_MAX));
    buildDST(tree, a, 0, 7);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 8; j++) {
            std::cout << tree[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << get(tree, 0, 3);
}
