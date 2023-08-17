#include <iostream>
#include <vector>
#include "DST.h"

int main() {
    std::vector<int> a = {3, 8, 7, -2, 9, 5, 18, 11};
    DST<int> dst = DST(a, INT32_MAX);
    std::cout << dst.queryDST(4, 5) << "\n";
    std::vector<double> b = {3.4, -4.9, 10.2, 4.5};
    DST<double> dst1 = DST(b, (double) INT32_MAX);
    std::cout << dst1.queryDST(0, 1);
}
