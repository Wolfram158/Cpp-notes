#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

std::vector<long long> h;
std::vector<long long> p;
long long number = 1000000007ll;

long long calc_hash(int n, int l, int r) {
	return (h[r + 1] - h[l]) * p[n - l];
}

void init_hash(std::string& string) {
	p.resize(string.size() + 1, 1);
	for (int i = 1; i < string.size() + 1; i++) {
		p[i] = p[i - 1] * number;
	}
	h.resize(string.size() + 1, 0);
	for (int i = 0; i < string.size(); i++) {
		long long phi = string[i];
		h[i + 1] = h[i] + p[i] * phi;
	}
}

void build_suffix_array(std::string& s, std::vector<int>& suffix_array) {
	std::sort(suffix_array.begin(), suffix_array.end(), [&](int a, int b) {
		int l = 0;
		int r = std::min(s.size() - a, s.size() - b);
		int max = r;
		while (l < r - 1) {
			int mid = (l + r) / 2;
			if (calc_hash(s.size(), a, a + mid) == calc_hash(s.size(), b, b + mid)) {
				l = mid;
			}
			else {
				r = mid;
			}
		}
		if (l == 0 && s[a] != s[b]) {
			return s[a] < s[b];
		}
		return s[a + l + 1] < s[b + l + 1];
	});
}

int main() {
	std::string s = "abbacaabbab$";
	init_hash(s);
	std::vector<int> sa(s.size());
	for (int i = 0; i < s.size(); i++) {
		sa[i] = i;
	}
	build_suffix_array(s, sa);
	for (int i = 0; i < s.size(); i++) {
		std::cout << sa[i] << "\n";
	}
}
