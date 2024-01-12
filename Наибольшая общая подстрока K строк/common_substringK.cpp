#include <vector>
#include <iostream>
#include <set>

int k = 0;
size_t max_len = 0;
size_t min_len = 100000;
std::vector<std::vector<long long>> h;
std::vector<long long> p;

long long calc_hash(int n, int l, int r, int i) {
	return (h[i][r+1] - h[i][l]) * p[n - l];
}

void init_hash(std::vector<std::string> & strings) {
	p.resize(max_len + 1);
	p[0] = 1;
	for (int i = 1; i <= max_len; i++) {
		p[i] = (p[i - 1] * 1000000007ll);
	}
	h.resize(k + 1, std::vector<long long>(max_len + 1));
	for (int j = 0; j < k; j++) {
		h[j][0] = 0;
		for (int i = 0; i < strings[j].size(); i++) {
			long long phi = strings[j][i];
			h[j][i + 1] = (h[j][i] + p[i] * phi);
		}
	}
}

int main() {
	std::cin >> k;
	std::vector<std::string> strings;
	for (int i = 0; i < k; i++) {
		std::string s;
		std::cin >> s;
		max_len = std::max(max_len, s.size());
		min_len = std::min(min_len, s.size());
		strings.push_back(s);
	}
	init_hash(strings);
	int l = 0;
	int r = min_len;
	int curl = 0, curr = 0;
	std::vector<std::set<long long>> hs(k);
	while (l < r - 1) {
		int mid = (l + r) / 2;
		for (int i = 0; i < k; i++) {
			hs[i].clear();
		}
		for (int j = mid; j < strings[0].size(); j++) {
			hs[0].insert(calc_hash(max_len, j - mid, j, 0));
		}
		for (int i = 1; i < k; i++) {
			for (int j = mid; j < strings[i].size(); j++) {
				long long t = calc_hash(max_len, j - mid, j, i);
				if (hs[i - 1].find(t) != hs[i - 1].end()) {
					if (i == k - 1) {
						curl = j - mid;
						curr = j;
					}
					hs[i].insert(t);
				}
			}
		}
		if (hs[k - 1].empty()) {
			r = mid;
		}
		else {
			l = mid;
		}
	}
	std::cout << strings[k - 1].substr(curl, curr - curl + 1) << "\n";
}
