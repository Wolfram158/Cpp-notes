#include <iostream>
#include <vector>

const constexpr int modulo = 998244353;
int n, m;

void print_p_plus_q(std::vector<int>& p, std::vector<int>& q) {
	std::cout << std::max(m, n) << "\n";
	for (int i = 0; i < std::min(m, n) + 1; i++) {
		std::cout << (p[i] + q[i]) % modulo << " ";
	}
	if (n < m) {
		for (int i = n + 1; i < m + 1; i++) {
			std::cout << q[i] << " ";
		}
	}
	else {
		for (int i = m + 1; i < n + 1; i++) {
			std::cout << p[i] << " ";
		}
	}
	std::cout << "\n";
}

void print_p_mul_q(std::vector<int>& p, std::vector<int>& q) {
	std::vector<int> r(m + n + 1, 0);
	std::cout << m + n << "\n";
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < m + 1; j++) {
			r[i + j] = (r[i + j] + (p[i] * q[j]) % modulo) % modulo;
		}
	}
	for (int i = 0; i < m + n + 1; i++) {
		std::cout << r[i] << " ";
	}
	std::cout << "\n";
}

void print_p_div_q(std::vector<int>& p, std::vector<int>& q) {
	std::vector<int> r(1000, 0);
	r[0] = p[0] / q[0];
	std::cout << r[0] << " ";
	for (int i = 1; i < 1000; i++) {
		int sum = 0;
		for (int j = 1; j < i + 1; j++) {
			sum = (sum + (q[j] * r[i - j]) % modulo) % modulo;
		}
		int res = (p[i] - sum) / q[0];
		r[i] = res;
		if (res < 0) {
			res += modulo;
		}
		std::cout << res << " ";
	}
	std::cout << "\n";
}

int main() {
	std::cin >> n >> m;
	std::vector<int> p(1000, 0);
	std::vector<int> q(1000, 0);
	for (int i = 0; i < n + 1; i++) {
		std::cin >> p[i];
	}
	for (int i = 0; i < m + 1; i++) {
		std::cin >> q[i];
	}
	print_p_plus_q(p, q);
	print_p_mul_q(p, q);
	print_p_div_q(p, q);
}