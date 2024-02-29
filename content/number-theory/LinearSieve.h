/**
 * Author: Nisiyama_Suzune (Codeforces)
 * Date: 2017-08-26
 * Source: https://codeforces.com/blog/entry/54090
 * Description: Finds prime powers and computes a multiplicative function ($f(ab)=f(a)f(b)$ for coprime $a,b$) on $[1,n)$.
 * Customize with $f(p)$ and $f(ip)$ where $p$ is $i$'s smallest prime factor. Currently computes phi.
 * Time: O(n)
 */
#pragma once

vi P; // primes
bool C[MAXN]; // composite?
int F[MAXN], K[MAXN]; // function value, power of min pf

void sieve(int n) {
	fill_n(C, n, 0); F[1] = 1; P.clear();
	rep(i, 2, n) {
		if (!C[i]) {
			P.push_back(i);
			K[i] = 1;
			F[i] = i - 1; // f(i) for i prime (custom)
		}
		for (int j = 0, x; j < sz(P) && (x = i * P[j]) < n; ++j) {
			C[x] = K[x] = 1;
			if (i % P[j]) F[x] = F[i] * F[P[j]]; // f(i*p), p doesn't divide i
			else {
				F[x] = F[i] * P[j]; // f(ip), p = i's min pf (custom)
				K[x] += K[i];
				break;
			}}}}
