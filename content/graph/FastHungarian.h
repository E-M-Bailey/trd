/**
 * Author: Evan Bailey
 * Date: 2020-04-04
 * Description: Like WeightedMatching.h, but faster for dense graphs. costs should be nxn.
 * Time: O(N^3)
 * Status: Tested on kattis
 */

#pragma once

typedef double T;
const T EPS = 1e-8, INF = numeric_limits<T>::infinity();
struct Hungarian {
	vector<vector<T>> C;
	int n;
	vector<T> lw, lj, slj;
	vi slwj, mw, mj, p;
	vector<bool> cw;
	Hungarian(vector<vector<T>> &costs) :
		C(costs), n(sz(C)), lw(n), lj(n), slj(n), slwj(n),
		mw(n, -1), mj(n, -1), p(n), cw(n) {}
	void reduce() {
		rep(w, 0, n) {
			T m = INF;
			rep(j, 0, n) m = min(C[w][j], m);
			rep(j, 0, n) C[w][j] -= m;
		}
		rep(j, 0, n) {
			T m = INF;
			rep(w, 0, n) m = min(m, C[w][j]);
			rep(w, 0, n) C[w][j] -= m;
		}}
	void init() {
		rep(j, 0, n) lj[j] = INF;
		rep(w, 0, n) rep(j, 0, n) lj[j] = min(lj[j], C[w][j]);
	}
	void match() {
		rep(w, 0, n) rep(j, 0, n)
			if (mw[w] < 0 && mj[j] < 0 && abs(C[w][j] - lw[w] - lj[j]) <= EPS) {
				mw[w] = j; mj[j] = w;
			}}
	int next() {
		rep(w, 0, n) if (mw[w] < 0) return w;
		return -1;
	}
	void init_phase(int w) {
		fill(all(cw), false); fill(all(p), -1); cw[w] = true;
		rep(j, 0, n) { slj[j] = C[w][j] - lw[w] - lj[j]; slwj[j] = w; }}
	void upd(T s) {
		rep(w, 0, n) lw[w] += s * cw[w];
		rep(j, 0, n) (p[j] < 0 ? slj : lj)[j] -= s;
	}
	void phase() {
		for (;;) {
			int sw = -1, sj = -1; T s = INF;
			rep(j, 0, n)
				if (p[j] < 0 && slj[j] < s) { s = slj[j]; sw = slwj[j]; sj = j; }
			if (abs(s) > EPS) upd(s);
			p[sj] = sw;
			if (mj[sj] < 0) {
				int j = sj, w = p[sj];
				for (;;) {
					mj[j] = w; swap(mw[w], j); if (j < 0) break;
					w = p[j];
				} return;
			}
			int w = mj[sj]; cw[w] = 1;
			rep(j, 0, n)
				if (p[j] < 0) {
					T sl = C[w][j] - lw[w] - lj[j];
					if (slj[j] > sl) { slj[j] = sl; slwj[j] = w; }}}}
	vi operator()() {
		reduce(); init(); match();
		for (int w = next(); w >= 0; w = next()) {
			init_phase(w); phase();
		} return mw;
	}};
