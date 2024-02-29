/**
 * Author: Evan Bailey
 * Date: 2022-11-01
 * Description: A persistent treap with lazy propagation.
 * Use a bump allocator for better performance.
 * Time: $O(\log N)$
 * Status: tested
 */
#pragma once

ll rnd() {
	static uniform_int_distribution<ll> d(MIN<ll>, MAX<ll>);
	return d(randy);
}
struct Node; typedef array<Node *, 2> ann;
int T = 0; // current version. T++ => next op is persistent.
struct Node {
	ll p;      // priority
	ann c{};   // children
	int t;     // timestamp
	int n;     // size
	int v;     // value
	ll s;      // sum
	int l = 0; // lazy value
	Node(int val, ll p) : p(p), c{}, t(T), v(val), s(val) {}
	Node(int val) : Node(val, rnd()) {}
	//~Node() { rep(i, 0, 2) delete c[i]; } // breaks persistence
	void upd(int x) { l += x; v += x; s += (ll)n * x; } // custom
};
// If persistence not needed, just return n.
// Use bump allocator (w/ free list?) for efficiency.
void cp(Node *&n) {
	if (n && n->t != T) (n = new Node(*n))->t = T;
}
// If lazy prop not needed, just call cp(n) and return n.
Node *push(Node *&n) { // custom; recurse for ST beats
	cp(n); if (n && n->l) { // if (has lazy value)
		rep(i, 0, 2) if (cp(n->c[i]), n->c[i]) n->c[i]->upd(n->l);
		n->l = 0; // reset lazy value
	} return n;
}
void fix(Node *n) { // custom
	n->n = 1; n->s = n->v; // reset sum
	rep(i, 0, 2) if (n->c[i]) {
		n->n += n->c[i]->n;
		n->s += n->c[i]->s; // update n's sum
	}}
// p determines which side a node goes on (1 -> right)
template<class P> void split(Node *n, ann &A, P &&p) {
	if (n) {
		bool d = p(push(n)); split(n->c[!d], A, p);
		n->c[!d] = A[d]; fix(A[d] = n);
	}}
template<class P> auto split(Node *n, P &&p) {
	ann A{}; split(n, A, p); return make_pair(A[0], A[1]);
}
// 1st node where p is true or nullptr
// Pd: push-down every node
template<bool Pd, class P> Node *upper_bound(Node *n, P &&p) {
	if (n) {
		bool d = p(Pd ? push(n) : n);
		Node *m = upper_bound(n->c[!d], p);
		if (m || !d) n = m; // m || d for last where !p
	} return n;
}
// p for implicit split/upper_bound (splits < i and >= i)
struct impl {
	int i;
	impl(int i) : i(i) {}
	bool operator()(Node *n) {
		int l = n->c[0] ? n->c[0]->n : 0, d = l >= i;
		if (!d) i -= l + 1; return d;
	}};
Node *join(ann A) {
	rep(i, 0, 2) if (!A[i]) return A[!i];
	bool d = A[0]->p < A[1]->p;
	Node *n = push(A[!d]); A[!d] = n->c[d];
	n->c[d] = join(A); fix(n); return n;
}
Node *join(Node *a, Node *b) { return join({ a, b }); }
template<bool Pd, class F> void each(Node *n, const F &&f) {
	if (Pd ? push(n) : n) {
		each(n->c[0], f); f(n); each(n->c[1], f);
	}}
Node *tree(vector<Node *> A) { // A's sorted, singleton, nonlazy
	static stack<Node *> S({ 0 }); Node *r;
	for (Node *a : A) {
		cp(a); for (r = 0; S.top() && a->p < S.top()->p; S.pop())
			fix(r = S.top());
		a->c[0] = r; if (S.top()) S.top()->c[1] = a; S.push(a);
	} for (r = 0; S.top(); S.pop()) fix(r = S.top()); return r;
}