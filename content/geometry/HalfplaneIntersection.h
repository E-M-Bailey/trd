/**
 * Author: Mike Nguyen
 * Date: 2023-05-31
 * Source: MIT ICPC Notebook
 * Description: Compute intersection of a set of halfplanes. A point is in a halfplane if it lies on the left of the plane boundary. 
 * Status: tested on kattis:birthdaycake
*/
#include "../geometry/Point.h"

using P = Point<double>;
constexpr double EPS = 1e-10;
int lineIntersection(const P& s1, const P& e1, const P& s2, const P& e2, P& r) {
	if ((e1-s1).cross(e2-s2)) { //if not parallell
		r = s2-(e2-s2)*(e1-s1).cross(s2-s1)/(e1-s1).cross(e2-s2);
		return 1;
	}
    return -((e1-s1).cross(s2-s1)==0 || s2==e2);
}
struct HalfPlane {
    P p1, p2;
    explicit HalfPlane(P a = P(), P b = P()) : p1(a), p2(b) {}
    P intpo(HalfPlane y) {
        P r; 
        assert(lineIntersection(p1, p2, y.p1, y.p2, r) == 1);
        return r;
    }
    P dir() { return p2-p1; }
    bool contains(P x) { return (p2-p1).cross(x-p1) < EPS; }
    bool out(P x) { return !contains(x); }};
vector<P> halfPlaneIntersection(vector<HalfPlane> b) {
    sort(begin(b), end(b), [](HalfPlane al, HalfPlane bl) {
        auto a = al.dir(), bb = bl.dir();
        if (a.x * bb.x < 0) return a.x < 0;
        if (abs(a.x) < EPS) {
            if (abs(bb.x) < EPS) return a.y > 0 && bb.y < 0;
            if (bb.x < 0) return a.y > 0;
            if (bb.x > 0) return true;
        } 
        if (abs(bb.x) < EPS) {
            if (a.x < 0) return bb.y < 0;
            if (a.x > 0) return false;
        }
        return a.cross(bb) > 0;
    });
    int n = sz(b), q = 1, h = 0, i;
    vector<HalfPlane> c(sz(b) + 10);
    for (i=0; i<n; i++) {
        while (q < h && b[i].out(c[h].intpo(c[h-1]))) h--;
        while (q < h && b[i].out(c[q].intpo(c[q+1]))) q++;
        c[++h] = b[i];
        if (q < h && abs(c[h].dir().cross(c[h-1].dir())) < EPS) {
            if (c[h].dir().dot(c[h-1].dir()) > 0)
                if (b[i].out(c[--h].p1)) c[h] = b[i];
            else return {};
        }}
    while (q < h-1 && c[q].out(c[h].intpo(c[h-1]))) h--;
    while (q < h-1 && c[h].out(c[q].intpo(c[q+1]))) q++;
    if (h-q <= 1) return {};
    c[h+1] = c[q];
    vector<P> s;
    for (i = q; i<=h; i++) s.push_back(c[i].intpo(c[i+1]));
    s.push_back(s[0]);
    return s;
} 