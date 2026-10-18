#pragma once
// fitness.h — objective function, its hard/soft decomposition, the per-gene
// local objective and the ScheduleIndex bucket structure.
//
// Objective (Eq. 4 of the article):
//   f(sigma) = sum_{i=1..5} beta_i * Pi_i(sigma)
//   Pi_1 = sum_{lecturers l} c_l^alpha,  Pi_2 = sum_{groups g} c_g^alpha,
//   Pi_3 = sum_{rooms a} c_a^alpha,      (c_e = weighted number of clashing
//                                          ordered pairs of classes of entity e)
//   Pi_4 = sum_{l,d} w_{l,d}^alpha,      Pi_5 = sum_{g,d} w_{g,d}^alpha,
//                                         (w_{e,d} = number of idle timeslots of
//                                          entity e on day d)
//   alpha = 2.
// Hard part  H(sigma) = beta1*Pi1 + beta2*Pi2 + beta3*Pi3  (H = 0 <=> feasible)
// Soft part  W(sigma) = beta4*Pi4 + beta5*Pi5
#include "types.h"
#include <cmath>
#include <vector>
#include <algorithm>

static constexpr double B1 = 150, B2 = 100, B3 = 50, B4 = 5, B5 = 20;
static constexpr double ALPHA = 2;

inline double powa(double x) { return ALPHA == 2 ? x * x : std::pow(x, ALPHA); }

inline double freq_weight(Freq f) { return f == Freq::Weekly ? 1.0 : 0.5; }

// Two classes held in the same (day, slot) clash unless one runs in numerator
// weeks and the other in denominator weeks.
inline bool freq_clash(Freq a, Freq b) {
    if (a == Freq::Weekly || b == Freq::Weekly) return true;
    return a == b;
}
inline bool time_conflict(const Gene& a, const Gene& b) {
    if (a.day != b.day || a.slot != b.slot) return false;
    return freq_clash(a.freq, b.freq);
}

// ── ScheduleIndex: bucket lists for O(b) conflict lookups ─────────────────────
// bucket[e][tau] = indices of genes of entity e placed in timeslot tau = d*S+s.
// Bucket occupancy b is small (typically 0..2), so per-gene conflict evaluation
// costs O(b * |G(r)|) instead of the O(V) of a pairwise scan.
struct ScheduleIndex {
    const Problem& p;
    int DS;
    std::vector<std::vector<std::vector<int>>> lec_slot;   // [L][DS]
    std::vector<std::vector<std::vector<int>>> grp_slot;   // [G][DS]
    std::vector<std::vector<std::vector<int>>> room_slot;  // [A][DS]

    explicit ScheduleIndex(const Problem& prob) : p(prob), DS(prob.days * prob.slots) {
        lec_slot .assign(p.L, std::vector<std::vector<int>>(DS));
        grp_slot .assign(p.G, std::vector<std::vector<int>>(DS));
        room_slot.assign(p.A, std::vector<std::vector<int>>(DS));
    }

    void build(const Schedule& s) {
        for (auto& v : lec_slot)  for (auto& b : v) b.clear();
        for (auto& v : grp_slot)  for (auto& b : v) b.clear();
        for (auto& v : room_slot) for (auto& b : v) b.clear();
        for (int i = 0; i < p.V; ++i) insert(i, s[i]);
    }

    void insert(int i, const Gene& g) {
        const int ts = g.day * p.slots + g.slot;
        lec_slot[p.req_lecturer[g.req]][ts].push_back(i);
        room_slot[g.room][ts].push_back(i);
        for (int gr : p.req_groups[g.req]) grp_slot[gr][ts].push_back(i);
    }
    void remove(int i, const Gene& g) {
        const int ts = g.day * p.slots + g.slot;
        auto rm = [i](std::vector<int>& v) { v.erase(std::find(v.begin(), v.end(), i)); };
        rm(lec_slot[p.req_lecturer[g.req]][ts]);
        rm(room_slot[g.room][ts]);
        for (int gr : p.req_groups[g.req]) rm(grp_slot[gr][ts]);
    }
    // Move gene i (currently stored in s[i]) to (d, sl, a) and update the buckets.
    void move(Schedule& s, int i, int d, int sl, int a) {
        remove(i, s[i]);
        s[i].day = d; s[i].slot = sl; s[i].room = a;
        insert(i, s[i]);
    }
};

// Number of idle timeslots ("windows") of one entity on day d.
//   eb       : the entity's bucket row (eb[tau] = genes in timeslot tau)
//   exclude  : a gene to ignore (-1 = none)
//   add_slot : a slot to treat as occupied in addition (-1 = none)
// Windows = empty slots strictly between the first and the last occupied slot
// of the day (a slot is occupied if at least one class is held there).
inline int day_windows(const std::vector<std::vector<int>>& eb, int d, int S,
                       int exclude = -1, int add_slot = -1) {
    int first = -1, last = -1, occupied = 0;
    for (int t = 0; t < S; ++t) {
        const auto& b = eb[d * S + t];
        int c = (int)b.size();
        if (exclude >= 0 && std::find(b.begin(), b.end(), exclude) != b.end()) --c;
        if (t == add_slot) ++c;
        if (c > 0) {
            if (first < 0) first = t;
            last = t; ++occupied;
        }
    }
    if (first < 0) return 0;
    return (last - first + 1) - occupied;
}

// ── Evaluation with decomposition ─────────────────────────────────────────────
struct Eval {
    double total = 0, hard = 0, soft = 0;
    double pi[5] = {0, 0, 0, 0, 0};
};

inline Eval evaluate(const Schedule& s, const Problem& p, ScheduleIndex& idx) {
    idx.build(s);
    Eval e;
    auto bucket_pairs = [&](const std::vector<int>& b) -> double {
        double cnt = 0;
        for (size_t x = 0; x < b.size(); ++x)
            for (size_t y = x + 1; y < b.size(); ++y)
                if (time_conflict(s[b[x]], s[b[y]]))
                    cnt += freq_weight(s[b[x]].freq) * freq_weight(s[b[y]].freq) * 2;
        return cnt;
    };
    for (int l = 0; l < p.L; ++l) {
        double c = 0;
        for (int ts = 0; ts < idx.DS; ++ts) c += bucket_pairs(idx.lec_slot[l][ts]);
        e.pi[0] += powa(c);
        for (int d = 0; d < p.days; ++d) e.pi[3] += powa(day_windows(idx.lec_slot[l], d, p.slots));
    }
    for (int g = 0; g < p.G; ++g) {
        double c = 0;
        for (int ts = 0; ts < idx.DS; ++ts) c += bucket_pairs(idx.grp_slot[g][ts]);
        e.pi[1] += powa(c);
        for (int d = 0; d < p.days; ++d) e.pi[4] += powa(day_windows(idx.grp_slot[g], d, p.slots));
    }
    for (int a = 0; a < p.A; ++a) {
        double c = 0;
        for (int ts = 0; ts < idx.DS; ++ts) c += bucket_pairs(idx.room_slot[a][ts]);
        e.pi[2] += powa(c);
    }
    e.hard  = B1 * e.pi[0] + B2 * e.pi[1] + B3 * e.pi[2];
    e.soft  = B4 * e.pi[3] + B5 * e.pi[4];
    e.total = e.hard + e.soft;
    return e;
}

inline double objective(const Schedule& s, const Problem& p, ScheduleIndex& idx) {
    return evaluate(s, p, idx).total;
}

// ── Local objective (LOF) of gene i — used to rank genes in crossover/mutation ─
inline double slot_conflicts(const std::vector<int>& bucket, const Schedule& s, int i) {
    double cnt = 0;
    for (int j : bucket) {
        if (j == i) continue;
        if (time_conflict(s[i], s[j])) cnt += freq_weight(s[j].freq);
    }
    return cnt;
}

inline double local_objective(const Schedule& s, int i, const Problem& p, const ScheduleIndex& idx) {
    const Gene& gi = s[i];
    const int l  = p.req_lecturer[gi.req];
    const int ts = gi.day * p.slots + gi.slot;
    const double psi1 = slot_conflicts(idx.lec_slot[l][ts], s, i);
    double psi2 = 0, psi5 = 0;
    for (int g : p.req_groups[gi.req]) {
        psi2 += slot_conflicts(idx.grp_slot[g][ts], s, i);
        psi5 += day_windows(idx.grp_slot[g], gi.day, p.slots);
    }
    const double psi3 = slot_conflicts(idx.room_slot[gi.room][ts], s, i);
    const double psi4 = day_windows(idx.lec_slot[l], gi.day, p.slots);
    return B1 * powa(psi1) + B2 * powa(psi2) + B3 * powa(psi3)
         + B4 * powa(psi4) + B5 * powa(psi5);
}
