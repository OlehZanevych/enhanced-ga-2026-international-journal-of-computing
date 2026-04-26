#pragma once
#include "types.h"
#include <cmath>
#include <vector>
#include <algorithm>

static constexpr double B1 = 150, B2 = 100, B3 = 50, B4 = 5, B5 = 20;
static constexpr double A1 = 2,   A2 = 2,   A3 = 2,  A4 = 2, A5 = 2;

inline double freq_weight(Freq f) { return f == Freq::Weekly ? 1.0 : 0.5; }

inline bool time_conflict(const Gene& a, const Gene& b) {
    if (a.day != b.day || a.slot != b.slot) return false;
    if (a.freq == Freq::Weekly || b.freq == Freq::Weekly) return true;
    return a.freq == b.freq;
}

// ── Index structure for O(1) bucket lookups ───────────────────────────────────
// Replaces O(V) scans inside LOF with direct bucket access.
struct ScheduleIndex {
    const Problem& p;
    int DS; // days * slots

    // lec_slot[l][day*slots+slot]  -> list of gene indices
    std::vector<std::vector<std::vector<int>>> lec_slot;
    // grp_slot[g][day*slots+slot]  -> list of gene indices
    std::vector<std::vector<std::vector<int>>> grp_slot;
    // room_slot[r][day*slots+slot] -> list of gene indices
    std::vector<std::vector<std::vector<int>>> room_slot;
    // lec_day[l][day]  -> sorted slot list
    std::vector<std::vector<std::vector<int>>> lec_day;
    // grp_day[g][day]  -> sorted slot list
    std::vector<std::vector<std::vector<int>>> grp_day;

    explicit ScheduleIndex(const Problem& prob) : p(prob), DS(prob.days * prob.slots) {
        lec_slot .assign(p.L, std::vector<std::vector<int>>(DS));
        grp_slot .assign(p.G, std::vector<std::vector<int>>(DS));
        room_slot.assign(p.A, std::vector<std::vector<int>>(DS));
        lec_day  .assign(p.L, std::vector<std::vector<int>>(p.days));
        grp_day  .assign(p.G, std::vector<std::vector<int>>(p.days));
    }

    void build(const Schedule& s) {
        // clear
        for (auto& v : lec_slot)  for (auto& b : v) b.clear();
        for (auto& v : grp_slot)  for (auto& b : v) b.clear();
        for (auto& v : room_slot) for (auto& b : v) b.clear();
        for (auto& v : lec_day)   for (auto& b : v) b.clear();
        for (auto& v : grp_day)   for (auto& b : v) b.clear();

        for (int i = 0; i < p.V; ++i) {
            const Gene& g = s[i];
            int ts = g.day * p.slots + g.slot;
            int l  = p.req_lecturer[g.req];
            lec_slot[l][ts].push_back(i);
            lec_day [l][g.day].push_back(g.slot);
            room_slot[g.room][ts].push_back(i);
            for (int gr : p.req_groups[g.req]) {
                grp_slot[gr][ts].push_back(i);
                grp_day [gr][g.day].push_back(g.slot);
            }
        }
        // sort day-slot lists for window computation
        for (auto& v : lec_day) for (auto& b : v) std::sort(b.begin(), b.end());
        for (auto& v : grp_day) for (auto& b : v) std::sort(b.begin(), b.end());
    }

    // Count conflicts in a bucket for gene i (same time slot)
    // Returns weighted overlap count excluding gene i itself
    double slot_conflicts(const std::vector<int>& bucket,
                          const Schedule& s, int i) const {
        double cnt = 0;
        for (int j : bucket) {
            if (j == i) continue;
            if (time_conflict(s[i], s[j]))
                cnt += freq_weight(s[j].freq);
        }
        return cnt;
    }

    static double windows(const std::vector<int>& sorted_slots) {
        double w = 0;
        for (size_t k = 1; k < sorted_slots.size(); ++k)
            w += sorted_slots[k] - sorted_slots[k-1] - 1;
        return w;
    }
};

// ── Global objective — O(V) with index ───────────────────────────────────────
double objective(const Schedule& s, const Problem& p, ScheduleIndex& idx) {
    idx.build(s);
    double pi1=0, pi2=0, pi3=0, pi4=0, pi5=0;

    // overlaps: for each bucket with >1 gene, count pairs
    auto overlap_score = [&](const std::vector<int>& bucket) -> double {
        double cnt = 0;
        for (int ii = 0; ii < (int)bucket.size(); ++ii)
            for (int jj = ii+1; jj < (int)bucket.size(); ++jj)
                if (time_conflict(s[bucket[ii]], s[bucket[jj]]))
                    cnt += freq_weight(s[bucket[ii]].freq)
                         * freq_weight(s[bucket[jj]].freq) * 2;
        return cnt;
    };

    for (int l = 0; l < p.L; ++l) {
        double cnt = 0;
        for (int ts = 0; ts < idx.DS; ++ts) cnt += overlap_score(idx.lec_slot[l][ts]);
        pi1 += std::pow(cnt, A1);
        for (int d = 0; d < p.days; ++d)
            pi4 += std::pow(ScheduleIndex::windows(idx.lec_day[l][d]), A4);
    }
    for (int g = 0; g < p.G; ++g) {
        double cnt = 0;
        for (int ts = 0; ts < idx.DS; ++ts) cnt += overlap_score(idx.grp_slot[g][ts]);
        pi2 += std::pow(cnt, A2);
        for (int d = 0; d < p.days; ++d)
            pi5 += std::pow(ScheduleIndex::windows(idx.grp_day[g][d]), A5);
    }
    for (int r = 0; r < p.A; ++r) {
        double cnt = 0;
        for (int ts = 0; ts < idx.DS; ++ts) cnt += overlap_score(idx.room_slot[r][ts]);
        pi3 += std::pow(cnt, A3);
    }

    return B1*pi1 + B2*pi2 + B3*pi3 + B4*pi4 + B5*pi5;
}

// ── Local objective — O(bucket_size) per gene ────────────────────────────────
double local_objective(const Schedule& s, int i,
                       const Problem& p, const ScheduleIndex& idx) {
    const Gene& gi = s[i];
    int l   = p.req_lecturer[gi.req];
    int ts  = gi.day * p.slots + gi.slot;

    double psi1 = idx.slot_conflicts(idx.lec_slot[l][ts], s, i);

    double psi2 = 0;
    for (int g : p.req_groups[gi.req])
        psi2 += idx.slot_conflicts(idx.grp_slot[g][ts], s, i);

    double psi3 = idx.slot_conflicts(idx.room_slot[gi.room][ts], s, i);

    double psi4 = ScheduleIndex::windows(idx.lec_day[l][gi.day]);

    double psi5 = 0;
    for (int g : p.req_groups[gi.req])
        psi5 += ScheduleIndex::windows(idx.grp_day[g][gi.day]);

    return B1*std::pow(psi1,A1) + B2*std::pow(psi2,A2) + B3*std::pow(psi3,A3)
         + B4*std::pow(psi4,A4) + B5*std::pow(psi5,A5);
}
