#pragma once
// ga.h — enhanced genetic algorithm with optional conflict-driven local search.
//
// Local-search modes (GAParams::ls_mode):
//   None            "GA"       : no local search.
//   Overlap         "GA+LS"    : overlap-only conflict-driven repair (basic variant).
//   WindowTiebreak  "GA+LS-W"  : ablation — the same repair, but window count is
//                                used as a tiebreaker among conflict-free placements.
//   Adaptive        "GA+AMLS"  : proposed — adaptive multi-neighbourhood LS:
//                                overlap-only relocation + ejection swaps for the
//                                hard constraints, followed by a separate,
//                                feasibility-preserving window-compaction
//                                neighbourhood with exact delta evaluation; the
//                                number of genes treated grows with stagnation.
#include "types.h"
#include "fitness.h"
#include "rng.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <numeric>
#include <vector>

enum class LSMode : int { None = 0, Overlap = 1, WindowTiebreak = 2, Adaptive = 3 };

inline const char* ls_mode_name(LSMode m) {
    switch (m) {
        case LSMode::None:           return "no_ls";
        case LSMode::Overlap:        return "ls";
        case LSMode::WindowTiebreak: return "ls_w";
        case LSMode::Adaptive:       return "amls";
    }
    return "?";
}

struct GAParams {
    int    S            = 20;     // population size (N in the article)
    double C            = 0.5;    // crossover offspring per generation = C*S
    double CR_min       = 0.05;   // crossover replaces r in [CR_min*V, CR_max*V] genes
    double CR_max       = 0.20;
    double M            = 0.10;   // mutation acts on the M*V worst genes
    double MR_day       = 0.25;   // per-operator mutation probabilities
    double MR_slot      = 0.25;
    double MR_room      = 0.25;
    double MR_freq      = 0.25;
    int    elites       = 2;
    int    tournament   = 3;
    int    stagnation   = 30;     // T_stag
    double restart_frac = 0.5;
    int    max_iter     = 5000;   // T_max
    double time_limit_sec = 0;    // > 0: additional wall-clock budget (equal-time study)
    LSMode ls_mode      = LSMode::None;
    double ls_fraction  = 0.10;   // LS candidate genes per mutant: floor(ls_fraction*V)
    int    swap_candidates = 8;   // AMLS: ejection-swap candidates examined per gene
};

struct LSStats {
    long long calls = 0, relocations = 0, swaps = 0, soft_moves = 0;
    double time_sec = 0;
    // Sampled every 10th LS call: number of conflicted genes (genes involved in
    // at least one hard-constraint clash) in the mutant before and after LS.
    long long sampled = 0, conflicted_before = 0, conflicted_after = 0;
};

// ── helpers ───────────────────────────────────────────────────────────────────
inline Gene random_gene(int req, const Problem& p, Rng& rng) {
    Gene g;
    g.req  = req;
    g.day  = rand_int(rng, 0, p.days  - 1);
    g.slot = rand_int(rng, 0, p.slots - 1);
    g.room = rand_int(rng, 0, p.A     - 1);
    g.freq = p.req_freq[req];
    return g;
}

inline Schedule random_schedule(const Problem& p, Rng& rng) {
    Schedule s(p.V);
    for (int i = 0; i < p.V; ++i) s[i] = random_gene(i, p, rng);
    return s;
}

// Top-k gene indices by local objective (descending; ties -> lower index first).
inline std::vector<int> worst_genes(const Schedule& s, int k, const Problem& p,
                                    const ScheduleIndex& idx) {
    std::vector<double> lof(s.size());
    for (int i = 0; i < (int)s.size(); ++i) lof[i] = local_objective(s, i, p, idx);
    std::vector<int> order(s.size());
    std::iota(order.begin(), order.end(), 0);
    std::partial_sort(order.begin(), order.begin() + k, order.end(), [&](int a, int b) {
        return lof[a] > lof[b] || (lof[a] == lof[b] && a < b);
    });
    order.resize(k);
    return order;
}

// ── conflict cost of a gene at a (possibly hypothetical) placement ────────────
// Weighted, frequency-aware clash degree of gene i against the genes already in
// a bucket (gene i itself is ignored), consistent with Pi_1..Pi_3 of f.
inline double bucket_clash(const std::vector<int>& b, const Schedule& s, int i) {
    const Freq fi = s[i].freq;
    double c = 0;
    for (int j : b)
        if (j != i && freq_clash(fi, s[j].freq)) c += freq_weight(fi) * freq_weight(s[j].freq);
    return c;
}
// Lecturer + group part at timeslot ts.
inline double lg_cost(const Schedule& s, int i, int ts, const Problem& p, const ScheduleIndex& idx) {
    const int l = p.req_lecturer[s[i].req];
    double c = B1 * bucket_clash(idx.lec_slot[l][ts], s, i);
    for (int g : p.req_groups[s[i].req]) c += B2 * bucket_clash(idx.grp_slot[g][ts], s, i);
    return c;
}
inline double room_cost(const Schedule& s, int i, int ts, int a, const ScheduleIndex& idx) {
    return B3 * bucket_clash(idx.room_slot[a][ts], s, i);
}
inline double gene_cost(const Schedule& s, int i, const Problem& p, const ScheduleIndex& idx) {
    const int ts = s[i].day * p.slots + s[i].slot;
    return lg_cost(s, i, ts, p, idx) + room_cost(s, i, ts, s[i].room, idx);
}

inline int count_conflicted(const Schedule& s, const Problem& p, const ScheduleIndex& idx) {
    int n = 0;
    for (int i = 0; i < p.V; ++i) n += gene_cost(s, i, p, idx) > 0;
    return n;
}

// Top-k conflicted genes (cost > 0), sorted by cost descending, ties -> lower index.
inline std::vector<int> conflicted_genes(const Schedule& s, int k, const Problem& p,
                                         const ScheduleIndex& idx) {
    std::vector<std::pair<double, int>> sc;
    sc.reserve(p.V);
    for (int i = 0; i < p.V; ++i) {
        const double c = gene_cost(s, i, p, idx);
        if (c > 0) sc.push_back({c, i});
    }
    const int n = std::min(k, (int)sc.size());
    std::partial_sort(sc.begin(), sc.begin() + n, sc.end(), [](const auto& a, const auto& b) {
        return a.first > b.first || (a.first == b.first && a.second < b.second);
    });
    std::vector<int> out(n);
    for (int t = 0; t < n; ++t) out[t] = sc[t].second;
    return out;
}

// ── LS mode 1: overlap-only conflict-driven repair (basic variant) ───────────
// For each of the k most conflicted genes, scan all (d, s, a) placements in
// lexicographic order and move the gene to the first conflict-free placement;
// if none exists, to the placement of minimum conflict cost. Window penalties
// are ignored. Buckets are updated incrementally after every move.
inline void ls_overlap(Schedule& s, const Problem& p, ScheduleIndex& idx, int k, LSStats& st) {
    for (int i : conflicted_genes(s, k, p, idx)) {
        const double cur = gene_cost(s, i, p, idx);
        if (cur == 0) continue;                           // fixed by an earlier move
        double best = cur;
        int bd = s[i].day, bsl = s[i].slot, br = s[i].room;
        for (int d = 0; d < p.days && best > 0; ++d)
            for (int sl = 0; sl < p.slots && best > 0; ++sl) {
                const int ts = d * p.slots + sl;
                const double lg = lg_cost(s, i, ts, p, idx);
                if (lg >= best) continue;
                for (int a = 0; a < p.A; ++a) {
                    const double tot = lg + room_cost(s, i, ts, a, idx);
                    if (tot < best) { best = tot; bd = d; bsl = sl; br = a; if (best == 0) break; }
                }
            }
        if (bd == s[i].day && bsl == s[i].slot && br == s[i].room) continue;
        idx.move(s, i, bd, bsl, br);
        ++st.relocations;
    }
}

// Windows of gene i's lecturer and groups on day d if i were placed in slot sl.
inline int local_windows_if(const Schedule& s, int i, int d, int sl, const Problem& p,
                            const ScheduleIndex& idx) {
    const int l = p.req_lecturer[s[i].req];
    int w = day_windows(idx.lec_slot[l], d, p.slots, i, sl);
    for (int g : p.req_groups[s[i].req]) w += day_windows(idx.grp_slot[g], d, p.slots, i, sl);
    return w;
}

// ── LS mode 2 (ablation): repair with window tiebreaker ───────────────────────
// Same candidates and conflict cost as mode 1, but the whole neighbourhood is
// scanned and, among conflict-free placements, the one minimising the windows
// of the repaired gene's own lecturer/groups on the target day is chosen. The
// windows this move may create on the day the gene leaves are not considered.
inline void ls_window_tiebreak(Schedule& s, const Problem& p, ScheduleIndex& idx, int k, LSStats& st) {
    for (int i : conflicted_genes(s, k, p, idx)) {
        const double cur = gene_cost(s, i, p, idx);
        if (cur == 0) continue;
        double best = cur;
        int bestW = 1 << 30;
        int bd = s[i].day, bsl = s[i].slot, br = s[i].room;
        for (int d = 0; d < p.days; ++d)
            for (int sl = 0; sl < p.slots; ++sl) {
                const int ts = d * p.slots + sl;
                const double lg = lg_cost(s, i, ts, p, idx);
                if (lg > best) continue;
                int w = -1;                                // lazily computed
                for (int a = 0; a < p.A; ++a) {
                    const double tot = lg + room_cost(s, i, ts, a, idx);
                    if (tot < best) {
                        best = tot; bd = d; bsl = sl; br = a;
                        if (tot == 0) { if (w < 0) w = local_windows_if(s, i, d, sl, p, idx); bestW = w; }
                    } else if (tot == 0 && best == 0) {
                        if (w < 0) w = local_windows_if(s, i, d, sl, p, idx);
                        if (w < bestW) { bestW = w; bd = d; bsl = sl; br = a; }
                    }
                }
            }
        if (bd == s[i].day && bsl == s[i].slot && br == s[i].room) continue;
        idx.move(s, i, bd, bsl, br);
        ++st.relocations;
    }
}

// Exact change of the soft penalty (beta4*Pi4 + beta5*Pi5) if gene i moves
// from (d0, s0) to (d, sl). Only the lecturer/group day-rows of the two days
// involved can change.
inline double soft_delta(const Schedule& s, int i, int d, int sl, const Problem& p,
                         const ScheduleIndex& idx) {
    const int d0 = s[i].day, S = p.slots;
    auto one = [&](const std::vector<std::vector<int>>& eb, double beta) {
        if (d == d0)
            return beta * (powa(day_windows(eb, d0, S, i, sl)) - powa(day_windows(eb, d0, S)));
        return beta * (powa(day_windows(eb, d0, S, i, -1)) + powa(day_windows(eb, d, S, -1, sl))
                     - powa(day_windows(eb, d0, S)) - powa(day_windows(eb, d, S)));
    };
    double delta = one(idx.lec_slot[p.req_lecturer[s[i].req]], B4);
    for (int g : p.req_groups[s[i].req]) delta += one(idx.grp_slot[g], B5);
    return delta;
}

// Genes involved in a clash with gene i at (ts, a). Returns the single blocking
// gene if exactly one other gene causes all clashes, otherwise -1.
inline int single_blocker(const Schedule& s, int i, int ts, int a, const Problem& p,
                          const ScheduleIndex& idx) {
    int blocker = -1;
    auto scan = [&](const std::vector<int>& b) {
        for (int j : b) {
            if (j == i || !freq_clash(s[i].freq, s[j].freq)) continue;
            if (blocker == -1) blocker = j;
            else if (blocker != j) return false;
        }
        return true;
    };
    if (!scan(idx.lec_slot[p.req_lecturer[s[i].req]][ts])) return -1;
    for (int g : p.req_groups[s[i].req]) if (!scan(idx.grp_slot[g][ts])) return -1;
    if (!scan(idx.room_slot[a][ts])) return -1;
    return blocker;
}

// ── LS mode 3 (proposed): adaptive multi-neighbourhood local search ──────────
// N1 relocation (as mode 1) and N2 ejection swap for genes that cannot be
// relocated without clashes, followed by N3 window compaction: a conflict-free
// gene is moved to another conflict-free placement only if the EXACT change of
// the soft penalty is negative, so the hard part is unchanged and f strictly
// decreases (no window inflation by construction). k grows with stagnation.
inline void ls_adaptive(Schedule& s, const Problem& p, ScheduleIndex& idx, int k,
                        int swap_candidates, LSStats& st) {
    // Phase A — hard constraints: N1 relocation, N2 ejection swap.
    struct Cand { int j, d, sl, a; };
    std::vector<Cand> cands;
    for (int i : conflicted_genes(s, k, p, idx)) {
        const double cur = gene_cost(s, i, p, idx);
        if (cur == 0) continue;
        double best = cur;
        int bd = s[i].day, bsl = s[i].slot, br = s[i].room;
        cands.clear();
        for (int d = 0; d < p.days && best > 0; ++d)
            for (int sl = 0; sl < p.slots && best > 0; ++sl) {
                const int ts = d * p.slots + sl;
                const double lg = lg_cost(s, i, ts, p, idx);
                for (int a = 0; a < p.A; ++a) {
                    if (d == s[i].day && sl == s[i].slot && a == s[i].room) continue;
                    const double tot = lg + room_cost(s, i, ts, a, idx);
                    if (tot < best) { best = tot; bd = d; bsl = sl; br = a; if (best == 0) break; }
                    if (tot > 0 && (int)cands.size() < swap_candidates) {
                        const int j = single_blocker(s, i, ts, a, p, idx);
                        if (j >= 0) cands.push_back({j, d, sl, a});
                    }
                }
            }
        if (best == 0) { idx.move(s, i, bd, bsl, br); ++st.relocations; continue; }
        // N2: i takes the blocked placement, the blocker j takes i's old placement.
        bool swapped = false;
        const int d0 = s[i].day, s0 = s[i].slot, a0 = s[i].room;
        for (const Cand& c : cands) {
            const int j = c.j;
            const int dj = s[j].day, sj = s[j].slot, aj = s[j].room;
            const double cj = gene_cost(s, j, p, idx);
            idx.move(s, i, c.d, c.sl, c.a);
            idx.move(s, j, d0, s0, a0);
            const double after = gene_cost(s, i, p, idx) + gene_cost(s, j, p, idx);
            // accept if the pair ends up better than with the best relocation of i
            if (after < best + cj) { swapped = true; ++st.swaps; break; }
            idx.move(s, j, dj, sj, aj);                    // revert
            idx.move(s, i, d0, s0, a0);
        }
        if (!swapped && best < cur) { idx.move(s, i, bd, bsl, br); ++st.relocations; }
    }

    // Phase B — soft constraints: N3 feasibility-preserving window compaction.
    std::vector<std::pair<double, int>> wsc;
    for (int i = 0; i < p.V; ++i) {
        if (gene_cost(s, i, p, idx) > 0) continue;
        const int l = p.req_lecturer[s[i].req];
        double w = B4 * day_windows(idx.lec_slot[l], s[i].day, p.slots);
        for (int g : p.req_groups[s[i].req]) w += B5 * day_windows(idx.grp_slot[g], s[i].day, p.slots);
        if (w > 0) wsc.push_back({w, i});
    }
    const int n = std::min(k, (int)wsc.size());
    std::partial_sort(wsc.begin(), wsc.begin() + n, wsc.end(), [](const auto& a, const auto& b) {
        return a.first > b.first || (a.first == b.first && a.second < b.second);
    });
    for (int t = 0; t < n; ++t) {
        const int i = wsc[t].second;
        if (gene_cost(s, i, p, idx) > 0) continue;
        double best = -1e-9;
        int bd = -1, bsl = -1, br = -1;
        for (int d = 0; d < p.days; ++d)
            for (int sl = 0; sl < p.slots; ++sl) {
                if (d == s[i].day && sl == s[i].slot) continue;
                const int ts = d * p.slots + sl;
                if (lg_cost(s, i, ts, p, idx) > 0) continue;
                int a = -1;
                if (room_cost(s, i, ts, s[i].room, idx) == 0) a = s[i].room;
                else for (int r = 0; r < p.A; ++r) if (room_cost(s, i, ts, r, idx) == 0) { a = r; break; }
                if (a < 0) continue;
                const double delta = soft_delta(s, i, d, sl, p, idx);
                if (delta < best) { best = delta; bd = d; bsl = sl; br = a; }
            }
        if (bd < 0) continue;
        idx.move(s, i, bd, bsl, br);
        ++st.soft_moves;
    }
}

// ── 3-parent conflict-aware crossover ─────────────────────────────────────────
// The child starts as a copy of P1; each of the r genes of P1 with the highest
// local objective is replaced by the version (from P1, P2 or P3) that has the
// lowest local objective in its own parent (Eq. 7 of the article).
inline Schedule crossover(const Schedule& s1, const Schedule& s2, const Schedule& s3,
                          int r, const Problem& p,
                          ScheduleIndex& idx1, ScheduleIndex& idx2, ScheduleIndex& idx3) {
    Schedule child = s1;
    for (int i : worst_genes(s1, r, p, idx1)) {
        const double l2 = local_objective(s2, i, p, idx2);
        const double l3 = local_objective(s3, i, p, idx3);
        const double l1 = local_objective(s1, i, p, idx1);
        if (l2 < l1 && l2 <= l3)      child[i] = s2[i];
        else if (l3 < l1 && l3 < l2)  child[i] = s3[i];
    }
    return child;
}

// ── Adaptive mutation ─────────────────────────────────────────────────────────
// Acts on the M*V genes with the highest local objective; each of the four
// operators fires independently with probability MR * mr_scale, where
// mr_scale = 1 + t_stag / T_stag grows linearly while the search stagnates.
inline Schedule mutate(const Schedule& s, const Problem& p, Rng& rng,
                       const GAParams& prm, double mr_scale, ScheduleIndex& idx) {
    Schedule child = s;
    const int k = std::max(1, (int)(prm.M * p.V));
    auto flip = [&](double base) { return rand_real(rng) < base * mr_scale; };
    for (int i : worst_genes(s, k, p, idx)) {
        Gene& g = child[i];
        if (flip(prm.MR_day))  g.day  = rand_int(rng, 0, p.days  - 1);
        if (flip(prm.MR_slot)) g.slot = rand_int(rng, 0, p.slots - 1);
        if (flip(prm.MR_room)) g.room = rand_int(rng, 0, p.A     - 1);
        if (flip(prm.MR_freq) && p.req_freq[g.req] != Freq::Weekly)
            g.freq = (g.freq == Freq::Numerator) ? Freq::Denominator : Freq::Numerator;
    }
    return child;
}

// ── Survivor selection: elitism + tournament over the (mu + lambda) pool ──────
inline std::vector<int> select_indices(const std::vector<double>& scores, int S, int elites,
                                       int tournament, Rng& rng) {
    std::vector<int> order(scores.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return scores[a] < scores[b]; });
    std::vector<int> next;
    next.reserve(S);
    for (int i = 0; i < elites && i < (int)scores.size(); ++i) next.push_back(order[i]);
    const int n = (int)scores.size();
    while ((int)next.size() < S) {
        int best = rand_int(rng, 0, n - 1);
        for (int t = 1; t < tournament; ++t) {
            const int c = rand_int(rng, 0, n - 1);
            if (scores[c] < scores[best]) best = c;
        }
        next.push_back(best);
    }
    return next;
}

// ── Main GA loop ──────────────────────────────────────────────────────────────
struct IterRecord {
    int    iteration;
    double best_score, best_hard, best_soft;   // best-so-far and its decomposition
    double iter_score;                         // best of this generation's pool
    double iter_time_sec, cumul_time_sec;
    bool   restart;
};

struct GAResult {
    Schedule best;
    Eval     best_eval;
    int      iterations = 0;
    double   elapsed_sec = 0;
    int      first_feasible_iter = -1;         // first iteration with H = 0 (0 = initial pop)
    double   first_feasible_sec  = -1;
    std::vector<IterRecord> records;
    LSStats  ls;
};

inline GAResult run_ga(const Problem& p, const GAParams& prm, Rng& rng, bool verbose = false) {
    using clk = std::chrono::steady_clock;
    ScheduleIndex idx(p), idx2(p), idx3(p);
    GAResult R;

    std::vector<Schedule> pop(prm.S);
    for (auto& s : pop) s = random_schedule(p, rng);
    std::vector<Eval> pev(prm.S);
    for (int i = 0; i < prm.S; ++i) pev[i] = evaluate(pop[i], p, idx);

    int bi = 0;
    for (int i = 1; i < prm.S; ++i) if (pev[i].total < pev[bi].total) bi = i;
    R.best = pop[bi]; R.best_eval = pev[bi];

    const auto t_start = clk::now();
    auto elapsed = [&] { return std::chrono::duration<double>(clk::now() - t_start).count(); };
    if (R.best_eval.hard == 0) { R.first_feasible_iter = 0; R.first_feasible_sec = 0; }

    int stag = 0;
    const int n_cross = std::max(1, (int)(prm.C * prm.S));

    for (int iter = 1; iter <= prm.max_iter; ++iter) {
        if (R.best_eval.total == 0.0) break;
        if (prm.time_limit_sec > 0 && elapsed() >= prm.time_limit_sec) break;
        const auto t_iter = clk::now();
        const double mr_scale = 1.0 + (double)stag / prm.stagnation;

        std::vector<Schedule> offspring;
        offspring.reserve(n_cross + prm.S);

        // Crossover offspring (parents drawn uniformly at random).
        for (int c = 0; c < n_cross; ++c) {
            const int i1 = rand_int(rng, 0, prm.S - 1);
            int i2, i3;
            do { i2 = rand_int(rng, 0, prm.S - 1); } while (i2 == i1);
            do { i3 = rand_int(rng, 0, prm.S - 1); } while (i3 == i1 || i3 == i2);
            int r = (int)(prm.CR_min * p.V) + (int)(rand_real(rng) * (prm.CR_max - prm.CR_min) * p.V);
            r = std::max(1, std::min(r, p.V - 1));
            idx.build(pop[i1]); idx2.build(pop[i2]); idx3.build(pop[i3]);
            offspring.push_back(crossover(pop[i1], pop[i2], pop[i3], r, p, idx, idx2, idx3));
        }

        // One mutant per population member, optionally repaired by local search.
        for (const auto& s : pop) {
            idx.build(s);
            Schedule m = mutate(s, p, rng, prm, mr_scale, idx);
            if (prm.ls_mode != LSMode::None) {
                const auto t_ls = clk::now();
                idx.build(m);
                const bool sample = (R.ls.calls % 10 == 0);
                if (sample) R.ls.conflicted_before += count_conflicted(m, p, idx);
                const int k = std::max(1, (int)(prm.ls_fraction * p.V + 1e-9));
                switch (prm.ls_mode) {
                    case LSMode::Overlap:
                        ls_overlap(m, p, idx, k, R.ls);
                        break;
                    case LSMode::WindowTiebreak:
                        ls_window_tiebreak(m, p, idx, k, R.ls);
                        break;
                    case LSMode::Adaptive: {
                        const int ka = std::max(1, (int)(prm.ls_fraction * p.V * mr_scale + 1e-9));
                        ls_adaptive(m, p, idx, ka, prm.swap_candidates, R.ls);
                        break;
                    }
                    default: break;
                }
                if (sample) { R.ls.conflicted_after += count_conflicted(m, p, idx); ++R.ls.sampled; }
                ++R.ls.calls;
                R.ls.time_sec += std::chrono::duration<double>(clk::now() - t_ls).count();
            }
            offspring.push_back(std::move(m));
        }

        // (mu + lambda) pool: parents keep their cached evaluations.
        std::vector<Schedule> pool;
        std::vector<Eval> pool_ev;
        pool.reserve(prm.S + offspring.size());
        pool_ev.reserve(prm.S + offspring.size());
        for (int i = 0; i < prm.S; ++i) { pool.push_back(std::move(pop[i])); pool_ev.push_back(pev[i]); }
        for (auto& s : offspring) { pool_ev.push_back(evaluate(s, p, idx)); pool.push_back(std::move(s)); }

        double iter_best = pool_ev[0].total;
        for (const auto& e : pool_ev) iter_best = std::min(iter_best, e.total);
        for (int i = 0; i < (int)pool.size(); ++i)
            if (pool_ev[i].total < R.best_eval.total) {
                R.best_eval = pool_ev[i]; R.best = pool[i]; stag = 0;
            }

        std::vector<double> totals(pool_ev.size());
        for (size_t i = 0; i < pool_ev.size(); ++i) totals[i] = pool_ev[i].total;
        const std::vector<int> sel = select_indices(totals, prm.S, prm.elites, prm.tournament, rng);
        for (int i = 0; i < prm.S; ++i) { pop[i] = pool[sel[i]]; pev[i] = pool_ev[sel[i]]; }
        ++stag;

        bool did_restart = false;
        if (stag >= prm.stagnation) {
            const int n_restart = (int)(prm.restart_frac * prm.S);
            std::vector<int> ord(prm.S);
            std::iota(ord.begin(), ord.end(), 0);
            std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return pev[a].total > pev[b].total; });
            for (int t = 0; t < n_restart; ++t) {
                pop[ord[t]] = random_schedule(p, rng);
                pev[ord[t]] = evaluate(pop[ord[t]], p, idx);
            }
            stag = 0;
            did_restart = true;
        }

        const double iter_sec  = std::chrono::duration<double>(clk::now() - t_iter).count();
        const double cumul_sec = elapsed();
        if (R.first_feasible_iter < 0 && R.best_eval.hard == 0) {
            R.first_feasible_iter = iter; R.first_feasible_sec = cumul_sec;
        }
        R.records.push_back({iter, R.best_eval.total, R.best_eval.hard, R.best_eval.soft,
                             iter_best, iter_sec, cumul_sec, did_restart});
        if (verbose) {
            printf("  iter %5d | best=%10.4f (H=%.2f W=%.2f) | iter_best=%10.4f | total=%.2fs%s\n",
                   iter, R.best_eval.total, R.best_eval.hard, R.best_eval.soft, iter_best, cumul_sec,
                   did_restart ? " [RESTART]" : "");
            fflush(stdout);
        }
    }
    R.iterations  = (int)R.records.size();
    R.elapsed_sec = elapsed();
    return R;
}
