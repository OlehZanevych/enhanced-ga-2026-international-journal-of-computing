#pragma once
#include "types.h"
#include "fitness.h"
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include <chrono>

struct GAParams {
    int    S            = 20;
    double C            = 0.5;
    double CR_min       = 0.05;
    double CR_max       = 0.20;
    double M            = 0.10;
    double MR_day       = 0.25;
    double MR_slot      = 0.25;
    double MR_room      = 0.25;
    double MR_freq      = 0.25;
    int    elites       = 2;
    int    tournament   = 3;
    int    stagnation   = 30;
    double restart_frac = 0.5;
    int    max_iter     = 500;
    bool   use_local_search = false; // toggle improvement procedure
    double ls_fraction = 0.75;       // fraction of worst mutants to apply LS to (0=none, 1=all)
};

// ── helpers ───────────────────────────────────────────────────────────────────

inline Gene random_gene(int req, const Problem& p, std::mt19937& rng) {
    Gene g;
    g.req  = req;
    g.day  = std::uniform_int_distribution<>(0, p.days  - 1)(rng);
    g.slot = std::uniform_int_distribution<>(0, p.slots - 1)(rng);
    g.room = std::uniform_int_distribution<>(0, p.A     - 1)(rng);
    g.freq = p.req_freq[req];
    return g;
}

inline Schedule random_schedule(const Problem& p, std::mt19937& rng) {
    Schedule s(p.V);
    for (int i = 0; i < p.V; ++i) s[i] = random_gene(i, p, rng);
    return s;
}

// Top-k gene indices by LOF descending — index must be pre-built for s
inline std::vector<int> worst_genes(const Schedule& s, int k,
                                    const Problem& p, const ScheduleIndex& idx) {
    std::vector<double> lof(s.size());
    for (int i = 0; i < (int)s.size(); ++i)
        lof[i] = local_objective(s, i, p, idx);
    std::vector<int> order(s.size());
    std::iota(order.begin(), order.end(), 0);
    std::partial_sort(order.begin(), order.begin() + k, order.end(),
                      [&](int a, int b){ return lof[a] > lof[b]; });
    order.resize(k);
    return order;
}

// ── Local search: conflict-driven repair with window tiebreaker ───────────────
// Single pass over the top-k conflicted genes (sorted by cheap overlap count).
// For each gene, finds the (day,slot,room) minimising overlaps, with window
// count as tiebreaker among zero-overlap placements.
// Index is updated incrementally after each move — no full rebuilds.

static void idx_remove(ScheduleIndex& idx, int i, const Gene& g,
                        int l, const std::vector<int>& groups) {
    int ts = g.day * idx.p.slots + g.slot;
    auto rm = [](std::vector<int>& v, int x){
        v.erase(std::find(v.begin(), v.end(), x));
    };
    rm(idx.lec_slot[l][ts], i);
    rm(idx.room_slot[g.room][ts], i);
    for (int gr : groups) rm(idx.grp_slot[gr][ts], i);
    rm(idx.lec_day[l][g.day], g.slot);
    for (int gr : groups) rm(idx.grp_day[gr][g.day], g.slot);
}

static void idx_insert(ScheduleIndex& idx, int i, const Gene& g,
                        int l, const std::vector<int>& groups) {
    int ts = g.day * idx.p.slots + g.slot;
    auto ins = [](std::vector<int>& v, int x){
        v.insert(std::upper_bound(v.begin(), v.end(), x), x);
    };
    idx.lec_slot[l][ts].push_back(i);
    idx.room_slot[g.room][ts].push_back(i);
    for (int gr : groups) idx.grp_slot[gr][ts].push_back(i);
    ins(idx.lec_day[l][g.day], g.slot);
    for (int gr : groups) ins(idx.grp_day[gr][g.day], g.slot);
}

// (window_delta kept for potential future use in objective analysis)

// ── Local search: overlap-only repair ────────────────────────────────────────
// Repairs only hard constraint violations (overlaps). Ignores window penalties
// entirely — the GA handles soft constraints through crossover/mutation.
// For each conflicted gene, finds the first (day,slot,room) with zero overlaps.
// If none exists, picks the placement with minimum overlap count.
// Index updated incrementally — no full rebuilds.
void local_search(Schedule& s, const Problem& p, ScheduleIndex& idx) {
    const int k = std::max(1, p.V / 10);

    std::vector<std::pair<int,int>> scored;
    scored.reserve(p.V);
    for (int i = 0; i < p.V; ++i) {
        int l  = p.req_lecturer[s[i].req];
        int ts = s[i].day * p.slots + s[i].slot;
        int ov = (int)idx.lec_slot[l][ts].size() - 1;
        for (int g : p.req_groups[s[i].req]) ov += (int)idx.grp_slot[g][ts].size() - 1;
        ov += (int)idx.room_slot[s[i].room][ts].size() - 1;
        if (ov > 0) scored.push_back({ov, i});
    }
    std::partial_sort(scored.begin(),
                      scored.begin() + std::min(k, (int)scored.size()),
                      scored.end(), std::greater<>());
    int n = std::min(k, (int)scored.size());

    for (int ki = 0; ki < n; ++ki) {
        int i = scored[ki].second;
        int l = p.req_lecturer[s[i].req];
        const auto& groups = p.req_groups[s[i].req];

        // Overlap count only — no window tiebreaker
        auto overlaps = [&](int d, int sl, int r) -> int {
            int ts = d * p.slots + sl;
            int ov = (int)idx.lec_slot[l][ts].size() - 1;
            for (int g : groups) ov += (int)idx.grp_slot[g][ts].size() - 1;
            ov += (int)idx.room_slot[r][ts].size() - 1;
            return ov;
        };

        int cur_ov = overlaps(s[i].day, s[i].slot, s[i].room);
        if (cur_ov == 0) continue; // fixed by previous move

        int best_ov = cur_ov;
        int bd = s[i].day, bsl = s[i].slot, br = s[i].room;

        for (int d = 0; d < p.days; ++d)
            for (int sl = 0; sl < p.slots; ++sl)
                for (int r = 0; r < p.A; ++r) {
                    int ov = overlaps(d, sl, r);
                    if (ov < best_ov) { best_ov = ov; bd=d; bsl=sl; br=r; }
                    if (best_ov == 0) goto found;
                }
        found:
        if (bd == s[i].day && bsl == s[i].slot && br == s[i].room) continue;

        idx_remove(idx, i, s[i], l, groups);
        s[i].day = bd; s[i].slot = bsl; s[i].room = br;
        idx_insert(idx, i, s[i], l, groups);
    }
}

// ── 3-parent crossover ────────────────────────────────────────────────────────
Schedule crossover(const Schedule& s1, const Schedule& s2, const Schedule& s3,
                   int r, const Problem& p,
                   ScheduleIndex& idx1, ScheduleIndex& idx2, ScheduleIndex& idx3) {
    Schedule child = s1;
    auto bad = worst_genes(s1, r, p, idx1);
    for (int i : bad) {
        double l2 = local_objective(s2, i, p, idx2);
        double l3 = local_objective(s3, i, p, idx3);
        child[i] = (l2 <= l3) ? s2[i] : s3[i];
    }
    return child;
}

// ── Adaptive mutation ─────────────────────────────────────────────────────────
Schedule mutate(const Schedule& s, const Problem& p, std::mt19937& rng,
                const GAParams& params, double mr_scale, ScheduleIndex& idx) {
    Schedule child = s;
    int k = std::max(1, (int)(params.M * p.V));
    auto bad = worst_genes(s, k, p, idx);

    auto flip = [&](double base) {
        return std::uniform_real_distribution<>(0,1)(rng) < base * mr_scale;
    };
    for (int i : bad) {
        Gene& g = child[i];
        if (flip(params.MR_day))  g.day  = std::uniform_int_distribution<>(0, p.days  - 1)(rng);
        if (flip(params.MR_slot)) g.slot = std::uniform_int_distribution<>(0, p.slots - 1)(rng);
        if (flip(params.MR_room)) g.room = std::uniform_int_distribution<>(0, p.A     - 1)(rng);
        if (flip(params.MR_freq) && p.req_freq[g.req] != Freq::Weekly)
            g.freq = (g.freq == Freq::Numerator) ? Freq::Denominator : Freq::Numerator;
    }
    return child;
}

// ── Tournament selection ──────────────────────────────────────────────────────
std::vector<Schedule> select(std::vector<Schedule>& pool,
                              const std::vector<double>& scores,
                              int S, int elites, int tournament,
                              std::mt19937& rng) {
    std::vector<int> order(pool.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b){ return scores[a] < scores[b]; });

    std::vector<Schedule> next;
    next.reserve(S);
    for (int i = 0; i < elites && i < (int)pool.size(); ++i)
        next.push_back(pool[order[i]]);

    std::uniform_int_distribution<> pick(0, (int)pool.size() - 1);
    while ((int)next.size() < S) {
        int best = pick(rng);
        for (int t = 1; t < tournament; ++t) {
            int c = pick(rng);
            if (scores[c] < scores[best]) best = c;
        }
        next.push_back(pool[best]);
    }
    return next;
}

// ── Main GA loop ──────────────────────────────────────────────────────────────
struct IterRecord {
    int    iteration;
    double best_score;       // best objective so far
    double iter_score;       // best in this iteration's pool
    double iter_time_sec;    // wall time for this iteration
    double cumul_time_sec;   // cumulative wall time
    bool   restart;          // was a restart triggered this iteration
};

struct GAResult {
    Schedule best;
    double   best_score;
    int      iterations;
    double   elapsed_sec;
    std::vector<IterRecord> records;
};

GAResult run_ga(const Problem& p, const GAParams& params,
                std::mt19937& rng, bool verbose = true) {
    using clk = std::chrono::steady_clock;
    ScheduleIndex idx(p), idx2(p), idx3(p);

    auto eval_one = [&](const Schedule& s) {
        idx.build(s); return objective(s, p, idx);
    };
    auto eval_all = [&](const std::vector<Schedule>& v) {
        std::vector<double> sc(v.size());
        for (int i = 0; i < (int)v.size(); ++i) sc[i] = eval_one(v[i]);
        return sc;
    };

    std::vector<Schedule> pop(params.S);
    for (auto& s : pop) s = random_schedule(p, rng);

    std::vector<double> scores = eval_all(pop);
    int bi = (int)(std::min_element(scores.begin(), scores.end()) - scores.begin());
    double best_score = scores[bi];
    Schedule best = pop[bi];

    // (ls_threshold retained in GAParams for reference, not used in main loop)

    std::vector<IterRecord> records;
    int stagnation_cnt = 0;
    std::uniform_int_distribution<> pop_pick(0, params.S - 1);
    std::uniform_real_distribution<> rdist(0, 1);

    auto t_start = clk::now();
    auto elapsed = [&]{ return std::chrono::duration<double>(clk::now() - t_start).count(); };

    for (int iter = 1; iter <= params.max_iter; ++iter) {
        if (best_score == 0.0) break;
        auto t_iter = clk::now();

        double mr_scale = 1.0 + (double)stagnation_cnt / params.stagnation;
        int n_cross = std::max(1, (int)(params.C * params.S));

        std::vector<Schedule> offspring;
        offspring.reserve(n_cross + params.S);

        for (int c = 0; c < n_cross; ++c) {
            int i1 = pop_pick(rng), i2, i3;
            do { i2 = pop_pick(rng); } while (i2 == i1);
            do { i3 = pop_pick(rng); } while (i3 == i1 || i3 == i2);
            int r = (int)(params.CR_min * p.V)
                  + (int)(rdist(rng) * (params.CR_max - params.CR_min) * p.V);
            r = std::max(1, std::min(r, p.V - 1));
            idx.build(pop[i1]); idx2.build(pop[i2]); idx3.build(pop[i3]);
            offspring.push_back(crossover(pop[i1], pop[i2], pop[i3], r, p,
                                          idx, idx2, idx3));
        }

        for (const auto& s : pop) {
            idx.build(s);
            Schedule m = mutate(s, p, rng, params, mr_scale, idx);
            if (params.use_local_search) { idx.build(m); local_search(m, p, idx); }
            offspring.push_back(std::move(m));
        }

        std::vector<Schedule> pool;
        pool.reserve(params.S + offspring.size());
        for (auto& s : pop)       pool.push_back(std::move(s));
        for (auto& s : offspring) pool.push_back(std::move(s));

        std::vector<double> pool_scores = eval_all(pool);
        double iter_best = *std::min_element(pool_scores.begin(), pool_scores.end());

        for (int i = 0; i < (int)pool.size(); ++i) {
            if (pool_scores[i] < best_score) {
                best_score = pool_scores[i]; best = pool[i]; stagnation_cnt = 0;
            }
        }

        pop = select(pool, pool_scores, params.S, params.elites, params.tournament, rng);
        scores = eval_all(pop);
        ++stagnation_cnt;

        bool did_restart = false;
        if (stagnation_cnt >= params.stagnation) {
            int n_restart = (int)(params.restart_frac * params.S);
            std::vector<int> ord(params.S);
            std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](int a, int b){ return scores[a] > scores[b]; });
            for (int k = 0; k < n_restart; ++k) pop[ord[k]] = random_schedule(p, rng);
            scores = eval_all(pop);
            stagnation_cnt = 0;
            did_restart = true;
        }

        double iter_sec  = std::chrono::duration<double>(clk::now() - t_iter).count();
        double cumul_sec = elapsed();
        records.push_back({ iter, best_score, iter_best, iter_sec, cumul_sec, did_restart });

        if (verbose) {
            printf("  iter %5d | best=%10.4f | iter_best=%10.4f | iter_t=%.3fs | total=%.1fs%s\n",
                   iter, best_score, iter_best, iter_sec, cumul_sec,
                   did_restart ? " [RESTART]" : "");
            fflush(stdout);
        }
    }

    return { best, best_score, (int)records.size(), elapsed(), records };
}
