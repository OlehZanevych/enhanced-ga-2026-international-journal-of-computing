// generate.cpp — generates the 5 synthetic UCTP instances used in the article.
// Fixed infrastructure (D = 5 days, S = 6 slots per day), growing academic load.
// Uses the portable RNG of rng.h, so the generated files are byte-identical on
// every platform (and identical to the files produced by the original build).
#include "types.h"
#include "rng.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static void write_variant(int variant, int V, int L, int G, int A,
                          int days = 5, int slots = 6) {
    std::string dir = "data/variant" + std::to_string(variant);
    if (system(("mkdir -p " + dir).c_str()) != 0) { fprintf(stderr, "mkdir failed\n"); exit(1); }

    Rng rng(42 + variant);
    const int total_ts = days * slots;   // 30

    // ── Lecturer assignment: load-balanced with random tie-breaking ──────────
    // Candidates: lecturers within 2 classes of the minimum load, below the
    // 70% cap of the available timeslots.
    std::vector<int> lec_load(L, 0), req_lec(V);
    const int max_lec_load = (int)(0.70 * total_ts);
    for (int i = 0; i < V; ++i) {
        std::vector<int> candidates;
        const int min_load = *std::min_element(lec_load.begin(), lec_load.end());
        for (int l = 0; l < L; ++l)
            if (lec_load[l] <= min_load + 2 && lec_load[l] < max_lec_load) candidates.push_back(l);
        if (candidates.empty())
            for (int l = 0; l < L; ++l)
                if (lec_load[l] < max_lec_load) candidates.push_back(l);
        const int l = candidates[rng() % candidates.size()];
        req_lec[i] = l;
        lec_load[l]++;
    }

    // ── Group assignment: load-balanced, 1-3 groups per requirement ──────────
    std::vector<int> grp_load(G, 0);
    const int max_grp_load = (int)(0.70 * total_ts);
    DiscreteDist ngroups_dist({55, 30, 15});
    std::vector<std::vector<int>> req_groups(V);
    for (int i = 0; i < V; ++i) {
        const int ng = 1 + ngroups_dist(rng);
        std::vector<int> gs;
        for (int k = 0; k < ng; ++k) {
            std::vector<int> candidates;
            const int min_load = *std::min_element(grp_load.begin(), grp_load.end());
            for (int g = 0; g < G; ++g) {
                if (std::find(gs.begin(), gs.end(), g) != gs.end()) continue;
                if (grp_load[g] <= min_load + 2 && grp_load[g] < max_grp_load) candidates.push_back(g);
            }
            if (candidates.empty()) break;
            const int g = candidates[rng() % candidates.size()];
            gs.push_back(g);
            grp_load[g]++;
        }
        req_groups[i] = gs;
    }

    // ── Frequency: ~65% weekly, ~20% numerator, ~15% denominator ─────────────
    DiscreteDist freq_dist({65, 20, 15});
    std::vector<int> req_freq(V);
    for (int i = 0; i < V; ++i) req_freq[i] = freq_dist(rng);

    {
        FILE* f = fopen((dir + "/requirements.csv").c_str(), "w");
        fprintf(f, "req_id,lecturer_id,freq,group_ids\n");
        for (int i = 0; i < V; ++i) {
            fprintf(f, "%d,%d,%d,", i, req_lec[i], req_freq[i]);
            for (int k = 0; k < (int)req_groups[i].size(); ++k)
                fprintf(f, k ? ";%d" : "%d", req_groups[i][k]);
            fprintf(f, "\n");
        }
        fclose(f);
    }
    {
        FILE* f = fopen((dir + "/meta.csv").c_str(), "w");
        fprintf(f, "V,L,G,A,days,slots\n");
        fprintf(f, "%d,%d,%d,%d,%d,%d\n", V, L, G, A, days, slots);
        fclose(f);
    }

    // Load statistics (weekly-equivalent load: biweekly classes count 0.5).
    std::vector<double> lec_w(L, 0), grp_w(G, 0);
    double weq = 0;
    for (int i = 0; i < V; ++i) {
        const double w = req_freq[i] == 0 ? 1.0 : 0.5;
        weq += w; lec_w[req_lec[i]] += w;
        for (int g : req_groups[i]) grp_w[g] += w;
    }
    const int max_ll = *std::max_element(lec_load.begin(), lec_load.end());
    const int max_gl = *std::max_element(grp_load.begin(), grp_load.end());
    int multi = 0;
    for (const auto& gs : req_groups) multi += gs.size() > 1;
    printf("Variant %d: V=%-4d L=%-3d G=%-3d A=%-3d  room_util=%.0f%% (weekly-eq %.0f%%)  "
           "max_lec=%d/%d(%.0f%%)  max_grp=%d/%d(%.0f%%)  max_lec_weq=%.1f  max_grp_weq=%.1f  "
           "multi_group=%d  biweekly=%d\n",
           variant, V, L, G, A, 100.0 * V / (days * slots * A), 100.0 * weq / (days * slots * A),
           max_ll, total_ts, 100.0 * max_ll / total_ts, max_gl, total_ts, 100.0 * max_gl / total_ts,
           *std::max_element(lec_w.begin(), lec_w.end()), *std::max_element(grp_w.begin(), grp_w.end()),
           multi, (int)std::count_if(req_freq.begin(), req_freq.end(), [](int f) { return f != 0; }));
}

int main() {
    struct Config { int V, L, G, A; };
    const Config configs[] = {
        { 30,  8,  6,  5},   // small department (util = 20%)
        { 60, 12, 10,  8},   // medium department (25%)
        {120, 20, 15, 10},   // large department (40%)
        {250, 35, 25, 15},   // faculty (56%)
        {500, 60, 40, 20},   // university (83%)
    };
    for (int k = 0; k < 5; ++k)
        write_variant(k, configs[k].V, configs[k].L, configs[k].G, configs[k].A);
    printf("Done. Generated 5 variants in data/\n");
    return 0;
}
