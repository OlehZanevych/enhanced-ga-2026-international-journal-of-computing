// generate.cpp — generates 5 realistic UCTP input variants
// Fixed infrastructure (rooms, days, slots), growing academic load
#include "types.h"
#include <cstdio>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

static void write_variant(int variant, int V, int L, int G, int A,
                          int days = 5, int slots = 6) {
    std::string dir = "data/variant" + std::to_string(variant);
    system(("mkdir -p " + dir).c_str());

    std::mt19937 rng(42 + variant);
    int total_ts = days * slots; // 30

    // ── Lecturer assignment: load-balanced with realistic variance ───────
    // Target: max lecturer load ≤ 70% of timeslots
    std::vector<int> lec_load(L, 0);
    std::vector<int> req_lec(V);
    int max_lec_load = (int)(0.70 * total_ts);

    for (int i = 0; i < V; ++i) {
        std::vector<int> candidates;
        int min_load = *std::min_element(lec_load.begin(), lec_load.end());
        for (int l = 0; l < L; ++l)
            if (lec_load[l] <= min_load + 2 && lec_load[l] < max_lec_load)
                candidates.push_back(l);
        if (candidates.empty())
            for (int l = 0; l < L; ++l)
                if (lec_load[l] < max_lec_load) candidates.push_back(l);
        int l = candidates[rng() % candidates.size()];
        req_lec[i] = l;
        lec_load[l]++;
    }

    // ── Group assignment: load-balanced, 1-3 groups per req ─────────────
    std::vector<int> grp_load(G, 0);
    int max_grp_load = (int)(0.70 * total_ts);
    std::discrete_distribution<> ngroups_dist({55, 30, 15});

    std::vector<std::vector<int>> req_groups(V);
    for (int i = 0; i < V; ++i) {
        int ng = 1 + ngroups_dist(rng);
        std::vector<int> gs;
        for (int k = 0; k < ng; ++k) {
            std::vector<int> candidates;
            int min_load = *std::min_element(grp_load.begin(), grp_load.end());
            for (int g = 0; g < G; ++g) {
                if (std::find(gs.begin(), gs.end(), g) != gs.end()) continue;
                if (grp_load[g] <= min_load + 2 && grp_load[g] < max_grp_load)
                    candidates.push_back(g);
            }
            if (candidates.empty()) break;
            int g = candidates[rng() % candidates.size()];
            gs.push_back(g);
            grp_load[g]++;
        }
        req_groups[i] = gs;
    }

    // ── Frequency: ~65% weekly, ~20% numerator, ~15% denominator ────────
    std::discrete_distribution<> freq_dist({65, 20, 15});
    std::vector<int> req_freq(V);
    for (int i = 0; i < V; ++i)
        req_freq[i] = freq_dist(rng);

    // --- requirements.csv ---
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

    // --- meta.csv ---
    {
        FILE* f = fopen((dir + "/meta.csv").c_str(), "w");
        fprintf(f, "V,L,G,A,days,slots\n");
        fprintf(f, "%d,%d,%d,%d,%d,%d\n", V, L, G, A, days, slots);
        fclose(f);
    }

    int max_ll = *std::max_element(lec_load.begin(), lec_load.end());
    int max_gl = *std::max_element(grp_load.begin(), grp_load.end());
    double util = 100.0 * V / (days * slots * A);
    printf("Variant %d: V=%-4d L=%-3d G=%-3d A=%-3d  "
           "room_util=%.0f%%  max_lec=%d/%d(%.0f%%)  max_grp=%d/%d(%.0f%%)\n",
           variant, V, L, G, A, util,
           max_ll, total_ts, 100.0*max_ll/total_ts,
           max_gl, total_ts, 100.0*max_gl/total_ts);
}

int main() {
    struct Config { int V, L, G, A; };
    Config configs[] = {
        { 30,  8,  6,  5},   // small dept   (util=20%)
        { 60, 12, 10,  8},   // medium dept  (util=25%)
        {120, 20, 15, 10},   // large dept   (util=40%)
        {250, 35, 25, 15},   // faculty      (util=56%)
        {500, 60, 40, 20},   // university   (util=83%)
    };

    for (int k = 0; k < 5; ++k)
        write_variant(k, configs[k].V, configs[k].L, configs[k].G, configs[k].A);

    printf("Done. Generated 5 variants in data/\n");
    return 0;
}
