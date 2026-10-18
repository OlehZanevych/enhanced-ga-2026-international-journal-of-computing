// schedule.cpp — runs one experiment (one instance, one algorithm, one seed).
//
// Usage:
//   ./schedule <variant 0-4> <mode 0-3> <seed> [options]
//     mode 0 = GA (no local search)            -> "no_ls"
//     mode 1 = GA+LS   (overlap-only repair)    -> "ls"     (basic variant)
//     mode 2 = GA+LS-W (window-tiebreak repair) -> "ls_w"   (ablation)
//     mode 3 = GA+AMLS (adaptive multi-neighbourhood LS) -> "amls" (proposed)
//   options:
//     --max-iter N     iteration budget T_max (default 5000)
//     --time-limit S   additional wall-clock budget in seconds (default: none)
//     --ls-frac F      LS candidate fraction (default 0.10)
//     --out DIR        output base directory (default results/main)
//     --verbose        print one line per iteration
//
// Output: <out>/v<variant>_<mode>/seed<seed>/
//   iterations.csv  per-iteration trace (best score and its hard/soft parts)
//   summary.csv     one row with all parameters and final metrics
//   schedule.csv    best timetable found
#include "types.h"
#include "fitness.h"
#include "ga.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

static bool read_line(FILE* f, char* buf, int n) { return fgets(buf, n, f) != nullptr; }

static Problem load_problem(const std::string& dir) {
    Problem p;
    {
        FILE* f = fopen((dir + "/meta.csv").c_str(), "r");
        if (!f) { fprintf(stderr, "Cannot open %s/meta.csv\n", dir.c_str()); exit(1); }
        char line[256];
        if (!read_line(f, line, sizeof(line)) || !read_line(f, line, sizeof(line))) exit(1);
        if (sscanf(line, "%d,%d,%d,%d,%d,%d", &p.V, &p.L, &p.G, &p.A, &p.days, &p.slots) != 6) exit(1);
        fclose(f);
    }
    p.req_lecturer.resize(p.V);
    p.req_groups.resize(p.V);
    p.req_freq.resize(p.V);
    {
        FILE* f = fopen((dir + "/requirements.csv").c_str(), "r");
        if (!f) { fprintf(stderr, "Cannot open %s/requirements.csv\n", dir.c_str()); exit(1); }
        char line[1024];
        if (!read_line(f, line, sizeof(line))) exit(1);
        while (read_line(f, line, sizeof(line))) {
            int req, lec, fq; char gs[512] = "";
            if (sscanf(line, "%d,%d,%d,%511s", &req, &lec, &fq, gs) < 3) continue;
            p.req_lecturer[req] = lec;
            p.req_freq[req] = fq == 0 ? Freq::Weekly : fq == 1 ? Freq::Numerator : Freq::Denominator;
            for (char* tok = strtok(gs, ";"); tok; tok = strtok(nullptr, ";"))
                p.req_groups[req].push_back(atoi(tok));
        }
        fclose(f);
    }
    return p;
}

static std::string timestamp() {
    time_t t = time(nullptr);
    char buf[32];
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime(&t));
    return buf;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <variant 0-4> <mode 0-3> <seed> [--max-iter N] [--time-limit S] "
                        "[--ls-frac F] [--out DIR] [--verbose]\n", argv[0]);
        return 1;
    }
    const int variant = atoi(argv[1]);
    const int mode    = atoi(argv[2]);
    const int seed    = atoi(argv[3]);
    if (variant < 0 || variant > 4 || mode < 0 || mode > 3) {
        fprintf(stderr, "variant must be 0-4, mode must be 0-3\n"); return 1;
    }

    GAParams prm;                      // defaults = values of Table 2 of the article
    prm.ls_mode = static_cast<LSMode>(mode);
    std::string out_base = "results/main";
    bool verbose = false;
    for (int a = 4; a < argc; ++a) {
        const std::string opt = argv[a];
        auto next = [&]() -> const char* {
            if (a + 1 >= argc) { fprintf(stderr, "missing value for %s\n", opt.c_str()); exit(1); }
            return argv[++a];
        };
        if      (opt == "--max-iter")   prm.max_iter = atoi(next());
        else if (opt == "--time-limit") prm.time_limit_sec = atof(next());
        else if (opt == "--ls-frac")    prm.ls_fraction = atof(next());
        else if (opt == "--out")        out_base = next();
        else if (opt == "--verbose")    verbose = true;
        else { fprintf(stderr, "unknown option %s\n", opt.c_str()); return 1; }
    }

    const std::string mode_str = ls_mode_name(prm.ls_mode);
    const std::string data_dir = "data/variant" + std::to_string(variant);
    const std::string out_dir  = out_base + "/v" + std::to_string(variant) + "_" + mode_str
                               + "/seed" + std::to_string(seed);
    if (system(("mkdir -p '" + out_dir + "'").c_str()) != 0) { fprintf(stderr, "mkdir failed\n"); return 1; }

    const Problem p = load_problem(data_dir);

    printf("[%s] START variant=%d mode=%s seed=%d V=%d L=%d G=%d A=%d max_iter=%d time_limit=%.1f\n",
           timestamp().c_str(), variant, mode_str.c_str(), seed, p.V, p.L, p.G, p.A,
           prm.max_iter, prm.time_limit_sec);
    fflush(stdout);

    Rng rng(seed);
    const GAResult res = run_ga(p, prm, rng, verbose);

    printf("[%s] DONE  variant=%d mode=%s seed=%d best=%.4f hard=%.4f soft=%.4f iters=%d time=%.2fs\n",
           timestamp().c_str(), variant, mode_str.c_str(), seed, res.best_eval.total,
           res.best_eval.hard, res.best_eval.soft, res.iterations, res.elapsed_sec);
    fflush(stdout);

    // ── iterations.csv ────────────────────────────────────────────────────────
    {
        FILE* f = fopen((out_dir + "/iterations.csv").c_str(), "w");
        // (variant, mode and seed are given by the directory name and summary.csv)
        fprintf(f, "iteration,best_score,best_hard,best_soft,iter_best_score,cumul_time_sec,restart\n");
        for (const auto& r : res.records)
            fprintf(f, "%d,%g,%g,%g,%g,%.5f,%d\n", r.iteration, r.best_score, r.best_hard,
                    r.best_soft, r.iter_score, r.cumul_time_sec, r.restart ? 1 : 0);
        fclose(f);
    }

    // ── summary.csv ───────────────────────────────────────────────────────────
    {
        int restarts = 0;
        for (const auto& r : res.records) restarts += r.restart ? 1 : 0;
        const double avg_iter = res.iterations > 0 ? res.elapsed_sec / res.iterations : 0;
        const Eval& e = res.best_eval;
        FILE* f = fopen((out_dir + "/summary.csv").c_str(), "w");
        fprintf(f, "variant,mode,seed,V,L,G,A,days,slots,"
                   "pop_size,crossover_frac,CR_min,CR_max,mutation_frac,MR_day,MR_slot,MR_room,MR_freq,"
                   "elites,tournament,stagnation_limit,restart_frac,ls_mode,ls_fraction,"
                   "max_iter,time_limit_sec,iterations_run,converged,feasible,"
                   "final_best_score,hard_penalty,soft_penalty,pi1,pi2,pi3,pi4,pi5,"
                   "first_feasible_iter,first_feasible_sec,total_time_sec,avg_iter_time_sec,restarts_count,"
                   "ls_calls,ls_relocations,ls_swaps,ls_soft_moves,ls_time_sec,"
                   "ls_sampled,ls_conflicted_before,ls_conflicted_after\n");
        fprintf(f, "%d,%s,%d,%d,%d,%d,%d,%d,%d,"
                   "%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,"
                   "%d,%d,%d,%.2f,%d,%.3f,"
                   "%d,%.3f,%d,%d,%d,"
                   "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,"
                   "%d,%.6f,%.6f,%.6f,%d,"
                   "%lld,%lld,%lld,%lld,%.6f,%lld,%lld,%lld\n",
                variant, mode_str.c_str(), seed, p.V, p.L, p.G, p.A, p.days, p.slots,
                prm.S, prm.C, prm.CR_min, prm.CR_max, prm.M, prm.MR_day, prm.MR_slot, prm.MR_room, prm.MR_freq,
                prm.elites, prm.tournament, prm.stagnation, prm.restart_frac, mode, prm.ls_fraction,
                prm.max_iter, prm.time_limit_sec, res.iterations, e.total == 0.0 ? 1 : 0, e.hard == 0.0 ? 1 : 0,
                e.total, e.hard, e.soft, e.pi[0], e.pi[1], e.pi[2], e.pi[3], e.pi[4],
                res.first_feasible_iter, res.first_feasible_sec, res.elapsed_sec, avg_iter, restarts,
                res.ls.calls, res.ls.relocations, res.ls.swaps, res.ls.soft_moves, res.ls.time_sec,
                res.ls.sampled, res.ls.conflicted_before, res.ls.conflicted_after);
        fclose(f);
    }

    // ── schedule.csv ──────────────────────────────────────────────────────────
    {
        FILE* f = fopen((out_dir + "/schedule.csv").c_str(), "w");
        fprintf(f, "req_id,day,slot,room,freq\n");
        for (const auto& g : res.best) {
            const char* fs = g.freq == Freq::Weekly ? "weekly" : g.freq == Freq::Numerator ? "numerator" : "denominator";
            fprintf(f, "%d,%d,%d,%d,%s\n", g.req, g.day, g.slot, g.room, fs);
        }
        fclose(f);
    }
    return 0;
}
