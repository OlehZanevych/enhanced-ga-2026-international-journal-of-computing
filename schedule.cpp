// schedule.cpp
// Usage: ./schedule <variant 0-4> <mode 0=no_ls 1=ls> <seed>
// Outputs to results/v<k>_<mode>/seed<N>/
//   iterations.csv  — per-iteration data
//   summary.csv     — single-row run summary with all parameters
//   schedule.csv    — best schedule found
#include "types.h"
#include "fitness.h"
#include "ga.h"
#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

static Problem load_problem(const std::string& dir) {
    Problem p;
    {
        FILE* f = fopen((dir + "/meta.csv").c_str(), "r");
        if (!f) { fprintf(stderr, "Cannot open %s/meta.csv\n", dir.c_str()); exit(1); }
        char line[256];
        fgets(line, sizeof(line), f);
        fgets(line, sizeof(line), f);
        sscanf(line, "%d,%d,%d,%d,%d,%d", &p.V,&p.L,&p.G,&p.A,&p.days,&p.slots);
        fclose(f);
    }
    p.req_lecturer.resize(p.V);
    p.req_groups.resize(p.V);
    p.req_freq.resize(p.V);
    {
        FILE* f = fopen((dir + "/requirements.csv").c_str(), "r");
        if (!f) { fprintf(stderr, "Cannot open %s/requirements.csv\n", dir.c_str()); exit(1); }
        char line[1024];
        fgets(line, sizeof(line), f);
        while (fgets(line, sizeof(line), f)) {
            int req, lec, fq; char gs[512] = "";
            sscanf(line, "%d,%d,%d,%511s", &req, &lec, &fq, gs);
            p.req_lecturer[req] = lec;
            p.req_freq[req] = fq==0 ? Freq::Weekly : fq==1 ? Freq::Numerator : Freq::Denominator;
            char* tok = strtok(gs, ";");
            while (tok) { p.req_groups[req].push_back(atoi(tok)); tok = strtok(nullptr,";"); }
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
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <variant 0-4> <mode 0=no_ls 1=ls> <seed>\n", argv[0]);
        return 1;
    }
    int variant = atoi(argv[1]);
    int mode    = atoi(argv[2]);
    int seed    = atoi(argv[3]);
    if (variant < 0 || variant > 4 || (mode != 0 && mode != 1)) {
        fprintf(stderr, "variant must be 0-4, mode must be 0 or 1\n"); return 1;
    }

    std::string mode_str  = mode ? "ls" : "no_ls";
    std::string data_dir  = "data/variant" + std::to_string(variant);
    std::string out_dir   = "results/v" + std::to_string(variant) + "_" + mode_str
                          + "/seed" + std::to_string(seed);
    system(("mkdir -p " + out_dir).c_str());

    Problem p = load_problem(data_dir);

    GAParams params;
    params.S = 20; params.C = 0.5; params.CR_min = 0.05; params.CR_max = 0.20;
    params.M = 0.10;
    params.MR_day = params.MR_slot = params.MR_room = params.MR_freq = 0.25;
    params.elites = 2; params.tournament = 3;
    params.stagnation = 30; params.restart_frac = 0.5;
    params.max_iter = 5000;
    params.use_local_search = (mode == 1);

    printf("[%s] START variant=%d mode=%s seed=%d V=%d L=%d G=%d A=%d max_iter=%d\n",
           timestamp().c_str(), variant, mode_str.c_str(), seed,
           p.V, p.L, p.G, p.A, params.max_iter);
    fflush(stdout);

    std::mt19937 rng(seed);
    GAResult res = run_ga(p, params, rng, /*verbose=*/false);

    printf("[%s] DONE  variant=%d mode=%s seed=%d best=%.4f iters=%d time=%.2fs\n",
           timestamp().c_str(), variant, mode_str.c_str(), seed,
           res.best_score, res.iterations, res.elapsed_sec);
    fflush(stdout);

    // ── iterations.csv ────────────────────────────────────────────────────────
    {
        FILE* f = fopen((out_dir + "/iterations.csv").c_str(), "w");
        fprintf(f, "variant,mode,seed,V,L,G,A,days,slots,max_iter,"
                   "iteration,best_score,iter_best_score,"
                   "iter_time_sec,cumul_time_sec,restart\n");
        for (const auto& r : res.records) {
            fprintf(f, "%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,"
                       "%d,%.6f,%.6f,%.6f,%.6f,%d\n",
                    variant, mode_str.c_str(), seed,
                    p.V, p.L, p.G, p.A, p.days, p.slots, params.max_iter,
                    r.iteration, r.best_score, r.iter_score,
                    r.iter_time_sec, r.cumul_time_sec, r.restart ? 1 : 0);
        }
        fclose(f);
    }

    // ── summary.csv ───────────────────────────────────────────────────────────
    {
        FILE* f = fopen((out_dir + "/summary.csv").c_str(), "w");
        fprintf(f, "variant,mode,seed,V,L,G,A,days,slots,"
                   "pop_size,crossover_frac,CR_min,CR_max,"
                   "mutation_frac,MR_day,MR_slot,MR_room,MR_freq,"
                   "elites,tournament,stagnation_limit,restart_frac,"
                   "use_local_search,max_iter,"
                   "iterations_run,converged,"
                   "final_best_score,total_time_sec,"
                   "avg_iter_time_sec,restarts_count\n");
        int restarts = 0;
        for (const auto& r : res.records) restarts += r.restart ? 1 : 0;
        double avg_iter = res.iterations > 0 ? res.elapsed_sec / res.iterations : 0;
        fprintf(f, "%d,%s,%d,%d,%d,%d,%d,%d,%d,"
                   "%d,%.2f,%.2f,%.2f,"
                   "%.2f,%.2f,%.2f,%.2f,%.2f,"
                   "%d,%d,%d,%.2f,"
                   "%d,%d,"
                   "%d,%d,"
                   "%.6f,%.6f,"
                   "%.6f,%d\n",
                variant, mode_str.c_str(), seed,
                p.V, p.L, p.G, p.A, p.days, p.slots,
                params.S, params.C, params.CR_min, params.CR_max,
                params.M, params.MR_day, params.MR_slot, params.MR_room, params.MR_freq,
                params.elites, params.tournament, params.stagnation, params.restart_frac,
                params.use_local_search ? 1 : 0, params.max_iter,
                res.iterations, res.best_score == 0.0 ? 1 : 0,
                res.best_score, res.elapsed_sec,
                avg_iter, restarts);
        fclose(f);
    }

    // ── schedule.csv ──────────────────────────────────────────────────────────
    {
        FILE* f = fopen((out_dir + "/schedule.csv").c_str(), "w");
        fprintf(f, "req_id,day,slot,room,freq\n");
        for (const auto& g : res.best) {
            const char* fstr = g.freq==Freq::Weekly ? "weekly"
                             : g.freq==Freq::Numerator ? "numerator" : "denominator";
            fprintf(f, "%d,%d,%d,%d,%s\n", g.req, g.day, g.slot, g.room, fstr);
        }
        fclose(f);
    }

    return 0;
}
