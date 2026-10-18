# Enhanced GA for University Course Timetabling

This application implements an enhanced genetic algorithm (GA) with optional conflict-driven local search (LS) for the University Course Timetabling Problem (UCTP). It was developed to run the experiments for the article:

> **"Enhanced Genetic Algorithm for University Course Timetabling with Conflict-Driven Local Search"**
> by Oleh Zanevych and Vitaliy Kukharskyy
>
> Planned for publication in [*International Journal of Computing*](https://www.computingonline.net/computing)

## Algorithms

| Mode | Name in the article | Description |
|------|--------------------|-------------|
| `0` (`no_ls`) | GA | Enhanced GA: three-parent conflict-aware crossover, adaptive mutation of the worst genes, (μ+λ) survivor selection with elitism and tournament, stagnation-triggered partial restart. |
| `1` (`ls`) | GA+LS | Basic variant: GA + overlap-only conflict-driven repair of every mutant. |
| `2` (`ls_w`) | GA+LS-W | Ablation: the same repair, with window count used as a tiebreaker among conflict-free placements. |
| `3` (`amls`) | GA+AMLS | Proposed: adaptive multi-neighbourhood LS — relocation + ejection swaps for hard constraints, then a separate feasibility-preserving window-compaction move with exact delta evaluation; intensity grows with stagnation. |

## Project Structure

```
├── types.h              # Gene, Schedule, Problem, Freq
├── rng.h                # Portable RNG (libc++-compatible distributions, identical on every platform)
├── fitness.h            # Objective with hard/soft decomposition, local objective, ScheduleIndex
├── ga.h                 # GA engine and the three local-search variants
├── generate.cpp         # Instance generator (5 synthetic instances)
├── schedule.cpp         # Runs one experiment (instance, algorithm, seed)
├── run_experiments.sh   # Complete experimental protocol (build + generate + all runs)
├── CMakeLists.txt
├── analyze.py           # Statistics, LaTeX tables and pgfplots figures of the article
├── data/variant{0..4}/  # meta.csv, requirements.csv
└── results/
    ├── system_info.txt, run.log, DONE.log, equal_time_budgets.txt
    ├── main/        v<k>_<mode>/seed<N>/{iterations,summary,schedule}.csv + all_summaries.csv
    └── equal_time/  (same layout)
```

## Running the experiments

```bash
./run_experiments.sh
```

The script deletes `results/`, builds both executables (`clang++` if available, otherwise `g++`; override with `CXX=...`), regenerates the instances (and checks that they are identical to the reference instances), and runs, sequentially:

| Experiment | Runs | Description |
|-----------|------|-------------|
| E1 `main` | 5 instances × 4 algorithms × 50 seeds = 1000 | Iteration budget T_max = 5000 |
| E2 `equal_time` | V ∈ {250, 500} × 4 algorithms × 50 seeds = 400 | Equal wall-clock budget = mean GA time of E1 for that instance |

Total runtime on an Apple M4 Max (MacBook Pro, macOS 15.1, Apple clang 16): about 70 minutes. On macOS the script keeps the computer awake with `caffeinate`. An interrupted run can be continued with `RESUME=1 ./run_experiments.sh`; `NUM_SEEDS=2 ./run_experiments.sh` performs a quick smoke test.

### Single runs

```bash
clang++ -std=c++17 -O3 -Wall -o generate generate.cpp
clang++ -std=c++17 -O3 -Wall -o schedule schedule.cpp
./generate
./schedule <variant 0-4> <mode 0-3> <seed> [--max-iter N] [--time-limit S] [--ls-frac F] [--out DIR] [--verbose]

./schedule 4 1 1042            # V=500, GA+LS, seed 1042 -> results/main/v4_ls/seed1042/
./schedule 3 3 1000 --verbose  # V=250, GA+AMLS, print every iteration
```

## Reproducibility

* `std::mt19937` is fully specified by the C++ standard, but the standard distributions are not. `rng.h` re-implements them with the arithmetic used by libc++, so a given seed produces bit-identical runs with GCC/libstdc++, Clang/libc++ and on x86-64 or ARM64. Ties in all sorting steps are broken by index. Only the wall-clock times are machine-dependent.

### Analysis

```bash
python3 analyze.py --results results --figures ../figures --tables analysis_out
```

Requires `numpy` and `scipy`. Writes the pgfplots figures and the LaTeX table bodies used in the article, `numbers.json` with every value quoted in the text, and a short `report.md`.

## Output format

`summary.csv` contains one row per run with all parameters plus: `iterations_run`, `converged` (f = 0), `feasible` (all hard constraints satisfied), `final_best_score`, `hard_penalty`, `soft_penalty`, `pi1`–`pi5`, `first_feasible_iter`, `first_feasible_sec`, `total_time_sec`, `avg_iter_time_sec`, `restarts_count`, and LS statistics: `ls_calls`, `ls_relocations`, `ls_swaps`, `ls_soft_moves`, `ls_time_sec`, and the number of conflicted genes in sampled mutants before and after LS (`ls_sampled`, `ls_conflicted_before`, `ls_conflicted_after`).

`iterations.csv`: `iteration, best_score, best_hard, best_soft, iter_best_score, cumul_time_sec, restart`.

`schedule.csv`: `req_id, day, slot, room, freq`.

## Requirements

C++17 compiler (tested with Apple Clang, Clang 18 and GCC 11/13), bash, no external libraries.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
