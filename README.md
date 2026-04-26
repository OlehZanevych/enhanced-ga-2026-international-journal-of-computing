# Enhanced GA for University Course Timetabling

This application implements an enhanced genetic algorithm (GA) with an optional conflict-driven local search (LS) procedure for the University Course Timetabling Problem (UCTP). It is developed for the purpose of conducting experiments for the article:

> **"Enhanced Genetic Algorithm for University Course Timetabling with Conflict-Driven Local Search"**
> by Oleh Zanevych and Vitaliy Kukharskyy
>
> Planned for publication in [*International Journal of Computing*](https://www.computingonline.net/computing)

## Project Structure

```
enhanced-ga/
├── types.h              # Core data types: Gene, Schedule, Problem, Freq enum
├── fitness.h            # Objective function, local objective, ScheduleIndex (O(1) bucket-index)
├── ga.h                 # GA engine: crossover, mutation, selection, local search, main loop
├── generate.cpp         # Test instance generator (5 realistic variants)
├── schedule.cpp         # Main experiment runner (reads data, runs GA, writes results)
├── CMakeLists.txt       # CMake build configuration
├── run_experiments.sh   # Automated experiment script (build + generate + run all)
├── data/                # Generated input data (5 variants)
│   ├── variant0/        # V=30,  L=8,  G=6,  A=5  (20% room utilisation)
│   ├── variant1/        # V=60,  L=12, G=10, A=8  (25% room utilisation)
│   ├── variant2/        # V=120, L=20, G=15, A=10 (40% room utilisation)
│   ├── variant3/        # V=250, L=35, G=25, A=15 (56% room utilisation)
│   └── variant4/        # V=500, L=60, G=40, A=20 (83% room utilisation)
│       ├── meta.csv           # Problem dimensions (V, L, G, A, days, slots)
│       └── requirements.csv   # Class requirements (lecturer, frequency, groups)
└── results/             # Experiment output (created by run_experiments.sh)
    ├── run.log          # Full execution log with system info and timing
    ├── DONE.log         # Completion report (written only on success)
    └── v<k>_<mode>/seed<N>/   # Per-run output
        ├── iterations.csv     # Per-iteration convergence trace
        ├── summary.csv        # Single-row summary with all parameters and metrics
        └── schedule.csv       # Best timetable found
```

## Source Files

| File | Description |
|------|-------------|
| `types.h` | Defines `Gene` (day, slot, room, frequency), `Schedule` (vector of genes), `Problem` (dimensions, lecturer/group assignments), and `Freq` enum (Weekly, Numerator, Denominator). |
| `fitness.h` | Implements the weighted penalty objective function (Eq. 1 in the article). Contains `ScheduleIndex` — a bucket-index structure providing O(1) conflict lookups for lecturers, groups, and rooms. Also provides `local_objective` for per-gene conflict evaluation used by crossover and mutation. |
| `ga.h` | The complete GA engine: random schedule initialisation, three-parent conflict-aware crossover, adaptive mutation with stagnation scaling, tournament selection with elite preservation, stagnation-triggered restart, and the conflict-driven local search procedure (overlap-only repair). The `run_ga()` function executes the main evolutionary loop and returns per-iteration records. |
| `generate.cpp` | Generates 5 test instances with load-balanced lecturer and group assignments. Lecturers are assigned using a min-load heuristic (70% max utilisation cap). Groups follow a weighted distribution (55% single, 30% double, 15% triple). Frequencies: 65% weekly, 20% numerator, 15% denominator. |
| `schedule.cpp` | Entry point for a single experiment run. Reads problem data from CSV, configures GA parameters, runs the algorithm with a given seed, and writes `iterations.csv`, `summary.csv`, and `schedule.csv` to the output directory. |
| `run_experiments.sh` | Orchestrates the full experiment pipeline: compiles both executables, generates input data, and runs all 500 experiments (5 variants × 2 modes × 50 seeds). Logs system info, per-configuration timing, and total wall time. |

## Algorithm Configurations

Two configurations are compared, differing only in the local search toggle:

- **GA** (`mode=0`): Enhanced GA without local search.
- **GA+LS** (`mode=1`): Enhanced GA with conflict-driven local search applied to every mutant offspring.

Both share identical GA parameters: population size 20, three-parent crossover (C=0.5), adaptive mutation (M=0.10), tournament selection (k=3), elite preservation (2), stagnation restart after 30 iterations, and a maximum of 5000 iterations.

## Requirements

- C++17 compiler (tested with Apple Clang 15.0)
- Bash shell
- No external libraries required

## Quick Start

### Run all experiments (recommended)

```bash
cd enhanced-ga
./run_experiments.sh
```

This will:
1. Build `generate` and `schedule` executables
2. Generate 5 test instances in `data/`
3. Run 500 experiments (5 variants × 2 modes × 50 seeds 1000–1049)
4. Write results to `results/`
5. Log everything to `results/run.log`

Total runtime: ~50 minutes on Apple M1 Pro.

### Manual execution

```bash
# Build
clang++ -std=c++17 -O3 -Wall -o generate generate.cpp
clang++ -std=c++17 -O3 -Wall -o schedule schedule.cpp

# Generate input data
./generate

# Run a single experiment
./schedule <variant 0-4> <mode 0=GA 1=GA+LS> <seed>

# Examples
./schedule 0 0 1000   # Variant 0 (V=30), GA only, seed 1000
./schedule 4 1 1042   # Variant 4 (V=500), GA+LS, seed 1042
```

### Build with CMake

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Output Format

### `summary.csv`

Contains one row per run with all experiment metadata and results:

```
variant, mode, seed, V, L, G, A, days, slots,
pop_size, crossover_frac, CR_min, CR_max,
mutation_frac, MR_day, MR_slot, MR_room, MR_freq,
elites, tournament, stagnation_limit, restart_frac,
use_local_search, max_iter,
iterations_run, converged,
final_best_score, total_time_sec,
avg_iter_time_sec, restarts_count
```

### `iterations.csv`

Per-iteration convergence trace:

```
variant, mode, seed, V, L, G, A, days, slots, max_iter,
iteration, best_score, iter_best_score,
iter_time_sec, cumul_time_sec, restart
```

### `schedule.csv`

The best timetable found:

```
req_id, day, slot, room, freq
```
