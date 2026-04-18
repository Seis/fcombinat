# Benchmark Suite — Combinatorial Structure Evaluation

A modular benchmark framework for comparing **counting**, **unranking**, and **ranking** performance across multiple combinatorial structures and algorithm backends.

## Structures

| Structure | Grammar Spec | Mode | OEIS |
|-----------|-------------|------|------|
| Permutations | `Set(Cycle(Z))` | Labeled | A000142 |
| Integer Partitions | `Set(Sequence(Z, card >= 1))` | Unlabeled | A000041 |
| Distinct Partitions | `PowerSet(Sequence(Z, card >= 1))` | Unlabeled | A000009 |
| Surjections | `Sequence(Set(Z, card >= 1))` | Labeled | A000670 |
| Cycles | `Cycle(Z)` | Labeled | A000142 |

## Quick Start

```bash
# Build fcombinat + run all benchmarks + generate plots
make all

# Run tests for a specific backend only
make all fcombinat
make maple

# Or step by step:
make build         # Build fcombinat binary
make run           # Run all backends × all structures
make plot          # Generate charts from results/

# Filter by backend or structure:
./run_suite.sh --backend fcombinat --structure permutations
```

## Directory Layout

```
benchmarks/
├── README.md                     # This file
├── Makefile
├── run_suite.sh                  # Main runner
├── structures/                   # Structure definitions
│   ├── permutations.conf
│   ├── integer_partitions.conf
│   ├── distinct_partitions.conf
│   ├── surjections.conf
│   └── cycles.conf
├── backends/                     # Backend plugins
│   ├── fcombinat/
│   │   ├── info.conf
│   │   └── run.sh
│   └── maple/
│       ├── info.conf
│       └── run.sh
├── results/                      # Output (auto-created)
│   ├── fcombinat_permutations.csv
│   ├── fcombinat_cycles.csv
│   └── *.png
└── plot_results.py
```

## Adding a New Structure

Create `structures/your_structure.conf`:

```bash
NAME="YourStructure"
SPEC='S = Set(Cycle(Z, card >= 2))'   # fcombinat grammar
LABELED=1                              # 1=labeled, 0=unlabeled
SYMBOL="S"                             # Target symbol in the grammar
N_VALUES="10 20 50 100 200 500"        # Space-separated N values
```

The suite will automatically pick up any `.conf` file in `structures/`.

---

## Adding a New Backend

This is the key extensibility point. To add a new algorithm (e.g., Maple, SageMath, Haskell):

### 1. Create the backend directory

```bash
mkdir -p backends/maple
```

### 2. Create `backends/maple/info.conf`

```bash
BACKEND_NAME="maple"
BACKEND_DESC="Maple CAS — Iterator and combstruct packages"
```

### 3. Create `backends/maple/run.sh`

This script must implement **3 commands** with a standard interface:

```bash
#!/bin/bash
# backends/maple/run.sh
# Standard interface for the benchmark suite.
#
# Arguments:
#   $1 = ACTION:  "count" | "unrank" | "rank"
#   $2 = SPEC:    Grammar specification string (e.g., 'P = Set(Cycle(Z))')
#   $3 = N:       Size parameter
#   $4 = LABELED: "1" (labeled) or "0" (unlabeled)
#   $5 = SYMBOL:  Target symbol (e.g., "P")
#   $6 = RANK or OBJECT (depending on action)
#
# Output: print result to stdout
# Exit:   0 = success, non-zero = failure

ACTION="$1"
SPEC="$2"
N="$3"
LABELED="$4"
SYMBOL="$5"

case "$ACTION" in
    count)
        # Return the number of objects of size N
        # Example: for Permutations n=5, print "120"
        maple -q -c "
            with(combstruct);
            spec := {$SPEC, Z = Atom};
            printf(\"%d\", count([$SYMBOL, spec, labeled], size=$N));
            quit;
        "
        ;;

    unrank)
        RANK="$6"
        # Return the object at the given rank
        # Must print a string representation that the rank command can parse back
        maple -q -c "
            with(combstruct);
            spec := {$SPEC, Z = Atom};
            obj := draw([$SYMBOL, spec, labeled], size=$N);
            # ... implement rank-based selection ...
            printf(\"%a\", obj);
            quit;
        "
        ;;

    rank)
        OBJECT="$6"
        # Return the integer rank of the given object
        maple -q -c "
            # ... implement ranking ...
            quit;
        "
        ;;
esac
```

### 4. Make it executable

```bash
chmod +x backends/maple/run.sh
```

### 5. Run it

```bash
# The suite auto-discovers backends in backends/
./run_suite.sh --backend maple

# Or run everything
./run_suite.sh
```

### Interface Contract

| Command | Input | Output (stdout) | Example |
|---------|-------|-----------------|---------|
| `count` | spec, n, labeled, symbol | Integer count | `120` |
| `unrank` | spec, n, labeled, symbol, rank | Object string | `Set(Cycle(3,1,2), Cycle(4,5))` |
| `rank` | spec, n, labeled, symbol, object | Integer rank | `42` |

**Important notes:**
- The **object string** returned by `unrank` must be parseable by `rank` in the same backend
- Cross-backend object comparison is not required (each backend outputs its own format)
- Exit code `0` = success, anything else = the runner records it as a failure and skips
- The runner adds timing externally — do **not** add timing inside `run.sh`
- If an operation is not supported, exit with code `1`

## Output Format

Each benchmark run produces `results/<backend>_<structure>.csv`:

```csv
n,count_time,count_value,unrank_time,rank_time,total_time
10,0.000100,3628800,0.000050,0.000040,0.000090
20,0.000300,2432902008176640000,0.000120,0.000100,0.000220
```

Charts are saved as PNG files in `results/`:
- `results/<structure>.png` — per-structure comparison across backends
- `results/<backend>_overview.png` — per-backend overview across structures

## Configuration

The runner supports these options:

```bash
./run_suite.sh --backend fcombinat    # Only run one backend
./run_suite.sh --structure cycles     # Only run one structure
./run_suite.sh --trials 5             # Number of trials (default: 3)
./run_suite.sh --timeout 600          # Timeout in seconds (default: 300)
```

## Prerequisites

- **fcombinat backend**: `gcc`, `make`, FLINT library
- **maple backend**: Maple CAS (`maple` in PATH)
- **Plotting**: `python3`, `matplotlib`
- **General**: `bc`, `date`, `timeout` (GNU coreutils)
