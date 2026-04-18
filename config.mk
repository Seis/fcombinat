# fcombinat configuration
#
# This is the single place to configure binary paths for all backends.
# A backend is automatically skipped if its binary is not found or not executable.

# fcombinat binary (built at repo root by 'make')
export FCOMBINAT_BIN ?= ./fcombinat

# Maple CAS — override if Maple is not in PATH, e.g.:
#   make bench MAPLE_BIN=/opt/maple2026/bin/maple
export MAPLE_BIN     ?= $(HOME)/maple2026/bin/maple

# Parallelism - number of concurrent test jobs
JOBS ?= $(shell nproc)

# Number of terms to verify in correctness tests
NTERMS ?= 10

# Max size for exhaustive uniqueness enumeration (--enum-all)
ENUM_DEPTH ?= 5
