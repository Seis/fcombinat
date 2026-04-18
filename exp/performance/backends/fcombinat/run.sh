#!/bin/bash
##############################################################################
# backends/fcombinat/run.sh — fcombinat backend for benchmark suite
#
# Standard interface:
#   ./run.sh count    <spec> <n> <labeled> <symbol>
#   ./run.sh unrank   <spec> <n> <labeled> <symbol> <rank>
#   ./run.sh rank     <spec> <n> <labeled> <symbol> <object>
#   ./run.sh draw     <spec> <n> <labeled> <symbol>
#   ./run.sh boltzmann <spec> <n> <labeled> <symbol>
#
# Prints result to stdout. Exit 0 = success.
##############################################################################

set -euo pipefail
export LC_NUMERIC=C

FCOMBINAT="${FCOMBINAT_BIN:-}"
if [ -z "$FCOMBINAT" ] || [ ! -x "$FCOMBINAT" ]; then
    echo "Error: fcombinat binary not found. Set FCOMBINAT_BIN in config.mk." >&2
    exit 1
fi

ACTION="$1"
SPEC="$2"
N="$3"
LABELED="$4"
SYMBOL="$5"

# Build labeled flag
LABELED_FLAG=""
if [ "$LABELED" = "1" ]; then
    LABELED_FLAG="--labeled"
fi

case "$ACTION" in
    count)
        # Run --terms and extract the coefficient at position N
        # fcombinat outputs one coefficient per line: a_0, a_1, ..., a_{N-1}
        # We need N+1 terms to get the coefficient at position N
        "$FCOMBINAT" -j "$SPEC" -n $((N + 1)) $LABELED_FLAG --terms --symbol "$SYMBOL" | sed -n "$((N + 1))p"
        ;;

    unrank)
        RANK="$6"
        "$FCOMBINAT" -j "$SPEC" -n "$N" $LABELED_FLAG --unrank "$RANK" --symbol "$SYMBOL"
        ;;

    rank)
        OBJECT="$6"
        "$FCOMBINAT" -j "$SPEC" -n "$N" $LABELED_FLAG --rank "$OBJECT" --symbol "$SYMBOL"
        ;;

    draw)
        "$FCOMBINAT" -j "$SPEC" -n "$N" $LABELED_FLAG --draw --symbol "$SYMBOL"
        ;;

    boltzmann)
        "$FCOMBINAT" -j "$SPEC" -n "$N" $LABELED_FLAG --boltzmann --symbol "$SYMBOL"
        ;;

    *)
        echo "Unknown action: $ACTION" >&2
        echo "Usage: $0 {count|unrank|rank|draw|boltzmann} <spec> <n> <labeled> <symbol> [rank|object]" >&2
        exit 1
        ;;
esac
