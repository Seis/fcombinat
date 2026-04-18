#!/bin/bash
##############################################################################
# backends/maple/run.sh — Maple CAS backend for benchmark suite
#
# Uses the combstruct package (grammar-based) for a fair comparison
# against fcombinat's grammar-based approach.
#
# Standard interface:
#   ./run.sh count  <spec> <n> <labeled> <symbol>
#   ./run.sh unrank <spec> <n> <labeled> <symbol> <rank>
#   ./run.sh rank   <spec> <n> <labeled> <symbol> <object>
#
# combstruct API:
#   count([Symbol, spec, labeled/unlabeled], size=n)
#   draw([Symbol, spec, labeled/unlabeled], size=n) → random object
#   allstructs([Symbol, spec, labeled/unlabeled], size=n) → all objects
#
# NOTE: combstruct does not natively support rank/unrank.
#       For small N, we enumerate via allstructs and index by position.
#       For large N this will naturally time out, which is the fair
#       comparison point: fcombinat supports rank/unrank natively while
#       Maple's grammar-based combstruct does not.
#
# Prints result to stdout. Exit 0 = success.
##############################################################################

set -euo pipefail
export LC_NUMERIC=C

# Maple binary path
MAPLE="${MAPLE_BIN:-}"
if [ -z "$MAPLE" ] || [ ! -x "$MAPLE" ]; then
    echo "Error: maple binary not found. Set MAPLE_BIN in config.mk." >&2
    exit 1
fi

ACTION="$1"
SPEC="$2"
N="$3"
LABELED="$4"
SYMBOL="$5"

# Determine labeled/unlabeled keyword for combstruct
if [ "$LABELED" = "1" ]; then
    LABEL_TYPE="labeled"
else
    LABEL_TYPE="unlabeled"
fi

# Generate Maple script and execute
TMPSCRIPT=$(mktemp /tmp/maple_bench_XXXXXX.mpl)
trap "rm -f $TMPSCRIPT" EXIT

case "$ACTION" in
    # ===================================================================
    # COUNT — number of objects of size N via combstruct grammar
    # This is the primary fair comparison: both fcombinat and Maple
    # use grammar-based counting via generating functions.
    # ===================================================================
    count)
        cat > "$TMPSCRIPT" <<MAPLE_EOF
with(combstruct):
spec := {$SPEC, Z = Atom}:
printf("%a", count([$SYMBOL, spec, $LABEL_TYPE], size=$N)):
quit:
MAPLE_EOF
        ;;

    # ===================================================================
    # UNRANK — return the object at a given rank
    # combstruct has no native unrank. We use allstructs to enumerate
    # all objects and index by rank. This is O(n!) for labeled structures
    # and will naturally time out for large N — which is precisely
    # the point of comparisons: fcombinat handles this natively.
    # ===================================================================
    unrank)
        RANK="$6"
        cat > "$TMPSCRIPT" <<MAPLE_EOF
with(combstruct):
spec := {$SPEC, Z = Atom}:
all_objs := allstructs([$SYMBOL, spec, $LABEL_TYPE], size=$N):
rk := $RANK:
if rk + 1 > nops(all_objs) then
    printf("ERROR"):
    quit:
fi:
printf("%a", all_objs[rk + 1]):
quit:
MAPLE_EOF
        ;;

    # ===================================================================
    # RANK — return the integer rank of a given object
    # Same limitation as unrank: enumerates all objects via allstructs
    # and searches linearly. O(n!) time and memory.
    # ===================================================================
    rank)
        OBJECT="$6"
        cat > "$TMPSCRIPT" <<MAPLE_EOF
with(combstruct):
spec := {$SPEC, Z = Atom}:
all_objs := allstructs([$SYMBOL, spec, $LABEL_TYPE], size=$N):
target := $OBJECT:
found := -1:
for i from 1 to nops(all_objs) do
    if all_objs[i] = target then
        found := i - 1:
        break:
    fi:
od:
printf("%d", found):
quit:
MAPLE_EOF
        ;;

    # ===================================================================
    # DRAW — uniformly random object of size N via combstruct draw()
    # ===================================================================
    draw)
        cat > "$TMPSCRIPT" <<MAPLE_EOF
with(combstruct):
spec := {$SPEC, Z = Atom}:
printf("%a", draw([$SYMBOL, spec, $LABEL_TYPE], size=$N)):
quit:
MAPLE_EOF
        ;;

    *)
        echo "Unknown action: $ACTION" >&2
        echo "Usage: $0 {count|unrank|rank|draw} <spec> <n> <labeled> <symbol> [rank|object]" >&2
        exit 1
        ;;
esac

# Execute the generated Maple script
"$MAPLE" -q "$TMPSCRIPT"
