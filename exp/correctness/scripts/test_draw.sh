#!/bin/bash
# scripts/test_draw.sh - Verify draw roundtrip: draw -> rank -> unrank
# Usage: test_draw.sh JSON_FILE BINARY METHOD N_TRIALS
#   METHOD: "draw" or "boltzmann"

JSON_FILE="$1"
BINARY="$2"
METHOD="${3:-draw}"
N_TRIALS="${4:-10}"

if [ ! -f "$JSON_FILE" ]; then
    echo "Error: JSON file $JSON_FILE not found." >&2
    exit 1
fi

SPEC=$(jq -r .specification "$JSON_FILE")
LABELED=$(jq -r 'if .labeled then "--labeled" else "" end' "$JSON_FILE")
SYMBOL=$(jq -r .symbol "$JSON_FILE")

# Only run for labeled grammars (draw needs labels for roundtrip)
if [ -z "$LABELED" ]; then
    exit 0
fi

# Find an N with non-zero objects
N=6
COUNT=$("$BINARY" -j "$SPEC" $LABELED --terms --symbol "$SYMBOL" -n $((N+1)) 2>/dev/null | sed -n "$((N+1))p" | tr -d '\r')
if [ "$COUNT" = "0" ] || [ -z "$COUNT" ]; then
    N=5
    COUNT=$("$BINARY" -j "$SPEC" $LABELED --terms --symbol "$SYMBOL" -n $((N+1)) 2>/dev/null | sed -n "$((N+1))p" | tr -d '\r')
    if [ "$COUNT" = "0" ] || [ -z "$COUNT" ]; then
        exit 0 # Skip if 6 and 5 are both zero
    fi
fi

failed=0

for i in $(seq 1 $N_TRIALS); do
    # Draw a random object
    OBJ=$("$BINARY" -j "$SPEC" $LABELED --"$METHOD" --symbol "$SYMBOL" -n "$N" 2>/dev/null)
    if [ $? -ne 0 ] || [ -z "$OBJ" ]; then
        echo "  FAIL: $METHOD returned error or empty (trial $i)" >&2
        failed=1
        continue
    fi

    # Rank the drawn object
    RANK=$("$BINARY" -j "$SPEC" $LABELED --rank "$OBJ" --symbol "$SYMBOL" -n "$N" 2>/dev/null)
    if [ $? -ne 0 ] || [ -z "$RANK" ]; then
        echo "  FAIL: rank returned error for obj=$OBJ (trial $i)" >&2
        failed=1
        continue
    fi

    # Unrank back
    OBJ2=$("$BINARY" -j "$SPEC" $LABELED --unrank "$RANK" --symbol "$SYMBOL" -n "$N" 2>/dev/null)
    if [ $? -ne 0 ] || [ -z "$OBJ2" ]; then
        echo "  FAIL: unrank returned error for rank=$RANK (trial $i)" >&2
        failed=1
        continue
    fi

    # Compare
    if [ "$OBJ" != "$OBJ2" ]; then
        echo "  FAIL: roundtrip mismatch (trial $i): draw=$OBJ unrank=$OBJ2" >&2
        failed=1
    fi
done

exit $failed
