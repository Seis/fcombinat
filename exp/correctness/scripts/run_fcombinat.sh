#!/bin/bash
# scripts/run_fcombinat.sh JSON_FILE BINARY [MODE] [NTERMS]

JSON_FILE="$1"
BINARY="$2"
MODE="${3:-verify}"
NTERMS="${4:-12}"

SPEC=$(jq -r .specification "$JSON_FILE")
LABELED=$(jq -r 'if .labeled then "--labeled" else "" end' "$JSON_FILE")
SYMBOL=$(jq -r .symbol "$JSON_FILE")

if [ "$MODE" = "bijection" ]; then
    "$BINARY" -j "$SPEC" $LABELED --bijection --symbol "$SYMBOL" -n "$NTERMS" > /dev/null 2>&1
    exit $?
else
    # Verify mode
    "$BINARY" -j "$SPEC" $LABELED --symbol "$SYMBOL" -n "$NTERMS" --terms 2>&1
fi
