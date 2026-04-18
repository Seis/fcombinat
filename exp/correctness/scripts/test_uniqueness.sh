#!/bin/bash
# tests/scripts/test_uniqueness.sh JSON_FILE BINARY N

JSON_FILE="$1"
BINARY="$2"
N="$3"

if [ ! -f "$JSON_FILE" ]; then
    echo "Error: JSON file $JSON_FILE not found."
    exit 1
fi

SPEC=$(jq -r .specification "$JSON_FILE")
LABELED=$(jq -r 'if .labeled then "--labeled" else "" end' "$JSON_FILE")
SYMBOL=$(jq -r .symbol "$JSON_FILE")

# Run enum-all and check for duplicates
# We use sort | uniq -d which outputs only duplicate lines.
DUPLICATES=$("$BINARY" -g "$SPEC" $LABELED --symbol "$SYMBOL" -n "$N" --enum-all | sort | uniq -d)

if [ -n "$DUPLICATES" ]; then
    echo "FAIL: Duplicates found for $JSON_FILE (N=$N):"
    echo "$DUPLICATES"
    exit 1
else
    exit 0
fi
