#!/bin/bash
# tests/scripts/setup_test.sh ID ECS_FILE N

ID="$1"
ECS_FILE="$2"
N="$3"

mkdir -p build
OUTPUT_JSON="build/test_${ID}.json"
OUTPUT_EXPECTED="build/test_${ID}.expected"

# Generate JSON test case (subset of ecs.json)
jq -r --arg id "$ID" '.[$id]' "$ECS_FILE" > "$OUTPUT_JSON"

# Generate Expected Output
# Format: t0\nt1\n... (one term per line)
jq -r --arg id "$ID" --argjson n "$N" '
    .[$id].terms[0:$n] | .[] | tostring
' "$ECS_FILE" > "$OUTPUT_EXPECTED"