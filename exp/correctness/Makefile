FCOMBINAT_BIN ?= ../../fcombinat
JOBS          ?= $(shell nproc)
NTERMS        ?= 10
ENUM_DEPTH    ?= 5

TARGET    = $(FCOMBINAT_BIN)
ECS_DATA  = ecs-data/ecs.json
TEST_IDS  := $(shell [ -f $(ECS_DATA) ] && jq -r 'keys[] | tonumber' $(ECS_DATA) | sort -n)
ERROR_LOG  = $(abspath ../../test_errors.log)

# Only phase-level targets are PHONY.
# Individual check-X / bijection-X / unique-X / draw-X / boltzmann-X targets
# are NOT listed here — doing so would shadow the pattern rules below.
.PHONY: all test test-bijection test-unique test-recursive_draw test-boltzmann generate-json clean clean-errors
TOTAL := $(words $(TEST_IDS))

all: clean-errors test test-bijection test-unique test-recursive_draw test-boltzmann clean
	@if [ -f "$(ERROR_LOG)" ]; then \
		NERR=$$(wc -l < "$(ERROR_LOG)"); \
		echo ""; \
		echo "=== correctness: $$NERR error(s) recorded in test_errors.log ==="; \
	else \
		echo ""; \
		echo "=== correctness: all $(TOTAL) tests x 5 phases passed ==="; \
	fi

generate-json:
	@$(MAKE) --no-print-directory -j$(JOBS) $(addprefix build/test_,$(addsuffix .json,$(TEST_IDS)))
	@echo "=== generated $(TOTAL) test JSON files in build/ ==="

clean-errors:
	@rm -f "$(ERROR_LOG)"

# Each phase invokes its targets in parallel.
# --no-print-directory suppresses "Entering/Leaving directory" noise.
test:
	@$(MAKE) --no-print-directory -j$(JOBS) $(addprefix check-,$(TEST_IDS))

test-bijection:
	@$(MAKE) --no-print-directory -j$(JOBS) $(addprefix bijection-,$(TEST_IDS))

test-unique:
	@$(MAKE) --no-print-directory -j$(JOBS) $(addprefix unique-,$(TEST_IDS))

test-recursive_draw:
	@$(MAKE) --no-print-directory -j$(JOBS) $(addprefix recursive_draw-,$(TEST_IDS))

test-boltzmann:
	@$(MAKE) --no-print-directory -j$(JOBS) $(addprefix boltzmann-,$(TEST_IDS))

# Verification
check-%:
	@scripts/setup_test.sh $* $(ECS_DATA) $(NTERMS)
	@scripts/run_fcombinat.sh build/test_$*.json $(TARGET) verify $(NTERMS) | cmp -s build/test_$*.expected - \
		&& echo "check $*: ok" || echo "check $*: FAIL" | tee -a "$(ERROR_LOG)"
	@rm -f build/test_$*.json build/test_$*.expected

# Bijection
bijection-%:
	@scripts/setup_test.sh $* $(ECS_DATA) $(NTERMS)
	@scripts/run_fcombinat.sh build/test_$*.json $(TARGET) bijection $(NTERMS) \
		&& echo "bijection $*: ok" || echo "bijection $*: FAIL" | tee -a "$(ERROR_LOG)"
	@rm -f build/test_$*.json build/test_$*.expected

# Uniqueness
unique-%:
	@scripts/setup_test.sh $* $(ECS_DATA) $(NTERMS)
	@scripts/test_uniqueness.sh build/test_$*.json $(TARGET) $(ENUM_DEPTH) \
		&& echo "unique $*: ok" || echo "unique $*: FAIL" | tee -a "$(ERROR_LOG)"
	@rm -f build/test_$*.json build/test_$*.expected

# Recursive draw roundtrip
recursive_draw-%:
	@scripts/setup_test.sh $* $(ECS_DATA) $(NTERMS)
	@scripts/test_draw.sh build/test_$*.json $(TARGET) draw 10 \
		&& echo "recursive_draw $*: ok" || echo "recursive_draw $*: FAIL" | tee -a "$(ERROR_LOG)"
	@rm -f build/test_$*.json build/test_$*.expected

# Boltzmann roundtrip
boltzmann-%:
	@scripts/setup_test.sh $* $(ECS_DATA) $(NTERMS)
	@scripts/test_draw.sh build/test_$*.json $(TARGET) boltzmann 5 \
		&& echo "boltzmann $*: ok" || echo "boltzmann $*: FAIL" | tee -a "$(ERROR_LOG)"
	@rm -f build/test_$*.json build/test_$*.expected

# Per-test JSON file rule (used by generate-json)
build/test_%.json:
	@scripts/setup_test.sh $* $(ECS_DATA) $(NTERMS)

clean:
	@rm -rf build/
