PERF_OPS_DIR := exp/perf-operations
RESULTS_DIR  := exp/performance/results

EXP_DIRS += $(PERF_OPS_DIR)

# Prevent make from trying to rebuild the templates themselves.
$(PERF_OPS_DIR)/bar_count_base.tex: ;
$(PERF_OPS_DIR)/bar_runr_base.tex: ;

.PHONY: perf-operations perf-obj

# --- count chart (fcombinat vs Maple) ---
$(PERF_OPS_DIR)/perf-obj-count.dat: $(PERF_OPS_DIR)/reshape_byobj.py
	@python3 $< count $(RESULTS_DIR) > $@ \
		|| { echo "No results for 'count' — run 'make bench' first"; rm -f $@; exit 1; }

$(PERF_OPS_DIR)/perf-obj-count.tex: $(PERF_OPS_DIR)/bar_count_base.tex $(PERF_OPS_DIR)/perf-obj-count.dat
	@cp $< $@

$(PERF_OPS_DIR)/perf-obj-count.pdf: $(PERF_OPS_DIR)/perf-obj-count.dat

# --- rank vs unrank chart (fcombinat only) ---
$(PERF_OPS_DIR)/perf-obj-rankvsunrank.dat: $(PERF_OPS_DIR)/reshape_rankvsunrank.py
	@python3 $< $(RESULTS_DIR) > $@ \
		|| { echo "No results — run 'make bench' first"; rm -f $@; exit 1; }

$(PERF_OPS_DIR)/perf-obj-rankvsunrank.tex: $(PERF_OPS_DIR)/bar_runr_base.tex $(PERF_OPS_DIR)/perf-obj-rankvsunrank.dat
	@cp $< $@

$(PERF_OPS_DIR)/perf-obj-rankvsunrank.pdf: $(PERF_OPS_DIR)/perf-obj-rankvsunrank.dat

perf-obj: $(PERF_OPS_DIR)/perf-obj-count.pdf $(PERF_OPS_DIR)/perf-obj-rankvsunrank.pdf

perf-operations: perf-obj
