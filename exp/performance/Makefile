# Makefile for benchmark suite
#
# Usage:
#   make all           — build + run suite + plot
#   make build         — build fcombinat
#   make run           — run full benchmark suite
#   make plot          — generate charts from results/
#   make run-fcombinat — run only fcombinat backend
#   make run-maple     — run only Maple backend
#   make run-boltzmann — run only boltzmann operation
#   make run-draw      — run only draw operation
#   make clean         — remove results
# FCOMBINAT_BIN is passed in from the root Makefile (defined in config.mk)

.PHONY: all build run plot clean run-fcombinat run-maple run-boltzmann run-draw

all: build run plot

build:
	@$(MAKE) -C $(patsubst %/,%,$(dir $(FCOMBINAT_BIN))) -j$$(nproc)

run: $(FCOMBINAT_BIN)
	@chmod +x run_suite.sh backends/*/run.sh
	@if [ -n "$(BACKEND)" ]; then \
		FCOMBINAT_BIN="$(FCOMBINAT_BIN)" MAPLE_BIN="$(MAPLE_BIN)" ./run_suite.sh --backend $(BACKEND); \
	else \
		FCOMBINAT_BIN="$(FCOMBINAT_BIN)" MAPLE_BIN="$(MAPLE_BIN)" ./run_suite.sh; \
	fi

run-boltzmann: $(FCOMBINAT_BIN)
	@chmod +x run_suite.sh backends/*/run.sh
	@if [ -n "$(BACKEND)" ]; then \
		FCOMBINAT_BIN="$(FCOMBINAT_BIN)" MAPLE_BIN="$(MAPLE_BIN)" ./run_suite.sh --backend $(BACKEND) --operation boltzmann; \
	else \
		FCOMBINAT_BIN="$(FCOMBINAT_BIN)" MAPLE_BIN="$(MAPLE_BIN)" ./run_suite.sh --operation boltzmann; \
	fi

run-draw: $(FCOMBINAT_BIN)
	@chmod +x run_suite.sh backends/*/run.sh
	@if [ -n "$(BACKEND)" ]; then \
		FCOMBINAT_BIN="$(FCOMBINAT_BIN)" MAPLE_BIN="$(MAPLE_BIN)" ./run_suite.sh --backend $(BACKEND) --operation draw; \
	else \
		FCOMBINAT_BIN="$(FCOMBINAT_BIN)" MAPLE_BIN="$(MAPLE_BIN)" ./run_suite.sh --operation draw; \
	fi

run-fcombinat: $(FCOMBINAT_BIN)
	@chmod +x run_suite.sh backends/fcombinat/run.sh
	@./run_suite.sh --backend fcombinat

run-maple:
	@chmod +x run_suite.sh backends/maple/run.sh
	@./run_suite.sh --backend maple

plot:
	@python3 plot_results.py

clean:
	@rm -rf results/
