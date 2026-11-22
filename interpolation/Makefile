FC = ifx
FFLAGS ?= -O2 -g
PRECISION ?= qp

CPPFLAGS := -cpp
ifeq ($(PRECISION),dp)
  PPFLAGS := -DPRECISION_DP
else
  PPFLAGS :=
endif

BUILD := build/$(PRECISION)
BIN := bin/$(PRECISION)

ifeq ($(findstring gfortran,$(FC)),gfortran)
  MODFLAG := -J $(BUILD)
else
  MODFLAG := -module $(BUILD)
endif

SRC_DIRS := src/interpolation src/test app
vpath %.f90 $(SRC_DIRS) tests

OBJS := $(BUILD)/precision_select.o $(BUILD)/scheme_diagnostics.o $(BUILD)/scheme_params.o $(BUILD)/interp_linear.o $(BUILD)/interp_deriv.o $(BUILD)/interp_qr.o $(BUILD)/interp_rbf.o $(BUILD)/interp_traditional.o $(BUILD)/interp_framework.o $(BUILD)/scheme_dispatch.o $(BUILD)/control.o $(BUILD)/standard_function.o $(BUILD)/diagnostics.o $(BUILD)/mesh.o $(BUILD)/solver.o $(BUILD)/io.o

TEST_UTIL := $(BUILD)/test_util.o
TEST_SRCS := $(filter-out tests/test_util.f90,$(wildcard tests/test_*.f90))
TEST_NAMES := $(notdir $(basename $(TEST_SRCS)))
TEST_BINS := $(addprefix $(BIN)/,$(TEST_NAMES))

all: $(BIN)/interpolation_test

$(BUILD) $(BIN):
	mkdir -p $@

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/main.o: $(OBJS)

$(BIN)/interpolation_test: $(BUILD)/main.o $(OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BUILD)/test_util.o: tests/test_util.f90 $(OBJS) | $(BUILD)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) -c $< -o $@

$(BIN)/test_%: tests/test_%.f90 $(OBJS) $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) $(OBJS) $(TEST_UTIL) $< -o $@

$(BIN)/invalid_choice_driver: tests/invalid_choice_driver.f90 $(OBJS) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) $(OBJS) $< -o $@

test: $(TEST_BINS)
	@set -e; \
	tmp=$$(mktemp -d); \
	mkdir -p $$tmp/results/convergence; \
	trap 'rm -rf $$tmp' EXIT; \
	for t in $(TEST_NAMES); do \
		echo "== $$t =="; \
		( cd $$tmp && "$(CURDIR)/$(BIN)/$$t" ); \
	done

run: all
	./$(BIN)/interpolation_test

smoke: all
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

results: all
	PRECISION=$(PRECISION) sh scripts/run_results.sh

clean:
	rm -rf build bin

.PHONY: all test clean run smoke results
