# 1D Euler test cases: unit tests, numerical solver, and exact-solution driver.
FC = ifx
PRECISION ?= dp
FFLAGS ?= -O2 -g

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

# Shared interpolation standards live in the interpolation subproject.
VPATH := src/common src/euler src/wcns app tests \
         ../interpolation/src/interpolation ../interpolation/tests

SHARED_OBJS := $(BUILD)/precision_select.o $(BUILD)/scheme_diagnostics.o \
               $(BUILD)/scheme_params.o $(BUILD)/interp_linear.o \
               $(BUILD)/interp_deriv.o $(BUILD)/interp_qr.o \
               $(BUILD)/interp_rbf.o $(BUILD)/interp_traditional.o \
               $(BUILD)/interp_framework.o $(BUILD)/scheme_dispatch.o
INIT_OBJS := $(BUILD)/e1d_globals.o $(BUILD)/e1d_initial.o
CTRL_OBJS := $(INIT_OBJS) $(BUILD)/e1d_control.o $(BUILD)/difference_operators.o
FLUX_OBJS := $(INIT_OBJS) $(BUILD)/e1d_flux.o $(BUILD)/difference_operators.o
SOLVER_OBJS := $(FLUX_OBJS) $(BUILD)/e1d_solver.o
EXACT_OBJS := $(INIT_OBJS) $(BUILD)/e1d_exact.o
IO_OBJS := $(BUILD)/e1d_io.o
APP_OBJS := $(SOLVER_OBJS) $(BUILD)/e1d_control.o $(IO_OBJS)
EXACT_APP_OBJS := $(EXACT_OBJS) $(BUILD)/e1d_flux.o $(BUILD)/e1d_control.o $(IO_OBJS)
TEST_UTIL := $(BUILD)/test_util.o

all: $(BIN)/euler1d $(BIN)/euler1d_exact

$(BUILD) $(BIN):
	mkdir -p $@

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) -c $< -o $@

# Fortran module dependencies (the pattern rule above does not infer them).
# Shared interpolation modules, listed in compile order.
$(BUILD)/scheme_diagnostics.o: $(BUILD)/precision_select.o
$(BUILD)/scheme_params.o: $(BUILD)/precision_select.o
$(BUILD)/interp_linear.o: $(BUILD)/precision_select.o
$(BUILD)/interp_deriv.o: $(BUILD)/precision_select.o
$(BUILD)/interp_qr.o: $(BUILD)/precision_select.o $(BUILD)/scheme_params.o \
                      $(BUILD)/scheme_diagnostics.o
$(BUILD)/interp_rbf.o: $(BUILD)/precision_select.o $(BUILD)/interp_linear.o \
                       $(BUILD)/interp_deriv.o $(BUILD)/scheme_params.o \
                       $(BUILD)/scheme_diagnostics.o
$(BUILD)/interp_traditional.o: $(BUILD)/precision_select.o $(BUILD)/interp_linear.o \
                               $(BUILD)/scheme_params.o $(BUILD)/scheme_diagnostics.o
$(BUILD)/interp_framework.o: $(BUILD)/precision_select.o $(BUILD)/interp_linear.o \
                             $(BUILD)/interp_rbf.o $(BUILD)/interp_qr.o \
                             $(BUILD)/scheme_params.o $(BUILD)/scheme_diagnostics.o
$(BUILD)/scheme_dispatch.o: $(BUILD)/precision_select.o $(BUILD)/scheme_params.o \
                            $(BUILD)/scheme_diagnostics.o $(BUILD)/interp_linear.o \
                            $(BUILD)/interp_qr.o $(BUILD)/interp_rbf.o \
                            $(BUILD)/interp_traditional.o $(BUILD)/interp_framework.o
# Local modules.
$(BUILD)/e1d_globals.o: $(BUILD)/precision_select.o
$(BUILD)/e1d_initial.o: $(BUILD)/precision_select.o $(BUILD)/e1d_globals.o
$(BUILD)/e1d_control.o: $(BUILD)/precision_select.o $(BUILD)/e1d_initial.o \
                        $(BUILD)/scheme_dispatch.o
$(BUILD)/difference_operators.o: $(BUILD)/precision_select.o
$(BUILD)/e1d_flux.o: $(BUILD)/precision_select.o $(BUILD)/e1d_globals.o \
                     $(BUILD)/scheme_dispatch.o
$(BUILD)/e1d_solver.o: $(BUILD)/precision_select.o $(BUILD)/e1d_globals.o \
                       $(BUILD)/e1d_flux.o $(BUILD)/difference_operators.o
$(BUILD)/e1d_exact.o: $(BUILD)/precision_select.o $(BUILD)/e1d_globals.o \
                      $(BUILD)/e1d_initial.o
$(BUILD)/e1d_io.o: $(BUILD)/precision_select.o $(BUILD)/e1d_globals.o \
                   $(BUILD)/e1d_flux.o
$(BUILD)/test_util.o: $(BUILD)/precision_select.o
$(BUILD)/test_initial.o: $(INIT_OBJS) $(BUILD)/test_util.o
$(BUILD)/test_control.o: $(CTRL_OBJS) $(SHARED_OBJS) $(BUILD)/test_util.o
$(BUILD)/test_flux.o: $(FLUX_OBJS) $(SHARED_OBJS) $(BUILD)/test_util.o
$(BUILD)/test_solver.o: $(SOLVER_OBJS) $(SHARED_OBJS) $(BUILD)/test_util.o
$(BUILD)/test_exact.o: $(EXACT_OBJS) $(BUILD)/test_util.o
$(BUILD)/euler1d.o: $(APP_OBJS) $(SHARED_OBJS)
$(BUILD)/euler1d_exact.o: $(EXACT_APP_OBJS) $(SHARED_OBJS)

$(BIN)/test_initial: $(BUILD)/test_initial.o $(INIT_OBJS) $(BUILD)/test_util.o | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_control: $(BUILD)/test_control.o $(CTRL_OBJS) $(SHARED_OBJS) \
                     $(BUILD)/test_util.o | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_flux: $(BUILD)/test_flux.o $(FLUX_OBJS) $(SHARED_OBJS) \
                  $(BUILD)/test_util.o | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_solver: $(BUILD)/test_solver.o $(SOLVER_OBJS) $(SHARED_OBJS) \
                    $(BUILD)/test_util.o | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_exact: $(BUILD)/test_exact.o $(EXACT_OBJS) $(BUILD)/test_util.o | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/euler1d: $(BUILD)/euler1d.o $(APP_OBJS) $(SHARED_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/euler1d_exact: $(BUILD)/euler1d_exact.o $(EXACT_APP_OBJS) $(SHARED_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

TEST_NAMES := test_initial test_control test_flux test_solver test_exact
$(TEST_NAMES): %: $(BIN)/%

test: $(addprefix $(BIN)/,$(TEST_NAMES)) all
	@set -e; \
	tmp=$$(mktemp -d); trap 'rm -rf $$tmp' EXIT; \
	for t in $(TEST_NAMES); do \
		echo "== $$t =="; \
		( cd $$tmp && "$(CURDIR)/$(BIN)/$$t" ); \
	done; \
	sh tests/test_docs.sh; \
	PRECISION=$(PRECISION) sh tests/test_smoke.sh; \
	sh tests/test_reference.sh

smoke: all
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

results: all
	PRECISION=$(PRECISION) sh scripts/run_results.sh

run: all
	./$(BIN)/euler1d

clean:
	rm -rf build bin

.SECONDARY:

.PHONY: all test smoke results run clean
