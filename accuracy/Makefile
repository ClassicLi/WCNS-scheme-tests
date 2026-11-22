# Accuracy tests: 1D convection and 2D Euler cases.
# Interpolation standards are shared with the interpolation subproject.
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
INTERP := ../interpolation/src/interpolation
VPATH := src/common src/wcns src/convection1d src/euler2d app tests \
         $(INTERP) ../interpolation/tests

SHARED_OBJS := $(BUILD)/precision_select.o $(BUILD)/scheme_diagnostics.o \
               $(BUILD)/scheme_params.o $(BUILD)/interp_linear.o \
               $(BUILD)/interp_deriv.o $(BUILD)/interp_qr.o \
               $(BUILD)/interp_rbf.o $(BUILD)/interp_traditional.o \
               $(BUILD)/interp_framework.o $(BUILD)/scheme_dispatch.o
ACC_OBJS := $(BUILD)/accuracy_control.o $(BUILD)/accuracy_order.o \
            $(BUILD)/difference_operators.o \
            $(BUILD)/c1d_initial.o $(BUILD)/c1d_exact.o $(BUILD)/c1d_solver.o \
            $(BUILD)/e2d_globals.o $(BUILD)/e2d_cases.o $(BUILD)/e2d_flux.o \
            $(BUILD)/e2d_solver.o $(BUILD)/e2d_io.o
TEST_UTIL := $(BUILD)/test_util.o

all: $(BIN)/accuracy_1d $(BIN)/accuracy_2d

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
# Local accuracy modules.
$(BUILD)/accuracy_control.o: $(BUILD)/precision_select.o $(BUILD)/scheme_params.o \
                             $(BUILD)/scheme_dispatch.o $(BUILD)/accuracy_order.o
$(BUILD)/accuracy_order.o: $(BUILD)/precision_select.o
$(BUILD)/difference_operators.o: $(BUILD)/precision_select.o
$(BUILD)/c1d_initial.o: $(BUILD)/precision_select.o
$(BUILD)/c1d_exact.o: $(BUILD)/precision_select.o $(BUILD)/c1d_initial.o
$(BUILD)/c1d_solver.o: $(BUILD)/precision_select.o $(BUILD)/scheme_dispatch.o \
                       $(BUILD)/difference_operators.o $(BUILD)/c1d_initial.o
$(BUILD)/e2d_globals.o: $(BUILD)/precision_select.o
$(BUILD)/e2d_cases.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o
$(BUILD)/e2d_flux.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                     $(BUILD)/scheme_dispatch.o
$(BUILD)/e2d_solver.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                       $(BUILD)/e2d_flux.o $(BUILD)/difference_operators.o
$(BUILD)/e2d_io.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o
$(BUILD)/test_util.o: $(BUILD)/precision_select.o

$(BUILD)/accuracy_1d.o: $(ACC_OBJS) $(SHARED_OBJS)
$(BUILD)/accuracy_2d.o: $(ACC_OBJS) $(SHARED_OBJS)

$(BIN)/accuracy_1d: $(BUILD)/accuracy_1d.o $(ACC_OBJS) $(SHARED_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/accuracy_2d: $(BUILD)/accuracy_2d.o $(ACC_OBJS) $(SHARED_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_%: tests/test_%.f90 $(ACC_OBJS) $(SHARED_OBJS) $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) \
		$(ACC_OBJS) $(SHARED_OBJS) $(TEST_UTIL) $< -o $@

TEST_NAMES := test_control test_order test_1d_core test_cases_2d test_2d_core
$(TEST_NAMES): %: $(BIN)/%

test: $(addprefix $(BIN)/,$(TEST_NAMES)) all
	@set -e; \
	tmp=$$(mktemp -d); trap 'rm -rf $$tmp' EXIT; \
	for t in $(TEST_NAMES); do \
		echo "== $$t =="; \
		( cd $$tmp && "$(CURDIR)/$(BIN)/$$t" ); \
	done
	sh tests/test_smoke.sh
	sh tests/test_docs.sh

docs:
	sh tests/test_docs.sh

run: all
	./$(BIN)/accuracy_1d

smoke: all
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

results: all
	PRECISION=$(PRECISION) sh scripts/run_results.sh

clean:
	rm -rf build bin

.SECONDARY:

.PHONY: all test docs run smoke results clean
