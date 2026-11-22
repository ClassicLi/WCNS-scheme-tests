# 2D Euler test cases: configuration and operator-pairing unit tests.
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
VPATH := src/common src/euler2d src/wcns app tests \
         ../interpolation/src/interpolation ../interpolation/tests

SHARED_OBJS := $(BUILD)/precision_select.o $(BUILD)/scheme_diagnostics.o \
               $(BUILD)/scheme_params.o $(BUILD)/interp_linear.o \
               $(BUILD)/interp_deriv.o $(BUILD)/interp_qr.o \
               $(BUILD)/interp_rbf.o $(BUILD)/interp_traditional.o \
               $(BUILD)/interp_framework.o $(BUILD)/scheme_dispatch.o
CTRL_OBJS := $(BUILD)/e2d_globals.o $(BUILD)/e2d_cases.o $(BUILD)/e2d_control.o
APP_OBJS := $(CTRL_OBJS) $(BUILD)/e2d_flux.o $(BUILD)/e2d_solver.o \
            $(BUILD)/e2d_io.o $(BUILD)/difference_operators.o $(SHARED_OBJS)
TEST_UTIL := $(BUILD)/test_util.o

all: $(BIN)/test_control $(BIN)/euler2d

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
$(BUILD)/difference_operators.o: $(BUILD)/precision_select.o
$(BUILD)/e2d_globals.o: $(BUILD)/precision_select.o
$(BUILD)/e2d_flux.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                     $(BUILD)/scheme_dispatch.o
$(BUILD)/e2d_cases.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o
$(BUILD)/e2d_solver.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                       $(BUILD)/e2d_flux.o $(BUILD)/e2d_cases.o \
                       $(BUILD)/difference_operators.o
$(BUILD)/e2d_io.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o
$(BUILD)/e2d_control.o: $(BUILD)/precision_select.o $(BUILD)/e2d_cases.o \
                        $(BUILD)/e2d_globals.o $(BUILD)/scheme_params.o \
                        $(BUILD)/scheme_dispatch.o
$(BUILD)/test_util.o: $(BUILD)/precision_select.o
$(BUILD)/test_control.o: $(CTRL_OBJS) $(SHARED_OBJS) $(TEST_UTIL)
$(BUILD)/test_operators.o: $(BUILD)/precision_select.o \
                           $(BUILD)/difference_operators.o $(TEST_UTIL)
$(BUILD)/test_flux.o: $(BUILD)/e2d_globals.o $(BUILD)/e2d_flux.o \
                      $(SHARED_OBJS) $(BUILD)/test_util.o
$(BUILD)/test_cases.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                        $(BUILD)/e2d_cases.o $(BUILD)/test_util.o
$(BUILD)/test_solver.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                        $(BUILD)/e2d_cases.o $(BUILD)/e2d_solver.o \
                        $(BUILD)/test_util.o
$(BUILD)/test_io.o: $(BUILD)/precision_select.o $(BUILD)/e2d_globals.o \
                    $(BUILD)/e2d_io.o $(BUILD)/test_util.o
$(BUILD)/euler2d.o: $(BUILD)/precision_select.o $(BUILD)/scheme_params.o \
                    $(APP_OBJS)

$(BIN)/test_control: $(BUILD)/test_control.o $(CTRL_OBJS) $(SHARED_OBJS) \
                     $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_operators: $(BUILD)/test_operators.o $(BUILD)/difference_operators.o \
                       $(BUILD)/precision_select.o $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_flux: $(BUILD)/test_flux.o $(BUILD)/e2d_globals.o $(BUILD)/e2d_flux.o \
                  $(SHARED_OBJS) $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_cases: $(BUILD)/test_cases.o $(BUILD)/e2d_globals.o $(BUILD)/e2d_cases.o \
                   $(BUILD)/precision_select.o $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_solver: $(BUILD)/test_solver.o $(BUILD)/e2d_solver.o \
                    $(BUILD)/e2d_globals.o $(BUILD)/e2d_flux.o \
                    $(BUILD)/e2d_cases.o $(BUILD)/difference_operators.o \
                    $(SHARED_OBJS) $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_io: $(BUILD)/e2d_io.o $(BUILD)/e2d_globals.o \
                $(BUILD)/precision_select.o $(BUILD)/test_io.o \
                $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/euler2d: $(BUILD)/euler2d.o $(APP_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

TEST_NAMES := test_control test_operators test_flux test_cases test_solver \
              test_io
$(TEST_NAMES): %: $(BIN)/%

test: $(addprefix $(BIN)/,$(TEST_NAMES)) all
	@set -e; \
	tmp=$$(mktemp -d); trap 'rm -rf $$tmp' EXIT; \
	for t in $(TEST_NAMES); do \
		echo "== $$t =="; \
		( cd $$tmp && "$(CURDIR)/$(BIN)/$$t" ); \
	done
	@sh tests/test_docs.sh
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

docs:
	@sh tests/test_docs.sh

smoke: all
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

clean:
	rm -rf build bin

.SECONDARY:

.PHONY: all test docs smoke clean
