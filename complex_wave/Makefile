# Complex-wave time-advancement test: WCNS interpolation plus selectable
# difference operators.  Default precision is double (dp).
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
VPATH := src/common src/wcns app tests $(INTERP) ../interpolation/tests

SHARED_OBJS := $(BUILD)/precision_select.o $(BUILD)/scheme_diagnostics.o \
               $(BUILD)/scheme_params.o $(BUILD)/interp_linear.o \
               $(BUILD)/interp_deriv.o $(BUILD)/interp_qr.o \
               $(BUILD)/interp_rbf.o $(BUILD)/interp_traditional.o \
               $(BUILD)/interp_framework.o $(BUILD)/scheme_dispatch.o
CW_OBJS := $(BUILD)/cw_control.o $(BUILD)/cw_periodic.o $(BUILD)/cw_initial.o \
           $(BUILD)/cw_time.o $(BUILD)/cw_io.o $(BUILD)/cw_params.o \
           $(BUILD)/difference_operators.o $(BUILD)/cw_wcns_rhs.o
TEST_UTIL := $(BUILD)/test_util.o

all: $(BIN)/complex_wave

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
# Local complex-wave modules.
$(BUILD)/cw_control.o: $(BUILD)/precision_select.o
$(BUILD)/cw_periodic.o: $(BUILD)/precision_select.o
$(BUILD)/cw_initial.o: $(BUILD)/precision_select.o
$(BUILD)/cw_time.o: $(BUILD)/precision_select.o $(BUILD)/cw_periodic.o
$(BUILD)/cw_io.o: $(BUILD)/precision_select.o
$(BUILD)/cw_params.o: $(BUILD)/precision_select.o
$(BUILD)/difference_operators.o: $(BUILD)/precision_select.o
$(BUILD)/cw_wcns_rhs.o: $(BUILD)/precision_select.o $(BUILD)/scheme_params.o \
                        $(BUILD)/scheme_dispatch.o $(BUILD)/cw_params.o \
                        $(BUILD)/difference_operators.o
$(BUILD)/test_util.o: $(BUILD)/precision_select.o
$(BUILD)/complex_wave.o: $(CW_OBJS) $(SHARED_OBJS)

$(BIN)/complex_wave: $(BUILD)/complex_wave.o $(CW_OBJS) $(SHARED_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_%: tests/test_%.f90 $(CW_OBJS) $(SHARED_OBJS) $(TEST_UTIL) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) \
		$(CW_OBJS) $(SHARED_OBJS) $(TEST_UTIL) $< -o $@

TEST_NAMES := test_params test_periodic test_initial test_difference_operators \
              test_rhs test_rk4
$(TEST_NAMES): %: $(BIN)/%

test: $(addprefix $(BIN)/,$(TEST_NAMES)) $(BIN)/complex_wave
	@set -e; \
	tmp=$$(mktemp -d); trap 'rm -rf $$tmp' EXIT; \
	for t in $(TEST_NAMES); do \
		echo "== $$t =="; \
		( cd $$tmp && "$(CURDIR)/$(BIN)/$$t" ); \
	done
	sh tests/test_params.sh
	sh tests/test_docs.sh

docs:
	sh tests/test_docs.sh

results: all
	PRECISION=$(PRECISION) sh scripts/run_results.sh

run: all
	./$(BIN)/complex_wave

smoke: all
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

clean:
	rm -rf build bin

.SECONDARY:

.PHONY: all test docs run smoke results clean
