# ADR spectral test: shared harness, WENO and WCNS drivers.
FC = ifx
PRECISION ?= qp
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

INTERP := ../interpolation/src/interpolation
VPATH := src/common src/weno src/wcns app tests $(INTERP) ../interpolation/tests

CORE_OBJS := $(BUILD)/precision_select.o $(BUILD)/adr_control.o $(BUILD)/adr_periodic.o \
             $(BUILD)/adr_dft.o $(BUILD)/adr_time.o $(BUILD)/adr_io.o
COMMON_OBJS := $(CORE_OBJS) $(BUILD)/test_util.o
WENO_OBJS := $(BUILD)/weno_params.o $(BUILD)/weno_schemes.o $(BUILD)/weno_dispatch.o
WCNS_PARAM_OBJS := $(BUILD)/wcns_params.o
DIFF_OBJS := $(BUILD)/difference_operators.o
WCNS_OBJS := $(BUILD)/scheme_diagnostics.o $(BUILD)/scheme_params.o \
             $(BUILD)/interp_linear.o $(BUILD)/interp_deriv.o $(BUILD)/interp_qr.o \
             $(BUILD)/interp_rbf.o $(BUILD)/interp_traditional.o $(BUILD)/interp_framework.o \
             $(BUILD)/scheme_dispatch.o $(BUILD)/wcns_params.o \
             $(BUILD)/difference_operators.o $(BUILD)/wcns_dispatch.o

all: $(BIN)/adr_weno $(BIN)/adr_wcns

$(BUILD) $(BIN):
	mkdir -p $@

$(BUILD)/%.o: %.f90 | $(BUILD)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) -c $< -o $@

$(BUILD)/adr_weno.o: $(CORE_OBJS) $(WENO_OBJS)
$(BUILD)/adr_wcns.o: $(CORE_OBJS) $(WCNS_OBJS)

$(BIN)/adr_weno: $(BUILD)/adr_weno.o $(CORE_OBJS) $(WENO_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/adr_wcns: $(BUILD)/adr_wcns.o $(CORE_OBJS) $(WCNS_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $^ -o $@

$(BIN)/test_%: tests/test_%.f90 $(COMMON_OBJS) $(WENO_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) $(COMMON_OBJS) $(WENO_OBJS) $< -o $@

$(BIN)/test_operators: tests/test_operators.f90 $(COMMON_OBJS) $(DIFF_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) $(COMMON_OBJS) $(DIFF_OBJS) $< -o $@

$(BIN)/test_params: tests/test_params.f90 $(COMMON_OBJS) $(WENO_OBJS) $(WCNS_PARAM_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) $(COMMON_OBJS) $(WENO_OBJS) $(WCNS_PARAM_OBJS) $< -o $@

$(BIN)/test_linear_adr_wcns: tests/test_linear_adr_wcns.f90 $(COMMON_OBJS) $(WCNS_OBJS) | $(BIN)
	$(FC) $(FFLAGS) $(CPPFLAGS) $(PPFLAGS) $(MODFLAG) $(COMMON_OBJS) $(WCNS_OBJS) $< -o $@

TEST_NAMES := test_periodic_dft test_rk3 test_params test_linear_adr_weno test_weno_schemes test_operators test_linear_adr_wcns

test: $(addprefix $(BIN)/,$(TEST_NAMES))
	@set -e; \
	tmp=$$(mktemp -d); trap 'rm -rf $$tmp' EXIT; \
	for t in $(TEST_NAMES); do \
		echo "== $$t =="; \
		(cd $$tmp && "$(CURDIR)/$(BIN)/$$t"); \
	done; \
	sh tests/test_docs.sh

docs:
	sh tests/test_docs.sh

smoke: $(BIN)/adr_weno $(BIN)/adr_wcns
	PRECISION=$(PRECISION) sh tests/test_smoke.sh

clean:
	rm -rf build bin

.SECONDARY:

.PHONY: all test docs smoke clean
