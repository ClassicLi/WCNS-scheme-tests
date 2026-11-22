.PHONY: all interpolation adr complex_wave accuracy euler1d euler2d test clean all-precisions

all: interpolation adr complex_wave accuracy euler1d euler2d

interpolation adr complex_wave accuracy euler1d euler2d:
	$(MAKE) -C $@

test:
	$(MAKE) -C interpolation test
	$(MAKE) -C adr test
	$(MAKE) -C complex_wave test
	$(MAKE) -C accuracy test
	$(MAKE) -C euler1d test
	$(MAKE) -C euler2d test

all-precisions:
	$(MAKE) -C interpolation PRECISION=qp
	$(MAKE) -C interpolation PRECISION=dp
	$(MAKE) -C adr PRECISION=qp
	$(MAKE) -C adr PRECISION=dp
	$(MAKE) -C complex_wave PRECISION=qp
	$(MAKE) -C complex_wave PRECISION=dp
	$(MAKE) -C accuracy PRECISION=qp
	$(MAKE) -C accuracy PRECISION=dp
	$(MAKE) -C euler1d PRECISION=qp
	$(MAKE) -C euler1d PRECISION=dp
	$(MAKE) -C euler2d PRECISION=qp
	$(MAKE) -C euler2d PRECISION=dp

clean:
	$(MAKE) -C interpolation clean
	$(MAKE) -C adr clean
	$(MAKE) -C complex_wave clean
	$(MAKE) -C accuracy clean
	$(MAKE) -C euler1d clean
	$(MAKE) -C euler2d clean
