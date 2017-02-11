TARGETS := multicore \
	multilore \
	multigore \
	hcore

MACHINE := $(shell uname -sm | sed 's/ /-/g')

ifeq "$(MACHINE)" "Linux-i686"
	FC := ifort
endif
ifeq "$(MACHINE)" "IRIX-IP32"
	FC := f77
endif

all: $(TARGETS)
clean:
	rm -f $(TARGETS)
