# Compile options
#	make cap d=WB			  # write a binary file with search data
#	make cap f=omp			  # run in parallel mode (openMP)
#	make cap d=WB f=omp		  # run in parallel and write a binary datafile 
#	make clean cap d=WB f=omp # Clean up and compile with all flags
# 
# UPDATE 2018-11-08. Ubuntu: CAP doesn't compile unless using flag `no-pie` 
#					 gcc version: (Ubuntu 7.3.0-27ubuntu1~18.04) 7.3.0.
# UPDATE 2022-04-25. Centos 7: Compiles fine without -no-pie
#
#FFLAGS = -O -no-pie	# enable for Ubuntu.
FC = gfortran
CFLAGS = ${FFLAGS}

ifeq (${d}, WB)
	CFLAGS += -D${d}
endif

ifeq (${f}, omp)
	FFLAGS += -fopenmp
	CFLAGS += -DOMP
endif

CAP  = cap cap_dir

SUBS = fft.o Complex.o radiats.o futterman.o sacio.o trap.o sub_tt2cmt.o sub_first_motion_misfit.o sub_fmp_print_params.o sub_misfit.o sub_inversion.o sub_uv2lune.o

all: $(CAP)

cap cap_dir: %:%.o $(SUBS) cap_sub.o
	$(LINK.f) -o $@ $^ -L$(SACHOME)/lib -lsac -lsacio

cap_dir.o: cap.c
	$(COMPILE.c) -DDIRECTIVITY -o $@ $<

clean:
	rm -f $(CAP) *.o
