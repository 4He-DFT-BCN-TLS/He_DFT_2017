#
# Makefile para montar el programa DFT4He3d 
#  (Density Functional Theory 4He 3dimensional program)
#
# (Version de compilacion para Pentium IV)
#
COMP = ifort
CFLAGS = -c  -static -module ./modules
CFLAGS = -c -O3 -tpp6 -axK -xK -tune pn3 -arch pn3 -static -reentrancy threaded -threads -module ./modules
CFLAGS = -c -O0 -tpp7 -static -reentrancy threaded -threads -module ./modules
CFLAGS = -c -g -check -static -module ./modules 
CFLAGS = -c -O3 -tpp7 -unroll2048 -reentrancy threaded -threads -module ./modules
CFLAGS = -c -O3 -tpp7 -unroll8192 -align -axN -xN -arch pn4 -static -module ./modules
CFLAGS = -c -O3 -static -module ./modules
#CFLAGS = -c -fast -fpp -fopenmp -parallel -opt-matmul -unroll -module ./modules
CFLAGS = -c -fast -fpp -fopenmp -parallel -opt-matmul -mkl -unroll -module ./modules -xAVX -axAVX -mavx -arch AVX
#CFLAGS = -c -g -check all -static -module ./modules 
#LD_FLAGS = -threads -parallel -opt-matmul
LD_FLAGS = -threads -parallel -opt-matmul -mkl -openmp
CFLAGS = -c  -O3  -xHost -align array64byte -qopenmp\
	          -qopt-matmul -unroll0 -module ./modules
LD_FLAGS = -threads -I${MKLROOT}/include/fftw -mkl=parallel -qopt-matmul

#   Libraries for the FFT
LIB1=fftw3
LIB2=fftw3_threads

#   Name of the program
PROGNAME=BCN4HeDFTv_bf


#   Fortran objects
objs=init_deriv_parallel.o	V_impur.o	modules.o	FT_V_spline.o	BCN4HeDFTv_bf.o	derden.o	dimen.o\
		energy.o	evolo.o\
		evolox.o	fforma.o	fft.o		initcg.o	ironing.o	ironingc.o\
		mates.o		poten.o		printout.o	r_cm.o		back_flow3D.o\
		readen.o	respar.o	term_alfa.o	timer.o		titols.o\
		vareps.o	varmu.o		tstgrid.o	s13adf.o	den_cil_average.o\
		newder.o	redef.o		steprkr.o	steppcr.o	spls3.o\
		diag.o		instates.o	rhoasin0.o	instates_external.o
#
.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(objs)
#	$(COMP)	-o $(PROGNAME) $(objs) -l$(LIB1) -l$(LIB2)  $(LD_FLAGS)
	$(COMP)	-o $(PROGNAME) $(objs)  $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
modules:	
	$(COMP) $(CFLAGS) $(mods)

clear:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
clean:
	rm -f *.o *.bak *.lst modules/*.mod $(PROGNAME);
