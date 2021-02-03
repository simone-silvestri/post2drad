
F77 = mpif90 -fdefault-real-8 -ffixed-line-length-none -mcmodel=large -O5
##F77 = ftn -132 -r8 -O2 -mcmodel=large

##CFLAGS = -g -Wall -Wno-unused-parameter -Wextra -Warray-temporaries
##CFLAGS+= -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 
##CFLAGS+= -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

DECOMP = $(HOME)/2decomp_fft
INC = -I$(DECOMP)/include
LIB = -L$(DECOMP)/lib -l2decomp_fft
RM = rm -f

PROGRAM = postOutput
OBJS    = params.o math.o write.o read.o comm.o fans.o budget.o vfft.o spectra.o post.o 

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(FLAGS) -o $(PROGRAM) $(OBJS) $(INC) $(LIB)

##decomp_2d.o : decomp_2d.f90 params.o
##	$(F77)  -DDOUBLE_PREC -DOVERWRITE -O2 -cpp  -c decomp_2d.f90
##io.o : io.f90 params.o
##	$(F77)  -DDOUBLE_PREC -DOVERWRITE -O2 -cpp  -c  io.f90
post.o: post.f90  write.o params.o fans.o comm.o
	$(F77) $(CFLAGS) -c post.f90 $(INC) $(LIB)
math.o: math.f
	$(F77) $(CFLAGS) -c math.f $(INC) $(LIB)
##pois.o: pois.f
##	$(F77) $(CFLAGS) -c pois.f $(INC) $(LIB)
vfft.o: vfft.f
	$(F77) -c vfft.f
write.o: write.f90 params.o
	$(F77) $(CFLAGS) -c write.f90 $(INC) $(LIB)
read.o: read.f90 params.o
	$(F77) $(CFLAGS) -c read.f90 $(INC) $(LIB)
fans.o:  fans.f90 params.o
	$(F77) $(CFLAGS) -c fans.f90 $(INC) $(LIB)
spectra.o: spectra.f90 params.o
	$(F77) $(CFLAGS) -c spectra.f90 $(INC) $(LIB)
budget.o: budget.f90 params.o comm.o
	$(F77) $(CFLAGS) -c budget.f90 $(INC) $(LIB)
##press.o: press.f90 params.o comm.o
##	$(F77) $(CFLAGS) -c press.f90 $(INC) $(LIB)
comm.o:  comm.f90 params.o
	$(F77) $(CFLAGS) -c comm.f90 $(INC) $(LIB)
params.o: params.f90
	$(F77) $(CFLAGS) -c params.f90 $(INC) $(LIB)


clean: 
	rm -rf *.mod *.o $(PROGRAM) *.vtk

