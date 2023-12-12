PROG = kraken

OBJDIR = obj
MODDIR = mod
SRCDIR = src_fortran

# N.B. If you change the compiler or compiler options, you need to also
# check the mexopts.sh file to be sure the compiler here is compatible
# with the mex file compiler.

# gnu compiler
# These options worked for me using the mex/mexopts.sh.gcc options file
# gcc and gfortran should be two corresponding C and Fortran compilers
FC=gfortran

## Changs done:
# -march=v8.6-A

# For MacOS
# FFLAGS= -cpp -fPIC -I $(MODDIR) -m64 -Wall -finline-functions -ffast-math -fno-strength-reduce -fomit-frame-pointer -falign-functions=2 -O3
# FFLAGS= -cpp -fPIC -I $(MODDIR) -m64 -Wall -finline-functions -fno-strength-reduce -fomit-frame-pointer -falign-functions=2 -O3
# For Linux
FFLAGS= -cpp -fPIC -I $(MODDIR) -m64 -mfpmath=sse -Wall -finline-functions -fno-strength-reduce -fomit-frame-pointer -falign-functions=2  -O3


all : $(PROG)
	@echo " "
	@echo "Mex KRAKEN built"
	@echo " "

clean : 
	rm -f $(OBJ) $(MOD)

$(MODDIR)/%.mod : $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< 
	mv -f $*.o $(OBJDIR)/ 
	mv -f $*.mod $(MODDIR)/

$(OBJDIR)/%.o : $(MOD) $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< 
	mv -f $*.o $(OBJDIR)/

MOD =				\
	$(MODDIR)/krakmod.mod	\
	$(MODDIR)/sdrdrmod.mod	\

OBJ =				\
	$(OBJDIR)/zsecx.o	\
	$(OBJDIR)/zbrentx.o	\
	$(OBJDIR)/bcimp.o	\
	$(OBJDIR)/kuping.o	\
	$(OBJDIR)/sinvitd.o	\
	$(OBJDIR)/mergev.o	\
	$(OBJDIR)/weight.o	\
	$(OBJDIR)/twersk.o	\
	$(OBJDIR)/subtab.o	\
	$(OBJDIR)/readin.o	\
	$(OBJDIR)/errout.o	\
	$(OBJDIR)/sorti.o	\
	$(OBJDIR)/splinec.o	\
	$(OBJDIR)/kraken.o	\
	$(OBJDIR)/krakmod.o	\
	$(OBJDIR)/sdrdrmod.o	

$(PROG) : $(MOD) $(OBJ)
	@echo " "
	@echo " "
	@echo "Building Kraken"
	@echo " "
	@echo " "
# Change the `-lgfortran` flag to -lSystem for MacOS
# 	gfortran obj/zsecx.o obj/zbrentx.o obj/bcimp.o obj/kuping.o obj/sinvitd.o obj/mergev.o obj/weight.o obj/twersk.o obj/subtab.o obj/readin.o obj/errout.o obj/sorti.o obj/splinec.o obj/krakmod.o obj/sdrdrmod.o obj/kraken.o -shared -o kraken.dylib
	gfortran $(OBJ) -shared -o kraken.dylib
#	strip mkrak.mexa64
#	make clean
