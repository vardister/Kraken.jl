PROG = kraken

OBJDIR = src_fortran/obj
MODDIR = src_fortran/mod
SRCDIR = src_fortran/source
FC=gfortran

# For MacOS and Linux
FFLAGS= -cpp -fPIC -I $(MODDIR) -m64 -Wall -finline-functions -fno-strength-reduce -fomit-frame-pointer -falign-functions=2  -O3

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

ifeq ($(shell uname), Linux)
	gfortran $(OBJ) -shared -o ./src/kraken.so
else ifeq ($(shell uname), Darwin)
	gfortran $(OBJ) -shared -o ./src/kraken.dylib
endif
	make clean