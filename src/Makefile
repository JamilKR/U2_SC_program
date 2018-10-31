#######################################################
#
# U2_SH_program compiler
#
######################################################
bin_path=../bin
#
test_program_src = par_def.f90 test.f90 mix_inter.f90
test_program_obj = $(test_program_src:.f90=.o)
#
.SUFFIXES: 
.SUFFIXES: .o .f90 .oOMP .f
#
fc = gfortran
#
fopt = -O3
prep = -cpp
#
libraries   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
#
global_opt = $(fopt)
global_opt += $(libraries)
#
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
#
LAPACK77 = -llapack
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
#
.f90.o:
	$(FC) -c $(prep) $(global_opt) "$<"
#
test: $(test_program_obj) Makefile
	$(info)
	$(fc) -o $(bin_path)/$@_$(fc) $(test_program_obj) $(global_opt) $(LIBS)
#
clean:
	@rm -f *.o *.mod
