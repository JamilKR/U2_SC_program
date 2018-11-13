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
ei_val_pgm_src = par_def.f90 mix_inter.f90 obs_def.f90 U2_SP_einv.f90
ei_val_pgm_obj = $(ei_val_pgm_src:.f90=.o)
#
min_sp_pgm_src= par_def.f90 mix_inter.f90 obs_def.f90 min_funct_mod.f90 minuit-cern.f U2_SP_min.f90
min_sp_pgm_obj= par_def.o mix_inter.o obs_def.o min_funct_mod.o minuit-cern.o U2_SP_min.o
#
.SUFFIXES: 
.SUFFIXES: .o .f90 .oOMP .f
#
fc = gfortran
#
fopt = -O3 -Wall
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
eigen_pgm:	$(ei_val_pgm_obj) Makefile
		$(info)
		$(fc) -o $(bin_path)/$@_$(fc) $(ei_val_pgm_obj) $(global_opt) $(LIBS)
#
min_SP_pgm:	$(min_sp_pgm_obj) Makefile
		$(info)
		$(fc) -o $(bin_path)/$@_$(fc) $(min_sp_pgm_obj) $(global_opt) $(LIBS)
#
clean:
	@rm -f *.o *.mod
