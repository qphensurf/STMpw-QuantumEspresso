
# Only line you need to update to point to the
# main directory of Quantum Espresso (QE)
# you might have to update the libraries called here
# if there are new libraries in QE

TOPDIR       = /path/to/qe-7.0/src

MPIF90        = mpiifort
#FC           = ifort
FC            = mpiifort
F77           = ifort

FLAGS =  -O 

MOD_FLAG      = -I
#BASEMOD_FLAGS= $(MOD_FLAG)$(TOPDIR)/iotk/src \
#               $(MOD_FLAG)$(TOPDIR)/upflib \
#               $(MOD_FLAG)$(TOPDIR)/Modules \
#               $(MOD_FLAG)$(TOPDIR)/FFTXlib \
#               $(MOD_FLAG)$(TOPDIR)/LAXlib \
#               $(MOD_FLAG)$(TOPDIR)/UtilXlib \
#               $(MOD_FLAG)$(TOPDIR)/PW/src/ \
#               $(MOD_FLAG)$(TOPDIR)/FoX/finclude

BASEMOD_FLAGS= $(MOD_FLAG)$(TOPDIR)/upflib \
               $(MOD_FLAG)$(TOPDIR)/XClib \
               $(MOD_FLAG)$(TOPDIR)/Modules \
               $(MOD_FLAG)$(TOPDIR)/FFTXlib \
	       $(MOD_FLAG)$(TOPDIR)/LAXlib \
	       $(MOD_FLAG)$(TOPDIR)/UtilXlib \
	       $(MOD_FLAG)$(TOPDIR)/MBD \
	       $(MOD_FLAG)$(TOPDIR)/KS_Solvers \
	       $(MOD_FLAG)$(TOPDIR)/FoX/finclude \
	       $(MOD_FLAG)$(TOPDIR)/dft-d3/ \
               $(MOD_FLAG)$(TOPDIR)/PW/src/ \
	       $(MOD_FLAG)$(TOPDIR)/external/devxlib/src 



PWOBJS = $(TOPDIR)/PW/src/libpw.a $(TOPDIR)/KS_Solvers/libks_solvers.a 
QEMODS = $(TOPDIR)/Modules/libqemod.a\
	 $(TOPDIR)/KS_Solvers/libks_solvers.a \
         $(TOPDIR)/FFTXlib/libqefft.a $(TOPDIR)/LAXlib/libqela.a $(TOPDIR)/UtilXlib/libutil.a  \
	 $(TOPDIR)/MBD/libmbd.a  $(TOPDIR)/XClib/xc_lib.a $(TOPDIR)/external/devxlib/src/libdevXlib.a
OBJS_GPU= \
$(TOPDIR)/upflib/dylmr2_gpu.o \
$(TOPDIR)/upflib/init_us_2_base_gpu.o \
$(TOPDIR)/upflib/interp_atwfc_gpu.o \
$(TOPDIR)/upflib/qvan2_gpu.o \
$(TOPDIR)/upflib/simpsn_gpu.o \
$(TOPDIR)/upflib/sph_bes_gpu.o \
$(TOPDIR)/upflib/ylmr2_gpu.o
OBJS_DEP= \
$(TOPDIR)/upflib/init_us_0.o \
$(TOPDIR)/upflib/init_us_b0.o \
$(TOPDIR)/upflib/init_us_1.o \
$(TOPDIR)/upflib/init_us_2_base.o \
$(TOPDIR)/upflib/init_tab_atwfc.o \
$(TOPDIR)/upflib/init_tab_beta.o \
$(TOPDIR)/upflib/init_tab_qrad.o
OBJS_NODEP= \
$(TOPDIR)/upflib/atom.o \
$(TOPDIR)/upflib/atomic_number.o \
$(TOPDIR)/upflib/dqvan2.o \
$(TOPDIR)/upflib/dylmr2.o \
$(TOPDIR)/upflib/gth.o \
$(TOPDIR)/upflib/interp_atwfc.o \
$(TOPDIR)/upflib/paw_variables.o \
$(TOPDIR)/upflib/pseudo_types.o \
$(TOPDIR)/upflib/qvan2.o \
$(TOPDIR)/upflib/radial_grids.o \
$(TOPDIR)/upflib/read_cpmd.o \
$(TOPDIR)/upflib/read_fhi.o \
$(TOPDIR)/upflib/read_ncpp.o \
$(TOPDIR)/upflib/read_ps.o \
$(TOPDIR)/upflib/read_upf_new.o \
$(TOPDIR)/upflib/read_upf_v1.o \
$(TOPDIR)/upflib/read_uspp.o \
$(TOPDIR)/upflib/spinor.o \
$(TOPDIR)/upflib/sph_ind.o \
$(TOPDIR)/upflib/sph_bes.o \
$(TOPDIR)/upflib/splinelib.o \
$(TOPDIR)/upflib/simpsn.o \
$(TOPDIR)/upflib/upf_auxtools.o \
$(TOPDIR)/upflib/upf_const.o \
$(TOPDIR)/upflib/upf_error.o \
$(TOPDIR)/upflib/upf_invmat.o \
$(TOPDIR)/upflib/upf_io.o \
$(TOPDIR)/upflib/upf_ions.o \
$(TOPDIR)/upflib/upf_kinds.o \
$(TOPDIR)/upflib/upf_params.o \
$(TOPDIR)/upflib/upf_parallel_include.o \
$(TOPDIR)/upflib/upf_spinorb.o \
$(TOPDIR)/upflib/upf_to_internal.o \
$(TOPDIR)/upflib/upf_utils.o \
$(TOPDIR)/upflib/uspp.o \
$(TOPDIR)/upflib/uspp_data.o \
$(TOPDIR)/upflib/write_upf_new.o \
$(TOPDIR)/upflib/xmltools.o \
$(TOPDIR)/upflib/ylmr2.o

MODULES = $(PWOBJS) $(QEMODS) $(OBJS_NODEP) $(OBJS_GPU) $(OBJS_DEP)

#LIBOBJS        = $(TOPDIR)/clib/clib.a  $(TOPDIR)/iotk/src/libiotk.a

LIBOBJS        =
#BLAS_LIBS      = -L$(MKLROOT)/lib/intel64/  -lmkl_intel_lp64  -lmkl_intel_thread -lmkl_core
BLAS_LIBS      =   -L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64  -lmkl_sequential -lmkl_core
LAPACK_LIBS    =  $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a
SCALAPACK_LIBS = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

FFT_LIBS       =  -L/opt/fftw3/fftw-3.3.10/lib -lfftw3
FOX_LIB  = -L$(TOPDIR)/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common\
            -lFoX_utils -lFoX_fsys
#BEEF_LIBS      = $(TOPDIR)/LIBBEEF/libbeef.a
BEEF_LIBS    = 


QELIBS         = -mkl $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FOX_LIB) $(FFT_LIBS) $(BLAS_LIBS) $(BEEF_LIBS) $(LD_LIBS)



OBJ =   declaration.o declfft.o  \
	volumen.o Fourier.o currentBardeen.o  \
	STMpw.o 

STMpw.out: $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(TOPDIR)/PP/src/libpp.a $(MODULES) $(LIBOBJS) $(QELIBS) -o STMpw.out

utils: Imagen_Bardeen_gnu.out Cond_gnu.out

all: utils STMpw.out

Imagen_Bardeen_gnu.out: 
	$(FC) $(FLAGS) Utils/Imagen_Bardeen_gnu.f90 -o Utils/Imagen_Bardeen_gnu.out

Cond_gnu.out:
	$(FC) $(FLAGS) Utils/Cond_gnu.f90 -o Utils/Cond_gnu.out

Fourier.o: Fourier.f90
	$(FC) Fourier.f90 $(FLAGS) -c

declaration.o: declaration.f90
	$(FC) declaration.f90 $(FLAGS) -c

declfft.o: declfft.f90
	$(FC) declfft.f90 $(FLAGS) -c

determinequantities.o: determinequantities.f90
	$(FC) determinequantities.f90 $(FLAGS) -c

currentBardeen.o: currentBardeen.f90
	$(FC) currentBardeen.f90 $(FLAGS) -c

volumen.o: volumen.f90
	$(FC) volumen.f90 $(FLAGS) -c

STMpw.o: STMpw.f90
	$(FC) STMpw.f90 $(BASEMOD_FLAGS) $(FLAGS) -c	

mappingscar_gen.o: mappingscar_gen.f90
	$(FC) mappingscar_gen.f90 $(FLAGS) -c	

norm_name.o: norm_name.f90
	$(FC) norm_name.f90 $(FLAGS) -c	

clean: 
	@echo "Cleaning objects..."
	rm -f *.o *.mod

cleanall: 
	@echo "Cleaning all..."
	rm -f *.o *.mod *.out */*.out
