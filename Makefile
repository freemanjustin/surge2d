###################################################################
#
# freeman.justin@gmail.com 
#
##################################################################

FC=	ifort

FSRC=	./src/

FFLAGS=	-O3 -g -module $(FSRC)

INC=	

LFLAGS= -lnetcdff

FOBJ=	\
	$(FSRC)mod_strings.o \
	$(FSRC)mod_kinds.o \
	$(FSRC)mod_param.o \
	$(FSRC)mod_iounits.o \
	$(FSRC)mod_scalars.o \
	$(FSRC)mod_parallel.o \
	$(FSRC)mod_stepping.o \
	$(FSRC)mod_grid.o \
	$(FSRC)mod_ocean.o \
	$(FSRC)diag.o \
	$(FSRC)mod_ncparam.o \
	$(FSRC)exchange_2d.o \
	$(FSRC)mod_boundary.o \
	$(FSRC)zetabc.o \
	$(FSRC)mod_clima.o \
	$(FSRC)mod_forces.o \
	$(FSRC)u2dbc_im.o \
	$(FSRC)v2dbc_im.o \
	$(FSRC)ini_fields.o \
	$(FSRC)bc_2d.o \
	$(FSRC)set_vbc.o \
	$(FSRC)mod_mixing.o \
	$(FSRC)mod_netcdf.o \
	$(FSRC)obc_volcons.o \
	$(FSRC)mod_coupling.o \
	$(FSRC)mod_sources.o \
	$(FSRC)wetdry.o \
	$(FSRC)step2d.o \
	$(FSRC)main2d.o \
	$(FSRC)ocean_control.o \
	$(FSRC)master.o \
	$(FSRC)get_data.o \
	$(FSRC)analytical.o \
	$(FSRC)set_2dfld.o \
	$(FSRC)set_data.o \
	$(FSRC)ntimestep.o \
	$(FSRC)mp_routines.o \
	$(FSRC)nf_fread2d.o \
	$(FSRC)nf_fread3d.o \
	$(FSRC)get_2dfld.o \
	$(FSRC)timers.o \
	$(FSRC)get_date.o \
	$(FSRC)get_ngfld.o \
	$(FSRC)set_ngfld.o \
	$(FSRC)interpolate.o \
	$(FSRC)regrid.o \
	$(FSRC)output.o \
	$(FSRC)inquire.o \
	$(FSRC)uv_rotate.o \
	$(FSRC)nf_fwrite2d.o \
	$(FSRC)wrt_his.o \
	$(FSRC)def_var.o \
	$(FSRC)def_his.o \
	$(FSRC)def_dim.o \
	$(FSRC)strings.o \
	$(FSRC)def_info.o \
	$(FSRC)lbc.o \
	$(FSRC)wrt_info.o \
	$(FSRC)get_cycle.o \
	$(FSRC)def_rst.o \
	$(FSRC)wrt_rst.o \
	$(FSRC)get_varcoords.o \
	$(FSRC)mod_arrays.o \
	$(FSRC)nrutil.o \
	$(FSRC)ran_state.o \
	$(FSRC)inp_par.o \
	$(FSRC)read_phypar.o \
	$(FSRC)get_bounds.o \
	$(FSRC)ini_hmixcoef.o \
	$(FSRC)metrics.o \
	$(FSRC)set_masks.o \
	$(FSRC)stiffness.o \
	$(FSRC)initial.o \
	$(FSRC)get_grid.o \
	$(FSRC)get_state.o \
	$(FSRC)get_wetdry.o \
	$(FSRC)check_multifile.o \
	$(FSRC)close_io.o \
	$(FSRC)checkvars.o \
	$(FSRC)checkdefs.o \
	$(FSRC)get_idata.o \
	$(FSRC)get_nudgcoef.o \
	$(FSRC)checkadj.o
	

OBJ=	$(FOBJ) 

EXEC=	./bin/oceanSerial	

$(EXEC):$(OBJ)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJ) $(LFLAGS)

$(FOBJ) : %.o : %.f90
	$(FC) $(INC) $(FFLAGS) -c $< -o $@

clean:
	rm $(FSRC)*.mod
	rm $(FOBJ)
	rm $(EXEC)
