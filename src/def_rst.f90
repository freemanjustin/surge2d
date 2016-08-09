      SUBROUTINE def_rst (ng)
!
!svn $Id: def_rst.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates restart NetCDF file, it defines its            !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, got_var(NV)
      integer, parameter :: Natt = 25
      integer :: i, j, nvd3, nvd4, nvd5
      integer :: recdim, status, varid
      integer :: DimIDs(32)
      integer :: r2dgrd(4), ru2dgrd(4), rv2dgrd(4)
      integer :: sp2dgrd(3), sr2dgrd(3), su2dgrd(3), sv2dgrd(3)
      integer :: sr3dgrd(4), su3dgrd(4), sv3dgrd(4)
      integer :: t2dgrd(4), u2dgrd(4), v2dgrd(4)
      integer :: Vsize(4)
      integer :: def_dim
      real(r8) :: Aval(6)
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_rst.F'
!
!=======================================================================
!  Create a new restart NetCDF file.
!=======================================================================
!
!  Activate creation of restart NetCDF file.  Create a new restart
!  file if during a restart run, the restart filename "RST(ng)%name"
!  is different than the initial filename "INI(ng)%name".
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=RST(ng)%name
      Ldefine=.FALSE.
      IF (((nrrec(ng).eq.0).and.(iic(ng).eq.ntstart(ng))).or.           &
     &    ((nrrec(ng).ne.0).and.                                        &
     &     (TRIM(ncname).ne.TRIM(INI(ng)%name)))) THEN
        Ldefine=.TRUE.
      END IF
!
      DEFINE : IF (Ldefine) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), RST(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN
        status=def_dim(ng, iNLM, RST(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN
        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
        nvd3=3
        nvd4=4
        nvd5=5
!
!  Define dimension vectors for staggered tracer type variables.
!
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        sr2dgrd(1)=DimIDs( 1)
        sr2dgrd(2)=DimIDs( 5)
        sr2dgrd(3)=DimIDs(12)
        t2dgrd(3)=DimIDs(12)
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        sp2dgrd(1)=DimIDs( 4)
        sp2dgrd(2)=DimIDs( 8)
        sp2dgrd(3)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momentum type variables.
!
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
!
!  Define dimension vectors for staggered v-momentum type variables.
!
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
!
!  Define dimension vector for sediment, radiation stress variables.
!
        su2dgrd(1)=DimIDs( 2)
        su2dgrd(2)=DimIDs( 6)
        su2dgrd(3)=DimIDs(12)
        sv2dgrd(1)=DimIDs( 3)
        sv2dgrd(2)=DimIDs( 7)
        sv2dgrd(3)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        RST(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, RST(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on PSI-points.
!
        Vinfo( 1)=Vname(1,idPwet)
        Vinfo( 2)=Vname(2,idPwet)
        Vinfo( 3)=Vname(3,idPwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idPwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idPwet,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idPwet),     &
     &                 NF_FOUT, nvd3, sp2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on RHO-points.
!
        Vinfo( 1)=Vname(1,idRwet)
        Vinfo( 2)=Vname(2,idRwet)
        Vinfo( 3)=Vname(3,idRwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idRwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRwet,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRwet),     &
     &                 NF_FOUT, nvd3, sr2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on U-points.
!
        Vinfo( 1)=Vname(1,idUwet)
        Vinfo( 2)=Vname(2,idUwet)
        Vinfo( 3)=Vname(3,idUwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idUwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUwet,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUwet),     &
     &                 NF_FOUT, nvd3, su2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define wet/dry mask on V-points.
!
        Vinfo( 1)=Vname(1,idVwet)
        Vinfo( 2)=Vname(2,idVwet)
        Vinfo( 3)=Vname(3,idVwet)
        Vinfo(14)=Vname(4,idVwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVwet,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVwet),     &
     &                 NF_FOUT, nvd3, sv2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),     &
     &                 NF_FRST, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)=Vname(3,idUbar)
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),     &
     &                 NF_FRST, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)=Vname(3,idVbar)
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),     &
     &                 NF_FRST, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, RST(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, RST(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing restart file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.Ldefine) THEN
        ncname=RST(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open restart file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, RST(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  restart variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            RST(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            RST(ng)%Vid(idRwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            RST(ng)%Vid(idUwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            RST(ng)%Vid(idVwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            RST(ng)%Vid(idPwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            RST(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            RST(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            RST(ng)%Vid(idVbar)=var_id(i)
          END IF
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPwet)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idPwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  Set unlimited time record dimension to current value.
!
        IF (LcycleRST(ng)) THEN
          RST(ng)%Rindex=0
        ELSE
          RST(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
  10  FORMAT (/,' DEF_RST - unable to create restart NetCDF file: ',a)
  20  FORMAT (1pe11.4,1x,'millimeter')
  30  FORMAT (/,' DEF_RST - unable to open restart NetCDF file: ',a)
  40  FORMAT (/,' DEF_RST - unable to find variable: ',a,2x,            &
     &        ' in restart NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_rst
