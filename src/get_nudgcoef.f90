      SUBROUTINE get_nudgcoef (ng, model)
!
!svn $Id: get_nudgcoef.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine reads grid time scales for nudging to climatology   !
!  fields.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_clima
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE exchange_2d_mod
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: tile, LBi, UBi, LBj, UBj
      integer :: gtype, i, ic, itrc
      integer :: nvatt, nvdim, status, vindex
      integer :: Vsize(4)
      real(r8) :: Fmax, Fmin, Fscl
      character (len=40 ) :: tunits
      character (len=256) :: ncname
!
      SourceFile='get_nudgcoef.F'
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=NUD(ng)%name
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
!  Nudging coefficients for 2D momentum.
!
      IF (LnudgeM2CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM2nc),            &
     &                       NUD(ng)%Vid(idM2nc))) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idM2nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Open grid NetCDF file for reading.
!
      IF (NUD(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, NUD(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read in nudging coefficients.
!-----------------------------------------------------------------------
!
!  Set 2D arrays bounds.
!
      tile=-1
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  If appropriate, read nudging coefficients for 2D momentum  (inverse
!  time scales, 1/day).
!
      IF (LnudgeM2CLM(ng)) THEN
        gtype=r2dvar
        vindex=idM2nc
        Fscl=1.0_r8/day2sec                    ! default units: 1/day
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       MyVarName = TRIM(Vname(1,vindex)),         &
     &                       VarID = NUD(ng)%Vid(vindex),               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (exit_flag.ne.NoError) RETURN
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_r8/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_r8
            END IF
          END IF
        END DO
!
        status=nf_fread2d(ng, model, ncname, NUD(ng)%ncid,              &
     &                    Vname(1,vindex), NUD(ng)%Vid(vindex),         &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fscl, Fmin, Fmax,                             &
     &                    GRID(ng) % rmask,                             &
     &                    CLIMA(ng) % M2nudgcof)
        IF (status.ne.nf90_noerr) THEN
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            CLIMA(ng) % M2nudgcof)
        END IF
      END IF
!
  10  FORMAT (/,' GET_NUDGCOEF - unable to find nudging variable: ',a,  &
     &        /,16x,'in NetCDF file: ',a)
  20  FORMAT (/,' GET_NUDGCOEF - unable to find nudging variable: ',a,  &
     &        /,16x,            '     or generc nudging variable: ',a,  &
     &        /,16x,'in nudging NetCDF file: ',a)
  30  FORMAT (/,' GET_NUDGCOEF - unable to open nudging NetCDF file: ',a)
      RETURN
      END SUBROUTINE get_nudgcoef
