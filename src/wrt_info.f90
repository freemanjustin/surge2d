      SUBROUTINE wrt_info (ng, model, ncid, ncname)
!
!svn $Id: wrt_info.F 806 2016-06-10 01:01:57Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer)                            !
!     model    Calling model identifier (integer)                      !
!     ncid     NetCDF file ID (integer)                                !
!     ncname   NetCDF file name (string)                               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_sources
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .TRUE.
      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, status, varid
      integer :: ifield = 0
      real(r8) :: scale
!
      SourceFile='wrt_info.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname, ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Time stepping parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'ntimes',                &
     &                      ntimes(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndtfast',               &
     &                      ndtfast(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dt',                    &
     &                      dt(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dtfast',                &
     &                      dtfast(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'dstart',                &
     &                      dstart, (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nHIS',                  &
     &                      nHIS(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'ndefHIS',               &
     &                      ndefHIS(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_ivar (ng, model, ncname, 'nRST',                  &
     &                      nRST(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Horizontal mixing coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'nl_visc2',              &
     &                      nl_visc2(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSponge',             &
     &                      LuvSponge(ng), (/0/), (/0/),                &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Drag coefficients.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg',                  &
     &                      rdrg(ng), (/0/), (/0/),                     &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'rdrg2',                 &
     &                      rdrg2(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Nudging inverse time scales used in various tasks.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'Znudg',                 &
     &                      Znudg(ng)/sec2day, (/0/), (/0/),            &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'M2nudg',                &
     &                      M2nudg(ng)/sec2day, (/0/), (/0/),           &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_in',            &
     &                        FSobc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'FSobc_out',           &
     &                        FSobc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_in',            &
     &                        M2obc_in(ng,:), (/1/), (/4/),             &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
        CALL netcdf_put_fvar (ng, model, ncname, 'M2obc_out',           &
     &                        M2obc_out(ng,:), (/1/), (/4/),            &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Equation of State parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'rho0',                  &
     &                      rho0, (/0/), (/0/),                         &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Slipperiness parameters.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'gamma2',                &
     &                      gamma2(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LuvSrc',                &
     &                      LuvSrc(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'LwSrc',                 &
     &                      LwSrc(ng), (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Logical switches to process climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LsshCLM',               &
     &                      LsshCLM(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_lvar (ng, model, ncname, 'Lm2CLM',                &
     &                      Lm2CLM(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Logical switches for nudging climatology fields.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'LnudgeM2CLM',           &
     &                      LnudgeM2CLM(ng), (/0/), (/0/),              &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch. Writing characters in parallel I/O is extremely
!  inefficient.  It is better to write this as an integer switch:
!  0=Cartesian, 1=spherical.
!
      CALL netcdf_put_lvar (ng, model, ncname, 'spherical',             &
     &                      spherical, (/0/), (/0/),                    &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Domain Length.
!
      CALL netcdf_put_fvar (ng, model, ncname, 'xl',                    &
     &                      xl(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      CALL netcdf_put_fvar (ng, model, ncname, 'el',                    &
     &                      el(ng), (/0/), (/0/),                       &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        CALL netcdf_put_fvar (ng, model, ncname, 'user',                &
     &                        user, (/1/), (/Nuser/),                   &
     &                        ncid = ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      GRID_VARS : IF (ncid.ne.FLT(ng)%ncid) THEN
!
!  Bathymetry.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_r8
          IF (ncid.ne.STA(ng)%ncid) THEN
            IF (find_string(var_name, n_var, 'h', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % h,                          &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'h', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'h', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Coriolis parameter.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'f', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % f,                          &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'f', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'f', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'pm', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pm,                         &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'pm', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pm', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'pn', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % pn,                         &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'pn', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'pn', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'lon_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % lonr,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'lat_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % latr,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'x_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % xr,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            scale=1.0_r8
            IF (ncid.ne.STA(ng)%ncid) THEN
              IF (find_string(var_name, n_var, 'y_rho', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % rmask,                    &
     &                             GRID(ng) % yr,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_rho', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_rho', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % lonu,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % latu,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % xu,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_u', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % umask,                    &
     &                             GRID(ng) % yu,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_u', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_u', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % lonv,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'lon_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % latv,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % xv,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,10) 'x_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_v', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % yv,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_v', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_v', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lon_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % lonp,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lon_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lon_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'lat_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % latp,                     &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'lat_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'lat_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
        IF (.not.spherical) THEN
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'x_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % pmask,                    &
     &                             GRID(ng) % xp,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'x_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'x_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
          IF (exit_flag.eq.NoError) THEN
            IF (ncid.ne.STA(ng)%ncid) THEN
              scale=1.0_r8
              IF (find_string(var_name, n_var, 'y_psi', varid)) THEN
                status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,   &
     &                             LBi, UBi, LBj, UBj, scale,           &
     &                             GRID(ng) % vmask,                    &
     &                             GRID(ng) % yp,                       &
     &                             SetFillVal = .FALSE.)
                IF (status.ne.nf90_noerr) THEN
                  IF (Master) WRITE (stdout,10) 'y_p', TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                IF (Master) WRITE (stdout,20) 'y_p', TRIM(ncname)
                exit_flag=3
                ioerror=nf90_enotvar
              END IF
            END IF
          END IF
        END IF
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (exit_flag.eq.NoError) THEN
          scale=1.0_r8
          IF (ncid.ne.STA(ng)%ncid) THEN
            IF (find_string(var_name, n_var, 'angle', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % angler,                     &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'angle', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'angle', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_rho', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, r2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % rmask,                      &
     &                           GRID(ng) % rmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_rho', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_rho', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_u', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, u2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % umask,                      &
     &                           GRID(ng) % umask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_u', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_u', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_v', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, v2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % vmask,                      &
     &                           GRID(ng) % vmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_v', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_v', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
        IF (exit_flag.eq.NoError) THEN
          IF (ncid.ne.STA(ng)%ncid) THEN
            scale=1.0_r8
            IF (find_string(var_name, n_var, 'mask_psi', varid)) THEN
              status=nf_fwrite2d(ng, model, ncid, varid, 0, p2dvar,     &
     &                           LBi, UBi, LBj, UBj, scale,             &
     &                           GRID(ng) % pmask,                      &
     &                           GRID(ng) % pmask,                      &
     &                           SetFillVal = .FALSE.)
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,10) 'mask_psi', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            ELSE
              IF (Master) WRITE (stdout,20) 'mask_psi', TRIM(ncname)
              exit_flag=3
              ioerror=nf90_enotvar
            END IF
          END IF
        END IF
      END IF GRID_VARS
!
!-----------------------------------------------------------------------
!  Synchronize NetCDF file to disk to allow other processes to access
!  data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (exit_flag.ne.NoError) RETURN
!
  10  FORMAT (/,' WRT_INFO - error while writing variable: ',a,/,       &
     &        12x,'into file: ',a)
  20  FORMAT (/,' WRT_INFO - error while inquiring ID for variable: ',  &
     &        a,/,12x,'in file: ',a)
  30  FORMAT (/,' WRT_INFO - unable to synchronize to disk file: ',     &
     &        /,12x,a)
      RETURN
      END SUBROUTINE wrt_info
