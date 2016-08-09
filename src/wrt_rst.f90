      SUBROUTINE wrt_rst (ng)
!
!svn $Id: wrt_rst.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes fields into restart NetCDF file.                !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, i, itrc, status, varid
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile='wrt_rst.F'
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      RST(ng)%Rindex=RST(ng)%Rindex+1
      Fcount=RST(ng)%Fcount
      RST(ng)%Nrec(Fcount)=RST(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (LcycleRST(ng)) THEN
        RST(ng)%Rindex=MOD(RST(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, RST(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/RST(ng)%Rindex/), (/1/),                  &
     &                      ncid = RST(ng)%ncid,                        &
     &                      varid = RST(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_r8
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idPwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % pmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idPwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % vmask,                              &
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out free-surface (m).
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % rmask,                              &
     &                   OCEAN(ng) % zeta(:,:,knew(ng)),                &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % umask_full,                         &
     &                   OCEAN(ng) % ubar(:,:,knew(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % vmask_full,                         &
     &                   OCEAN(ng) % vbar(:,:,knew(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, RST(ng)%name, RST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) knew(ng), RST(ng)%Rindex
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
      RETURN
      END SUBROUTINE wrt_rst
