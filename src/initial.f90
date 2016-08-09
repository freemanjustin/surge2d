      SUBROUTINE initial
!
!svn $Id: initial.F 798 2016-05-21 17:04:58Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes all model variables.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE analytical_mod
      USE ini_hmixcoef_mod, ONLY : ini_hmixcoef
      USE metrics_mod, ONLY : metrics
      USE set_masks_mod, ONLY : set_masks
      USE stiffness_mod, ONLY : stiffness
      USE wetdry_mod, ONLY : wetdry
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: First = .TRUE.
      logical :: update = .FALSE.
      integer :: Fcount
      integer :: ng, thread, tile
      integer, dimension(Ngrids) :: IniRec, Tindex
!
!=======================================================================
!   Initialize model variables.
!=======================================================================
!
!$OMP MASTER
      IF (Master) THEN
        WRITE (stdout,20) 'INITIAL: Configuring and initializing ',     &
     &                    'forward nonlinear model ...'
 20     FORMAT (/,1x,a,a,/,1x,'*******')
      END IF
!$OMP END MASTER
!
!-----------------------------------------------------------------------
!  Initialize time stepping indices and counters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        iif(ng)=1
        indx1(ng)=1
        kstp(ng)=1
        krhs(ng)=1
        knew(ng)=1
        PREDICTOR_2D_STEP(ng)=.FALSE.
!
        iic(ng)=0
        nstp(ng)=1
        nrhs(ng)=1
        nnew(ng)=1
!
        IniRec(ng)=nrrec(ng)
        Tindex(ng)=1
!
        synchro_flag(ng)=.TRUE.
        first_time(ng)=0
        tdays(ng)=dstart
        time(ng)=tdays(ng)*day2sec
!$OMP MASTER
        ntstart(ng)=INT((time(ng)-dstart*day2sec)/dt(ng))+1
        ntend(ng)=ntimes(ng)
        ntfirst(ng)=ntstart(ng)
!$OMP END MASTER
!$OMP BARRIER
        step_counter(ng)=0
        CALL time_string (time(ng), time_code(ng))
      END DO
!
!=======================================================================
!  On first pass of ensemble/perturbation/iteration loop, initialize
!  model configuration.
!=======================================================================
!
      IF (Nrun.eq.ERstr) THEN
!
!-----------------------------------------------------------------------
!  Set horizontal grid, bathymetry, and Land/Sea masking (if any).
!  Use analytical functions or read in from a grid NetCDF.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP MASTER
          CALL get_grid (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  Compute various metric term combinations.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL metrics (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
!
!-----------------------------------------------------------------------
!  If appropriate, set spatially varying nudging coefficients time
!  scales.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
          IF (Lnudging(ng)) THEN
!$OMP MASTER
            CALL get_nudgcoef (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Initialize horizontal mixing coefficients. If applicable, scale
!  mixing coefficients according to the grid size (smallest area).
!  Also increase their values in sponge areas using the "visc_factor"
!  and/or "diff_factor" read from input Grid NetCDF file.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL ini_hmixcoef (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
!
!=======================================================================
!  Initialize model state variables and forcing.  This part is
!  executed for each ensemble/perturbation/iteration run.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Set primitive variables initial conditions.
!-----------------------------------------------------------------------
!
!  Analytical initial conditions for momentum and active tracers.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).eq.0) THEN
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ana_initial (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END IF
      END DO
!
!  If restart, read in initial conditions restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
          CALL get_state (ng, 0, 1, INI(ng)%name,                       &

     &                    IniRec(ng), Tindex(ng))
!$OMP END MASTER
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
          time(ng)=io_time                   ! needed for shared-memory
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Process initial wet/dry masks.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!
!  If restart, read in wet/dry masks.
!
        IF (nrrec(ng).ne.0) THEN
!$OMP MASTER
          CALL get_wetdry (ng, iNLM, IniRec(ng))
!$OMP END MASTER
!$OMP BARRIER
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL wetdry (ng, tile, Tindex(ng), .TRUE.)
          END DO
!$OMP BARRIER
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Read in initial forcing, climatology and assimilation data from
!  input NetCDF files.  It loads the first relevant data record for
!  the time-interpolation between snapshots.
!-----------------------------------------------------------------------
!
!  If applicable, check input multi-files where the time records for
!  a particular field are split into several NetCDF files.  Initialize
!  several parameters in the file structure so the appropriate input
!  file is selected during initialization/restart.
!
      DO ng=1,Ngrids
!$OMP MASTER
        CALL check_multifile (ng, iNLM)
!$OMP END MASTER
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  If applicable, read in input data.
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iNLM)
!$OMP MASTER
        CALL get_idata (ng)
        CALL get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Set internal I/O mask arrays.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL set_masks (ng, tile, iNLM)
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Compute grid stiffness.
!-----------------------------------------------------------------------
!
      IF (Lstiffness) THEN
        Lstiffness=.FALSE.
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL stiffness (ng, tile, iNLM)
          END DO
!$OMP BARRIER
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time-stepping counter and clock.
!-----------------------------------------------------------------------
!
!  Subsract one time unit to avoid special case due to initialization
!  in the main time-stepping routine.
!
      DO ng=1,Ngrids
        iic(ng)=ntstart(ng)-1
        time(ng)=time(ng)-dt(ng)
      END DO
      RETURN
      END SUBROUTINE initial
