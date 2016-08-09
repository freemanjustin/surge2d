      SUBROUTINE main2d (RunInterval)
!
!svn $Id: main2d.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine is the main driver for nonlinear ROMS/TOMS when     !
!  configurated as shallow water (barotropic) ocean model only. It     !
!  advances forward  the vertically integrated primitive equations     !
!  for all nested grids,  if any,  by the specified  time interval     !
!  (seconds), RunInterval.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
      USE diag_mod, ONLY : diag
      USE ini_fields_mod, ONLY : ini_fields, ini_zeta
      USE set_vbc_mod, ONLY: set_vbc
      USE step2d_mod, ONLY : step2d
!
      implicit none
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: RunInterval
!
!  Local variable declarations.
!
      logical :: DoNestLayer, Time_Step
      integer :: Nsteps, Rsteps
      integer :: ig, il, istep, ng, nl, tile
      integer :: next_indx1
!
!=======================================================================
!  Time-step vertically integrated equations.
!=======================================================================
!
!  Time-step the 3D kernel for the specified time interval (seconds),
!  RunInterval.
!
      Time_Step=.TRUE.
      DoNestLayer=.TRUE.
!
      KERNEL_LOOP : DO WHILE (Time_Step)
!
!  In nesting applications, the number of nesting layers (NestLayers) is
!  used to facilitate refinement grids and composite/refinament grids
!  combinations. Otherwise, the solution it is looped once for a single
!  grid application (NestLayers = 1).
!
        nl=0
!
        NEST_LAYER : DO WHILE (DoNestLayer)
!
!  Determine number of time steps to compute in each nested grid layer
!  based on the specified time interval (seconds), RunInterval. Non
!  nesting applications have NestLayers=1. Notice that RunInterval is
!  set in the calling driver. Its value may span the full period of the
!  simulation, or a multi-model coupling interval, or just a single
!  step.
!
          CALL ntimesteps (iNLM, RunInterval, nl, Nsteps, Rsteps)
          IF (exit_flag.ne.NoError) RETURN
          IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
!
!  Time-step governing equations for Nsteps.
!
          STEP_LOOP : DO istep=1,Nsteps
!
!  Set time indices and time clock.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iic(ng)=iic(ng)+1
              time(ng)=time(ng)+dt(ng)
              tdays(ng)=time(ng)*sec2day
              CALL time_string (time(ng), time_code(ng))
              IF (step_counter(ng).eq.Rsteps) Time_Step=.FALSE.
            END DO
!
!-----------------------------------------------------------------------
!  Read in required data, if any, from input NetCDF files.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
!$OMP MASTER
              CALL get_data (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF (exit_flag.ne.NoError) RETURN
            END DO
!
!-----------------------------------------------------------------------
!  If applicable, process input data: time interpolate between data
!  snapshots.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL set_data (ng, tile)
              END DO
!$OMP BARRIER
            END DO
            IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Initialize all time levels and compute other initial fields.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Initialize free-surface.
!
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL ini_zeta (ng, tile, iNLM)
                END DO
!$OMP BARRIER
!
!  Initialize other state variables.
!
                DO tile=last_tile(ng),first_tile(ng),-1
                  CALL ini_fields (ng, tile, iNLM)
                END DO
!$OMP BARRIER
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Compute and report diagnostics. If appropriate, accumulate time-
!  averaged output data which needs a irreversible loop in shared-memory
!  jobs.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1     ! irreversible
                CALL diag (ng, tile)
              END DO
!$OMP BARRIER
            END DO
            IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set vertical boundary conditions. Process tidal forcing.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL set_vbc (ng, tile)
              END DO
!$OMP BARRIER
            END DO
!
!-----------------------------------------------------------------------
!  If appropriate, write out fields into output NetCDF files.  Notice
!  that IO data is written in delayed and serial mode.  Exit if last
!  time step.
!-----------------------------------------------------------------------
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
!$OMP MASTER
              CALL output (ng)
!$OMP END MASTER
!$OMP BARRIER
              IF ((exit_flag.ne.NoError).or.                            &
     &            ((iic(ng).eq.(ntend(ng)+1)).and.(ng.eq.Ngrids))) THEN
                RETURN
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Solve the vertically integrated primitive equations for the
!  free-surface and momentum components.
!-----------------------------------------------------------------------
!
!  Set time indices for predictor step. The PREDICTOR_2D_STEP switch
!  it is assumed to be false before the first time-step.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              iif(ng)=1
              nfast(ng)=1
              next_indx1=3-indx1(ng)
              IF (.not.PREDICTOR_2D_STEP(ng)) THEN
                PREDICTOR_2D_STEP(ng)=.TRUE.
                IF (iic(ng).eq.ntfirst(ng)) THEN
                  kstp(ng)=indx1(ng)
                ELSE
                  kstp(ng)=3-indx1(ng)
                END IF
                knew(ng)=3
                krhs(ng)=indx1(ng)
              END IF
!
!  Predictor step - Advance barotropic equations using 2D time-step
!  ==============   predictor scheme.
!
              DO tile=last_tile(ng),first_tile(ng),-1
                CALL step2d (ng, tile)
              END DO
!$OMP BARRIER
            END DO
!
!  Set time indices for corrector step.
!
            DO ig=1,GridsInLayer(nl)
              ng=GridNumber(ig,nl)
              IF (PREDICTOR_2D_STEP(ng)) THEN
                PREDICTOR_2D_STEP(ng)=.FALSE.
                knew(ng)=next_indx1
                kstp(ng)=3-knew(ng)
                krhs(ng)=3
                IF (iif(ng).lt.(nfast(ng)+1)) indx1(ng)=next_indx1
              END IF
!
!  Corrector step - Apply 2D time-step corrector scheme.  Notice that
!  ==============   there is not need for a corrector step during the
!  auxiliary (nfast+1) time-step.
!
              IF (iif(ng).lt.(nfast(ng)+1)) THEN
                DO tile=first_tile(ng),last_tile(ng),+1
                  CALL step2d (ng, tile)
                END DO
!$OMP BARRIER
              END IF
            END DO
          END DO STEP_LOOP
        END DO NEST_LAYER
      END DO KERNEL_LOOP
      RETURN
      END SUBROUTINE main2d
