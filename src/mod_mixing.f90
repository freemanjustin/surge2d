      MODULE mod_mixing
!
!svn $Id: mod_mixing.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Horizontal and vertical mixing coefficients:                        !
!                                                                      !
!  Akt          Vertical mixing coefficient (m2/s) for tracers.        !
!  Akv          Vertical mixing coefficient (m2/s) for momentum.       !
!  visc2_r      Horizontal, time invariant harmonic viscosity          !
!                 coefficient (m2/s) at RHO-points.                    !
!  visc2_p      Horizontal, time invariant harmonic viscosity          !
!                 coefficient (m2/s) at PSI-points.                    !
!  visc4_r      Horizontal, time invariant harmonic viscosity          !
!                 coefficient SQRT(m4/s) at RHO-points.                !
!  visc4_p      Horizontal, time invariant harmonic viscosity          !
!                 coefficient SQRT(m4/s) at RHO-points.                !
!  visc_factor  Horizontal viscosity factor (nondimensional) for       !
!                 increasing mixing in sponge areas at RHO-points.     !
!                                                                      !
!  Variables associated with the equation of state:                    !
!                                                                      !
!  alpha        Surface thermal expansion coefficient (1/Celsius).     !
!  beta         Surface saline contraction coefficient (1/PSU).        !
!  bvf          Brunt-Vaisala frequency squared (1/s2).                !
!  neutral      Coefficient to convert "in situ" density to neutral    !
!                 surface.                                             !
!                                                                      !
!  tke          Turbulent energy squared (m2/s2) at horizontal         !
!                 at W-points.                                         !
!  gls          Turbulent energy squared times turbulent length        !
!                 scale (m3/s2) at W-points.                           !
!                                                                      !
!  Large/McWilliams/Doney interior vertical mixing variables:          !
!                                                                      !
!  alfaobeta    Ratio of thermal expansion and saline contraction      !
!                 coefficients (Celsius/PSU) used in double            !
!                 diffusion.                                           !
!                                                                      !
!  Water clarity parameters:                                           !
!                                                                      !
!  Jwtype       Water clarity (Jerlov water type classification).      !
!                                                                      !
!  Large/McWilliams/Doney oceanic boundary layer variables:            !
!                                                                      !
!  ghats        Boundary layer nonlocal transport (T units/m).         !
!  hbbl         Depth of bottom oceanic boundary layer (m).            !
!  hsbl         Depth of surface oceanic boundary layer (m).           !
!  kbbl         Index of grid level above bottom  boundary layer.      !
!  ksbl         Index of grid level below surface boundary layer.      !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_MIXING
!
!  Nonlinear model state.
!
          real(r8), pointer :: visc_factor(:,:)
          real(r8), pointer :: visc2_p(:,:)
          real(r8), pointer :: visc2_r(:,:)
        END TYPE T_MIXING
        TYPE (T_MIXING), allocatable :: MIXING(:)
      CONTAINS
      SUBROUTINE allocate_mixing (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( MIXING(Ngrids) )
!
!  Nonlinear model state.
!
      IF (LuvSponge(ng)) THEN
        allocate ( MIXING(ng) % visc_factor(LBi:UBi,LBj:UBj) )
      END IF
      allocate ( MIXING(ng) % visc2_p(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc2_r(LBi:UBi,LBj:UBj) )
      RETURN
      END SUBROUTINE allocate_mixing
      SUBROUTINE initialize_mixing (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_mixing" for all nested grids.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j
      real(r8), parameter :: IniVal = 0.0_r8
      real(r8) :: cff1, cff2, cff3, cff4
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!  Set array initialization range.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          IF (LuvSponge(ng)) THEN
            DO i=Imin,Imax
              MIXING(ng) % visc_factor(i,j) = IniVal
            END DO
          END IF
          DO i=Imin,Imax
            MIXING(ng) % visc2_p(i,j) = IniVal
            MIXING(ng) % visc2_r(i,j) = IniVal
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_mixing
      END MODULE mod_mixing
