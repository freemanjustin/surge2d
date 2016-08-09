      MODULE mod_boundary
!
!svn $Id: mod_boundary.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Open boundary conditions arrays:                                    !
!                                                                      !
!  zeta_west      Free-surface (m) western boundary conditions.        !
!  zeta_east      Free-surface (m) eastern boundary conditions.        !
!  zeta_south     Free-surface (m) southern boundary conditions.       !
!  zeta_north     Free-surface (m) northern boundary conditions.       !
!  ubar_west      2D u-momentum (m/s) western boundary conditions.     !
!  vbar_west      2D v-momentum (m/s) western boundary conditions.     !
!  ubar_east      2D u-momentum (m/s) eastern boundary conditions.     !
!  vbar_east      2D v-momentum (m/s) eastern boundary conditions.     !
!  ubar_south     2D u-momentum (m/s) southern boundary conditions.    !
!  vbar_south     2D v-momentum (m/s) southern boundary conditions.    !
!  ubar_north     2D u-momentum (m/s) northern boundary conditions.    !
!  vbar_north     2D v-momentum (m/s) northern boundary conditions.    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
!
!-----------------------------------------------------------------------
!  Lateral boundary condition apply switches.
!-----------------------------------------------------------------------
!
!  The following switches are used to control which grid points are
!  processed by the lateral boundary conditions. These switches are
!  set to TRUE by default.  However in composite grids, the points
!  processed by nesting are set to FALSE to allow mixed boundary
!  conditions along the grid edges.
!
        TYPE T_APPLY
          logical, pointer :: west(:)
          logical, pointer :: east(:)
          logical, pointer :: south(:)
          logical, pointer :: north(:)
        END TYPE
        TYPE (T_APPLY), allocatable :: LBC_apply(:)
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions structure.
!-----------------------------------------------------------------------
!
        TYPE T_BOUNDARY
!
!  Nonlinear model state.
!
          real(r8), pointer :: zeta_west(:)
          real(r8), pointer :: zeta_east(:)
          real(r8), pointer :: zeta_south(:)
          real(r8), pointer :: zeta_north(:)
          real(r8), pointer :: ubar_west(:)
          real(r8), pointer :: vbar_west(:)
          real(r8), pointer :: ubar_east(:)
          real(r8), pointer :: vbar_east(:)
          real(r8), pointer :: ubar_south(:)
          real(r8), pointer :: vbar_south(:)
          real(r8), pointer :: ubar_north(:)
          real(r8), pointer :: vbar_north(:)
        END TYPE T_BOUNDARY
        TYPE (T_BOUNDARY), allocatable ::BOUNDARY(:)
      CONTAINS
      SUBROUTINE allocate_boundary (ng)
!
!=======================================================================
!                                                                      !
!  This routine initializes all variables in the module for all nested !
!  grids.  Currently, there is not parallel tiling in boundary arrays. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: my_tile
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Se dimension ranges. Notice that the boundary arrays are dimensioned
!  with the global dimensions of grid. That is, no tiling ranges in
!  distributed-memory. This is done to facilitate processing.
!
      my_tile=-1                           ! for global values
      LBi=BOUNDS(ng)%LBi(my_tile)
      UBi=BOUNDS(ng)%UBi(my_tile)
      LBj=BOUNDS(ng)%LBj(my_tile)
      UBj=BOUNDS(ng)%UBj(my_tile)
!
!  Allocate structures.
!
      IF (ng.eq.1) THEN
        allocate ( LBC_apply(Ngrids) )
        allocate ( BOUNDARY(Ngrids) )
      END IF
!
!  Lateral boundary conditions apply switches.  These switches need to
!  be initilized to TRUE here because 'initialize_boundary' is called
!  several times in adjoint-based application to clear state arrays.
!  These switches are part of the application grid and will be set to
!  FALSE elsewhere, if the boundary point is assigned by a nested grid.
!
      allocate ( LBC_apply(ng) % west(LBj:UBj) )
      LBC_apply(ng) % west = .TRUE.
      allocate ( LBC_apply(ng) % east(LBj:UBj) )
      LBC_apply(ng) % east = .TRUE.
      allocate ( LBC_apply(ng) % south(LBi:UBi) )
      LBC_apply(ng) % south = .TRUE.
      allocate ( LBC_apply(ng) % north(LBi:UBi) )
      LBC_apply(ng) % north = .TRUE.
!
!  Nonlinear model state.
!
      IF (LBC(iwest,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_west(LBj:UBj) )
      END IF
      IF (LBC(ieast,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_east(LBj:UBj) )
      END IF
      IF (LBC(isouth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_south(LBi:UBi) )
      END IF
      IF (LBC(inorth,isFsur,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % zeta_north(LBi:UBi) )
      END IF
!
      IF (LBC(iwest,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_west(LBj:UBj) )
      END IF
      IF (LBC(ieast,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_east(LBj:UBj) )
      END IF
      IF (LBC(isouth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_south(LBi:UBi) )
      END IF
      IF (LBC(inorth,isUbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % ubar_north(LBi:UBi) )
      END IF
!
      IF (LBC(iwest,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_west(LBj:UBj) )
      END IF
      IF (LBC(ieast,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_east(LBj:UBj) )
      END IF
      IF (LBC(isouth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_south(LBi:UBi) )
      END IF
      IF (LBC(inorth,isVbar,ng)%acquire) THEN
        allocate ( BOUNDARY(ng) % vbar_north(LBi:UBi) )
      END IF
      RETURN
      END SUBROUTINE allocate_boundary
      SUBROUTINE initialize_boundary (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      real(r8), parameter :: IniVal = 0.0_r8
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
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isFsur,ng)%acquire) THEN
          BOUNDARY(ng) % zeta_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isUbar,ng)%acquire) THEN
          BOUNDARY(ng) % ubar_north = IniVal
        END IF
!
        IF (DOMAIN(ng)%NorthWest_Test(tile).and.                        &
     &      LBC(iwest,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_west = IniVal
        END IF
        IF (DOMAIN(ng)%SouthEast_Test(tile).and.                        &
     &      LBC(ieast,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_east = IniVal
        END IF
        IF (DOMAIN(ng)%SouthWest_Test(tile).and.                        &
     &      LBC(isouth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_south = IniVal
        END IF
        IF (DOMAIN(ng)%NorthEast_Test(tile).and.                        &
     &      LBC(inorth,isVbar,ng)%acquire) THEN
          BOUNDARY(ng) % vbar_north = IniVal
        END IF
      END IF
      RETURN
      END SUBROUTINE initialize_boundary
      END MODULE mod_boundary
