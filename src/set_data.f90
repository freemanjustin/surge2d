      SUBROUTINE set_data (ng, tile)
!
!svn $Id: set_data.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine processes forcing, boundary, climatology, and       !
!  other input data. It time-interpolates between snapshots.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL set_data_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS)
      RETURN
      END SUBROUTINE set_data
!
!***********************************************************************
      SUBROUTINE set_data_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
      USE mod_scalars
      USE mod_sources
!
      USE analytical_mod
      USE exchange_2d_mod
      USE set_2dfld_mod
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      logical :: SetBC
      logical :: update = .FALSE.
      logical :: bry_update
      integer :: ILB, IUB, JLB, JUB
      integer :: i, ic, itrc, j, k, my_tile
      real(r8) :: cff, cff1, cff2
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
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (m2/s2).
!-----------------------------------------------------------------------
!
      CALL set_2dfld_tile (ng, tile, iNLM, idUsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%sustrG,                           &
     &                     FORCES(ng)%sustr,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
      CALL set_2dfld_tile (ng, tile, iNLM, idVsms,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%svstrG,                           &
     &                     FORCES(ng)%svstr,                            &
     &                     update)
      IF (exit_flag.ne.NoError) RETURN
!
!  If input point wind stress, rotate to curvilinear grid.  Notice
!  that rotation is done at RHO-points. It does not matter.
!
      IF (.not.Linfo(1,idUsms,ng).or.                                   &
     &    (Iinfo(5,idUsms,ng).ne.Lm(ng)+1).or.                          &
     &    (Iinfo(6,idUsms,ng).ne.Mm(ng)+2)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            cff1=FORCES(ng)%sustr(i,j)*GRID(ng)%CosAngler(i,j)+         &
     &           FORCES(ng)%svstr(i,j)*GRID(ng)%SinAngler(i,j)
            cff2=FORCES(ng)%svstr(i,j)*GRID(ng)%CosAngler(i,j)-         &
     &           FORCES(ng)%sustr(i,j)*GRID(ng)%SinAngler(i,j)
            FORCES(ng)%sustr(i,j)=cff1
            FORCES(ng)%svstr(i,j)=cff2
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%sustr)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            FORCES(ng)%svstr)
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set surface air pressure (mb).
!-----------------------------------------------------------------------
!
      SetBC=.TRUE.
!     SetBC=.FALSE.
      CALL set_2dfld_tile (ng, tile, iNLM, idPair,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     FORCES(ng)%PairG,                            &
     &                     FORCES(ng)%Pair,                             &
     &                     update, SetBC)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set point Sources/Sinks (river runoff).
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        IF (LuvSrc(ng).or.LwSrc(ng)) THEN
          CALL set_ngfld (ng, iNLM, idRtra, 1, Nsrc(ng), 1,             &
     &                    1, Nsrc(ng), 1,                               &
     &                    SOURCES(ng) % QbarG,                          &
     &                    SOURCES(ng) % Qbar,                           &
     &                    update)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Set open boundary conditions fields.  In grid refinement, only the
!  coarser grid (RefineScale(ng)=0) open boundary conditions data is
!  processed and needed.
!-----------------------------------------------------------------------
!
!  Free-surface
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        CALL ana_fsobc (ng, tile, iNLM)
!
! Ensure that water level on boundary cells is above bed elevation.
!
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrR,JendR
              cff=Dcrit(ng)-GRID(ng)%h(0,j)
              IF (BOUNDARY(ng)%zeta_west(j).le.cff) THEN
                BOUNDARY(ng)%zeta_west(j)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
        END IF
        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrR,JendR
              cff=Dcrit(ng)-GRID(ng)%h(Lm(ng)+1,j)
              IF (BOUNDARY(ng)%zeta_east(j).le.cff) THEN
                BOUNDARY(ng)%zeta_east(j)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
        END IF
        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrR,IendR
              cff=Dcrit(ng)-GRID(ng)%h(i,0)
              IF (BOUNDARY(ng)%zeta_south(i).le.cff) THEN
                BOUNDARY(ng)%zeta_south(i)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
        END IF
        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          bry_update=.FALSE.
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrR,IendR
              cff=Dcrit(ng)-GRID(ng)%h(i,Mm(ng)+1)
              IF (BOUNDARY(ng)%zeta_north(i).le.cff) THEN
                BOUNDARY(ng)%zeta_north(i)=cff
              END IF
            END DO
            bry_update=.TRUE.
          END IF
        END IF
      END IF
!
!  2D momentum.
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        CALL ana_m2obc (ng, tile, iNLM)
      END IF
!
!-----------------------------------------------------------------------
!  Set sea surface height climatology (m).
!-----------------------------------------------------------------------
!
      IF (LsshCLM(ng)) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idSSHc,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%sshG,                            &
     &                       CLIMA(ng)%ssh,                             &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Set 2D momentum climatology (m/s).
!-----------------------------------------------------------------------
!
      IF (Lm2CLM(ng)) THEN
        CALL set_2dfld_tile (ng, tile, iNLM, idUbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%ubarclmG,                        &
     &                       CLIMA(ng)%ubarclm,                         &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL set_2dfld_tile (ng, tile, iNLM, idVbcl,                    &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       CLIMA(ng)%vbarclmG,                        &
     &                       CLIMA(ng)%vbarclm,                         &
     &                       update)
        IF (exit_flag.ne.NoError) RETURN
      END IF
      RETURN
      END SUBROUTINE set_data_tile
