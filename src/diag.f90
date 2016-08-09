      MODULE diag_mod
!
!svn $Id: diag.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes various diagnostic fields.                    !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: diag
      CONTAINS
!
!***********************************************************************
      SUBROUTINE diag (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
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
      CALL diag_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                nstp(ng), krhs(ng),                               &
     &                GRID(ng) % h,                                     &
     &                GRID(ng) % pm,                                    &
     &                GRID(ng) % pn,                                    &
     &                GRID(ng) % omn,                                   &
     &                OCEAN(ng) % ubar,                                 &
     &                OCEAN(ng) % vbar,                                 &
     &                OCEAN(ng) % zeta)
      RETURN
      END SUBROUTINE diag
!
!***********************************************************************
      SUBROUTINE diag_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp, krhs,                                 &
     &                      h, pm, pn, omn,                             &
     &                      ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, krhs
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: NSUB, i, ispace, j, k, trd
      integer :: my_max_Ci, my_max_Cj, my_max_Ck
      integer :: my_threadnum
      real(r8) :: cff, my_avgke, my_avgpe, my_volume
      real(r8) :: my_C , my_max_C
      real(r8) :: my_Cu, my_max_Cu
      real(r8) :: my_Cv, my_max_Cv
      real(r8) :: my_maxspeed, u2v2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ke2d
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: pe2d
      character (len=8 ) :: kechar, pechar
      character (len=60) :: frmt
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
!  Compute and report out volume averaged kinetic, potential
!  total energy, volume, Courant numbers.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng)-1,ninfo(ng)).eq.0) THEN
        my_max_C =0.0_r8
        my_max_Cu=0.0_r8
        my_max_Cv=0.0_r8
        my_max_Ci=0
        my_max_Cj=0
        my_max_Ck=0
        my_maxspeed=0.0_r8
        DO j=Jstr,Jend
          cff=0.5_r8*g
          DO i=Istr,Iend
            u2v2=ubar(i  ,j,krhs)*ubar(i  ,j,krhs)+                     &
     &           ubar(i+1,j,krhs)*ubar(i+1,j,krhs)+                     &
     &           vbar(i,j  ,krhs)*vbar(i,j  ,krhs)+                     &
     &           vbar(i,j+1,krhs)*vbar(i,j+1,krhs)
            ke2d(i,j)=(zeta(i,j,krhs)+h(i,j))*0.25_r8*u2v2
            pe2d(i,j)=cff*zeta(i,j,krhs)*zeta(i,j,krhs)
            my_Cu=0.5_r8*ABS(ubar(i,j,krhs)+ubar(i+1,j,krhs))*          &
     &            dt(ng)*pn(i,j)
            my_Cv=0.5_r8*ABS(vbar(i,j,krhs)+vbar(i,j+1,krhs))*          &
     &            dt(ng)*pn(i,j)
            my_C=my_Cu+my_Cv
            IF (my_C.gt.my_max_C) THEN
              my_max_C =my_C
              my_max_Cu=my_Cu
              my_max_Cv=my_Cv
              my_max_Ci=i
              my_max_Cj=j
            END IF
            my_maxspeed=MAX(my_maxspeed,SQRT(0.5_r8*u2v2))
          END DO
        END DO
!
!  Integrate horizontally within one tile. In order to reduce the
!  round-off errors, the summation is performed in two stages. First,
!  the index j is collapsed and then the accumulation is carried out
!  along index i. In this order, the partial sums consist on much
!  fewer number of terms than in a straightforward two-dimensional
!  summation. Thus, adding numbers which are orders of magnitude
!  apart is avoided.
!
        DO i=Istr,Iend
          pe2d(i,Jend+1)=0.0_r8
          pe2d(i,Jstr-1)=0.0_r8
          ke2d(i,Jstr-1)=0.0_r8
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend
            pe2d(i,Jend+1)=pe2d(i,Jend+1)+                              &
     &                     omn(i,j)*(zeta(i,j,krhs)+h(i,j))
            pe2d(i,Jstr-1)=pe2d(i,Jstr-1)+omn(i,j)*pe2d(i,j)
            ke2d(i,Jstr-1)=ke2d(i,Jstr-1)+omn(i,j)*ke2d(i,j)
          END DO
        END DO
        my_volume=0.0_r8
        my_avgpe=0.0_r8
        my_avgke=0.0_r8
        DO i=Istr,Iend
          my_volume=my_volume+pe2d(i,Jend+1)
          my_avgpe =my_avgpe +pe2d(i,Jstr-1)
          my_avgke =my_avgke +ke2d(i,Jstr-1)
        END DO
!
!  Perform global summation: whoever gets first to the critical region
!  resets global sums before global summation starts; after the global
!  summation is completed, thread, which is the last one to enter the
!  critical region, finalizes the computation of diagnostics and prints
!  them out.
!
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
!$OMP CRITICAL (NL_DIAGNOSTICS)
        IF (tile_count.eq.0) THEN
          volume=my_volume
          avgke=my_avgke
          avgpe=my_avgpe
          maxspeed(ng)=my_maxspeed
          max_C =my_max_C
          max_Cu=my_max_Cu
          max_Cv=my_max_Cv
          max_Ci=my_max_Ci
          max_Cj=my_max_Cj
        ELSE
          volume=volume+my_volume
          avgke=avgke+my_avgke
          avgpe=avgpe+my_avgpe
          maxspeed(ng)=MAX(maxspeed(ng),my_maxspeed)
          IF (my_max_C.eq.max_C) THEN
            max_Ci=MIN(max_Ci,my_max_Ci)
            max_Cj=MIN(max_Cj,my_max_Cj)
          ELSE IF (my_max_C.gt.max_C) THEN
            max_C =my_max_C
            max_Cu=my_max_Cu
            max_Cv=my_max_Cv
            max_Ci=my_max_Ci
            max_Cj=my_max_Cj
          END IF
        END IF
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
          trd=my_threadnum()
          avgke=avgke/volume
          avgpe=avgpe/volume
          avgkp=avgke+avgpe
          IF (first_time(ng).eq.0) THEN
            first_time(ng)=1
            IF (Master.and.(ng.eq.1)) THEN
              WRITE (stdout,10) 'STEP', 'Day HH:MM:SS',                 &
     &                          'KINETIC_ENRG', 'POTEN_ENRG',           &
     &                          'TOTAL_ENRG', 'NET_VOLUME'
              WRITE (stdout,20) '  C => (i,j)', 'Cu', 'Cv',             &
     &                          ' C Max', 'Max Speed'
 10           FORMAT (/,3x,a,3x,a,2x,a,3x,a,4x,a,4x,a)
 20           FORMAT (10x,a,7x,a,12x,a,10x,a,7x,a,/)
            END IF
          END IF
          IF (Master) THEN
            WRITE (stdout,30) iic(ng)-1, time_code(ng),                 &
     &                        avgke, avgpe, avgkp, volume
            ispace=24-(5+Idigits(ng)+Jdigits(ng))
            WRITE (frmt,40) ispace,                                     &
     &                      '"("', Idigits(ng), Idigits(ng),            &
     &                      '","', Jdigits(ng), Jdigits(ng), '")"'
            WRITE (stdout,frmt) max_Ci, max_Cj,                         &
     &                          max_Cu, max_Cv, max_C,                  &
     &                          maxspeed(ng)
            CALL my_flush (stdout)
 30         FORMAT (i7,1x,a,4(1pe14.6))
 40         FORMAT ('(',i2.2,'x,',a,',i',i1,'.',i1,',',                 &
     &                            a,',i',i1,'.',i1,',',                 &
     &                            a,',t24,4(1pe13.6,1x))')
          END IF
!
!  If blowing-up, set exit_flag to stop computations.
!
          WRITE (kechar,'(1pe8.1)') avgke
          WRITE (pechar,'(1pe8.1)') avgpe
          DO i=1,8
            IF ((kechar(i:i).eq.'N').or.(pechar(i:i).eq.'N').or.        &
     &          (kechar(i:i).eq.'n').or.(pechar(i:i).eq.'n').or.        &
     &          (kechar(i:i).eq.'*').or.(pechar(i:i).eq.'*')) THEN
              exit_flag=1
            END IF
          END DO
!
!  Stop computations if exceeding maximum speed allowed.  This will be
!  useful during debugging to avoid NaNs in output NetCDF files.
!
          IF (maxspeed(ng).gt.max_speed) THEN
            exit_flag=1
          END IF
        END IF
!$OMP END CRITICAL (NL_DIAGNOSTICS)
      END IF
      RETURN
      END SUBROUTINE diag_tile
      END MODULE diag_mod
