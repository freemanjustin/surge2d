      MODULE extract_sta_mod
!
!svn $Id: extract_sta.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine extracts field at the requested (Xpos,Ypos,Zpos)       !
!  positions.  The extraction is done using linear interpolation.      !
!  The (Xpos,Ypos) positions are in fractional grid coordinates.       !
!  Zpos is in fractional grid coordinates (Zpos >= 0) or actual        !
!  depths (Zpos < 0), if applicable.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Cgrid      Switch to interpolate at native C-grid (TRUE) or to   !
!                  interpolate at RHO-points (FALSE).                  !
!     ifield     Field ID.                                             !
!     gtype      Grid type.                                            !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     Ascl       Factor to scale field after extraction.               !
!     A          Tile array (2D or 3D) to process.                     !
!     Npos       Number of values to extract.                          !
!     Xpos       X-extraction positions (grid coordinates).            !
!     Ypos       Y-extraction positions (grid coordinates).            !
!     Zpos       Z-extraction positions (grid coordinates or depth).   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Apos       Extracted values.                                     !
!                                                                      !
!  Note:                                                               !
!                                                                      !
!  Starting F95 zero values can be signed (-0 or +0) following the     !
!  IEEE 754 floating point standard. This can be advantageous in       !
!  some computations but not here when "Ascl" is negative and "Apos"   !
!  is zero.  This will produce different output files in serial        !
!  and distributed memory applications. Since comparing serial and     !
!  parallel output is essential for tracking parallel partition        !
!  bugs, "positive zero" is enforced.                                  !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE extract_sta2d (ng, model, Cgrid, ifield, gtype,        &
     &                          LBi, UBi, LBj, UBj, Ascl, A,            &
     &                          Npos, Xpos, Ypos, Apos)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      logical, intent(in) :: Cgrid
      integer, intent(in) :: ng, model, ifield, gtype, Npos
      integer, intent(in) :: LBi, UBi, LBj, UBj
      real(r8), intent(in) :: Ascl
      real(r8), intent(in) :: A(LBi:,LBj:)
      real(r8), intent(in) :: Xpos(:), Ypos(:)
      real(r8), intent(out) :: Apos(Npos)
!
!  Local variable declarations.
!
      integer :: i1, i2, j1, j2, np
      real(r8), parameter :: Aspv = 0.0_r8
      real(r8) :: Xmin, Xmax, Ymin, Ymax
      real(r8) :: Xgrd, Xoff, Ygrd, Yoff
      real(r8) :: p1, p2, q1, q2, r1, r2, wsum
      real(r8) :: w111, w211, w121, w221
      real(r8), dimension(Npos) :: bounded
!
!-----------------------------------------------------------------------
!  Interpolate from 2D field at RHO-points.
!-----------------------------------------------------------------------
!
      IF (gtype.eq.r2dvar) THEN
        Xmin=rXmin(ng)
        Xmax=rXmax(ng)
        Ymin=rYmin(ng)
        Ymax=rYmax(ng)
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w111=w111*GRID(ng)%rmask(i1,j1)
            w211=w211*GRID(ng)%rmask(i2,j1)
            w121=w121*GRID(ng)%rmask(i1,j2)
            w221=w221*GRID(ng)%rmask(i2,j2)
            wsum=w111+w211+w121+w221
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
            ELSE
              bounded(np)=0.0_r8
            ENDIF
            Apos(np)=Ascl*(w111*A(i1,j1)+                               &
     &                     w211*A(i2,j1)+                               &
     &                     w121*A(i1,j2)+                               &
     &                     w221*A(i2,j2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 2D field at U-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.u2dvar) THEN
        IF (Cgrid) THEN
          Xmin=uXmin(ng)+0.5_r8
          Xmax=uXmax(ng)+0.5_r8
          Ymin=uYmin(ng)
          Ymax=uYmax(ng)
          Xoff=0.0_r8
        ELSE
          Xmin=rXmin(ng)
          Xmax=rXmax(ng)
          Ymin=rYmin(ng)
          Ymax=rYmax(ng)
          Xoff=0.5_r8
        END IF
        DO np=1,Npos
          Xgrd=Xpos(np)+Xoff
          Ygrd=Ypos(np)
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w111=w111*GRID(ng)%umask(i1,j1)
            w211=w211*GRID(ng)%umask(i2,j1)
            w121=w121*GRID(ng)%umask(i1,j2)
            w221=w221*GRID(ng)%umask(i2,j2)
            wsum=w111+w211+w121+w221
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
            Apos(np)=Ascl*(w111*A(i1,j1)+                               &
     &                     w211*A(i2,j1)+                               &
     &                     w121*A(i1,j2)+                               &
     &                     w221*A(i2,j2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Interpolate from 2D field at V-points.
!-----------------------------------------------------------------------
!
      ELSE IF (gtype.eq.v2dvar) THEN
        IF (Cgrid) THEN
          Xmin=vXmin(ng)
          Xmax=vXmax(ng)
          Ymin=vYmin(ng)+0.5_r8
          Ymax=vYmax(ng)+0.5_r8
          Yoff=0.0_r8
        ELSE
          Xmin=rXmin(ng)
          Xmax=rXmax(ng)
          Ymin=rYmin(ng)
          Ymax=rYmax(ng)
          Yoff=0.5_r8
        END IF
        DO np=1,Npos
          Xgrd=Xpos(np)
          Ygrd=Ypos(np)+Yoff
          bounded(np)=0.0_r8
          IF (((Xmin.le.Xgrd).and.(Xgrd.lt.Xmax)).and.                  &
     &        ((Ymin.le.Ygrd).and.(Ygrd.lt.Ymax))) THEN
            i1=INT(Xgrd)
            j1=INT(Ygrd)
            i2=i1+1
            j2=j1+1
            IF (i2.gt.Lm(ng)+1) THEN
              i2=i1                   ! station at the eastern boundary
            END IF
            IF (j2.gt.Mm(ng)+1) THEN
              j2=j1                   ! station at the northern boundary
            END IF
            bounded(np)=1.0_r8
            p2=REAL(i2-i1,r8)*(Xgrd-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Ygrd-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            w111=p1*q1
            w211=p2*q1
            w121=p1*q2
            w221=p2*q2
            w111=w111*GRID(ng)%vmask(i1,j1)
            w211=w211*GRID(ng)%vmask(i2,j1)
            w121=w121*GRID(ng)%vmask(i1,j2)
            w221=w221*GRID(ng)%vmask(i2,j2)
            wsum=w111+w211+w121+w221
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              w111=w111*wsum
              w211=w211*wsum
              w121=w121*wsum
              w221=w221*wsum
            ELSE
              bounded(np)=0.0_r8
            END IF
            Apos(np)=Ascl*(w111*A(i1,j1)+                               &
     &                     w211*A(i2,j1)+                               &
     &                     w121*A(i1,j2)+                               &
     &                     w221*A(i2,j2))
            IF (ABS(Apos(np)).eq.0.0_r8) Apos(np)=0.0_r8      ! positive
          ELSE                                                ! zero
            Apos(np)=Aspv
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Set unbounded data to special value.
!-----------------------------------------------------------------------
!
      DO np=1,Npos
        IF (bounded(np).lt.1.0_r8) THEN
          Apos(np)=spval
        END IF
      END DO
      RETURN
      END SUBROUTINE extract_sta2d
      END MODULE extract_sta_mod
