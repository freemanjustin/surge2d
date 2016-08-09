      SUBROUTINE wclock_on (ng, model, region)
!
!svn $Id: timers.F 795 2016-05-11 01:42:43Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine turns on wall clock to meassure the elapsed time in    !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: iregion, MyModel, NSUB
      integer :: my_getpid
      integer :: my_threadnum
      real(r8), dimension(2) :: wtime
      real(r8) :: my_wtime
!
!-----------------------------------------------------------------------
! Initialize timing for all threads.
!-----------------------------------------------------------------------
!
!  Set number of subdivisions, same as for global reductions.
!
      NSUB=numthreads
!
!  Insure that MyModel is not zero.  Notice that zero value is used to
!  indicate restart of the nonlinear model.
!
      MyModel=MAX(1,model)
      Cstr(region,MyModel,ng)=my_wtime(wtime)
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.0)) THEN
        DO iregion=1,Nregion
          Cend(iregion,MyModel,ng)=0.0_r8
          Csum(iregion,MyModel,ng)=0.0_r8
        END DO
        proc(1,MyModel,ng)=1
        proc(0,MyModel,ng)=my_getpid()
!$OMP CRITICAL (START_WCLOCK)
        WRITE (stdout,10) ' Thread #', MyThread,                        &
     &                    ' (pid=',proc(0,MyModel,ng),') is active.'
 10     FORMAT (a,i3,a,i8,a)
        thread_count=thread_count+1
        IF (thread_count.eq.NSUB) thread_count=0
!$OMP END CRITICAL (START_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_on
      SUBROUTINE wclock_off (ng, model, region)
!
!=======================================================================
!                                                                      !
!  This routine turns off wall clock to meassure the elapsed time in   !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: imodel, iregion, MyModel, NSUB
      integer :: my_threadnum
      real(r8) :: percent, sumcpu, sumper, total
      real(r8), dimension(2) :: wtime
      real(r8) :: my_wtime
      character (len=14), dimension(4) :: label
!
!-----------------------------------------------------------------------
!  Compute elapsed wall time for all threads.
!-----------------------------------------------------------------------
!
!  Set number of subdivisions, same as for global reductions.
!
      NSUB=numthreads
!
!  Insure that MyModel is not zero.  Notice that zero value is used to
!  indicate restart of the nonlinear model.
!
      MyModel=MAX(1,model)
      IF (region.ne.0) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
      END IF
!
!  Report elapsed wall time.
!
      IF ((region.eq.0).and.(proc(1,MyModel,ng).eq.1)) THEN
        Cend(region,MyModel,ng)=Cend(region,MyModel,ng)+                &
     &                          (my_wtime(wtime)-                       &
     &                           Cstr(region,MyModel,ng))
        DO imodel=1,4
          proc(1,imodel,ng)=0
        END DO
!$OMP CRITICAL (FINALIZE_WCLOCK)
!
!  Report total elapsed time (seconds) for each CPU.  We get the same
!  time for all nested grids.
!
        IF (ng.eq.1) THEN
         WRITE (stdout,10) ' Thread #', MyThread, ' CPU:',              &
     &                     Cend(region,MyModel,ng)
 10      FORMAT (a,i3,a,f12.3)
       END IF
!
! Report elapsed time profile for each region of the code.
!
        thread_count=thread_count+1
        DO imodel=1,4
          Csum(region,imodel,ng)=Csum(region,imodel,ng)+                &
     &                           Cend(region,imodel,ng)
          DO iregion=1,Nregion
            Csum(iregion,imodel,ng)=Csum(iregion,imodel,ng)+            &
     &                              Cend(iregion,imodel,ng)
          END DO
        END DO
        DO imodel=1,4
          IF (imodel.ne.MyModel) THEN
            DO iregion=1,Nregion
              Csum(region,imodel,ng)=Csum( region,imodel,ng)+           &
     &                               Csum(iregion,imodel,ng)
            END DO
          END IF
        END DO
        IF (thread_count.eq.NSUB) THEN
          thread_count=0
          IF (Master.and.(ng.eq.1)) THEN
            total=0.0_r8
            DO imodel=1,4
              total=total+Csum(region,imodel,ng)
            END DO
            WRITE (stdout,20) ' Total:', total
 20         FORMAT (a,8x,f14.3)
          END IF
        END IF
!$OMP END CRITICAL (FINALIZE_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_off
