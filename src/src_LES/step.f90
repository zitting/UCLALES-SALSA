!----------------------------------------------------------------------------
! This file is part of UCLALES.
!testi
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE step

   IMPLICIT NONE

   INTEGER :: istpfl = 1
   REAL    :: timmax = 18000.
   LOGICAL :: corflg = .FALSE.

   REAL    :: frqhis =  9000.
   REAL    :: frqanl =  3600.
   REAL    :: radfrq =  0.

   REAL    :: time   =  0.
   REAL    :: strtim =  0.0
   REAL    :: cntlat =  31.5 ! 30.0
   LOGICAL :: outflg = .TRUE.

   ! Nudging options (nudge_*: 0=disabled, 1=soft, 2=hard), total nudging time (s), altitude range, and time constants (tau [s])
   INTEGER :: nudge_theta = 0, & ! (liquid water) potential temperature, depending on the microphysical level
              nudge_rv = 0, & ! Water vapor mixing ratio (maintain total water)
              nudge_u = 0, nudge_v = 0, & ! Horizontal winds
              nudge_ccn = 0 ! Sectional aerosol for levels 4 and 5 (maintain aerosol+cloud+ice)
   REAL    :: nudge_time = 3600., nudge_zmin = -1.e10, nudge_zmax = 1.e10
   REAL    :: tau_theta = 300., tau_rv = 300., tau_u = 300., tau_v = 300., tau_ccn = 300.
   REAL, SAVE, ALLOCATABLE :: theta_ref(:), rv_ref(:), u_ref(:), v_ref(:), aero_ref(:,:)

CONTAINS
   !
   ! ----------------------------------------------------------------------
   ! Subroutine model:  This is the main driver for the model's time
   ! integration.  It calls the routine tstep, which steps through the
   ! physical processes active on a time-step and updates variables.  It
   ! then checks to see whether or not different output options are
   ! satisfied.
   SUBROUTINE stepper

      USE mpi_interface, ONLY : myid, double_scalar_par_max

      USE grid, ONLY : dtl, dzt, zt, zm, nzp, dn0, u0, v0, a_up, a_vp, a_wp, &
                       a_uc, a_vc, a_wc, write_hist, write_anal, close_anal, dtlt,  &
                       dtlv, dtlong, nzp, nyp, nxp, level,                          &
                       ! For mass budged
                       a_rp, a_rc, a_srp, a_dn, calc_anal

      USE stat, ONLY : sflg, savg_intvl, ssam_intvl, write_ps, close_stat, mcflg, acc_massbudged,  &
                       write_massbudged
      USE thrm, ONLY : thermo

      LOGICAL, PARAMETER :: StopOnCFLViolation = .FALSE.
      REAL, PARAMETER :: cfl_upper = 0.50, cfl_lower = 0.30

      REAL         :: t1,t2,tplsdt,begtime
      REAL(kind=8) :: cflmax,gcflmax
      INTEGER      :: istp, iret
      LOGICAL :: cflflg
      !
      ! Timestep loop for program
      !
      begtime = time
      istp = 0

      CALL cpu_time(t1)

      DO WHILE (time + 0.1*dtl < timmax)

         istp = istp+1
         tplsdt = time + dtl + 0.1*dtl
         sflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dtl  &
            .OR. tplsdt >= timmax  .OR. tplsdt < 2.*dtl)

         CALL t_step(cflflg,cflmax)

         time = time + dtl

         CALL double_scalar_par_max(cflmax,gcflmax)
         cflmax = gcflmax

         IF (cflmax > cfl_upper .OR. cflmax < cfl_lower) THEN
            CALL tstep_reset(nzp,nxp,nyp,a_up,a_vp,a_wp,a_uc,a_vc,a_wc,     &
                             dtl,dtlong,cflmax,cfl_upper,cfl_lower)
            dtlv = 2.*dtl
            dtlt = dtl
         END IF

         !
         ! output control
         !
         IF (mod(tplsdt,savg_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_ps(nzp,dn0,u0,v0,zm,zt,time)

         IF ((mod(tplsdt,frqhis) < dtl .OR. time >= timmax) .AND. outflg)   &
            CALL write_hist(2, time)
         IF (mod(tplsdt,savg_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_hist(1, time)

         IF ((mod(tplsdt,frqanl) < dtl .OR. time >= timmax) .AND. outflg) THEN
            CALL thermo(level)
            !CALL calc_anal()
            CALL write_anal(time)
         END IF

         IF (cflflg) THEN
            cflflg = .FALSE.
            IF (StopOnCFLViolation) CALL write_hist(-1,time)
         END IF

         IF(myid == 0) THEN
            CALL cpu_time(t2) ! t2-t1 is the actual time from output
            IF (mod(istp,istpfl) == 0 ) THEN
               PRINT "('   Timestep # ',i6," //     &
                  "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3)",     &
                  istp, time, t2-t1
               CALL cpu_time(t1)
            END IF
         END IF

      END DO

      IF (mcflg) THEN
         !
         ! Juha:
         ! Get the final statistics of atmospheric water for mass budged
         CALL acc_massbudged(nzp,nxp,nyp,1,dtlt,dzt,a_dn%data,    &
                             rv=a_rp%data,rc=a_rc%data,prc=a_srp%data)

         CALL write_massbudged

      END IF ! mcflg

      CALL write_hist(1, time)
      iret = close_anal()
      iret = close_stat()

   END SUBROUTINE stepper
   !
   !----------------------------------------------------------------------
   ! Subroutine tstep_reset: Called to adjust current velocity and reset
   ! timestep based on cfl limits
   !
   SUBROUTINE tstep_reset(n1,n2,n3,up,vp,wp,uc,vc,wc,dtl,dtmx,cfl,c1,c2)

      INTEGER, INTENT (in)      :: n1,n2,n3
      REAL, INTENT (in)         :: up(n1,n2,n3),vp(n1,n2,n3),wp(n1,n2,n3),dtmx,c1,c2
      REAL(kind=8), INTENT (IN) :: cfl
      REAL, INTENT (inout)      :: uc(n1,n2,n3),vc(n1,n2,n3),wc(n1,n2,n3),dtl

      INTEGER :: i,j,k
      REAL    :: cbar, dtl_old

      cbar = (c1+c2)*0.5
      dtl_old = dtl

      IF (cfl > c1) dtl = min(dtmx,dtl*cbar/c1)
      IF (cfl < c2) dtl = min(dtmx,dtl*cbar/c2)

      DO j = 1, n3
         DO i = 1, n2
            DO k = 1, n1
               uc(k,i,j) = up(k,i,j) + (uc(k,i,j)-up(k,i,j))*dtl/dtl_old
               vc(k,i,j) = vp(k,i,j) + (vc(k,i,j)-vp(k,i,j))*dtl/dtl_old
               wc(k,i,j) = wp(k,i,j) + (wc(k,i,j)-wp(k,i,j))*dtl/dtl_old
            END DO
         END DO
      END DO

   END SUBROUTINE tstep_reset

   !
   !----------------------------------------------------------------------
   ! Subroutine t_step: Called by driver to timestep through the LES
   ! routines.  Within many subroutines, data is accumulated during
   ! the course of a timestep for the purposes of statistical analysis.
   !
   SUBROUTINE t_step(cflflg,cflmax)

      USE grid, ONLY : level, dtl, dtlt, Tspinup,                                         &
                       ! Added parameters for interfacing with SALSA
                       nxp, nyp, nzp, a_press, a_temp, a_rsl,                             &
                       a_rc, a_wp, a_rp, a_rt, a_rh,                                      &
                       a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                       a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                       a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                       a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                       a_gaerop, a_gaerot, a_dn,  a_nactd,  a_vactd,   prtcl,             &
                       sst, a_rsi, a_temp0


      USE stat, ONLY : sflg, statistics
      USE sgsm, ONLY : diffuse
      USE srfc, ONLY : surface
      USE thrm, ONLY : thermo
      USE mcrp, ONLY : micro
      USE prss, ONLY : poisson
      USE advf, ONLY : fadvect, newdroplet
      USE advl, ONLY : ladvect
      USE forc, ONLY : forcings
      USE util, ONLY : maskactiv !Juha: Included for SALSA
      USE nudg, ONLY : nudge_theta, nudge_rv, nudge_u, nudge_v, nudge_ccn, nudge_time, nudge_zmin, nudge_zmax, &
                       tau_theta, tau_rv, tau_u, tau_v, tau_ccn, theta_ref, rv_ref, u_ref, v_ref, aero_ref, &
                       nudging, nudge_any, nudge_any_2d, usenudge

      USE mo_salsa_driver, ONLY : run_SALSA
      USE class_ComponentIndex, ONLY : GetNcomp

      LOGICAL, INTENT (out)      :: cflflg
      REAL(KIND=8), INTENT (out) :: cflmax

      REAL :: xtime

      LOGICAL :: zactmask(nzp,nxp,nyp)
      REAL    :: zwp(nzp,nxp,nyp), &  !! FOR SINGLE-COLUMN RUNS
                 ztkt(nzp,nxp,nyp)
      INTEGER :: zrm

      INTEGER :: n4

      zwp = 0.5

      xtime = time/86400. + strtim
      cflflg = .FALSE.

      ! The runmode parameter zrm is used by SALSA only
      zrm = 3
      IF ( time < Tspinup ) zrm = 2


      ! Reset ALL tendencies here.
      !----------------------------------------------------------------
      ! "Scalar" timestep
      CALL tend0(.FALSE.)

      ! Put the newly activated to zero
      IF (level >= 4) THEN
         a_vactd = 0.
         a_nactd = 0.
      END IF

      IF (level >= 4 .AND. time < 1.) THEN
         CALL thermo(level)
         CALL SALSA_diagnostics
      END IF

      CALL surface(sst)

      CALL diffuse

      CALL sponge(0)

      IF (level >= 1) THEN

         CALL thermo(level)

         CALL forcings(xtime,cntlat,sst)

         IF (level >= 4) THEN

            n4 = GetNcomp(prtcl) + 1 ! Aerosol components + water

            CALL tend_constrain(n4)
            CALL update_sclrs
            CALL tend0(.TRUE.)

            ! Rate of change in absolute temperature (for some ice processes)
            IF (time >= 1.) THEN
               ztkt = a_temp-a_temp0
               a_temp0 = a_temp
            ELSE IF (time == 0.) THEN
               a_temp0 = a_temp
               ztkt = 0.
            END IF

            IF ( nxp == 5 .AND. nyp == 5 ) THEN
               ! 1D -runs
               CALL run_SALSA(nxp, nyp, nzp, n4, a_press, a_temp, ztkt,                   &
                              a_rp, a_rt, a_rsl, a_rsi, zwp, a_dn,    &
                              a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                              a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                              a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                              a_nicep,   a_nicet,   a_micep,   a_micet,    &
                              a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                              a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                              zrm, prtcl, dtlt, level   )

            ELSE
               !! for 2D or 3D runs
               CALL run_SALSA(nxp, nyp, nzp, n4, a_press, a_temp, ztkt,                   &
                              a_rp, a_rt, a_rsl, a_rsi, a_wp, a_dn,    &
                              a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                              a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                              a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                              a_nicep,   a_nicet,   a_micep,   a_micet,    &
                              a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                              a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                              zrm, prtcl, dtlt, level   )

            END IF !nxp==5 and nyp == 5
          
            CALL tend_constrain(n4)

         END IF

      END IF ! level

      CALL update_sclrs

      !-------------------------------------------
      ! "Deposition" timestep
      ! -- Reset only scalar tendencies
      CALL tend0(.TRUE.)

      ! Dont perform sedimentation or level 3 autoconversion during spinup
      IF (zrm == 3) CALL micro(level)

      IF (level >= 4) CALL tend_constrain(n4)
      CALL update_sclrs

      !-------------------------------------------
      ! "Advection" timestep
      ! -- Reset only scalar tendencies
      CALL tend0(.TRUE.)

      ! Mask for cloud base activation
      IF (level >= 4) CALL maskactiv(zactmask,nxp,nyp,nzp,2,a_rh%data,rc=a_rc%data,w=a_wp)
      ! Get tendencies from cloud base activation
      IF (level >= 4) CALL newdroplet(zactmask)

      CALL fadvect

      IF (level >= 4)  &
         CALL tend_constrain(n4)

      CALL update_sclrs

      CALL thermo(level)

      ! Nudging time is independent of the spinup!
      IF (usenudge) THEN
         IF (time <= nudge_time .AND. (nudge_theta /= 0 .OR. nudge_rv /= 0 .OR. &
            (level > 3 .AND. nudge_ccn /= 0) ) ) THEN

            ! Reset tendencies
            CALL tend0(.TRUE.)

            CALL nudging(time)

            CALL update_sclrs

            CALL thermo(level)
         END IF
      END IF

      IF (level >= 4)  THEN
         CALL SALSA_diagnostics
         CALL thermo(level)
      END IF

      CALL corlos

      CALL ladvect

      CALL buoyancy

      CALL sponge(1)

      CALL poisson

      CALL cfl (cflflg, cflmax)

      CALL thermo(level)

      IF (sflg) THEN
         CALL statistics (time+dtl)
      END IF

   END SUBROUTINE t_step

   !
   !----------------------------------------------------------------------
   ! Subroutine tend0: sets all tendency arrays to zero
   !
   SUBROUTINE tend0(sclONLY)

      USE grid, ONLY : a_ut, a_vt, a_wt, nscl, a_st, newsclr

      LOGICAL, INTENT(in) :: sclONLY ! If true, only put scalar tendencies to zero

      INTEGER :: n

      IF( .NOT. sclONLY) THEN
         a_ut = 0.; a_vt = 0.; a_wt = 0.
      END IF
      DO n = 1, nscl
         CALL newsclr(n)
         a_st = 0.
      END DO

   END SUBROUTINE tend0
   !
   !----------------------------------------------------------------------
   ! In case of negative tendencies to SALSA arrays, put some constrains
   ! in order to avoid concentrations going negative. This will possibly
   ! slightly affect the conservation of mass - needs testing/revision
   ! Juha Tonttila, FMI, 2014
   !
   SUBROUTINE tend_constrain(nn)

      USE grid, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                       a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                       a_nicep,  a_nicet, a_nsnowp, a_nsnowt,                            & ! ice'n'snow
                       a_micep,  a_micet, a_msnowp, a_msnowt,                            & ! ice'n'snow
                       dtlt, nxp,nyp,nzp
      USE mo_submctl, ONLY : nbins, ncld, nprc, &
                             nice, nsnw          !ice'n'snow

      INTEGER, INTENT(in) :: nn

      INTEGER :: cc, ii,jj,kk,ni

      DO jj = 3, nyp-2

         DO ii = 3, nxp-2

            DO kk = 1, nzp

               ! Aerosols
               DO cc = 1, nbins

                  IF ( a_naerop%data(kk,ii,jj,cc)+a_naerot%data(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_naerot%data(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_naerop%data(kk,ii,jj,cc))/dtlt,a_naerot%data(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_maerot%data(kk,ii,jj,(ni-1)*nbins+cc) = &
                                    MAX( ((1.e-10-1.0)*a_maerop%data(kk,ii,jj,(ni-1)*nbins+cc))/dtlt,  &
                                        a_maerot%data(kk,ii,jj,(ni-1)*nbins+cc) )
                     END DO

                  END IF

               END DO

               ! Cloud droplets
               DO cc = 1, ncld

                  IF ( a_ncloudp%data(kk,ii,jj,cc)+a_ncloudt%data(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_ncloudt%data(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_ncloudp%data(kk,ii,jj,cc))/dtlt,a_ncloudt%data(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_mcloudt%data(kk,ii,jj,(ni-1)*ncld+cc) = &
                                       MAX( ((1.e-10-1.0)*a_mcloudp%data(kk,ii,jj,(ni-1)*ncld+cc))/dtlt,  &
                                           a_mcloudt%data(kk,ii,jj,(ni-1)*ncld+cc) )
                     END DO

                  END IF

               END DO

               ! Precipitation
               DO cc = 1, nprc

                  IF ( a_nprecpp%data(kk,ii,jj,cc)+a_nprecpt%data(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nprecpt%data(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nprecpp%data(kk,ii,jj,cc))/dtlt,a_nprecpt%data(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_mprecpt%data(kk,ii,jj,(ni-1)*nprc+cc) = &
                                       MAX( ((1.e-10-1.0)*a_mprecpp%data(kk,ii,jj,(ni-1)*nprc+cc))/dtlt,  &
                                           a_mprecpt%data(kk,ii,jj,(ni-1)*nprc+cc) )
                     END DO

                  END IF

               END DO

               ! ice particles
               DO cc = 1, nice

                  IF ( a_nicep%data(kk,ii,jj,cc)+a_nicet%data(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nicet%data(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nicep%data(kk,ii,jj,cc))/dtlt,a_nicet%data(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_micet%data(kk,ii,jj,(ni-1)*ncld+cc) = &
                                     MAX( ((1.e-10-1.0)*a_micep%data(kk,ii,jj,(ni-1)*nice+cc))/dtlt,  &
                                         a_micet%data(kk,ii,jj,(ni-1)*nice+cc) )
                     END DO

                  END IF

               END DO

               ! Snow
               DO cc = 1, nsnw

                  IF ( a_nsnowp%data(kk,ii,jj,cc)+a_nsnowt%data(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nsnowt%data(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nsnowp%data(kk,ii,jj,cc))/dtlt,a_nsnowt%data(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_msnowt%data(kk,ii,jj,(ni-1)*nprc+cc) = &
                                      MAX( ((1.e-10-1.0)*a_msnowp%data(kk,ii,jj,(ni-1)*nsnw+cc))/dtlt,  &
                                          a_msnowt%data(kk,ii,jj,(ni-1)*nsnw+cc) )
                     END DO

                  END IF

               END DO

            END DO ! kk

         END DO ! ii

      END DO ! jj

   END SUBROUTINE tend_constrain
   !
   !----------------------------------------------------------------------
   ! Subroutine cfl: Driver for calling CFL computation subroutine
   !
   SUBROUTINE cfl(cflflg,cflmax)

      USE grid, ONLY : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dtlt
      USE stat, ONLY : fill_scalar

      LOGICAL, INTENT(out) :: cflflg
      REAL(KIND=8), INTENT (out)   :: cflmax
      REAL, PARAMETER :: cflnum = 0.95

      cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dtlt)

      cflflg = (cflmax > cflnum)
      IF (cflflg) PRINT *, 'Warning CFL Violation :', cflmax
      CALL fill_scalar(1,REAL(cflmax))

   END SUBROUTINE cfl
   !
   !----------------------------------------------------------------------
   ! Subroutine cfll: Checks CFL criteria, brings down the model if the
   ! maximum thershold is exceeded
   !
   REAL(KIND=8) FUNCTION cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

      INTEGER, INTENT (in) :: n1, n2, n3
      REAL, DIMENSION (n1,n2,n3), INTENT (in) :: u, v, w
      REAL, INTENT (in)    :: dxi,dyi,dzt(n1),dtlt

      INTEGER :: i, j, k
      cfll = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               cfll = max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                      abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
            END DO
         END DO
      END DO

   END FUNCTION cfll
   !
   !----------------------------------------------------------------------
   ! Subroutine update_sclrs:  Updates scalars by applying tendency and
   ! boundary conditions
   !
   SUBROUTINE update_sclrs

      USE grid, ONLY : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
                       dtlt, newsclr, isgstyp
      USE sgsm, ONLY : tkeinit
      USE util, ONLY : sclrset

      INTEGER :: n

      DO n = 1, nscl
         CALL newsclr(n)
         CALL update(nzp,nxp,nyp,a_sp,a_st,dtlt)
         CALL sclrset('mixd',nzp,nxp,nyp,a_sp,dzt)
      END DO

      IF (isgstyp == 2) THEN
         CALL tkeinit(nxyzp,a_qp%data)
      END IF

   END SUBROUTINE update_sclrs
   !
   ! ----------------------------------------------------------------------
   ! Subroutine update:
   !
   SUBROUTINE update(n1,n2,n3,a,fa,dt)

      INTEGER, INTENT(in)   :: n1, n2, n3
      REAL, INTENT (in)     :: fa(n1,n2,n3),dt
      REAL, INTENT (inout)  :: a(n1,n2,n3)
      INTEGER :: i, j, k

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               a(k,i,j) = a(k,i,j) + fa(k,i,j)*dt
            END DO
         END DO
      END DO

   END SUBROUTINE update
   !
   ! ----------------------------------------------------------------------
   ! Subroutine buoyancy:
   !
   SUBROUTINE buoyancy

      USE grid, ONLY : a_uc, a_vc, a_wc, a_wt, a_rv, a_rc, a_theta, &
                       a_rp, nxp, nyp, nzp, dzm, th00, level, pi1
      USE stat, ONLY : sflg, comp_tke
      USE util, ONLY : ae1mm
      USE thrm, ONLY : update_pi1

      REAL, DIMENSION (nzp) :: awtbar, a_tmp1(nzp,nxp,nyp)

      IF (level < 4) THEN
         CALL boyanc(nzp,nxp,nyp,level,a_wt,a_theta%data,a_rp%data,th00,a_tmp1,a_rv%data)
      ELSE IF (level >= 4) THEN
         CALL boyanc(nzp,nxp,nyp,level,a_wt,a_theta%data,a_rp%data,th00,a_tmp1,a_rc%data)
      END IF

      CALL ae1mm(nzp,nxp,nyp,a_wt,awtbar)
      CALL update_pi1(nzp,awtbar,pi1)

      IF (sflg)  CALL comp_tke(nzp,nxp,nyp,dzm,th00,a_uc,a_vc,a_wc,a_tmp1)

   END SUBROUTINE buoyancy
   !
   ! ----------------------------------------------------------------------
   ! Subroutine boyanc:
   !
   SUBROUTINE boyanc(n1,n2,n3,level,wt,th,rt,th00,scr,rx)

      USE defs, ONLY : g, ep2

      INTEGER, INTENT(in) :: n1,n2,n3,level
      REAL, INTENT(in)    :: th00,th(n1,n2,n3),  &
                             rt(n1,n2,n3)  ! This is total water mix rat for level < 4
                                           ! and water vapour mix rat for level = 4
      REAL, INTENT(in)    :: rx(n1,n2,n3)  ! This should be water vapour mix rat for level < 4
                                           ! and cloud liquid water mix rat for level = 4 (including rain??)
      REAL, INTENT(inout) :: wt(n1,n2,n3)
      REAL, INTENT(out)   :: scr(n1,n2,n3)

      INTEGER :: k, i, j
      REAL    :: gover2

      gover2 = 0.5*g

      DO j = 3, n3-2
         DO i = 3, n2-2
            IF (level >= 2 .AND. level < 4) THEN
               DO k = 1, n1
                  scr(k,i,j) = gover2*((th(k,i,j)*(1.+ep2*rx(k,i,j))-th00)       &
                                      /th00-(rt(k,i,j)-rx(k,i,j)))
               END DO
            ELSE IF (level >= 4) THEN
               DO k = 1, n1
                  scr(k,i,j) = gover2*((th(k,i,j)*(1.+ep2*rt(k,i,j))-th00)       &
                                      /th00-(rx(k,i,j)))
               END DO
            ELSE
               DO k = 1, n1
                  scr(k,i,j) = gover2*(th(k,i,j)/th00-1.)
               END DO
            END IF

            DO k = 2, n1-2
               wt(k,i,j) = wt(k,i,j)+scr(k,i,j)+scr(k+1,i,j)
            END DO
         END DO
      END DO

   END SUBROUTINE boyanc
   !
   ! ----------------------------------------------------------------------
   ! Subroutine corlos:  This is the coriolis driver, its purpose is to
   ! from the coriolis accelerations for u and v and add them into the
   ! accumulated tendency arrays of ut and vt.
   !
   SUBROUTINE corlos

      USE defs, ONLY : omega
      USE grid, ONLY : a_uc, a_vc, a_ut, a_vt, nxp, nyp, nzp, u0, v0

      LOGICAL, SAVE :: initialized = .FALSE.
      REAL, SAVE    :: fcor

      INTEGER :: i, j, k

      IF (corflg) THEN
         IF (.NOT. initialized) fcor = 2.*omega*sin(cntlat*0.01745329)
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = 2, nzp
                  a_ut(k,i,j) = a_ut(k,i,j) - fcor*(v0(k)-0.25*                   &
                                (a_vc(k,i,j)+a_vc(k,i+1,j)+a_vc(k,i,j-1)+a_vc(k,i+1,j-1)))
                  a_vt(k,i,j) = a_vt(k,i,j) + fcor*(u0(k)-0.25*                   &
                                (a_uc(k,i,j)+a_uc(k,i-1,j)+a_uc(k,i,j+1)+a_uc(k,i-1,j+1)))
               END DO
            END DO
         END DO
         initialized = .TRUE.
      END IF

   END SUBROUTINE corlos
   !
   ! ----------------------------------------------------------------------
   ! Subroutine sponge: does the rayleigh friction for the momentum terms,
   ! and newtonian damping of thermal term the damping is accumulated with the
   ! other tendencies
   !
   SUBROUTINE sponge (isponge)

      USE grid, ONLY : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
                       nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00

      INTEGER, INTENT (in) :: isponge

      INTEGER :: i, j, k, kk

      IF (maxval(spng_tfct) > epsilon(1.) .AND. nfpt > 1) THEN
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = nzp-nfpt, nzp-1
                  kk = k+1-(nzp-nfpt)
                  IF (isponge == 0) THEN
                     a_tt%data(k,i,j) = a_tt%data(k,i,j) - spng_tfct(kk)*                   &
                                        (a_tp%data(k,i,j)-th0(k)+th00)
                  ELSE
                     a_ut(k,i,j) = a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-u0(k))
                     a_vt(k,i,j) = a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-v0(k))
                     a_wt(k,i,j) = a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                  END IF
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE sponge

   !
   ! ---------------------------------------------------------------------
   ! SALSA_diagnostics: Update properties for the current timestep:
   !                    E.g. if enough water has evaporated from droplets,
   !                    deplete the cloud droplet bins and move CCN material
   !                    back to the aerosol regime.
   !                    In addition, update the diagnostic scalars for total grid-cell
   !                    liquid water contents.
   !
   ! Juha Tonttila, FMI, 2014
   ! Tomi Raatikainen, FMI, 2016

   SUBROUTINE SALSA_diagnostics
      USE grid, ONLY : nxp,nyp,nzp,    &
                       a_naerop,a_maerop,a_ncloudp,a_mcloudp,a_nprecpp,a_mprecpp,a_gaerop, &
                       a_rc, a_srp,a_snrp, binMixrat, prtcl,   &
                       a_rh, a_temp, a_ri,a_srs,a_snrs,a_rhi,                                      &
                       a_nicep,a_micep,a_nsnowp,a_msnowp
      USE mo_submctl, ONLY : nbins,ncld,nprc,ica,fca,icb,fcb,ira,fra,              &
                             in1a,fn2a,fn2b,                        &
                             nice,nsnw,iia,fia,iib,fib,isa,fsa,                    &
                             rhosu,rhowa,rhoic,rhosn,      &
                             msu,moc,mno,mnh,mss,mwa,avog,pi6,                     &
                             surfw0,surfi0, rg, nlim, prlim, pi, &
                             lscndgas
      USE class_ComponentIndex, ONLY : GetIndex, GetNcomp, IsUsed


      IMPLICIT NONE

      INTEGER :: i,j,k,bc,ba,s,sc,sa,str,end,nc,c,nn,iba

      REAL :: zvol,zvola,zvolnew
      REAL, PARAMETER :: rempty = 1.e-10
      REAL :: zdh2o,zddry
      REAL :: ns, bb, aa ! Number of moles, Raoult effect, Kelvin effect; For calculating the critical radius
      REAL :: cdcld(nzp,nxp,nyp,ncld),cdprc(nzp,nxp,nyp,nprc),  & ! Critical diameter for cloud droplets and precipitation
         cdice(nzp,nxp,nyp,nice),cdsnw(nzp,nxp,nyp,nsnw)   ! Critical diameter for cloud droplets and precipitation
      REAL :: vsum

      ! Remove negative values
      a_naerop%data  = MAX(0.,a_naerop%data)
      a_ncloudp%data = MAX(0.,a_ncloudp%data)
      a_nprecpp%data = MAX(0.,a_nprecpp%data)
      a_maerop%data  = MAX(0.,a_maerop%data)
      a_mcloudp%data = MAX(0.,a_mcloudp%data)
      a_mprecpp%data = MAX(0.,a_mprecpp%data)

      a_nicep%data  = MAX(0.,a_nicep%data)
      a_nsnowp%data = MAX(0.,a_nsnowp%data)
      a_micep%data  = MAX(0.,a_micep%data)
      a_msnowp%data = MAX(0.,a_msnowp%data)

      nn = GetNcomp(prtcl)+1 ! total number of species

      ! Critical radius for cloud droplets and precipitation
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 1, nzp
               ! Aerosols
               DO c = 1, nbins
                  vsum = 0.
                  DO s = 1, nn
                     vsum = vsum + a_maerop%data(k,i,j,(s-1)*nbins+c)
                  END DO
                  IF (a_naerop%data(k,i,j,c) > 0. .AND. vsum == 0.) THEN
                     a_naerop%data(k,i,j,c) = 0.
                     DO s = 1, nn
                        a_maerop%data(k,i,j,(s-1)*nbins+c) = 0.
                     END DO
                  END IF
               END DO

               ! Clouds
               DO c = 1, ncld
                  vsum = 0.
                  DO s = 1, nn
                     vsum = vsum + a_mcloudp%data(k,i,j,(s-1)*ncld+c)
                  END DO
                  IF (a_ncloudp%data(k,i,j,c) > 0. .AND. vsum == 0.) THEN
                     a_ncloudp%data(k,i,j,c) = 0.
                     DO s = 1, nn
                        a_mcloudp%data(k,i,j,(s-1)*ncld+c) = 0.
                     END DO
                  END IF

                  ! Critical radius -----------------
                  IF (a_ncloudp%data(k,i,j,c) > nlim .AND. a_rh%data(k,i,j) < 0.999) THEN
                     ! Moles of solute
                     ns = 0.
                     IF (IsUsed(prtcl,'SO4')) THEN
                        s = GetIndex(prtcl,'SO4')
                        str = (s-1)*ncld + c
                        ns = ns + 3.*a_mcloudp%data(k,i,j,str)/msu
                     END IF
                     IF (IsUsed(prtcl,'OC')) THEN
                        s = GetIndex(prtcl,'OC')
                        str = (s-1)*ncld + c
                        ns = ns + a_mcloudp%data(k,i,j,str)/moc
                     END IF
                     IF (IsUsed(prtcl,'NO')) THEN
                        s = GetIndex(prtcl,'NO')
                        str = (s-1)*ncld + c
                        ns = ns + a_mcloudp%data(k,i,j,str)/mno
                     END IF
                     IF (IsUsed(prtcl,'NH')) THEN
                        s = GetIndex(prtcl,'NH')
                        str = (s-1)*ncld + c
                        ns = ns + a_mcloudp%data(k,i,j,str)/mnh
                     END IF
                     IF (IsUsed(prtcl,'SS')) THEN
                        s = GetIndex(prtcl,'SS')
                        str = (s-1)*ncld + c
                        ns = ns + 2.*a_mcloudp%data(k,i,j,str)/mss
                     END IF
                     ns = ns/a_ncloudp%data(k,i,j,c)

                     bb = 6.*mwa*ns/(pi*rhowa)
                     aa = 4.*mwa*surfw0/(rg*rhowa*a_temp(k,i,j))
                     cdcld(k,i,j,c) = SQRT(3.*bb/aa)
                  ELSE
                     cdcld(k,i,j,c) = rempty
                  END IF ! nlim
                  ! -----------------------------------

               END DO ! ncld

               ! Precipitation
               DO c = 1, nprc
                  IF (a_nprecpp%data(k,i,j,c) > 0. .AND. a_mprecpp%data(k,i,j,(nn-1)*nprc+c) == 0.) THEN
                     a_nprecpp%data(k,i,j,c) = 0.
                     DO s = 1, nn
                        a_mprecpp%data(k,i,j,(s-1)*nprc+c) = 0.
                     END DO
                  END IF

                  ! Critical radius -----------------
                  IF (a_nprecpp%data(k,i,j,c) > prlim .AND. a_rh%data(k,i,j) < 0.999) THEN
                     ! Moles of solute
                     ns = 0.
                     IF (IsUsed(prtcl,'SO4')) THEN
                        s = GetIndex(prtcl,'SO4')
                        str = (s-1)*nprc + c
                        ns = ns + 3.*a_mprecpp%data(k,i,j,str)/msu
                     END IF
                     IF (IsUsed(prtcl,'OC')) THEN
                        s = GetIndex(prtcl,'OC')
                        str = (s-1)*nprc + c
                        ns = ns + a_mprecpp%data(k,i,j,str)/moc
                     END IF
                     IF (IsUsed(prtcl,'NO')) THEN
                        s = GetIndex(prtcl,'NO')
                        str = (s-1)*nprc + c
                        ns = ns + a_mprecpp%data(k,i,j,str)/mno
                     END IF
                     IF (IsUsed(prtcl,'NH')) THEN
                        s = GetIndex(prtcl,'NH')
                        str = (s-1)*nprc + c
                        ns = ns + a_mprecpp%data(k,i,j,str)/mnh
                     END IF
                     IF (IsUsed(prtcl,'SS')) THEN
                        s = GetIndex(prtcl,'SS')
                        str = (s-1)*nprc + c
                        ns = ns + a_mprecpp%data(k,i,j,str)/mss
                     END IF
                     ns = ns/a_nprecpp%data(k,i,j,c)

                     bb = 6.*mwa*ns/(pi*rhowa)
                     aa = 4.*mwa*surfw0/(rg*rhowa*a_temp(k,i,j))
                     cdprc(k,i,j,c) = SQRT(3.*bb/aa)
                  ELSE
                     cdprc(k,i,j,c) = rempty
                  END IF !prlim
                  ! -----------------------------------

               END DO ! nprc

               ! Ice
               DO c = 1, nice
                  vsum = 0.
                  DO s = 1, nn
                     vsum = vsum + a_micep%data(k,i,j,(s-1)*nice+c)
                  END DO
                  IF (a_nicep%data(k,i,j,c) > 0. .AND. vsum == 0.) THEN
                     a_nicep%data(k,i,j,c) = 0.
                     DO s = 1, nn
                        a_micep%data(k,i,j,(s-1)*nice+c) = 0.
                     END DO
                  END IF

                  ! Critical radius -----------------
                  IF (a_nicep%data(k,i,j,c) > prlim .AND. a_rhi%data(k,i,j) < 0.999) THEN
                     ! Moles of solute
                     ns = 0.
                     IF (IsUsed(prtcl,'SO4')) THEN
                        s = GetIndex(prtcl,'SO4')
                        str = (s-1)*nice + c
                        ns = ns + 3.*a_micep%data(k,i,j,str)/msu
                     END IF
                     IF (IsUsed(prtcl,'OC')) THEN
                        s = GetIndex(prtcl,'OC')
                        str = (s-1)*nice + c
                        ns = ns + a_micep%data(k,i,j,str)/moc
                     END IF
                     IF (IsUsed(prtcl,'NO')) THEN
                        s = GetIndex(prtcl,'NO')
                        str = (s-1)*nice + c
                        ns = ns + a_micep%data(k,i,j,str)/mno
                     END IF
                     IF (IsUsed(prtcl,'NH')) THEN
                        s = GetIndex(prtcl,'NH')
                        str = (s-1)*nice + c
                        ns = ns + a_micep%data(k,i,j,str)/mnh
                     END IF
                     IF (IsUsed(prtcl,'SS')) THEN
                        s = GetIndex(prtcl,'SS')
                        str = (s-1)*nice + c
                        ns = ns + 2.*a_micep%data(k,i,j,str)/mss
                     END IF
                     ns = ns/a_nicep%data(k,i,j,c)

                     bb = 3.*mwa*ns/(4.*pi*rhoic)
                     aa = 4.*mwa*surfi0/(rg*rhoic*a_temp(k,i,j))
                     cdice(k,i,j,c) = max(rempty,SQRT(3.*bb/aa))
                  ELSE
                     cdice(k,i,j,c) = rempty
                  END IF ! nlim
                  ! -----------------------------------

               END DO ! nice

               ! Snow
               DO c = 1, nsnw
                  IF (a_nsnowp%data(k,i,j,c) > 0. .AND. a_msnowp%data(k,i,j,(nn-1)*nsnw+c) == 0.) THEN
                     a_nsnowp%data(k,i,j,c) = 0.
                     DO s = 1, nn
                        a_msnowp%data(k,i,j,(s-1)*nsnw+c) = 0.
                     END DO
                  END IF

                  ! Critical radius -----------------
                  IF (a_nsnowp%data(k,i,j,c) > prlim .AND. a_rhi%data(k,i,j) < 0.999) THEN
                     ! Moles of solute
                     ns = 0.
                     IF (IsUsed(prtcl,'SO4')) THEN
                        s = GetIndex(prtcl,'SO4')
                        str = (s-1)*nsnw + c
                        ns = ns + 3.*a_msnowp%data(k,i,j,str)/msu
                     END IF
                     IF (IsUsed(prtcl,'OC')) THEN
                        s = GetIndex(prtcl,'OC')
                        str = (s-1)*nsnw + c
                        ns = ns + a_msnowp%data(k,i,j,str)/moc
                     END IF
                     IF (IsUsed(prtcl,'NO')) THEN
                        s = GetIndex(prtcl,'NO')
                        str = (s-1)*nsnw + c
                        ns = ns + a_msnowp%data(k,i,j,str)/mno
                     END IF
                     IF (IsUsed(prtcl,'NH')) THEN
                        s = GetIndex(prtcl,'NH')
                        str = (s-1)*nsnw + c
                        ns = ns + a_msnowp%data(k,i,j,str)/mnh
                     END IF
                     IF (IsUsed(prtcl,'SS')) THEN
                        s = GetIndex(prtcl,'SS')
                        str = (s-1)*nsnw + c
                        ns = ns + a_msnowp%data(k,i,j,str)/mss
                     END IF
                     ns = ns/a_nsnowp%data(k,i,j,c)

                     bb = 3.*mwa*ns/(4.*pi*rhoic)
                     aa = 4.*mwa*surfi0/(rg*rhoic*a_temp(k,i,j))
                     cdsnw(k,i,j,c) = max(rempty,SQRT(3.*bb/aa))
                  ELSE
                     cdsnw(k,i,j,c) = rempty
                  END IF !prlim
                  ! -----------------------------------

               END DO ! nsnw

            END DO !k
         END DO !i
      END DO !j

      ! Ghost species, particle radiae and diagnostic stracers
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 1, nzp

               !!!!!!!!!!!!!!!!!!!!!!!
               ! Ghost species
               !!!!!!!!!!!!!!!!!!!!!!!

               ! Loop over cloud droplet bins
               DO bc = ica%cur, fcb%cur

                  IF ( a_ncloudp%data(k,i,j,bc) > nlim .AND. a_rh%data(k,i,j) < 0.999) THEN

                     CALL binMixrat('cloud','wet',bc,i,j,k,zvol)
                     zvol = zvol/rhowa
                     zdh2o = (zvol/a_ncloudp%data(k,i,j,bc)/pi6)**((1./3.))

                     ! Loose the droplets if smaller than the critical size
                     IF ( zdh2o < MAX(0.2*cdcld(k,i,j,bc),2.e-6) ) THEN
                        IF (bc <= fca%cur) THEN
                           ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                        ELSE
                           ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                        END IF
                        ! Move the number of particles from cloud to aerosol bins
                        a_naerop%data(k,i,j,ba) = a_naerop%data(k,i,j,ba) + a_ncloudp%data(k,i,j,bc)
                        a_ncloudp%data(k,i,j,bc) = 0.

                        ! Move ccn material back to aerosol regime (including water)
                        DO s = 1, nn
                           sc = (s-1)*ncld + bc
                           sa = (s-1)*nbins + ba
                           a_maerop%data(k,i,j,sa) = a_maerop%data(k,i,j,sa) + a_mcloudp%data(k,i,j,sc)
                           a_mcloudp%data(k,i,j,sc) = 0.
                        END DO

                     END IF ! critical radius

                  END IF  ! blim

               END DO ! bc

               ! Loop over precipitation bins
               DO bc = ira, fra

                  IF ( a_nprecpp%data(k,i,j,bc) > prlim .AND. a_rh%data(k,i,j) < 0.999 ) THEN

                     CALL binMixrat('precp','wet',bc,i,j,k,zvol)
                     zvol = zvol/rhowa
                     zdh2o = (zvol/a_nprecpp%data(k,i,j,bc)/pi6)**((1./3.))

                     ! Loose the droplets if smaller than critical radius
                     IF ( zdh2o < MAX(0.02*cdprc(k,i,j,bc),2.e-6)  ) THEN

                        ! Move evaporating rain drops to a soluble aerosol bin with
                        ! the closest match in dry particle mass. Ain't perfect but
                        ! the bin update subroutine in SALSA will take care of the rest.
                        CALL binMixrat('precp','dry',bc,i,j,k,zvol)

                        zvol = zvol/a_nprecpp%data(k,i,j,bc)

                        ba = 0
                        zvola = -1.
                        DO iba = in1a, fn2a
                           IF (a_naerop%data(k,i,j,iba) > nlim) THEN
                              CALL binMixrat('aerosol','dry',iba,i,j,k,zvolnew)
                              zvolnew = zvolnew/a_naerop%data(k,i,j,iba)
                              IF (abs(zvolnew-zvol) < abs(zvola-zvol)) THEN
                                 ! New closest match
                                 ba = iba
                                 zvola = zvolnew
                              END IF
                           END IF
                        END DO
                        IF (ba == 0) STOP 'FAIL: no sink for evaporating rain drops'

                        ! Move the number of particles from cloud to aerosol bins
                        a_naerop%data(k,i,j,ba) = a_naerop%data(k,i,j,ba) + a_nprecpp%data(k,i,j,bc)
                        a_nprecpp%data(k,i,j,bc) = 0.

                        ! Move ccn material back to aerosol regime (including water)
                        DO s = 1, nn
                           sc = (s-1)*nprc + bc
                           sa = (s-1)*nbins + ba
                           a_maerop%data(k,i,j,sa) = a_maerop%data(k,i,j,sa) + a_mprecpp%data(k,i,j,sc)
                           a_mprecpp%data(k,i,j,sc) = 0.
                        END DO

                     END IF ! Critical radius

                  END IF ! prlim

               END DO ! bc

               ! Loop over ice bins
               DO bc = iia%cur, fib%cur

                  IF ( a_nicep%data(k,i,j,bc) > prlim .AND. a_rhi%data(k,i,j) < 0.999 ) THEN

                     CALL binMixrat('ice','wet',bc,i,j,k,zvol)
                     zvol = zvol/rhoic
                     zdh2o = (zvol/a_nicep%data(k,i,j,bc)/pi6)**((1./3.))

                     ! Loose the droplets if smaller than the critical size !! ice'n'snow
                     IF ( zdh2o < MAX(0.2*cdice(k,i,j,bc),2.e-6) ) THEN
                        IF (bc <= fia%cur) THEN
                           ba = iia%par + (bc-iia%cur) ! Index for parallel aerosol bin
                        ELSE
                           ba = iib%par + (bc-iib%cur) ! Index for parallel aerosol bin
                        END IF

                        ! Move the number of particles from cloud to aerosol bins
                        a_naerop%data(k,i,j,ba) = a_naerop%data(k,i,j,ba) + a_nicep%data(k,i,j,bc)
                        a_nicep%data(k,i,j,bc) = 0.

                        ! Move ccn material back to aerosol regime (including water)
                        DO s = 1, nn
                           sc = (s-1)*nice + bc
                           sa = (s-1)*nbins + ba
                           a_maerop%data(k,i,j,sa) = a_maerop%data(k,i,j,sa) + a_micep%data(k,i,j,sc)
                           a_micep%data(k,i,j,sc) = 0.
                        END DO

                     END IF ! critical radius

                  END IF  ! prlim

               END DO ! bc

               ! Loop over snow bins
               DO bc = isa, fsa

                  IF ( a_nsnowp%data(k,i,j,bc) > prlim .AND. a_rhi%data(k,i,j) < 0.999 ) THEN

                     CALL binMixrat('snow','wet',bc,i,j,k,zvol)
                     zvol = zvol/rhosn
                     zdh2o = (zvol/a_nsnowp%data(k,i,j,bc)/pi6)**((1./3.))

                     ! Loose the droplets if smaller than critical radius !! a_rhi ice'n'snow
                     IF ( zdh2o < MAX(0.02*cdsnw(k,i,j,bc),2.e-6) ) THEN

                        ! Move evaporating snow drops to a soluble aerosol bin with
                        ! the closest match in dry particle mass. Ain't perfect but
                        ! the bin update subroutine in SALSA will take care of the rest.
                        CALL binMixrat('snow','dry',bc,i,j,k,zvol)
                        zvol = zvol/a_nsnowp%data(k,i,j,bc)

                        ba = 0
                        zvola = -1.
                        DO iba = in1a, fn2a
                           IF (a_naerop%data(k,i,j,iba) > nlim) THEN
                              CALL binMixrat('aerosol','dry',iba,i,j,k,zvolnew)
                              zvolnew = zvolnew/a_naerop%data(k,i,j,iba)
                              IF (abs(zvolnew-zvol) < abs(zvola-zvol)) THEN
                                 ! New closest match
                                 ba = iba
                                 zvola = zvolnew
                              END IF
                           END IF
                        END DO
                        IF (ba == 0) STOP 'FAIL: no sink for evaporating snow'

                        ! Move the number of particles from cloud to aerosol bins
                        a_naerop%data(k,i,j,ba) = a_naerop%data(k,i,j,ba) + a_nsnowp%data(k,i,j,bc)
                        a_nsnowp%data(k,i,j,bc) = 0.

                        ! Move ccn material back to aerosol regime (including water)
                        DO s = 1, nn
                           sc = (s-1)*nsnw + bc
                           sa = (s-1)*nbins + ba
                           a_maerop%data(k,i,j,sa) = a_maerop%data(k,i,j,sa) + a_msnowp%data(k,i,j,sc)
                           a_msnowp%data(k,i,j,sc) = 0.
                        END DO

                     END IF ! Critical radius

                  END IF ! prlim

               END DO ! bc

               ! Loop over aerosol bins
               DO ba = 1, nbins
                  IF (a_naerop%data(k,i,j,ba) > nlim) THEN
                     CALL binMixrat('aerosol','dry',ba,i,j,k,zvol)
                     zvol = zvol/rhosu

                     ! Particles smaller than 0.1 nm diameter are set to zero
                     zddry = (zvol/a_naerop%data(k,i,j,ba)/pi6)**((1./3.))
                     IF ( zddry < 1.e-10 ) THEN
                        ! Volatile species to the gas phase
                        IF (IsUsed(prtcl,'SO4') .AND. lscndgas) THEN
                           nc = GetIndex(prtcl,'SO4')
                           s = (nc-1)*nbins + ba
                           a_gaerop%data(k,i,j,1) = a_gaerop%data(k,i,j,1) + a_maerop%data(k,i,j,s) / msu * avog
                        END IF
                        IF (IsUsed(prtcl,'OC') .AND. lscndgas) THEN
                           nc = GetIndex(prtcl,'OC')
                           s = (nc-1)*nbins + ba
                           a_gaerop%data(k,i,j,5) = a_gaerop%data(k,i,j,5) + a_maerop%data(k,i,j,s) / moc * avog
                        END IF
                        IF (IsUsed(prtcl,'NO') .AND. lscndgas) THEN
                           nc = GetIndex(prtcl,'NO')
                           s = (nc-1)*nbins + ba
                           a_gaerop%data(k,i,j,2) = a_gaerop%data(k,i,j,2) + a_maerop%data(k,i,j,s) / mno * avog
                        END IF
                        IF (IsUsed(prtcl,'NH') .AND. lscndgas) THEN
                           nc = GetIndex(prtcl,'NH')
                           s = (nc-1)*nbins + ba
                           a_gaerop%data(k,i,j,3) = a_gaerop%data(k,i,j,3) + a_maerop%data(k,i,j,s) / mnh * avog
                        END IF

                        ! Mass and number to zero (insolube species and water are lost)
                        a_maerop%data(k,i,j,ba:(nn-1)*nbins+ba:nbins) = 0.
                        a_naerop%data(k,i,j,ba) = 0.
                        zvol = 0.
                     END IF
                  END IF
               END DO

            END DO   ! k
         END DO   ! i
      END DO   ! j

      !!!!!!!!!!!!!!!!!!!!!!!
      ! Update diagnostic tracers
      !!!!!!!!!!!!!!!!!!!!!!!

      ! Liquid water content
      nc = GetIndex(prtcl,'H2O')
      ! Aerosols, regimes a and b
      str = (nc-1)*nbins + in1a
      end = (nc-1)*nbins + fn2b
      a_rc%data = SUM(a_maerop%data(:,:,:,str:end),DIM=4)
      ! Clouds, regime a and b
      str = (nc-1)*ncld+ica%cur
      end = (nc-1)*ncld+fcb%cur
      a_rc%data = a_rc%data + SUM(a_mcloudp%data(:,:,:,str:end),DIM=4)
      ! Precipitation
      str = (nc-1)*nprc+ira
      end = (nc-1)*nprc+fra
      a_srp%data = SUM(a_mprecpp%data(:,:,:,str:end),DIM=4)
      a_snrp%data = SUM(a_nprecpp%data(:,:,:,ira:fra),DIM=4)

      ! ice, regimes a and b
      str = (nc-1)*nice+iia%cur
      end = (nc-1)*nice+fib%cur
      a_ri%data = SUM(a_micep%data(:,:,:,str:end),DIM=4)
      ! Snow
      str = (nc-1)*nsnw+isa
      end = (nc-1)*nsnw+fsa
      a_srs%data = SUM(a_msnowp%data(:,:,:,str:end),DIM=4)
      a_snrs%data = SUM(a_nsnowp%data(:,:,:,isa:fsa),DIM=4)

   END SUBROUTINE SALSA_diagnostics


END MODULE step
