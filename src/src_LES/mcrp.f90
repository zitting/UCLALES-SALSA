!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
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
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE mcrp

   USE defs, ONLY : alvl, alvi, rowt, pi, Rm, cp, kb, g, vonk
   USE grid, ONLY : dtlt, dzt, nxp, nyp, nzp,a_pexnr, a_rp, a_tp, th00, CCN,                 &
                    dn0, pi0, a_rt, a_tt, a_rpp, a_rpt, a_npp, a_npt, a_rv, a_rc, a_theta,   &
                    a_press, a_temp, a_rsl, precip, a_dn, a_ustar,                           &
                    a_naerop,  a_naerot,  a_maerop,  a_maerot,                               &
                    a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,                              &
                    a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt, mc_ApVdom,                   &
                    a_nicep,   a_nicet,   a_micep,   a_micet,                                &
                    a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,                               &
                    snowin,    prtcl, calc_wet_radius
   USE thrm, ONLY : thermo
   USE stat, ONLY : sflg, updtst, acc_removal, mcflg, acc_massbudged, cs_rem_set
   IMPLICIT NONE

   LOGICAL, PARAMETER :: droplet_sedim = .FALSE., khairoutdinov = .FALSE.

   LOGICAL :: sed_aero = .TRUE.,  &
              sed_cloud = .TRUE., &
              sed_precp = .TRUE., &
              sed_ice = .TRUE.,   &
              sed_snow = .TRUE.
   !
   ! drop sizes definition is based on vanZanten (2005)
   ! cloud droplets' diameter: 2-50 e-6 m
   ! drizzle drops' diameter: 50-1000 e-6 m
   !
   REAL, PARAMETER :: eps0 = 1e-20       ! small number
   REAL, PARAMETER :: eps1 = 1e-9        ! small number
   REAL, PARAMETER :: rho_0 = 1.21       ! air density at surface

   REAL, PARAMETER :: D_min = 2.e-6      ! minimum diameter of cloud droplets
   REAL, PARAMETER :: D_bnd = 80.e-6     ! precip/cloud diameter boundary
   REAL, PARAMETER :: D_max = 1.e-3      ! maximum diameter of prcp drops

   REAL, PARAMETER :: X_min = (D_min**3)*rowt*pi/6. ! min cld mass
   REAL, PARAMETER :: X_bnd = (D_bnd**3)*rowt*pi/6. ! precip/cld bound mass
   REAL, PARAMETER :: X_max = (D_max**3)*rowt*pi/6. ! max prcp mass

   REAL, PARAMETER :: prw = pi * rowt / 6.

CONTAINS

   !
   ! ---------------------------------------------------------------------
   ! MICRO: sets up call to microphysics
   !
   SUBROUTINE micro(level)
      USE class_componentIndex, ONLY : GetNcomp
      INTEGER, INTENT (in) :: level
      INTEGER :: nn

      SELECT CASE (level)
         CASE(2)
            IF (droplet_sedim)  &
               CALL sedim_cd(nzp,nxp,nyp,a_theta,a_temp,a_rc,precip,a_rt,a_tt)
         CASE(3)
            CALL mcrph(nzp,nxp,nyp,dn0,a_theta,a_temp,a_rv,a_rsl,a_rc,a_rpp,   &
                       a_npp,precip,a_rt,a_tt,a_rpt,a_npt)
         CASE(4,5)
            IF (level < 5) THEN
               sed_ice = .FALSE.; sed_snow = .FALSE.
            END IF
            nn = GetNcomp(prtcl)+1
            CALL sedim_SALSA(nzp,nxp,nyp,nn,dtlt, a_temp, a_theta,                &
                             a_naerop,  a_naerot,  a_maerop,  a_maerot,           &
                             a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,          &
                             a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,          &
                             a_nicep,   a_nicet,   a_micep,   a_micet,            &
                             a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,           &
                             a_ustar, precip, snowin, a_tt                )
      END SELECT

   END SUBROUTINE micro

   !
   ! ---------------------------------------------------------------------
   ! MCRPH: calls microphysical parameterization
   !
   SUBROUTINE mcrph(n1,n2,n3,dn0,th,tk,rv,rs,rc,rp,np,rrate,         &
                    rtt,tlt,rpt,npt)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, DIMENSION(n1,n2,n3), INTENT (in)    :: th, tk, rv, rs
      REAL, DIMENSION(n1)      , INTENT (in)    :: dn0
      REAL, DIMENSION(n1,n2,n3), INTENT (inout) :: rc, rtt, tlt, rpt, npt, np, rp
      REAL, INTENT (out)                        :: rrate(n1,n2,n3)

      INTEGER :: i, j, k

      !
      ! Microphysics following Seifert Beheng (2001, 2005)
      ! note that the order below is important as the rc array is
      ! redefined in cld_dgn below and is assumed to be cloud water
      ! after that and total condensate priort to that
      !

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               rp(k,i,j) = max(0., rp(k,i,j))
               np(k,i,j) = max(min(rp(k,i,j)/X_bnd,np(k,i,j)),rp(k,i,j)/X_max)
            END DO
         END DO
      END DO

      CALL wtr_dff_SB(n1,n2,n3,dn0,rp,np,rc,rs,rv,tk,rpt,npt)

      CALL auto_SB(n1,n2,n3,dn0,rc,rp,rpt,npt)

      CALL accr_SB(n1,n2,n3,dn0,rc,rp,np,rpt,npt)

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               rp(k,i,j)  = rp(k,i,j) + max(-rp(k,i,j)/dtlt,rpt(k,i,j))*dtlt
               np(k,i,j)  = np(k,i,j) + max(-np(k,i,j)/dtlt,npt(k,i,j))*dtlt
               rpt(k,i,j) = 0.
               npt(k,i,j) = 0.
               rp(k,i,j)  = max(0., rp(k,i,j))
               np(k,i,j)  = max(min(rp(k,i,j)/X_bnd,np(k,i,j)),rp(k,i,j)/X_max)
            END DO
         END DO
      END DO

      CALL sedim_rd(n1,n2,n3,dtlt,dn0,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)

      IF (droplet_sedim) CALL sedim_cd(n1,n2,n3,th,tk,rc,rrate,rtt,tlt)

   END SUBROUTINE mcrph
   !
   ! ---------------------------------------------------------------------
   ! WTR_DFF_SB: calculates the evolution of the both number- and
   ! mass mixing ratio large drops due to evaporation in the absence of
   ! cloud water.
   !
   SUBROUTINE wtr_dff_SB(n1,n2,n3,dn0,rp,np,rl,rs,rv,tk,rpt,npt)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: tk(n1,n2,n3), np(n1,n2,n3), rp(n1,n2,n3),        &
                              rs(n1,n2,n3),rv(n1,n2,n3), dn0(n1)
      REAL, INTENT (inout) :: rpt(n1,n2,n3), npt(n1,n2,n3), rl(n1,n2,n3)

      REAL, PARAMETER    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
      REAL, PARAMETER    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]

      INTEGER             :: i, j, k
      REAL                :: Xp, Dp, G, S, cerpt, cenpt, xnpts
      REAL, DIMENSION(n1) :: v1

      IF(sflg) THEN
         xnpts = 1./((n3-4)*(n2-4))
         DO k = 1, n1
            v1(k) = 0.
         END DO
      END IF

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1
               IF (rp(k,i,j) > rl(k,i,j)) THEN
                  Xp = rp(k,i,j)/ (np(k,i,j)+eps0)
                  Xp = MIN(MAX(Xp,X_bnd),X_max)
                  Dp = ( Xp / prw )**(1./3.)

                  G = 1. / (1. / (dn0(k)*rs(k,i,j)*Dv) + &
                      alvl*(alvl/(Rm*tk(k,i,j))-1.) / (Kt*tk(k,i,j)))
                  S = rv(k,i,j)/rs(k,i,j) - 1.

                  IF (S < 0) THEN
                     cerpt = 2. * pi * Dp * G * S * np(k,i,j)
                     cenpt = cerpt * np(k,i,j) / rp(k,i,j)
                     rpt(k,i,j) = rpt(k,i,j) + cerpt
                     npt(k,i,j) = npt(k,i,j) + cenpt
                     IF (sflg) v1(k) = v1(k) + cerpt * xnpts
                  END IF
               END IF
               rl(k,i,j) = max(0.,rl(k,i,j) - rp(k,i,j))
            END DO
         END DO
      END DO

      IF (sflg) CALL updtst(n1,'prc',2,v1,1)

   END SUBROUTINE wtr_dff_SB
   !
   ! ---------------------------------------------------------------------
   ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
   ! drizzle drops due to autoconversion. The autoconversion rate assumes
   ! f(x)=A*x**(nu_c)*exp(-Bx), an exponential in drop MASS x. It can
   ! be reformulated for f(x)=A*x**(nu_c)*exp(-Bx**(mu)), where formu=1/3
   ! one would get a gamma dist in drop diam -> faster rain formation.
   !
   SUBROUTINE auto_SB(n1,n2,n3,dn0,rc,rp,rpt,npt)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: dn0(n1), rc(n1,n2,n3), rp(n1,n2,n3)
      REAL, INTENT (inout) :: rpt(n1,n2,n3), npt(n1,n2,n3)

      REAL, PARAMETER :: nu_c = 0            ! width parameter of cloud DSD
      REAL, PARAMETER :: k_c  = 9.44e+9      ! Long-Kernel
      REAL, PARAMETER :: k_1  = 6.e+2        ! Parameter for phi function
      REAL, PARAMETER :: k_2  = 0.68         ! Parameter for phi function

      REAL, PARAMETER :: Cau = 4.1e-15 ! autoconv. coefficient in KK param.
      REAL, PARAMETER :: Eau = 5.67    ! autoconv. exponent in KK param.
      REAL, PARAMETER :: mmt = 1.e+6   ! transformation from m to \mu m

      INTEGER :: i, j, k
      REAL    :: k_au, Xc, Dc, au, tau, phi

      k_au = k_c / (20.*X_bnd) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               Xc = rc(k,i,j)/(CCN+eps0)
               IF (Xc > 0.) THEN
                  Xc = MIN(MAX(Xc,X_min),X_bnd)
                  au = k_au * dn0(k) * rc(k,i,j)**2 * Xc**2
                  !
                  ! small threshold that should not influence the result
                  !
                  IF (rc(k,i,j) > 1.e-6) THEN
                     tau = 1.0-rc(k,i,j)/(rc(k,i,j)+rp(k,i,j)+eps0)
                     tau = MIN(MAX(tau,eps0),0.9)
                     phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
                     au  = au * (1.0 + phi/(1.0 - tau)**2)
                  END IF
                  !
                  ! Khairoutdinov and Kogan
                  !
                  IF (khairoutdinov) THEN
                     Dc = ( Xc / prw )**(1./3.)
                     au = Cau * (Dc * mmt / 2.)**Eau
                  END IF

                  rpt(k,i,j) = rpt(k,i,j) + au
                  npt(k,i,j) = npt(k,i,j) + au/X_bnd
                  !
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE auto_SB
   !
   ! ---------------------------------------------------------------------
   ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
   ! and self collection following Seifert & Beheng (2001).  Included is
   ! an alternative formulation for accretion only, following
   ! Khairoutdinov and Kogan
   !
   SUBROUTINE accr_SB(n1,n2,n3,dn0,rc,rp,np,rpt,npt)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: rc(n1,n2,n3), rp(n1,n2,n3), np(n1,n2,n3), dn0(n1)
      REAL, INTENT (inout) :: rpt(n1,n2,n3),npt(n1,n2,n3)

      REAL, PARAMETER :: k_r = 5.78
      REAL, PARAMETER :: k_1 = 5.e-4
      REAL, PARAMETER :: Cac = 67.     ! accretion coefficient in KK param.
      REAL, PARAMETER :: Eac = 1.15    ! accretion exponent in KK param.

      INTEGER :: i, j, k
      REAL    :: tau, phi, ac, sc

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               IF (rc(k,i,j) > 0. .AND. rp(k,i,j) > 0.) THEN
                  tau = 1.0-rc(k,i,j)/(rc(k,i,j)+rp(k,i,j)+eps0)
                  tau = MIN(MAX(tau,eps0),1.)
                  phi = (tau/(tau+k_1))**4
                  ac  = k_r * rc(k,i,j) * rp(k,i,j) * phi * sqrt(rho_0*dn0(k))
                  !
                  ! Khairoutdinov and Kogan
                  !
                  !ac = Cac * (rc(k,i,j) * rp(k,i,j))**Eac
                  !
                  rpt(k,i,j) = rpt(k,i,j) + ac

               END IF
               sc = k_r * np(k,i,j) * rp(k,i,j) * sqrt(rho_0*dn0(k))
               npt(k,i,j) = npt(k,i,j) - sc
            END DO
         END DO
      END DO

   END SUBROUTINE accr_SB
   !
   ! ---------------------------------------------------------------------
   ! SEDIM_RD: calculates the sedimentation of the rain drops and its
   ! effect on the evolution of theta_l and r_t.  This is expressed in
   ! terms of Dp the mean diameter, not the mass weighted mean diameter
   ! as is used elsewhere.  This is just 1/lambda in the exponential
   ! distribution
   !
   SUBROUTINE sedim_rd(n1,n2,n3,dt,dn0,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)

      INTEGER, INTENT (in)                      :: n1,n2,n3
      REAL, INTENT (in)                         :: dt
      REAL, INTENT (in),    DIMENSION(n1)       :: dn0
      REAL, INTENT (in),    DIMENSION(n1,n2,n3) :: rp, np, th, tk
      REAL, INTENT (out),   DIMENSION(n1,n2,n3) :: rrate
      REAL, INTENT (inout), DIMENSION(n1,n2,n3) :: rtt, tlt, rpt, npt

      REAL, PARAMETER :: a2 = 9.65       ! in SI [m/s]
      REAL, PARAMETER :: c2 = 6e2        ! in SI [1/m]
      REAL, PARAMETER :: Dv = 25.0e-6    ! in SI [m/s]
      REAL, PARAMETER :: cmur1 = 10.0    ! mu-Dm-relation for rain following
      REAL, PARAMETER :: cmur2 = 1.20e+3 ! Milbrandt&Yau 2005, JAS, but with
      REAL, PARAMETER :: cmur3 = 1.5e-3  ! revised constants
      REAL, PARAMETER :: aq = 6.0e3
      REAL, PARAMETER :: bq = -0.2
      REAL, PARAMETER :: an = 3.5e3
      REAL, PARAMETER :: bn = -0.1


      INTEGER :: i, j, k, kp1, kk, km1
      REAL    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz, xnpts
      REAL, DIMENSION(n1) :: nslope,rslope,dn,dr, rfl, nfl, vn, vr, cn, cr, v1

      IF(sflg) THEN
         xnpts = 1./((n3-4)*(n2-4))
         DO k = 1, n1
            v1(k) = 0.
         END DO
      END IF

      b2 = a2*exp(c2*Dv)

      DO j = 3, n3-2
         DO i = 3, n2-2

            nfl(n1) = 0.
            rfl(n1) = 0.

            DO k = n1-1, 2, -1
               Xp = rp(k,i,j) / (np(k,i,j)+eps0)
               Xp = MIN(MAX(Xp,X_bnd),X_max)
               !
               ! Adjust Dm and mu-Dm and Dp=1/lambda following Milbrandt & Yau
               !
               Dm = ( 6. / (rowt*pi) * Xp )**(1./3.)
               mu = cmur1*(1.+tanh(cmur2*(Dm-cmur3)))
               Dp = (Dm**3/((mu+3.)*(mu+2.)*(mu+1.)))**(1./3.)

               vn(k) = sqrt(dn0(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(1.+mu)))
               vr(k) = sqrt(dn0(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(4.+mu)))
               !
               ! Set fall speeds following Khairoutdinov and Kogan

               IF (khairoutdinov) THEN
                  vn(k) = max(0.,an * Dp + bn)
                  vr(k) = max(0.,aq * Dp + bq)
               END IF

            END DO

            DO k = 2, n1-1
               kp1 = min(k+1,n1-1)
               km1 = max(k,2)
               cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzt(k)*dt
               cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzt(k)*dt
            END DO

            !...piecewise linear method: get slopes
            DO k = n1-1, 2, -1
               dn(k) = np(k+1,i,j)-np(k,i,j)
               dr(k) = rp(k+1,i,j)-rp(k,i,j)
            END DO
            dn(1)  = dn(2)
            dn(n1) = dn(n1-1)
            dr(1)  = dr(2)
            dr(n1) = dr(n1-1)

            DO k = n1-1, 2, -1
               !...slope with monotone limiter for np
               sk = 0.5 * (dn(k-1) + dn(k))
               mini = min(np(k-1,i,j),np(k,i,j),np(k+1,i,j))
               maxi = max(np(k-1,i,j),np(k,i,j),np(k+1,i,j))
               nslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(np(k,i,j)-mini), &
                  &                              2.*(maxi-np(k,i,j)))
               !...slope with monotone limiter for rp
               sk = 0.5 * (dr(k-1) + dr(k))
               mini = min(rp(k-1,i,j),rp(k,i,j),rp(k+1,i,j))
               maxi = max(rp(k-1,i,j),rp(k,i,j),rp(k+1,i,j))
               rslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(rp(k,i,j)-mini), &
                  &                              2.*(maxi-rp(k,i,j)))
            END DO

            rfl(n1-1) = 0.
            nfl(n1-1) = 0.
            DO k = n1-2, 2, -1

               kk = k
               tot = 0.0
               zz  = 0.0
               cc  = min(1.,cn(k))
               DO WHILE (cc > 0 .AND. kk <= n1-1)
                  tot = tot + dn0(kk)*(np(kk,i,j)+nslope(kk)*(1.-cc))*cc/dzt(kk)
                  zz  = zz + 1./dzt(kk)
                  kk  = kk + 1
                  cc  = min(1.,cn(kk) - zz*dzt(kk))
               END DO
               nfl(k) = -tot /dt

               kk = k
               tot = 0.0
               zz  = 0.0
               cc  = min(1.,cr(k))
               DO WHILE (cc > 0 .AND. kk <= n1-1)
                  tot = tot + dn0(kk)*(rp(kk,i,j)+rslope(kk)*(1.-cc))*cc/dzt(kk)
                  zz  = zz + 1./dzt(kk)
                  kk  = kk + 1
                  cc  = min(1.,cr(kk) - zz*dzt(kk))
               END DO
               rfl(k) = -tot /dt

               kp1 = k+1
               flxdiv = (rfl(kp1)-rfl(k))*dzt(k)/dn0(k)
               rpt(k,i,j) =rpt(k,i,j)-flxdiv
               rtt(k,i,j) =rtt(k,i,j)-flxdiv
               tlt(k,i,j) =tlt(k,i,j)+flxdiv*(alvl/cp)*th(k,i,j)/tk(k,i,j)

               npt(k,i,j) = npt(k,i,j)-(nfl(kp1)-nfl(k))*dzt(k)/dn0(k)

               rrate(k,i,j) = -rfl(k)/dn0(k) * alvl*0.5*(dn0(k)+dn0(kp1))

               IF (sflg) v1(k) = v1(k) + rrate(k,i,j)*xnpts

            END DO
         END DO
      END DO

      IF (sflg) CALL updtst(n1,'prc',1,v1,1)

   END SUBROUTINE sedim_rd
   !
   ! ---------------------------------------------------------------------
   ! SEDIM_CD: calculates the cloud-droplet sedimentation flux and its effect
   ! on the evolution of r_t and theta_l assuming a log-normal distribution
   !
   SUBROUTINE sedim_cd(n1,n2,n3,th,tk,rc,rrate,rtt,tlt)


      INTEGER, INTENT (in):: n1,n2,n3
      REAL, INTENT (in),   DIMENSION(n1,n2,n3) :: th,tk,rc
      REAL, INTENT (out),  DIMENSION(n1,n2,n3) :: rrate
      REAL, INTENT (inout),DIMENSION(n1,n2,n3) :: rtt,tlt

      REAL, PARAMETER :: c = 1.19e8 ! Stokes fall velocity coef [m^-1 s^-1]
      REAL, PARAMETER :: sgg = 1.2  ! geometric standard dev of cloud droplets

      INTEGER :: i, j, k, kp1
      REAL    :: Dc, Xc, vc, flxdiv
      REAL    :: rfl(n1)

      !
      ! calculate the precipitation flux and its effect on r_t and theta_l
      !
      DO j = 3, n3-2
         DO i = 3, n2-2
            rfl(n1) = 0.
            DO k = n1-1, 2, -1
               Xc = rc(k,i,j) / (CCN+eps0)
               Dc = ( Xc / prw )**(1./3.)
               Dc = MIN(MAX(Dc,D_min),D_bnd)
               vc = min(c*(Dc*0.5)**2 * exp(4.5*(log(sgg))**2),1./(dzt(k)*dtlt))
               rfl(k) = - rc(k,i,j) * vc
               !
               kp1 = k+1
               flxdiv = (rfl(kp1)-rfl(k))*dzt(k)
               rtt(k,i,j) = rtt(k,i,j)-flxdiv
               tlt(k,i,j) = tlt(k,i,j)+flxdiv*(alvl/cp)*th(k,i,j)/tk(k,i,j)
               rrate(k,i,j) = -rfl(k)
            END DO
         END DO
      END DO

   END SUBROUTINE sedim_cd

   ! ---------------------------------------------------------------------
   ! SEDIM_AERO: calculates the salsa particles sedimentation and dry deposition flux  (.. Zubair) !
   !
   ! Juha: The code below is a modified version of the original one by Zubair
   !
   ! Juha: Rain is now treated completely separately (20151013)
   !
   ! Jaakko: Modified for the use of ice and snow bins

   SUBROUTINE sedim_SALSA(n1,n2,n3,n4,tstep,tk,th,          &
                          naerop, naerot, maerop, maerot,   &
                          ncloudp,ncloudt,mcloudp,mcloudt,  &
                          nprecpp,nprecpt,mprecpp,mprecpt,  &
                          nicep,  nicet,  micep,  micet,    &
                          nsnowp, nsnowt, msnowp, msnowt,   &
                          ustar,rrate, srate, tlt          )

      USE mo_submctl, ONLY : nbins, ncld, nprc,           &
                             nice,  nsnw,                 &
                             nlim,prlim,  &
                             rhowa,rhoic
      USE class_ComponentIndex, ONLY : GetIndex, GetNcomp
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3,n4
      REAL, INTENT(in) :: tstep,                     &
                          tk(n1,n2,n3),              &
                          th(n1,n2,n3),              &
                          ustar(n2,n3),              &
                          naerop(n1,n2,n3,nbins),    &
                          maerop(n1,n2,n3,n4*nbins), &
                          ncloudp(n1,n2,n3,ncld),    &
                          mcloudp(n1,n2,n3,n4*ncld), &
                          nprecpp(n1,n2,n3,nprc),    &
                          mprecpp(n1,n2,n3,n4*nprc), &
                          nicep(n1,n2,n3,nice),      &
                          micep(n1,n2,n3,n4*nice),   &
                          nsnowp(n1,n2,n3,nsnw),     &
                          msnowp(n1,n2,n3,n4*nsnw)

      REAL, INTENT(inout) :: naerot(n1,n2,n3,nbins),    &
                             maerot(n1,n2,n3,n4*nbins), &
                             ncloudt(n1,n2,n3,ncld),    &
                             mcloudt(n1,n2,n3,n4*ncld), &
                             nprecpt(n1,n2,n3,nprc),    &
                             mprecpt(n1,n2,n3,n4*nprc), &
                             nicet(n1,n2,n3,nice),      &
                             micet(n1,n2,n3,n4*nice),   &
                             nsnowt(n1,n2,n3,nsnw),     &
                             msnowt(n1,n2,n3,n4*nsnw),  &
                             tlt(n1,n2,n3),             & ! Liquid water pot temp tendency
                             rrate(n1,n2,n3),           &
                             srate(n1,n2,n3)              ! snowing rate

      INTEGER :: ss
      INTEGER :: i,j,k,nc,istr,iend
      REAL :: xnpts
      REAL :: v1(n1)

      REAL :: prnt(n1,n2,n3,nprc), prvt(n1,n2,n3,n4*nprc)     ! Number and mass tendencies due to fallout

      REAL :: srnt(n1,n2,n3,nsnw), srvt(n1,n2,n3,n4*nsnw)     ! Number and mass tendencies due to fallout

      ! PArticle removal arrays, given in kg/(m2 s)
      REAL :: remaer(n2,n3,n4*nbins),   &
              remcld(n2,n3,n4*ncld),    &
              remprc(n2,n3,n4*nprc),    &
              remice(n2,n3,n4*nice),    &
              remsnw(n2,n3,n4*nsnw)

      REAL :: mctmp(n2,n3) ! Helper for mass conservation calculations

      ! Divergence fields
      REAL :: amdiv(n1,n2,n3,n4*nbins),    &
              cmdiv(n1,n2,n3,n4*ncld),     &
              imdiv(n1,n2,n3,n4*nice)
      REAL :: andiv(n1,n2,n3,nbins),       &
              cndiv(n1,n2,n3,ncld),        &
              indiv(n1,n2,n3,nice)
    
      ! Deposition fields (given as particle divergence)
      REAL :: amdep(n2,n3,n4*nbins),  &
              cmdep(n2,n3,n4*ncld),   &
              imdep(n2,n3,n4*nice)
      REAL :: andep(n2,n3,nbins),     &
              cndep(n2,n3,ncld),      &
              indep(n2,n3,nice)


      remaer = 0.; remcld = 0.; remprc = 0.; remice = 0.; remsnw = 0.
      amdiv = 0.; amdep = 0.; cmdiv = 0.; cmdep = 0.; imdiv = 0.; imdep = 0.
      andiv = 0.; andep = 0.; cndiv = 0.; cndep = 0.; indiv = 0.; indep = 0.
      prnt = 0.; prvt = 0.; srnt = 0.; srvt = 0.

      IF(sflg) THEN
         xnpts = 1./((n3-4)*(n2-4))
         DO k = 1, n1
            v1(k) = 0.
         END DO
      END IF

      ! Sedimentation for slow (non-precipitating) particles
      !-------------------------------------------------------
      IF (sed_aero) THEN

         CALL NumMassDivergence(n1,n2,n3,n4,nbins,tk,a_dn,1500.,naerop,maerop,dzt,nlim,andiv,amdiv)
       
         CALL NumMassDepositionSlow(n1,n2,n3,n4,nbins,a_dn,1500.,tk,ustar,naerop,maerop,nlim,andep,amdep)
       
         naerot = naerot - andiv
         naerot(2,:,:,:) = naerot(2,:,:,:) - andep
         maerot = maerot - amdiv
         maerot(2,:,:,:) = maerot(2,:,:,:) - amdep

         remaer(:,:,:) = amdep(:,:,:)

         ! Account for changes in in liquid water pot temperature
         nc = GetIndex(prtcl,'H2O')
         istr = (nc-1)*nbins+1
         iend = nc*nbins
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1
                  tlt(k,i,j) = tlt(k,i,j) + SUM(amdiv(k,i,j,istr:iend))*(alvl/cp)*th(k,i,j)/tk(k,i,j)
               END DO
            END DO
         END DO

      END IF ! sed_aero

      IF (sed_cloud) THEN

         CALL NumMassDivergence(n1,n2,n3,n4,ncld,tk,a_dn,rhowa,ncloudp,mcloudp,dzt,nlim,cndiv,cmdiv)
    
         CALL NumMassDepositionSlow(n1,n2,n3,n4,ncld,a_dn,rhowa,tk,ustar,ncloudp,mcloudp,nlim,cndep,cmdep)

         ncloudt = ncloudt - cndiv
         ncloudt(2,:,:,:) = ncloudt(2,:,:,:) - cndep
         mcloudt = mcloudt - cmdiv
         mcloudt(2,:,:,:) = mcloudt(2,:,:,:) - cmdep

         remcld(:,:,:) = cmdep(:,:,:)

         ! Account for changes in in liquid water pot temperature
         nc = GetIndex(prtcl,'H2O')
         istr = (nc-1)*ncld+1
         iend = nc*ncld
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1
                  tlt(k,i,j) = tlt(k,i,j) + SUM(cmdiv(k,i,j,istr:iend))*(alvl/cp)*th(k,i,j)/tk(k,i,j)
               END DO
            END DO
         END DO

      END IF ! sed_cloud

      IF (sed_ice) THEN

         CALL NumMassDivergence(n1,n2,n3,n4,nice,tk,a_dn,rhoic,nicep,micep,dzt,prlim,indiv,imdiv)

         CALL NumMassDepositionSlow(n1,n2,n3,n4,nice,a_dn,rhoic,tk,ustar,nicep,micep,prlim,indep,imdep)

         nicet = nicet - indiv
         nicet(2,:,:,:) = nicet(2,:,:,:) - indep
         micet = micet - imdiv
         micet(2,:,:,:) = micet(2,:,:,:) - imdep

         remice(:,:,:) = imdep(:,:,:)

         ! Account for changes in in liquid water pot temperature
         nc = GetIndex(prtcl,'H2O')
         istr = (nc-1)*nice+1
         iend = nc*nice
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1
                  tlt(k,i,j) = tlt(k,i,j) + SUM(imdiv(k,i,j,istr:iend))*(alvi/cp)*th(k,i,j)/tk(k,i,j)
               END DO
            END DO
         END DO

      END IF ! sed_ice

      ! ---------------------------------------------------------
      ! SEDIMENTATION/DEPOSITION OF FAST PRECIPITATING PARTICLES
      IF (sed_precp) THEN
         CALL DepositionFast(n1,n2,n3,n4,nprc,tk,a_dn,rowt,nprecpp,mprecpp,tstep,dzt,prnt,prvt,remprc,prlim,rrate)
       
         nprecpt(:,:,:,:) = nprecpt(:,:,:,:) + prnt(:,:,:,:)/tstep
         mprecpt(:,:,:,:) = mprecpt(:,:,:,:) + prvt(:,:,:,:)/tstep

         ! Convert mass flux to heat flux (W/m^2)
         rrate(:,:,:) = rrate(:,:,:)*alvl

         IF (sflg) v1(:) = v1(:) + SUM(SUM(rrate(:,:,:),DIM=3),DIM=2)*xnpts

         nc = GetIndex(prtcl,'H2O')
         istr = (nc-1)*nprc + 1
         iend = nc*nprc

         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1-1

                  tlt(k,i,j) = tlt(k,i,j) - SUM(prvt(k,i,j,istr:iend))/tstep*(alvl/cp)*th(k,i,j)/tk(k,i,j)

               END DO
            END DO
         END DO
      END IF
    
      IF (sed_snow) THEN
         CALL DepositionFast(n1,n2,n3,n4,nsnw,tk,a_dn,rowt,nsnowp,msnowp,tstep,dzt,srnt,srvt,remsnw,prlim,srate)
           
         nsnowt(:,:,:,:) = nsnowt(:,:,:,:) + srnt(:,:,:,:)/tstep
         msnowt(:,:,:,:) = msnowt(:,:,:,:) + srvt(:,:,:,:)/tstep

         ! Convert mass flux to heat flux (W/m^2)
         rrate(:,:,:) = rrate(:,:,:)*alvi

         IF (sflg) v1(:) = v1(:) + SUM(SUM(srate(:,:,:),DIM=3),DIM=2)*xnpts

         nc = GetIndex(prtcl,'H2O')
         istr = (nc-1)*nsnw + 1
         iend = nc*nsnw
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1-1

                  tlt(k,i,j) = tlt(k,i,j) - SUM(srvt(k,i,j,istr:iend))/tstep*(alvi/cp)*th(k,i,j)/tk(k,i,j)

               END DO
            END DO
         END DO
      END IF

      IF (mcflg) THEN
         ! For mass conservation statistics
         mctmp(:,:) = 0.
         ss = getIndex(prtcl,'H2O')
         istr = (ss-1)*nbins; iend = ss*nbins
         mctmp(:,:) = mctmp(:,:) + SUM(remaer(:,:,istr:iend),dim=3)
         istr = (ss-1)*ncld; iend = ss*ncld
         mctmp(:,:) = mctmp(:,:) + SUM(remcld(:,:,istr:iend),dim=3)
         istr = (ss-1)*nprc; iend = ss*nprc
         mctmp(:,:) = mctmp(:,:) + SUM(remprc(:,:,istr:iend),dim=3)
         istr = (ss-1)*nice; iend = ss*nice
         mctmp(:,:) = mctmp(:,:) + SUM(remice(:,:,istr:iend),dim=3)
         istr = (ss-1)*nsnw; iend = ss*nsnw
         mctmp(:,:) = mctmp(:,:) + SUM(remsnw(:,:,istr:iend),dim=3)
         CALL acc_massbudged(n1,n2,n3,3,tstep,dzt,a_dn,rdep=mctmp,ApVdom=mc_ApVdom)
      END IF !mcflg
      IF (sflg) CALL updtst(n1,'prc',1,v1,1)
      ! Aerosol removal statistics
      IF (sflg) CALL acc_removal(n2,n3,n4,remaer,remcld,remprc,remice,remsnw)
      IF (sflg) CALL cs_rem_set(n2,n3,n4,remaer,remcld,remprc,remice,remsnw)

   END SUBROUTINE !sedim_SALSA

   ! -----------------------------------------------------------------


   SUBROUTINE NumMassDivergence(n1,n2,n3,n4,nn,tk,adn,pdn,numc,mass,dzt,clim,flxdivn,flxdivm)
  
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3,n4       ! Grid numbers, number of chemical species
      INTEGER, INTENT(in) :: nn                ! Number of bins
      REAL, INTENT(in) :: tk(n1,n2,n3)         ! Absolute temprature
      REAL, INTENT(in) :: adn(n1,n2,n3)        ! Air density
      REAL, INTENT(in) :: pdn                  ! Particle density
      REAL, INTENT(in) :: numc(n1,n2,n3,nn)    ! Particle number concentration
      REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4) ! Particle mass mixing ratio
      REAL, INTENT(in) :: dzt(n1)              ! Inverse of grid level thickness
      REAL, INTENT(IN) :: clim                ! Concentration limit
      REAL, INTENT(OUT) :: flxdivm(n1,n2,n3,nn*n4)
      REAL, INTENT(OUT) :: flxdivn(n1,n2,n3,nn)

      INTEGER :: i,j,k,kp1
      INTEGER :: bin,ss,bs

      REAL, PARAMETER :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                        ! modelling [kg molec-1]
      REAL, PARAMETER :: A = 1.249      ! fundamentals of atm. modelling pg509
      REAL, PARAMETER :: B = 0.42
      REAL, PARAMETER :: C = 0.87

      REAL :: avis, kvis           ! Viscosity of air, kinematic viscosity
      REAL :: lambda              ! Mean free path
      REAL :: va                    ! Thermal speed of air molecule
      REAL :: Kn, GG               ! Knudsen number, slip correction
      REAL :: vc                    ! Critical fall speed
      REAL :: rflm(n1,nn,n4), rfln(n1,nn), prvolc(n4), rwet
      flxdivm = 0.
      flxdivn = 0.

      DO j = 3, n3-2

         DO i = 3, n2-2

            rflm = 0.
            rfln = 0.
          
            DO k = n1-1, 2, -1
               kp1 = k+1

               ! atm modelling Eq.4.54
               avis = 1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
               kvis = avis/adn(k,i,j) !actual density ???
               va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
               lambda = 2*avis/(adn(k,i,j)*va) !mean free path

               ! Fluxes
               !------------------
               ! -- Calculate the *corrections* for small particles
               DO bin = 1, nn
                  IF (numc(k,i,j,bin) < clim) CYCLE

                  ! Calculate wet size
                  !   n4 = number of active species
                  !   bin = size bin
                  prvolc(:) = mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                  rwet = calc_wet_radius(n4,numc(k,i,j,bin),prvolc)

                  ! Terminal velocity
                  Kn = lambda/rwet
                  GG = 1.+ Kn*(A+B*exp(-C/Kn))
                  vc = terminal_vel(rwet,pdn,adn(k,i,j),avis,GG)

                  ! Flux for the particle mass
                  IF ( k > 2 ) THEN
                     DO ss = 1, n4
                        bs = (ss-1)*nn + bin
                        rflm(k,bin,ss) = -mass(k,i,j,bs)*vc
                     END DO
                  END IF
                  ! Flux for the particle number
                  IF ( k > 2 ) rfln(k,bin) = -numc(k,i,j,bin)*vc
               END DO

               DO bin = 1, nn
                  DO ss = 1, n4
                     bs = (ss-1)*nn + bin
                     flxdivm(k,i,j,bs) = (rflm(kp1,bin,ss)-rflm(k,bin,ss))*dzt(k)
                  END DO
               END DO

               flxdivn(k,i,j,:) = (rfln(kp1,:)-rfln(k,:))*dzt(k)

            END DO ! k
         END DO ! i
      END DO ! j

   END SUBROUTINE NumMassDivergence


   SUBROUTINE NumMassDepositionSlow(n1,n2,n3,n4,nn,adn,pdn,tk,ustar,numc,mass,clim,depflxn,depflxm)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3,n4,nn
      REAL, INTENT(in) :: adn(n1,n2,n3)
      REAL, INTENT(in) :: pdn
      REAL, INTENT(in) :: tk(n1,n2,n3)
      REAL, INTENT(in) :: ustar(n2,n3)
      REAL, INTENT(in) :: numc(n1,n2,n3,nn)
      REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4)
      REAL, INTENT(IN) :: clim                ! Concentration limit
      REAL, INTENT(OUT) :: depflxn(n2,n3,nn), depflxm(n2,n3,nn*n4)

      INTEGER :: i,j,k,bin
      INTEGER :: ss,bs

      REAL, PARAMETER :: A = 1.249 ! fundamentals of atm. modelling pg509
      REAL, PARAMETER :: B = 0.42
      REAL, PARAMETER :: C = 0.87
      REAL, PARAMETER :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                        ! modelling [kg molec-1]

      REAL :: GG
      REAL :: St, Sc, Kn
      REAL :: lambda      ! Mean free path
      REAL :: mdiff       ! Particle diffusivity
      REAL :: avis,kvis   ! Air viscosity, kinematic viscosity
      REAL :: va          ! Thermal speed of air molecule
      REAL :: vc,vd       ! Particle fall speed, deposition velocity
      REAL :: rt

      REAL :: rflm(nn,n4), rfln(nn), prvolc(n4), rwet

      depflxm = 0.
      depflxn = 0.

      ! Fix level to lowest above ground level
      k = 2

      DO j = 3, n3-2
         DO i = 3, n2-2

            ! atm modelling Eq.4.54
            avis = 1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
            kvis = avis/adn(k,i,j) !actual density ???
            va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
            lambda = 2*avis/(adn(k,i,j)*va) !mean free path

            rflm = 0.
            rfln = 0.

            DO bin = 1, nn
               IF (numc(k,i,j,bin) < clim) CYCLE

               ! Calculate wet size
               !   n4 = number of active species
               !   bin = size bin
               prvolc(:) = mass(k,i,j,bin:(n4-1)*nn+bin:nn)
               rwet = calc_wet_radius(n4,numc(k,i,j,bin),prvolc)

               Kn = lambda/rwet
               GG = 1.+ Kn*(A+B*exp(-C/Kn))

               ! Terminal velocity
               vc = terminal_vel(rwet,pdn,adn(k,i,j),avis,GG)

               ! Particle diffusitivity  (15.29) in jacobson book
               mdiff = (kb*tk(k,i,j)*GG)/(6.0*pi*avis*rwet)
             
               Sc = kvis/mdiff
               St = vc*ustar(i,j)**2/g*kvis !changed exponent
               IF (St < 0.01) St = 0.01
               rt = 1.0/MAX(epsilon(1.0),(ustar(i,j)*(Sc**(-2.0/3.0)+10**(-3.0/St)))) ! atm chem&phy eq19.18

               vd = (1./rt) + vc
    
               DO ss = 1, n4
                  bs = (ss-1)*nn + bin
                  rflm(bin,ss) = -mass(k,i,j,bs)*vd
                  depflxm(i,j,bs) = -rflm(bin,ss)*dzt(k)
               END DO

               rfln(bin) = -numc(k,i,j,bin)*vd

            END DO ! bin

            depflxn(i,j,:) = -rfln(:)*dzt(k)

         END DO ! i
      END DO ! j

   END SUBROUTINE NumMassDepositionSlow
  
   !------------------------------------------------------------------

   SUBROUTINE DepositionFast(n1,n2,n3,n4,nn,tk,adn,pdn,numc,mass,tstep,dzt,prnt,prvt,remprc,clim,rate)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3,n4,nn
      REAL, INTENT(in)    :: tk(n1,n2,n3)
      REAL, INTENT(in)    :: adn(n1,n2,n3)
      REAL, INTENT(in)    :: pdn
      REAL, INTENT(in)    :: numc(n1,n2,n3,nn)
      REAL, INTENT(in)    :: mass(n1,n2,n3,nn*n4)
      REAL, INTENT(in)    :: tstep
      REAL, INTENT(in)    :: dzt(n1)
      REAL, INTENT(IN)    :: clim                ! Concentration limit

      REAL, INTENT(out) :: prnt(n1,n2,n3,nn), prvt(n1,n2,n3,nn*n4)     ! Number and mass tendencies due to fallout
      REAL, INTENT(out) :: remprc(n2,n3,nn*n4)
      REAL, INTENT(out) :: rate(n1,n2,n3) ! Rain rate (kg/s/m^2)

      INTEGER :: k,i,j,bin
      INTEGER :: istr,iend

      REAL, PARAMETER :: A = 1.249 ! fundamentals of atm. modelling pg509
      REAL, PARAMETER :: B = 0.42
      REAL, PARAMETER :: C = 0.87
      REAL, PARAMETER :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                        ! modelling [kg molec-1]

      REAL :: GG
      REAL :: Kn
      REAL :: lambda      ! Mean free path
      REAL :: avis,kvis   ! Air viscosity, kinematic viscosity
      REAL :: va          ! Thermal speed of air molecule
      REAL :: vc
      REAL :: rwet

      ! For precipitation:
      REAL :: fd,fdmax,fdos ! Fall distance for rain drops, max fall distance, overshoot from nearest grid level
      REAL :: prnchg(n1,nn), prvchg(n1,nn,n4) ! Instantaneous changes in precipitation number and mass (volume)
 
      REAL    :: prnumc, prvolc(n4)  ! Instantaneous source number and volume
      INTEGER :: kf, ni,fi
      LOGICAL :: prcdep  ! Deposition flag

      remprc(:,:,:) = 0.
      rate(:,:,:) = 0.
      prnt(:,:,:,:) = 0.
      prvt(:,:,:,:) = 0.

      DO j = 3, n3-2
       
         DO i = 3, n2-2

            prnchg = 0.
            prvchg = 0.
          
            DO k = n1-1, 2, -1
          
               ! atm modelling Eq.4.54
               avis = 1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
               kvis = avis/adn(k,i,j) !actual density ???
               va = sqrt(8.*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
               lambda = 2.*avis/(adn(k,i,j)*va) !mean free path

               ! Precipitation bin loop
               DO bin = 1, nn
                  IF (numc(k,i,j,bin) < clim) CYCLE
             
                  ! Calculate wet size
                  !   n4 = number of active species
                  !   bin = size bin
                  prvolc(:) = mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                  rwet = calc_wet_radius(n4,numc(k,i,j,bin),prvolc)

                  ! Terminal velocity
                  Kn = lambda/rwet
                  GG = 1.+ Kn*(A+B*exp(-C/Kn))
                  vc = terminal_vel(rwet,pdn,adn(k,i,j),avis,GG)

                  ! Rain rate statistics: removal of water from the current bin is accounted for
                  ! Water is the last (n4) species and rain rate is given here kg/s/m^2
                  rate(k,i,j) = rate(k,i,j)+mass(k,i,j,(n4-1)*nn+bin)*adn(k,i,j)*vc

                  ! Determine output flux for current level: Find the closest level to which the
                  ! current drop parcel can fall within 1 timestep. If the lowest atmospheric level
                  ! is reached, the drops are sedimented.
                
                  ! Maximum fall distance:
                  fdmax = tstep*vc
             
                  fd = 0.
                  fi = 0
                  prcdep = .FALSE. ! deposition flag
                  DO WHILE ( fd < fdmax )
                     fd = fd + ( 1./dzt(k-fi) )
                     fi = fi + 1
                     ! Check if sedimentation occurs for current parcel
                     IF (k-fi <= 1) THEN
                        prcdep = .TRUE.
                        EXIT
                     END IF
                  END DO

                  ! Back up one level
                  fi = fi-1
                  kf = k - fi
                  fd = fd - ( 1./dzt(kf) )
             
                  ! How much the actual fall distance overshoots below the layer kf
                  fdos = MIN(MAX(fdmax-fd,0.),1./dzt(kf))
             
                  ! Remove the drops from the original level
                  prnumc = numc(k,i,j,bin)
                  prnchg(k,bin) = prnchg(k,bin) - prnumc
                  DO ni = 1, n4
                     prvolc(ni) = mass(k,i,j,(ni-1)*nn+bin)
                     prvchg(k,bin,ni) = prvchg(k,bin,ni) - prvolc(ni)
                  END DO ! ni
             
                  ! Removal statistics
                  IF (prcdep) THEN
                     DO ni = 1, n4
                        remprc(i,j,(ni-1)*nn+bin) = remprc(i,j,(ni-1)*nn+bin) +    &
                                                    prvolc(ni)*adn(k,i,j)*vc
                     END DO
                  END IF ! prcdep

                  ! Put the drops to new positions (may be partially the original grid cell as well!)
                  IF (fdos*dzt(kf) > 0.5) THEN  ! Reduce numerical diffusion
                     prnchg(kf-1,bin) = prnchg(kf-1,bin) + prnumc
                     DO ni = 1, n4
                        prvchg(kf-1,bin,ni) = prvchg(kf-1,bin,ni) + prvolc(ni)
                     END DO
                  ELSE
                     prnchg(kf-1,bin) = prnchg(kf-1,bin) + ( fdos*dzt(kf) )*prnumc
                     prnchg(kf,bin) = prnchg(kf,bin) + ( 1. - fdos*dzt(kf) )*prnumc
                     DO ni = 1, n4
                        prvchg(kf-1,bin,ni) = prvchg(kf-1,bin,ni) + ( fdos*dzt(kf) )*prvolc(ni)
                        prvchg(kf,bin,ni) = prvchg(kf,bin,ni) + ( 1. - fdos*dzt(kf) )*prvolc(ni)
                     END DO
                  END IF ! diffusion
             
               END DO !bin
          
            END DO ! k

            prnt(:,i,j,:) = prnt(:,i,j,:) + prnchg(:,:)
            DO ni = 1, n4
               istr = (ni-1)*nn+1
               iend = ni*nn
               prvt(:,i,j,istr:iend) = prvt(:,i,j,istr:iend) + prvchg(:,:,ni)
            END DO !ni

         END DO ! i

      END DO ! j

   END SUBROUTINE DepositionFast

   !********************************************************************
   ! Function for calculating terminal velocities for different particles size ranges.
   !     Tomi Raatikainen (2.5.2017)
   REAL FUNCTION terminal_vel(radius,rhop,rhoa,visc,beta)

      IMPLICIT NONE

      REAL, INTENT(in) :: radius, rhop ! Particle radius and density
      REAL, INTENT(in) :: rhoa, visc, beta ! Air density, viscocity and Cunningham correction factor

      ! Constants
      REAL, PARAMETER :: rhoa_ref = 1.225 ! reference air density (kg/m^3)

      IF (radius < 40.0e-6) THEN

         ! Stokes law with Cunningham slip correction factor
         terminal_vel = (4.*radius**2)*(rhop-rhoa)*g*beta/(18.*visc) ![m s-1]

      ELSE IF (radius < 0.6e-3) THEN

         ! Droplets from 40 um to 0.6 mm: linear dependence on particle radius and a correction for reduced pressure
         !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
         terminal_vel = 8.e3*radius*sqrt(rhoa_ref/rhoa)

      ELSE

         ! Droplets larger than 0.6 mm: square root dependence on particle radius and a correction for reduced pressure
         !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
         ! Note: this is valid up to 2 mm or 9 m/s (at 1000 mbar), where droplets start to break
         terminal_vel = 2.01e2*sqrt( min(radius,2.0e-3)*rhoa_ref/rhoa)

      END IF

   END FUNCTION terminal_vel
   !********************************************************************

END MODULE mcrp
