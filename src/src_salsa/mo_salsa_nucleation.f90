!****************************************************************
!*                                                              *
!*   MODULE MO_SALSA_NUCLEATION                             *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to new particle formation by nucleation                    *
!*                                                              *
!****************************************************************

MODULE mo_salsa_nucleation
    IMPLICIT NONE

CONTAINS

  ! fxm: currently only sulphuric acid grows particles from 1 to 3 nm
  !  (if asked from Markku, this is terribly wrong!!!
  !   nonvolatile OC should be added???)
  !********************************************************************
  !
  ! Subroutine NUCLEATION(kbdim,klev, &
  !       )
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates the particle number and volume increase,
  !  and gas-phase concentration decrease due to nucleation
  !  subsequent growth to detectable size of 3 nm
  !
  !
  ! Method:
  ! -------
  ! When the formed clusters grow by condensation (possibly also by
  !  self-coagulation), their number is reduced due to
  !  scavenging to pre-existing particles. Thus, the apparent
  !  nucleation rate at 3 nm is significantly lower than
  !  the real nucleation rate (at ~1 nm).
  !
  ! The formation rate of detectable particles at 3 nm is
  !  calculated with a parameterization presented in:
  !
  !  Kerminen, V.-M. and Kulmala, M. (2002) Analytical formulae
  !  connecting the 'real' and the 'apparent' nucleation rate
  !  and the nuclei number concentration for atmospheric
  !  nucleation events, J. Aerosol Sci., 33, 609-622.
  !
  !
  ! Interface:
  ! ----------
  ! Called from subroutine condensation
  !
  !
  ! Externals:
  ! ----------
  ! calls one of the following subroutines:
  !  - binnucl
  !  - ternucl
  !  - kinnucl
  !  - actnucl
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005
  ! Harri Kokkola (FMI) 2006
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  !
  !---------------------------------------------------------------------

  SUBROUTINE nucleation(kbdim,  klev,      &
                        paero,  ptemp,  prh,    ppres,  &
                        pcsa,   pcocnv, ptstep, pj3n3,  &
                        pxsa,   pxocnv, ppbl            )

    !USE mo_aero_mem_salsa, ONLY: d_jnuc,d_j3 !ALN002

    USE mo_submctl,   ONLY:  &
         t_section,              &
         act_coeff,              &
         nj3,                    &
         nsnucl,                 &
         fn2b,                   &
         pstand,                 &
         d_sa,                   &
         massacc,                &
         n3,                     &
         msu,                    &
         rhosu,                  &
         boltz,                  &
         avog,                   &
         rg,                     &
         pi,                     &
         reglim

    IMPLICIT NONE

    !-- Input and output variables -------------
    INTEGER, INTENT(IN) ::        &
         kbdim,                   & ! dimension for arrays
         klev                       ! number of vertical klev


    REAL, INTENT(IN) ::       &
         ptemp(kbdim,klev),       & ! ambient temperature [K]
         prh(kbdim, klev),        & ! relative humidity
         ppres(kbdim,klev),       & ! ambient pressure [Pa]
         ptstep,                  & ! timestep [s]
         pcsa(kbdim,klev),        & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev)         ! concentration of organic matter [#/m3]

    INTEGER, INTENT(IN) :: ppbl(kbdim) ! Planetary boundary layer top level

    REAL, INTENT(INOUT) ::    &
         pj3n3(kbdim,klev,2)        ! change in concentration [#/m3]

    TYPE(t_section), INTENT(inout) :: &
         paero(kbdim,klev,fn2b)


    REAL, INTENT(OUT) ::      &
         pxsa(kbdim,klev),        & ! ratio of sulphuric acid and organic vapor in 3nm particles [ ]
         pxocnv(kbdim,klev)

    !-- Local variables ------------------------
    INTEGER :: ii, jj, iteration
    INTEGER :: ppbl_bin(kbdim)      ! when nsnucl is chose to be binary nucleation, it is calculated at all vertical klev.
                                    ! with other nucleation types binary nucleation is calculated in planetary boundary layer and above
                                    ! according to the chosen nucleation type.

    REAL ::                   &
         zjnuc(kbdim,klev),       & ! nucleation rate at ~ 1 nm [#/m3s]
         zmixnh3(kbdim,klev),     & ! ammonia mixing ratio [ppt]
         zc_h2so4(kbdim,klev),    & ! sulphuric acid concentration (note units!) [#/cm3]
         zcsa_local(kbdim,klev),  & ! sulphuric acid concentration (note units!) [#/m3]
         zc_org(kbdim,klev),      & ! organic vapour concentration [#/cm3]
         zcocnv_local(kbdim,klev),& ! organic vapour concentration [#/m3]
         znsa(kbdim,klev),        & ! number of H2SO4 molecules in critical cluster [1]
         znoc(kbdim,klev),        & ! number of ORGANIC molecules in critical cluster [1]
         zksa(kbdim,klev),        & ! Lever: if k_sa = 1, h2so4 is involved in nucleation. [1]
         zkocnv(kbdim,klev) ,     & ! Lever: if k_ocnv = 1, organic compounds are involved in nucleation. [1]
         zdcrit(kbdim,klev),      & ! diameter of critical cluster [m]
         zknud(fn2b),             & ! particle Knudsen number [1]
         zbeta(fn2b),             & ! transitional correction factor [1]
         zGRclust,                & ! growth rate of formed clusters [nm/h]
         zdfvap,                  & ! air diffusion coefficient [m2/s]
         zmfp,                    & ! mean free path of condensing vapour [m]
         zcsink,                  & ! condensational sink [#/m2]
         zdmean,                  & ! mean diameter of existing particles [nm]
         zgamma,                  & ! proportionality factor [(nm2 m2)/h]
         zeta,                    & ! constant; proportional to ratio of CS/GR [m]
         zj3,                     & ! number conc. of formed 3 nm particles [#/m3]
         zdelta_vap,              & ! change in H2SO4 and Organic vapor concentration [#/m3]
         zNnuc,                   & ! Number of clusters/particles at the size range [d1,dx] (in /m3)
         zlambda,                 & ! A parameter for "adjusting" the growth rate due to self-coagulation
         zKeff,                   & ! "effective" coagulation coefficient between fresly-nucleated particles
         zGRtot,                  & ! Total growth rate
         zCoagStot,               & ! Total losses due to coagulation, includes condensation and self-coagulation

        !variables determined for the m-parameter
         zRc2(fn2b), &
         zRx2(fn2b), &
         zm_c, &
         zm_x, &
         zm_2(fn2b), &
         zcv_c, &       ! zcv = gas-phase velocity of clusters
         zcv_x, &
         zcv_2(fn2b), &
         zcv_c2(fn2b), &
         zcv_x2(fn2b), &
         zknud_c,&
         zknud_x,&
         zknud_2(fn2b),&
         zCc_c,&
         zCc_x,&
         zCc_2(fn2b),&
         zmyy, &
         zDc_c(fn2b), &
         zDc_x(fn2b), &
         zDc_2(fn2b), &
         zDc_c2(fn2b), &
         zDc_x2(fn2b), &
         zgammaF_c(fn2b), &
         zgammaF_x(fn2b), &
         zgammaF_2(fn2b), &
         zomega_c(fn2b), &
         zomega_x(fn2b), &
         zomega_2c(fn2b), &
         zomega_2x(fn2b), &
         zsigma_c2(fn2b), &
         zsigma_x2(fn2b), &
         zK_c2(fn2b), &
         zK_x2(fn2b), &
         zCoagS_c, &
         zCoagS_x, &
         zm_para
    !--------------------------------------------------------------------------------------

    !-- 1) Nucleation rate and diameter of critical cluster -------------------
    zjnuc = 0
    SELECT CASE (nsnucl)

    CASE(1) ! binary H2SO4-H2O nucleation

       zc_h2so4 = pcsa*1.e-6           ! [#/cm3]

       ppbl_bin = klev

       CALL binnucl(kbdim, klev, &
                    zc_h2so4, ptemp,  prh,     &
                    zjnuc,    znsa, znoc, zdcrit, &
                    ppbl_bin, zksa, zkocnv)

    CASE(2) ! activation type nucleation

       zc_h2so4 = pcsa*1.e-6           ! [#/cm3]

       CALL binnucl(kbdim, klev,     &
                    zc_h2so4, ptemp,  prh,         &
                    zjnuc,    znsa,   znoc, zdcrit, &
                    ppbl, zksa, zkocnv)

       CALL actnucl(kbdim, klev, &
                    pcsa,    zjnuc,  zdcrit, ppbl ,&
                    znsa,    znoc,   zksa,  zkocnv,act_coeff)

    CASE(3) ! kinetically limited nucleation of (NH4)HSO4 clusters

       zc_h2so4 = pcsa*1.e-6           ! [#/cm3]

       CALL binnucl( kbdim, klev,      &
                    zc_h2so4, ptemp,  prh,          &
                    zjnuc,    znsa,  znoc, zdcrit,   &
                    ppbl, zksa, zkocnv)

       CALL kinnucl(kbdim, klev, &
                    zc_h2so4,ptemp,          &
                    zjnuc,   zdcrit, ppbl,   &
                    znsa,    znoc,   zksa, zkocnv)

    CASE(4) ! ternary H2SO4-H2O-NH3 nucleation

       STOP 'mo_salsa_nucleation: Specify NH3 mixing ratio'

       zmixnh3 = 0.

       zc_h2so4 = pcsa*1.e-6           ! [#/cm3]

       CALL ternucl(kbdim,  klev,         &
                    zc_h2so4, zmixnh3, ptemp,   prh,     &
                    zjnuc,    znsa,    znoc,    zdcrit,  &
                    zksa,    zkocnv)

    CASE(5) ! organic nucleation, J~[ORG] or J~[ORG]**2

       zc_org = pcocnv*1.e-6           ! [#/cm3]
       zc_h2so4 = pcsa*1.e-6           ! [#/cm3]

       CALL binnucl(kbdim, klev, &
                    zc_h2so4, ptemp,  prh,     &
                    zjnuc,    znsa,   znoc, zdcrit, &
                    ppbl, zksa, zkocnv)

       CALL orgnucl(kbdim, klev, &
                    pcocnv,  zjnuc,  zdcrit, ppbl, &
                    znsa, znoc, zksa, zkocnv)

    CASE(6) ! sum of h2so4 and organic activation type nucleation
            ! J~[H2SO34]+[ORG]

       zc_h2so4 = pcsa*1.e-6

       CALL binnucl(kbdim, klev, &
                    zc_h2so4, ptemp,  prh,     &
                    zjnuc,    znsa,   znoc, zdcrit, &
                    ppbl, zksa, zkocnv)

       CALL sumnucl(kbdim, klev,  &
                    pcsa, pcocnv, zjnuc,&
                    zdcrit, ppbl, &
                    znsa, znoc, zksa, zkocnv)

    CASE(7) ! heteromolecular nucleation (J~[H2SO4]*[ORG])

       zc_h2so4 = pcsa*1.e-6
       zc_org = pcocnv*1.e-6

       CALL binnucl(kbdim, klev, &
                    zc_h2so4, ptemp,  prh,     &
                    zjnuc,    znsa,   znoc, zdcrit, &
                    ppbl, zksa, zkocnv)

       CALL hetnucl(kbdim, klev,  &
                    zc_h2so4, zc_org, zjnuc,   &
                    zdcrit, ppbl, &
                    znsa, znoc, zksa, zkocnv)

    CASE(8) ! homomolecular nucleation of H2SO4 and heteromolecular nucleation of
            ! H2SO4 and organic vapour (J~[H2SO4]**2 + [H2SO4]*[ORG]) !EUCAARI project

       zc_h2so4 = pcsa*1.e-6
       zc_org = pcocnv*1.e-6

       CALL binnucl(kbdim, klev, &
                    zc_h2so4, ptemp,  prh,     &
                    zjnuc,    znsa,   znoc, zdcrit, &
                    ppbl, zksa, zkocnv)

       CALL SAnucl(kbdim, klev,   &
                   zc_h2so4, zc_org, zjnuc, &
                   zdcrit, ppbl, &
                   znsa, znoc, zksa, zkocnv)

    CASE(9) ! homomolecular nucleation of H2SO4 and organic vapour and
            ! heteromolecular nucleation of H2SO4 and organic vapour
            ! (J~[H2SO4]**2 + [H2SO4]*[ORG]+[ORG]**2) ! EUCAARI project

       zc_h2so4 = pcsa*1.e-6
       zc_org = pcocnv*1.e-6

       CALL binnucl(kbdim, klev, &
                    zc_h2so4, ptemp,  prh,     &
                    zjnuc,    znsa,   znoc, zdcrit, &
                    ppbl, zksa, zkocnv)

       CALL SAORGnucl(kbdim, klev, &
                      zc_h2so4, zc_org, zjnuc,  zdcrit, ppbl, &
                      znsa, znoc, zksa, zkocnv)

    END SELECT

    zcsa_local = pcsa
    zcocnv_local = pcocnv

    !-- 2) Change of particle and gas concentrations ------------------------------------
    ! loops over
    DO jj = 1, klev !  vertical grid
       DO ii = 1, kbdim !  horizontal kbdim in the slab

!          !-- very small nucleation rates neglected
!          IF (zjnuc(ii,jj) < 1.e3) CYCLE
!
!          IF (pcsa(ii,jj)+pcocnv(ii,jj) < 1.) CYCLE


          !-- 2.1) Check that there is enough H2SO4 and Organic vapor to produce the nucleation -----

          IF (nsnucl <= 4) THEN
             !--------------------------------------------------------------------------------
             ! If the chosen nucleation scheme is 1-4, the nucleation occurs only due to H2SO4
             ! All of the total vapor concentration that is taking part to the nucleation is
             ! there for sulphuric acid (sa) and nonvolatile organic vapor is zero.

             pxsa(ii,jj) = 1.                ! ratio of sulphuric acid in 3nm particles
             pxocnv(ii,jj) = 0.              ! ratio of nonvolatile organic vapor in 3nm particles

          ELSE IF (nsnucl > 4) THEN
             !--------------------------------------------------------------------------------
             ! If the chosen nucleation scheme is 5-9, the nucleation occurs due to organic
             ! vapour or due to the combination of organic vapour and  H2SO4
             !
             ! The number of needed molecules depends on the chosen nucleation type and it has
             ! an effect also on the minimum ratio of the molecules present.

             IF (pcsa(ii,jj)*znsa(ii,jj)+pcocnv(ii,jj)*znoc(ii,jj) < 1.E-14) THEN
                pxsa(ii,jj) = 0.
                pxocnv(ii,jj) = 0.
             ELSE
                pxsa(ii,jj) =  pcsa(ii,jj)*znsa(ii,jj)/(pcsa(ii,jj)*znsa(ii,jj)+pcocnv(ii,jj)*znoc(ii,jj)) !pcsa(ii,jj)/(pcsa(ii,jj)+pcocnv(ii,jj))
                pxocnv(ii,jj) = pcocnv(ii,jj)*znoc(ii,jj)/(pcsa(ii,jj)*znsa(ii,jj)+pcocnv(ii,jj)*znoc(ii,jj)) ! 1.-pxsa
             END IF
          END IF

          !-- very small nucleation rates neglected
          IF (zjnuc(ii,jj) < 1.e3) CYCLE

          IF (pcsa(ii,jj)+pcocnv(ii,jj) < 1.) CYCLE


          ! The change in total vapor concentration is the sum of the concentrations of the
          ! vapors taking part to the nucleation (depends on the chosen nucleation scheme)

          zdelta_vap = MIN(zjnuc(ii,jj)*(znoc(ii,jj)+znsa(ii,jj)),(pcocnv(ii,jj)*zkocnv(ii,jj)+pcsa(ii,jj)*zksa(ii,jj))/ptstep)
          zjnuc(ii,jj) = zdelta_vap/(znoc(ii,jj)+znsa(ii,jj))
          zcsa_local(ii,jj) = MAX(1.,pcsa(ii,jj) - zdelta_vap*pxsa(ii,jj))  ! H2SO4 concentration after nucleation [#/m3]
          zcocnv_local(ii,jj) = MAX(1.,pcocnv(ii,jj) - zdelta_vap*pxocnv(ii,jj))  ! organic concentration after nucleation [#/m3]

          !-- 2.2) Formation rate of 3 nm particles (Kerminen & Kulmala, 2002) -----------------

          !--- 2.2.1) Growth rate of formed clusters by H2SO4
          !
          !  Kerminen & Kulmala (2002):
          !  GR = 3d-15/dens_clus*sum(molecspeed*molarmass*concent)
          !
          !   dens_clus = density of the clusters (here 1830 kg/m3)
          !   molarmass = molar mass of condensing species (here 98.08 g/mol)
          !   concent  = concentration of - " - [#/m3]
          !   molecspeed = molecular speed of - " - [m/s]
          !      = sqrt(8.*Rg*ptemp/(pi*molarmass))  (Seinfeld & Pandis, 1998)

          ! Growth rate by H2SO4 and Organic vapor
          zGRclust = 2.3623e-15*sqrt(ptemp(ii,jj))*(zcsa_local(ii,jj)+zcocnv_local(ii,jj)) ! [nm/h] (21) = GRcond

          !--- 2.2.2) Condensational sink of pre-existing particle population

          !--- mean free path of condensing vapour [m]
          !          for the formula, see subroutine CONDENSATION
          !
          zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj) ! diffusion coefficient [m2/s]
          zmfp   = 3.*zdfvap*sqrt(pi*msu/(8.*rg*ptemp(ii,jj)))       ! mean free path [m]

          !--- transitional regime correction factor
          zknud = 2.*zmfp/(paero(ii,jj,:)%dwet+d_sa)     ! Knudsen number
          ! correction factor according to
          ! Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
          ! Topics in current aerosol research, Pergamon.
          zbeta = (zknud + 1.)/ &
                  (0.377*zknud+1.+4./(3.*massacc)*(zknud+zknud**2))   !(4)
          !--- condensational sink [#/m2]
          zcsink = sum(paero(ii,jj,:)%dwet*zbeta*paero(ii,jj,:)%numc)                    !(3)

          !---------------------------------------------------------------------------------
          ! Parameterized formation rate of detectable 3 nm particles

          IF (nj3 == 1) THEN       ! Parametrization of J3 by Kerminen and Kulmala (2002)

             !--- 2.2.3) Parameterized formation rate of detectable 3 nm particles

             !--- Constants needed for the parameterization
             !
             ! The following values have been used to
             ! simplify the parameterization presented in
             ! Kerminen & Kulmala (2002):
             !
             ! dapp  = 3 nm
             ! dens_nuc = 1830 kg/m3
             !

             IF(zcsink < 1.e-30) THEN
                zeta = 0.
             ELSE
                zdmean = 1./sum(paero(ii,jj,:)%numc)* &       ! mean diameter of backgroud population [nm]
                  sum(paero(ii,jj,:)%numc*paero(ii,jj,:)%dwet)*1.e9

                zgamma = 0.23*(zdcrit(ii,jj)*1.e9)**0.2 &! fxm: can we use simple version of zgamma given in Kerminen et al.?
                         *(zdmean/150.)**0.048 &
                         *(ptemp(ii,jj)/293.)**(-0.75)*(rhosu/1000.)**(-0.33)
                ! [nm2*m2/h] (22)
                zeta = min(zgamma*zcsink/zGRclust, zdcrit(ii,jj)*1e11)! [nm] (11)
             END IF

             !-- Number conc. of clusters surviving to 3 nm in a time step

             zj3 = zjnuc(ii,jj)*exp(min(0.,zeta/3. - zeta/(zdcrit(ii,jj)*1.e9))) ! [#/m3] (14)

          ELSE IF (nj3 > 1) THEN

             ! ----------------------------------------------------------------------------------------
             ! Defining the value for zm_para
             ! The growth is investigated between [d1,reglim(1)] = [zdcrit,3nm]
             ! m = log(CoagS_dx/CoagX_zdcrit)/log(reglim/zdcrit) (Lehtinen et al. 2007)
             !
             ! The steps for the coagulation sink for reglim = 3nm and zdcrit ~= 1nm particles
             ! are explained in article of Kulmala et al. 2001.
             !
             ! The particles of diameter zdcrit ~1.14 nm  and reglim = 3nm are both in turn the
             ! "number 1" variables (according to the Kulmala et al. 2001)

             ! zRc = zdcrit(ii,jj)/2.       ! R1, [m]
             ! zRx = reglim(1)/2.          ! R1, [m]
             ! zR2 = pdwet(ii,jj,:)/2.      ! zR2  [m] a vector... (!)

             ! Sum of the radii, R12 = R1 + zR2         ! (m)
             zRc2 = zdcrit(ii,jj)/2. + paero(ii,jj,:)%dwet/2. ! [m]
             zRx2 = reglim(1)/2. + paero(ii,jj,:)%dwet/2.  ! [m]

             zm_c = 4./3.*pi*(zdcrit(ii,jj)/2.)**3*rhosu!(rhosu*pxsa+rhooc*pxocnv)    ! [kg] particle mass assumed to be only H2SO4
             zm_x = 4./3.*pi*(reglim(1)/2.)**3*rhosu!(rhosu*pxsa+rhooc*pxocnv)       ! [kg]
             zm_2 = 4./3.*pi*(paero(ii,jj,:)%dwet/2.)**3*rhosu!(rhosu*pxsa+rhooc*pxocnv)   ! [kg]

             zcv_c = SQRT(8.*boltz*ptemp(ii,jj)/(pi*zm_c))     ! [m/s]
             zcv_x = SQRT(8.*boltz*ptemp(ii,jj)/(pi*zm_x))     ! [m/s]
             zcv_2 = SQRT(8.*boltz*ptemp(ii,jj)/(pi*zm_2))     ! [m/s]

             zcv_c2 = SQRT(zcv_c**2+zcv_2**2)
             zcv_x2 = SQRT(zcv_x**2+zcv_2**2)

             zknud_c = 2.*zmfp/zdcrit(ii,jj) !Knudsen number
             zknud_x = 2.*zmfp/reglim(1)
             zknud_2 = 2.*zmfp/paero(ii,jj,:)%dwet

             zCc_c = 1.+zknud_c*(1.142+0.558*exp(-0.999/zknud_c)) ! Cunningham correction factor, [1]
             zCc_x = 1.+zknud_x*(1.142+0.558*exp(-0.999/zknud_x))
             zCc_2 = 1.+zknud_2*(1.142+0.558*exp(-0.999/zknud_2))

             zmyy = 1.81e-5*(ptemp(ii,jj)/293.)**(0.74) !
             ! gas dynamic viscosity, viscocity(air @20C) = 1.81e-5 N/m2 *s (Hinds, p. 25)

             zDc_c = boltz*ptemp(ii,jj)*zCc_c/(3.*pi*zmyy*zdcrit(ii,jj))  ! D1, Diffusion coefficient of zdcrit-particle, [m2/s]
             zDc_x = boltz*ptemp(ii,jj)*zCc_x/(3.*pi*zmyy*reglim(1))     ! D1, Diffusion coefficient of dx-particle
             zDc_2 = boltz*ptemp(ii,jj)*zCc_2/(3*pi*zmyy*paero(ii,jj,:)%dwet)
             zDc_c2 = zDc_c+zDc_2                                ! D12
             zDc_x2 = zDc_x+zDc_2                                ! D12

             ! zgammaF = 8*D/pi/zcv
             zgammaF_c = 8.*zDc_c/pi/zcv_c    ! [m]
             zgammaF_x = 8.*zDc_x/pi/zcv_x    ! [m]
             zgammaF_2 = 8.*zDc_2/pi/zcv_2    ! [m]

             zomega_c = ((zRc2+zgammaF_c)**3-(zRc2**2+zgammaF_c)**(3./2.))/&
                        (3.*zRc2*zgammaF_c)-zRc2 ! zomega1, [m]
             zomega_x = ((zRx2+zgammaF_x)**3-(zRx2**2+zgammaF_x)**(3./2.))/&
                        (3.*zRx2*zgammaF_x)-zRx2 ! zomega1

             zomega_2c = ((zRc2+zgammaF_2)**3-(zRc2**2+zgammaF_2)**(3./2.))/&
                         (3.*zRc2*zgammaF_2)-zRc2 ! zomega2
             zomega_2x = ((zRx2+zgammaF_2)**3-(zRx2**2+zgammaF_2)**(3./2.))/&
                         (3.*zRx2*zgammaF_2)-zRx2 ! zomega2

             zsigma_c2 = SQRT(zomega_c**2+zomega_2c**2) ! zigma12, [m]
             zsigma_x2 = SQRT(zomega_x**2+zomega_2x**2) ! zigma12

             zK_c2 = 4.*pi*zRc2*zDc_c2/(zRc2/(zRc2+zsigma_c2)+4.*zDc_c2/(zcv_c2*zRc2)) ![m*m2/s]
             zK_x2 = 4.*pi*zRx2*zDc_x2/(zRx2/(zRx2+zsigma_x2)+4.*zDc_x2/(zcv_x2*zRx2))

             zCoagS_c= MAX(1.e-20, sum(zK_c2*paero(ii,jj,:)%numc))  ![1/s]       Coagulation sink for critical cluster
             zCoagS_x= MAX(1.e-20, sum(zK_x2*paero(ii,jj,:)%numc))  !            Coagulation sink for 3 nm cluster

             zm_para = LOG(zCoagS_x/zCoagS_c)/LOG(reglim(1)/zdcrit(ii,jj))
             !-End of m-parater----------------------------------------------------------------------

             zgamma= (((reglim(1)/zdcrit(ii,jj))**(zm_para+1.))-1.)/(zm_para+1.) ! Anttila et al. 2010, eq. 5

             IF (nj3 == 2) THEN

                zj3 = zjnuc(ii,jj)*exp(min(0.,-zgamma*zdcrit(ii,jj)*zCoagS_c/(zGRclust*1.e-9/(60.**2)))) !J before iteration, in #/m3s

             ELSE IF (nj3 == 3) THEN

                ! If polluted air... then the self-coagulation becomes important
                !------------------------------------------------------------------------
                ! Self-coagulation of small particles < 3 nm
                !
                ! [Anttila, T., Kerminen, V-M, Lehtinen, K.E.J. Parametrizing the formation
                !     rate of new particles: the effect of nuclei self-coagulation, 2010]
                !
                !------------------------------------------------------------------------

                ! "effective" coagulation coefficient between freshly-nucleated particles:
                zKeff = 5.e-16              !cm3/s

                ! Use different values at sizes >10 nm:
                IF (reglim(1) >= 10.d-9) THEN
                   zlambda = 3.
                   zKeff = 5.e-17
                END IF

                ! zlambda parameter for "adjusting" the growth rate due to the self coagulation
                zlambda = 6.

                ! Initial values for coagulation sink and growth rate
                zCoagStot = zCoagS_c
                zGRtot = zGRclust*1.e-9/(60.**2)          ! zGRtot = GRcond = zGRclust (m/s)

                ! Number of clusters/particles at the size range [d1,dx] (in /m3):
                zNnuc = zjnuc(ii,jj)/zCoagStot                   ! initial guess [#/m3s*(1/s)] = [#/m3]

                ! Coagulation sink and growth rate due to self-coagulation are calculated as follows:
                ! CoagSscg = zKeff*zNnuc*1.e-6  ! [zKeff]=cm3/s, [zNnuc] = #/m3 = #/cm3*1e6, [1/s]
                ! GRscg = 1.5708e-6*zlambda*zdcrit(ii,jj)**3*(zNnuc*1.e-6)*zcv_c*avog*1.e-9/3600., [m/s]

                DO iteration = 1, 5

                   zCoagStot = zCoagS_c + zKeff*zNnuc*1.e-6  ! = CoagS + CoagSscg(d1) (1/s)

                   zGRtot = zGRclust*1.e-9/(3600.) + &
                            1.5708e-6*zlambda*zdcrit(ii,jj)**3*(zNnuc*1.e-6)*zcv_c*avog*1.e-9/3600.
                   ! = GRcond + GRscg, m/s

                   zeta = -zCoagStot/((zm_para+1.)*zGRtot*(zdcrit(ii,jj)**zm_para)) ! eq (7b)
                   zNnuc = zNnuc_tayl(zdcrit(ii,jj),reglim(1),zm_para,zjnuc(ii,jj),zeta,zGRtot)

                END DO

                ! calculate the final values with new zNnuc:
                zCoagStot = zCoagS_c + zKeff*zNnuc*1.e-6 ! = CoagS + CoagSscg(d1) (1/s)

                zGRtot = zGRclust*1.e-9/(3600) + &
                         1.5708e-6*zlambda*zdcrit(ii,jj)**3*(zNnuc*1.e-6)*zcv_c*avog*1.e-9/3600. ![m/s]

                zj3 = zjnuc(ii,jj)*exp(min(0.,-zgamma*zdcrit(ii,jj)*zCoagStot/zGRtot))  !(5a) [#/m3s]

             ELSE
                STOP 'mo_salsa_nucleation: Error in chosen nj3'
             END IF

          END IF

          !-- If J3 very small (< 1 #/cm3), neglect particle formation
          !  In real atmosphere this would mean that clusters form
          !  but coagulate to pre-existing particles who gain sulphate.
          !  Since CoagS ~ CS (4piD*CS'), we do *not* update H2SO4 concentration
          !  here but let condensation take care of it.

!          pcsa(ii,jj) = pcsa(ii,jj) + zdelta_vap*pxsa
!          pcocnv(ii,jj) = pcocnv(ii,jj) + zdelta_vap*pxocnv

          pj3n3(ii,jj,1) = zj3*n3*pxsa(ii,jj)
          pj3n3(ii,jj,2) = zj3*n3*pxocnv(ii,jj)

          !d_jnuc(ii,jj,krow)=d_jnuc(ii,jj,krow)+zjnuc(ii,jj) alaak
          !d_j3(ii,jj,krow)=d_j3(ii,jj,krow)+zj3 alaak

       END DO
    END DO

  END SUBROUTINE nucleation

  !********************************************************************
  !
  ! Subroutine BINNUCL
  ! Subroutine TERNUCL
  ! Subroutine KINNUCL
  ! Subroutine ACTNUCL
  !
  ! Subroutine ORGNUCL
  ! Subroutine SUMNUCL
  ! Subroutine HETNUCL
  ! Subroutine SANUCL
  ! Subroutine SAORGTNUCL
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculate the nucleation rate and the size of critical clusters
  !
  !
  ! Method:
  ! -------
  ! 1) Binary nucleation calculated according to
  !
  !  Vehkamaki et al. (2002): An improved parameterization for
  !  sulphuric acid/water nucleation rates for tropospheric
  !  and stratospheric conditions, JGR, 107, D22, 4622.
  !
  ! 2) Ternary nucleation calculated according to
  !
  !  Napari et al. (2002), An improved model for ternary nucleation of
  !  sulfuric acid - ammonia- water, J. Chem. Phys., 116, 4221-4227
  !
  !  Napari et al. (2002), Parameterization of ternary nucleation rates for
  !  H2SO4 - NH3 - H2O vapors, J. Geophys. Res., 107(D19), AAC 6-1
  !
  ! 3) Kinetic nucleation calculated assuming that each sulphuric acid
  !  molecule forms an (NH4)HSO4 molecule in the atmosphere and that
  !  two colliding (NH4)HSO4 molecules form a stable cluster
  !
  !
  ! Interface:
  ! ----------
  ! Called from subroutine nucleation
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI)      2005
  ! Hanna Vehkamaki (University of Helsinki) 2002
  ! Ismo Napari (University of Helsinki)  2002
  ! Sanna-Liisa Sihto (University of Helsinki) 2004
  !
  !---------------------------------------------------------------------

  SUBROUTINE binnucl(kbdim,      klev,                  &
                     pc_sa,     ptemp,      prh,                   &
                     pnuc_rate, pn_crit_sa, pn_crit_ocnv, pd_crit, &
                     ppbl,      pk_sa,      pk_ocnv)

    USE mo_submctl, ONLY : rg, pi, avog

    IMPLICIT NONE

    !-- Input variables -------------------
    INTEGER, INTENT(IN) :: &
         kbdim,  &      ! dimension for arrays
         klev           ! number of vertical klev

    REAL, INTENT(IN) ::   &
         pc_sa(kbdim,klev),   &  ! sulphuric acid concentration [#/cm3]
         ptemp(kbdim,klev),   &  ! ambient temperature [K]
         prh(kbdim,klev)         ! relative humidity [0-1]

    INTEGER :: ppbl(kbdim)            ! Boundary layer top level
    !-- Output variables ------------------
    REAL, INTENT(OUT) ::      &
         pnuc_rate(kbdim,klev),   &  ! nucleation rate [#/(m3 s)]
         pn_crit_sa(kbdim,klev),  &  ! number of sulphuric acid molecules in cluster [1]
         pn_crit_ocnv(kbdim,klev),&  ! number of organic molecules in cluster [1]
         pd_crit(kbdim,klev),     &  ! diameter of critical cluster [m]
         pk_sa(kbdim,klev),       &  ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
         pk_ocnv(kbdim,klev)         ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

    !-- Local variables -------------------
    INTEGER :: ii, jj   ! loop indices
    INTEGER :: zpbl(kbdim)

    REAL ::  &
         zx,     &   ! mole fraction of sulphate in critical cluster
         zntot,  &   ! number of molecules in critical cluster
         zt,     &   ! temperature
         zpcsa,  &   ! sulfuric acid concentration
         zrh,    &   ! relative humidity
         zma,zmw,zxmass, za, zb, zc, zroo, zm1, zm2, zv1, zv2, zcoll
    LOGICAL :: lnuctropo = .TRUE.
    !------------------------------------------------------------------------------------------------

    ! loops over

    pnuc_rate = 0.
    pd_crit   = 1.e-9

    zpbl(:) = ppbl(:)
    IF (lnuctropo) THEN

   DO ii = 1, kbdim !  horizontal kbdim in the slab
      DO jj = 1, zpbl(ii) !  vertical grid
             !-- 1) Checking that we are in the validity range of the parameterization -----------

             zt = max(ptemp(ii,jj), 190.15)
             zt = min(zt, 300.15)

             zpcsa = max(pc_sa(ii,jj), 1.e4)
             zpcsa = min(zpcsa, 1.e9) !TB maximum for atmospheric conditions

             zrh = max(prh(ii,jj), 0.0001)
             zrh = min(zrh, 1.)

             !-- 2) Mole fraction of sulphate in a critical cluster ------------------------------

             zx = 0.7409967177282139                            &
                  - 0.002663785665140117*zt                     &
                  + 0.002010478847383187*log(zrh)               &
                  - 0.0001832894131464668*zt*log(zrh)           &
                  + 0.001574072538464286*log(zrh)**2            &
                  - 0.00001790589121766952*zt*log(zrh)**2       &
                  + 0.0001844027436573778*log(zrh)**3           &
                  - 1.503452308794887e-6*zt*log(zrh)**3         &
                  - 0.003499978417957668*log(zpcsa)             &
                  + 0.0000504021689382576*zt*log(zpcsa)

             !-- 3) Nucleation rate -----------------------------------------------------------

             pnuc_rate(ii,jj) = 0.1430901615568665                      &
                  + 2.219563673425199*zt                                &
                  - 0.02739106114964264*zt**2                           &
                  + 0.00007228107239317088*zt**3                        &
                  + 5.91822263375044/zx                                 &
                  + 0.1174886643003278*log(zrh)                         &
                  + 0.4625315047693772*zt*log(zrh)                      &
                  - 0.01180591129059253*zt**2*log(zrh)                  &
                  + 0.0000404196487152575*zt**3*log(zrh)                &
                  + (15.79628615047088*log(zrh))/zx                     &
                  - 0.215553951893509*log(zrh)**2                       &
                  - 0.0810269192332194*zt*log(zrh)**2                   &
                  + 0.001435808434184642*zt**2*log(zrh)**2              &
                  - 4.775796947178588e-6*zt**3*log(zrh)**2              &
                  - (2.912974063702185*log(zrh)**2)/zx                  &
                  - 3.588557942822751*log(zrh)**3                       &
                  + 0.04950795302831703*zt*log(zrh)**3                  &
                  - 0.0002138195118737068*zt**2*log(zrh)**3             &
                  + 3.108005107949533e-7*zt**3*log(zrh)**3              &
                  - (0.02933332747098296*log(zrh)**3)/zx                &
                  + 1.145983818561277*log(zpcsa)                        &
                  - 0.6007956227856778*zt*log(zpcsa)                    &
                  + 0.00864244733283759*zt**2*log(zpcsa)                &
                  - 0.00002289467254710888*zt**3*log(zpcsa)             &
                  - (8.44984513869014*log(zpcsa))/zx                    &
                  + 2.158548369286559*log(zrh)*log(zpcsa)               &
                  + 0.0808121412840917*zt*log(zrh)*log(zpcsa)           &
                  - 0.0004073815255395214*zt**2*log(zrh)*log(zpcsa)     &
                  - 4.019572560156515e-7*zt**3*log(zrh)*log(zpcsa)      &
                  + (0.7213255852557236*log(zrh)*log(zpcsa))/zx         &
                  + 1.62409850488771*log(zrh)**2*log(zpcsa)             &
                  - 0.01601062035325362*zt*log(zrh)**2*log(zpcsa)       &
                  + 0.00003771238979714162*zt**2*log(zrh)**2*log(zpcsa) &
                  + 3.217942606371182e-8*zt**3*log(zrh)**2*log(zpcsa)   &
                  - (0.01132550810022116*log(zrh)**2*log(zpcsa))/zx     &
                  + 9.71681713056504*log(zpcsa)**2                      &
                  - 0.1150478558347306*zt*log(zpcsa)**2                 &
                  + 0.0001570982486038294*zt**2*log(zpcsa)**2           &
                  + 4.009144680125015e-7*zt**3*log(zpcsa)**2            &
                  + (0.7118597859976135*log(zpcsa)**2)/zx               &
                  - 1.056105824379897*log(zrh)*log(zpcsa)**2            &
                  + 0.00903377584628419*zt*log(zrh)*log(zpcsa)**2       &
                  - 0.00001984167387090606*zt**2*log(zrh)*log(zpcsa)**2 &
                  + 2.460478196482179e-8*zt**3*log(zrh)*log(zpcsa)**2   &
                  - (0.05790872906645181*log(zrh)*log(zpcsa)**2)/zx     &
                  - 0.1487119673397459*log(zpcsa)**3                    &
                  + 0.002835082097822667*zt*log(zpcsa)**3               &
                  - 9.24618825471694e-6*zt**2*log(zpcsa)**3             &
                  + 5.004267665960894e-9*zt**3*log(zpcsa)**3            &
                  - (0.01270805101481648*log(zpcsa)**3)/zx

             pnuc_rate(ii,jj) = exp(pnuc_rate(ii,jj)) ! [#/(cm3 s)]

             IF (pnuc_rate(ii,jj) < 1.e-7) THEN ! validity of parameterization
                pnuc_rate(ii,jj) = 0.
                pd_crit(ii,jj) = 1.e-9
             END IF

             !-- 4) Total number of molecules in the critical cluster -------------------------

             zntot = - 0.002954125078716302                                   &
                  - 0.0976834264241286*zt                                     &
                  + 0.001024847927067835*zt**2                                &
                  - 2.186459697726116e-6*zt**3                                &
                  - 0.1017165718716887/zx                                     &
                  - 0.002050640345231486*log(zrh)                             &
                  - 0.007585041382707174*zt*log(zrh)                          &
                  + 0.0001926539658089536*zt**2*log(zrh)                      &
                  - 6.70429719683894e-7*zt**3*log(zrh)                        &
                  - (0.2557744774673163*log(zrh))/zx                          &
                  + 0.003223076552477191*log(zrh)**2                          &
                  + 0.000852636632240633*zt*log(zrh)**2                       &
                  - 0.00001547571354871789*zt**2*log(zrh)**2                  &
                  + 5.666608424980593e-8*zt**3*log(zrh)**2                    &
                  + (0.03384437400744206*log(zrh)**2)/zx                      &
                  + 0.04743226764572505*log(zrh)**3                           &
                  - 0.0006251042204583412*zt*log(zrh)**3                      &
                  + 2.650663328519478e-6*zt**2*log(zrh)**3                    &
                  - 3.674710848763778e-9*zt**3*log(zrh)**3                    &
                  - (0.0002672510825259393*log(zrh)**3)/zx                    &
                  - 0.01252108546759328*log(zpcsa)                            &
                  + 0.005806550506277202*zt*log(zpcsa)                        &
                  - 0.0001016735312443444*zt**2*log(zpcsa)                    &
                  + 2.881946187214505e-7*zt**3*log(zpcsa)                     &
                  + (0.0942243379396279*log(zpcsa))/zx                        &
                  - 0.0385459592773097*log(zrh)*log(zpcsa)                    &
                  - 0.0006723156277391984*zt*log(zrh)*log(zpcsa)              &
                  + 2.602884877659698e-6*zt**2*log(zrh)*log(zpcsa)            &
                  + 1.194163699688297e-8*zt**3*log(zrh)*log(zpcsa)            &
                  - (0.00851515345806281*log(zrh)*log(zpcsa))/zx              &
                  - 0.01837488495738111*log(zrh)**2*log(zpcsa)                &
                  + 0.0001720723574407498*zt*log(zrh)**2*log(zpcsa)           &
                  - 3.717657974086814e-7*zt**2*log(zrh)**2*log(zpcsa)         &
                  - 5.148746022615196e-10*zt**3*log(zrh)**2*log(zpcsa)        &
                  + (0.0002686602132926594*log(zrh)**2*log(zpcsa))/zx         &
                  - 0.06199739728812199*log(zpcsa)**2                         &
                  + 0.000906958053583576*zt*log(zpcsa)**2                     &
                  - 9.11727926129757e-7*zt**2*log(zpcsa)**2                   &
                  - 5.367963396508457e-9*zt**3*log(zpcsa)**2                  &
                  - (0.007742343393937707*log(zpcsa)**2)/zx                   &
                  + 0.0121827103101659*log(zrh)*log(zpcsa)**2                 &
                  - 0.0001066499571188091*zt*log(zrh)*log(zpcsa)**2           &
                  + 2.534598655067518e-7*zt**2*log(zrh)*log(zpcsa)**2         &
                  - 3.635186504599571e-10*zt**3*log(zrh)*log(zpcsa)**2        &
                  + (0.0006100650851863252*log(zrh)*log(zpcsa)**2)/zx         &
                  + 0.0003201836700403512*log(zpcsa)**3                       &
                  - 0.0000174761713262546*zt*log(zpcsa)**3                    &
                  + 6.065037668052182e-8*zt**2*log(zpcsa)**3                  &
                  - 1.421771723004557e-11*zt**3*log(zpcsa)**3                 &
                  + (0.0001357509859501723*log(zpcsa)**3)/zx

             zntot = exp(zntot)

             !-- 5) Size of the critical cluster ----------------------------------
             pn_crit_sa(ii,jj) = zx*zntot
             pd_crit(ii,jj) = 2.e-9*exp(-1.6524245+0.42316402*zx+0.33466487*log(zntot)) ! [m]

             !-- 6) Organic compounds not involved when binary nucleation is assumed.
             pn_crit_ocnv(ii,jj) = 0.
             pk_sa(ii,jj) = 1.
             pk_ocnv(ii,jj) = 0.
             !--

             IF (pn_crit_sa(ii,jj) < 4.) THEN
                !set nucleation rate to collision rate

                !volumes of the colliding objects
                zma = 96.   ![g/mol]
                zmw = 18.   ![g/mol]
                zxmass = 1.
                !          zxmass = zxmole*zma/((1.0-zxmole)*zmw + zxmole*zma) !mass fraction of h2so4
                za = 0.7681724  +zxmass*(2.1847140   +zxmass*(7.1630022 + &
                                         zxmass*(-44.31447 + zxmass * (88.75606  + zxmass*(-75.73729+ &
                                         zxmass*23.43228)))))
                zb = 1.808225e-3+zxmass*(-9.294656e-3+zxmass*(-0.03742148+&
                                          zxmass*(0.2565321 + zxmass * (-0.5362872+ zxmass*(0.4857736- &
                                          zxmass*0.1629592)))))

                zc = -3.478524e-6+zxmass*(1.335867e-5+zxmass*(5.195706e-5+ &
                                          zxmass*(-3.717636e-4+zxmass * (7.990811e-4+zxmass*(-7.458060e-4+ &
                                          zxmass*2.58139e-4)))))

                zroo = za+zt*(zb+zc*zt) ! g/cm^3
                zroo = zroo*1.e+3 !kg/m^3

                zm1 = 0.098 ! [kg/mol]/ Na !kg
                zm2 = zm1
                zv1 = zm1/avog/zroo
                zv2 = zv1

                zcoll = zpcsa*zpcsa* &
                        (3.*pi/4.)**(1./6.) * &
                        SQRT(6.*rg*zt/zm1+6.*rg*zt/zm2) * &
                        (zv1**(1./3.) + zv2**(1./3.))**2 * &
                        1.e+6        ! m3 -> cm3

                zcoll = MIN(zcoll,1.e10)

                pnuc_rate(ii,jj) = zcoll

             ELSE

                pnuc_rate(ii,jj) = MIN(pnuc_rate(ii,jj),1.e10)

             END IF

             pnuc_rate(ii,jj) = pnuc_rate(ii,jj)*1.e6 ! [#/(m3 s)]

          END DO
       END DO
    ELSE
       pnuc_rate(:,:) = 0.
    END IF
  END SUBROUTINE binnucl

  !********************************************************************
  !********************************************************************


  SUBROUTINE ternucl(kbdim,      klev,                  &
                     pc_sa,     pc_nh3,     ptemp,        prh,     &
                     pnuc_rate, pn_crit_sa, pn_crit_ocnv, pd_crit, &
                     pk_sa,     pk_ocnv)

    IMPLICIT NONE
    !-- Input variables -------------------
    INTEGER, INTENT(IN) ::       &
         kbdim,                  & ! dimension for arrays
         klev                      ! number of vertical klev

    REAL, INTENT(IN) ::      &
         pc_sa(kbdim,klev),      & ! sulphuric acid concentration [#/cm3]
         pc_nh3(kbdim,klev),     & ! ammonia mixing ratio [ppt]
         ptemp(kbdim,klev),      & ! ambient temperature [K]
         prh(kbdim,klev)           ! relative humidity

    !-- Output variables -------------------
    REAL, INTENT(OUT) ::     &
         pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
         pn_crit_sa(kbdim,klev),  & ! number of H2SO4 molecules in cluster [1]
         pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
         pd_crit(kbdim,klev),     & ! diameter of critical cluster [m]
         pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
         pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

    !-- Local variables --------------------
    INTEGER :: ii, jj  ! loop indices

    REAL :: &
         zlnj     ! logarithm of nucleation rate

    !------------------------------------------------------------------------------------------------
    ! loops over
    DO jj = 1, klev !  vertical grid
        DO ii = 1, kbdim !  horizontal kbdim in the slab

          !-- 1) Checking that we are in the validity range of the parameterization -----------
          ! validity of parameterization : DO NOT REMOVE!
          IF (ptemp(ii,jj) < 240.)  STOP '  INVALID INPUT VALUE (ter. nucleation): ptemp(ii,jj)erature < 240 K'
          IF (ptemp(ii,jj) > 300.)  STOP '  INVALID INPUT VALUE (ter. nucleation) ptemp(ii,jj)erature > 300 K'
          IF (prh(ii,jj) < 0.05)    STOP '  INVALID INPUT VALUE (ter. nucleation) relative humidity < 5 %'
          IF (prh(ii,jj) > 0.95)    STOP '  INVALID INPUT VALUE (ter. nucleation) relative humidity > 95 %'
          IF (pc_sa(ii,jj) < 1.e4)  STOP '  INVALID INPUT VALUE (ter. nucleation) H2SO4 concentration < 10^4 1/cm3'
          IF (pc_sa(ii,jj) > 1.e9)  STOP '  INVALID INPUT VALUE (ter. nucleation) H2SO4 concentration > 10^9 1/cm3'
          IF (pc_nh3(ii,jj) < 0.1)  STOP '  INVALID INPUT VALUE (ter. nucleation) ammonia mixing ratio < 0.1 ppt'
          IF (pc_nh3(ii,jj) > 100.) STOP '  INVALID INPUT VALUE (ter. nucleation) ammonia mixing ratio > 100 ppt'

          !-- 2) Nucleation rate ---------------------------------------------------
          zlnj = -84.7551114741543 + 0.3117595133628944*prh(ii,jj)                               &
               + 1.640089605712946*prh(ii,jj)*ptemp(ii,jj)                                          &
               - 0.003438516933381083*prh(ii,jj)*ptemp(ii,jj)**2                                    &
               - 0.00001097530402419113*prh(ii,jj)*ptemp(ii,jj)**3                                  &
               - 0.3552967070274677/log(pc_sa(ii,jj))                                               &
               - (0.06651397829765026*prh(ii,jj))/log(pc_sa(ii,jj))                                 &
               - (33.84493989762471*ptemp(ii,jj))/log(pc_sa(ii,jj))                                 &
               - (7.823815852128623*prh(ii,jj)*ptemp(ii,jj))/log(pc_sa(ii,jj))                      &
               + (0.3453602302090915*ptemp(ii,jj)**2)/log(pc_sa(ii,jj))                             &
               + (0.01229375748100015*prh(ii,jj)*ptemp(ii,jj)**2)/log(pc_sa(ii,jj))                 &
               - (0.000824007160514956*ptemp(ii,jj)**3)/log(pc_sa(ii,jj))                           &
               + (0.00006185539100670249*prh(ii,jj)*ptemp(ii,jj)**3)/log(pc_sa(ii,jj))              &
               + 3.137345238574998*log(pc_sa(ii,jj))                                                &
               + 3.680240980277051*prh(ii,jj)*log(pc_sa(ii,jj))                                     &
               - 0.7728606202085936*ptemp(ii,jj)*log(pc_sa(ii,jj))                                  &
               - 0.204098217156962*prh(ii,jj)*ptemp(ii,jj)*log(pc_sa(ii,jj))                        &
               + 0.005612037586790018*ptemp(ii,jj)**2*log(pc_sa(ii,jj))                             &
               + 0.001062588391907444*prh(ii,jj)*ptemp(ii,jj)**2*log(pc_sa(ii,jj))                  &
               - 9.74575691760229e-6*ptemp(ii,jj)**3*log(pc_sa(ii,jj))                              &
               - 1.265595265137352e-6*prh(ii,jj)*ptemp(ii,jj)**3*log(pc_sa(ii,jj))                  &
               + 19.03593713032114*log(pc_sa(ii,jj))**2                                             &
               - 0.1709570721236754*ptemp(ii,jj)*log(pc_sa(ii,jj))**2                               &
               + 0.000479808018162089*ptemp(ii,jj)**2*log(pc_sa(ii,jj))**2                          &
               - 4.146989369117246e-7*ptemp(ii,jj)**3*log(pc_sa(ii,jj))**2                          &
               + 1.076046750412183*log(pc_nh3(ii,jj))                                               &
               + 0.6587399318567337*prh(ii,jj)*log(pc_nh3(ii,jj))                                   &
               + 1.48932164750748*ptemp(ii,jj)*log(pc_nh3(ii,jj))                                   &
               + 0.1905424394695381*prh(ii,jj)*ptemp(ii,jj)*log(pc_nh3(ii,jj))                      &
               - 0.007960522921316015*ptemp(ii,jj)**2*log(pc_nh3(ii,jj))                            &
               - 0.001657184248661241*prh(ii,jj)*ptemp(ii,jj)**2*log(pc_nh3(ii,jj))                 &
               + 7.612287245047392e-6*ptemp(ii,jj)**3*log(pc_nh3(ii,jj))                            &
               + 3.417436525881869e-6*prh(ii,jj)*ptemp(ii,jj)**3*log(pc_nh3(ii,jj))                 &
               + (0.1655358260404061*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))                          &
               + (0.05301667612522116*prh(ii,jj)*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))              &
               + (3.26622914116752*ptemp(ii,jj)*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))               &
               - (1.988145079742164*prh(ii,jj)*ptemp(ii,jj)*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))    &
               - (0.04897027401984064*ptemp(ii,jj)**2*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))          &
               + (0.01578269253599732*prh(ii,jj)*ptemp(ii,jj)**2*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))&
               + (0.0001469672236351303*ptemp(ii,jj)**3*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))          &
               - (0.00002935642836387197*prh(ii,jj)*ptemp(ii,jj)**3*log(pc_nh3(ii,jj)))/log(pc_sa(ii,jj))&
               + 6.526451177887659*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))                                &
               - 0.2580021816722099*ptemp(ii,jj)*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))                 &
               + 0.001434563104474292*ptemp(ii,jj)**2*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))           &
               - 2.020361939304473e-6*ptemp(ii,jj)**3*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))           &
               - 0.160335824596627*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))                          &
               + 0.00889880721460806*ptemp(ii,jj)*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))           &
               - 0.00005395139051155007*ptemp(ii,jj)**2*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))     &
               + 8.39521718689596e-8*ptemp(ii,jj)**3*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))        &
               + 6.091597586754857*log(pc_nh3(ii,jj))**2                                            &
               + 8.5786763679309*prh(ii,jj)*log(pc_nh3(ii,jj))**2                                   &
               - 1.253783854872055*ptemp(ii,jj)*log(pc_nh3(ii,jj))**2                               &
               - 0.1123577232346848*prh(ii,jj)*ptemp(ii,jj)*log(pc_nh3(ii,jj))**2                   &
               + 0.00939835595219825*ptemp(ii,jj)**2*log(pc_nh3(ii,jj))**2                          &
               + 0.0004726256283031513*prh(ii,jj)*ptemp(ii,jj)**2*log(pc_nh3(ii,jj))**2             &
               - 0.00001749269360523252*ptemp(ii,jj)**3*log(pc_nh3(ii,jj))**2                       &
               - 6.483647863710339e-7*prh(ii,jj)*ptemp(ii,jj)**3*log(pc_nh3(ii,jj))**2              &
               +(0.7284285726576598*log(pc_nh3(ii,jj))**2)/log(pc_sa(ii,jj))                        &
               + (3.647355600846383*ptemp(ii,jj)*log(pc_nh3(ii,jj))**2)/log(pc_sa(ii,jj))           &
               - (0.02742195276078021*ptemp(ii,jj)**2*log(pc_nh3(ii,jj))**2)/log(pc_sa(ii,jj))      &
               + (0.00004934777934047135*ptemp(ii,jj)**3*log(pc_nh3(ii,jj))**2)/log(pc_sa(ii,jj))   &
               + 41.30162491567873*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))**2                          &
               - 0.357520416800604*ptemp(ii,jj)*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))**2             &
               + 0.000904383005178356*ptemp(ii,jj)**2*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))**2       &
               - 5.737876676408978e-7*ptemp(ii,jj)**3*log(pc_sa(ii,jj))*log(pc_nh3(ii,jj))**2       &
               - 2.327363918851818*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))**2                       &
               + 0.02346464261919324*ptemp(ii,jj)*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))**2        &
               - 0.000076518969516405*ptemp(ii,jj)**2*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))**2    &
               + 8.04589834836395e-8*ptemp(ii,jj)**3*log(pc_sa(ii,jj))**2*log(pc_nh3(ii,jj))**2     &
               - 0.02007379204248076*log(prh(ii,jj))                                                &
               - 0.7521152446208771*ptemp(ii,jj)*log(prh(ii,jj))                                    &
               + 0.005258130151226247*ptemp(ii,jj)**2*log(prh(ii,jj))                               &
               - 8.98037634284419e-6*ptemp(ii,jj)**3*log(prh(ii,jj))                                &
               + (0.05993213079516759*log(prh(ii,jj)))/log(pc_sa(ii,jj))                            &
               + (5.964746463184173*ptemp(ii,jj)*log(prh(ii,jj)))/log(pc_sa(ii,jj))                 &
               - (0.03624322255690942*ptemp(ii,jj)**2*log(prh(ii,jj)))/log(pc_sa(ii,jj))            &
               + (0.00004933369382462509*ptemp(ii,jj)**3*log(prh(ii,jj)))/log(pc_sa(ii,jj))         &
               - 0.7327310805365114*log(pc_nh3(ii,jj))*log(prh(ii,jj))                              &
               - 0.01841792282958795*ptemp(ii,jj)*log(pc_nh3(ii,jj))*log(prh(ii,jj))                &
               + 0.0001471855981005184*ptemp(ii,jj)**2*log(pc_nh3(ii,jj))*log(prh(ii,jj))           &
               - 2.377113195631848e-7*ptemp(ii,jj)**3*log(pc_nh3(ii,jj))*log(prh(ii,jj))

          pnuc_rate(ii,jj) = exp(zlnj) ! [#/(cm3 s)]

          IF (pnuc_rate(ii,jj) < 1.e-5) THEN ! validity of parametrization
             pnuc_rate(ii,jj) = 0.
             pd_crit = 1.e-9
             CYCLE
          END IF
          ! validity of parametrization
          IF (pnuc_rate(ii,jj) > 1.e6) STOP '  INVALID OUTPUT VALUE for ternary nucleation: nucleation rate > 10^6 1/cm3s'

          pnuc_rate(ii,jj) = pnuc_rate(ii,jj)*1.e6 ! [#/(m3 s)]
          !-- 3) Number of H2SO4 molecules in a critical cluster -----------------------
          pn_crit_sa(ii,jj) = 38.16448247950508 + 0.7741058259731187*zlnj +                      &
                              0.002988789927230632*zlnj**2 - 0.3576046920535017*ptemp(ii,jj) -   &
                              0.003663583011953248*zlnj*ptemp(ii,jj) + 0.000855300153372776*ptemp(ii,jj)**2

          pn_crit_sa(ii,jj) = MAX(pn_crit_sa(ii,jj),2.e0) ! kinetic limit: at least 2 H2SO4 molecules in a cluster

          !-- 4) Size of the critical cluster ------------------------------------------
          pd_crit(ii,jj) = 0.1410271086638381 - 0.001226253898894878*zlnj -                      &
                           7.822111731550752e-6*zlnj**2 - 0.001567273351921166*ptemp(ii,jj) -    &
                           0.00003075996088273962*zlnj*ptemp(ii,jj) + 0.00001083754117202233*ptemp(ii,jj)**2  ! [nm]

          pd_crit(ii,jj) = pd_crit(ii,jj)*2.e-9  ! [m]

          !-- 5) Organic compounds not involved when ternary nucleation is assumed.
          pn_crit_ocnv(ii,jj) = 0.
          pk_sa(ii,jj) = 1.
          pk_ocnv(ii,jj) = 0.

       END DO
    END DO

  END SUBROUTINE ternucl

  !****************************************************************************
  !****************************************************************************

  !**************************************************************
  !
  ! Below the following assumption have been made:
  !
  !  nucrate = coagcoeff*zpcsa**2
  !  coagcoeff = 8*sqrt(3*boltz*ptemp*r_abs/dens_abs)
  !  r_abs = 0.315d-9 radius of bisulphate molecule [m]
  !  dens_abs = 1465  density of - " - [kg/m3]
  !
  !**************************************************************

  SUBROUTINE kinnucl(kbdim,        klev,          &
                     pc_sa,      ptemp,                       &
                     pnuc_rate,  pd_crit,      ppbl,          &
                     pn_crit_sa, pn_crit_ocnv, pk_sa, pk_ocnv)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) :: &
          kbdim,   &    ! dimension for arrays
          klev          ! number of vertical klev

     REAL, INTENT(IN) :: &
          pc_sa(kbdim,klev), & ! sulphuric acid concentration [#/m3]
          ptemp(kbdim,klev)    ! ambient temperature [K]

     REAL, INTENT(OUT) ::      &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! critical diameter of clusters [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

     INTEGER :: ppbl(kbdim)                      !boundary layer top

     !-- Local variables --------------------
     INTEGER :: ii, jj     ! loop indices
     INTEGER :: zpbl(kbdim)

     !-------------------------------------------------------------------------
     zpbl(:) = ppbl(:)

     ! loops ove
     DO ii = 1, kbdim !  horizontal kbdim in the slab
        DO jj = zpbl(ii), klev !  vertical grid

           !  pnuc_rate(ii,jj) = 1.19386e-17*sqrt(ptemp(ii,jj))*pc_sa(ii,jj)**2 ! [#/(m3 s)]
           pnuc_rate(ii,jj) = 5.0e-13*pc_sa(ii,jj)**2*1.E6 ! [#/(m3 s)]

           !-- Organic compounds not involved when kinetic nucleation is assumed.
           pn_crit_sa(ii,jj) = 2.
           pn_crit_ocnv(ii,jj) = 0.
           pk_sa(ii,jj) = 1.
           pk_ocnv(ii,jj) = 0.
        END DO
     END DO

     pd_crit = 7.9375e-10 ! [m]

   END SUBROUTINE kinnucl

   !********************************************************************
   !********************************************************************

   SUBROUTINE actnucl( kbdim,        klev,                   &
                       psa_conc,   pnuc_rate,    pd_crit, ppbl,          &
                       pn_crit_sa, pn_crit_ocnv, pk_sa,   pk_ocnv, activ)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) ::   &
          kbdim,              & ! dimension for arrays
          klev                  ! number of vertical klev

     REAL, INTENT(IN) ::  psa_conc(kbdim,klev) ! sulphuric acid concentration [#/m3]
     REAL, INTENT(IN) :: activ
     INTEGER :: ppbl(kbdim)                      !boundary layer top

     REAL, INTENT(OUT) ::                         &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! critical diameter of clusters [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

     !-- Local variables --------------------
     INTEGER :: ii, jj                                  ! loop indices
     INTEGER :: zpbl(kbdim)

     !-------------------------------------------------------------------------
     ! loops over
     !  previous loop structure
     !    DO ii = 1,kbdim !  horizontal kbdim in the slab
     !       DO jj = 1,klev !  vertical grid
     !          pnuc_rate(ii,jj) = 1.e-7*psa_conc(ii,jj) ! [#/(m3 s)]
     !       END DO
     !    END DO
     !    pd_crit = 7.9375e-10 ! [m]

     zpbl(:) = ppbl(:)

     DO ii = 1, kbdim !  horizontal kbdim in the slab
        DO jj = zpbl(ii), klev !  vertical grid
           !  gone through for boundary layer klev only!
           ! act_coeff 1e-7 by default, namelist controllable.
           pnuc_rate(ii,jj) = activ*psa_conc(ii,jj) ! [#/(m3 s)]

           !-- Organic compounds not involved when kinetic nucleation is assumed.
           pn_crit_sa(ii,jj) = 2.
           pn_crit_ocnv(ii,jj) = 0.
           pk_sa(ii,jj) = 1.
           pk_ocnv(ii,jj) = 0.
        END DO
     END DO

     pd_crit = 7.9375e-10 ! [m]

   END SUBROUTINE actnucl

   !********************************************************************
   !********************************************************************
   ! ---------------------------------------------------------------
   ! Paasonen et al. (2010) determined particle formation rates for
   ! 2 nm particles, j2, from different kind of combinations of
   ! sulphuric acid and organic matter concentration
   !
   ! ORGNUCL scheme conciders only the organic matter in nucleation
   ! ---------------------------------------------------------------

   SUBROUTINE orgnucl(kbdim,        klev,            &
                      pc_org,     pnuc_rate,    pd_crit, ppbl,   &
                      pn_crit_sa, pn_crit_ocnv, pk_sa,   pk_ocnv)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) :: &
          kbdim,   &     ! dimension for arrays
          klev           ! number of vertical klev

     REAL, INTENT(IN) :: &
          pc_org(kbdim,klev)    ! organic vapour  concentration [#/m3]

     !-- Output variables ------------------
     REAL, INTENT(OUT) ::      &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! diameter of critical cluster [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

     INTEGER :: ppbl(kbdim)             !boundary layer top

     !-- Local variables --------------------
     INTEGER :: ii, jj     ! loop indices

     INTEGER :: zpbl(kbdim)
     REAL    :: Aorg = 1.9e-7 ![1/s] [fxm] (Paasonen et al. Table 3.????)

     !-------------------------------------
     zpbl(:) = ppbl(:)
     ! loops over
     DO ii = 1, kbdim !  horizontal kbdim in the slab
        DO jj = zpbl(ii), klev !  vertical grid
           pnuc_rate(ii,jj) = Aorg*pc_org(ii,jj)
           !          pnuc_rate(ii,jj) = Korg*pc_org(ii,jj)**2       ! homomolecular nuleation - which one?

           !-- H2SO4 not involved when pure organic nucleation is assumed.
           pn_crit_sa(ii,jj) = 0.
           pn_crit_ocnv(ii,jj) = 1.
           pk_sa(ii,jj) = 0.
           pk_ocnv(ii,jj) = 1.
        END DO
     END DO

     pd_crit = 1.5E-9 ! [m]

   END SUBROUTINE orgnucl

   !********************************************************************
   !********************************************************************

   ! ---------------------------------------------------------------
   ! Nucleation by Paasonen et al. 2010
   !
   ! SUMNUCL scheme conciders both the organic vapor and H2SO4 in
   ! nucleation - activation type of nucleation
   ! ---------------------------------------------------------------

   SUBROUTINE sumnucl(kbdim,        klev,                     &
                      pc_sa,      pc_org,       pnuc_rate, pd_crit, ppbl, &
                      pn_crit_sa, pn_crit_ocnv, pk_sa,     pk_ocnv)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) :: &
          kbdim,  &      ! dimension for arrays
          klev           ! number of vertical levels

     REAL, INTENT(IN) ::  &
          pc_org(kbdim,klev), &  ! organic vapour concentration [#/m3]
          pc_sa(kbdim,klev)      ! sulphuric acid concentration [#/m3]

     !-- Output variables ------------------
     REAL, INTENT(OUT) :: &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! diameter of critical cluster [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.


     INTEGER :: ppbl(kbdim)             !boundary layer top

     !-- Local variables --------------------
     INTEGER :: ii, jj     ! loop indices
     INTEGER :: zpbl(kbdim)
     REAL    :: As1 = 7.3e-7, As2 = 0.61e-7 ![1/s] (Paasonen et al. Table 3.)

     !-------------------------------------
     zpbl(:) = ppbl(:)
    DO ii = 1, kbdim !  horizontal kbdim in the slab
        DO jj = zpbl(ii), klev !  vertical grid
           pnuc_rate(ii,jj) = As1*pc_sa(ii,jj)+As2*pc_org(ii,jj) ![#/m3/s]

           !-- Both Organic compounds and H2SO4 are involved when SUMnucleation is assumed.
           pn_crit_sa(ii,jj) = 1.
           pn_crit_ocnv(ii,jj) = 1.
           pk_sa(ii,jj) = 1.
           pk_ocnv(ii,jj) = 1.
        END DO
     END DO

     pd_crit = 1.5E-9

   END SUBROUTINE sumnucl

   !********************************************************************
   !********************************************************************

   ! ---------------------------------------------------------------
   ! Nucleation by Paasonen et al. 2010
   !
   ! HETNUCL scheme conciders both the organic vapor and H2SO4 in
   ! nucleation - heteromolecular nucleation
   ! ---------------------------------------------------------------

   SUBROUTINE hetnucl(kbdim,        klev,          &
                      pc_sa,      pc_org,                      &
                      pnuc_rate,  pd_crit,      ppbl,          &
                      pn_crit_sa, pn_crit_ocnv, pk_sa, pk_ocnv)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) :: &
          kbdim,  &      ! dimension for arrays
          klev           ! number of vertical klev

     REAL, INTENT(IN) ::    &
          pc_org(kbdim,klev),   & ! organic vapour concentration [#/cm3]
          pc_sa(kbdim,klev)       ! sulphuric acid concentration [#/cm3]

     !-- Output variables ------------------
     REAL, INTENT(OUT) ::      &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! diameter of critical cluster [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

     INTEGER :: ppbl(kbdim)             !boundary layer top

     !-- Local variables --------------------
     INTEGER :: ii, jj     ! loop indices
     INTEGER :: zpbl(kbdim)
     REAL    :: zKhet = 5.8e-14 ![cm3/s] (Paasonen et al. Table 3.)

     !-------------------------------------
     zpbl(:) = ppbl(:)
     DO ii = 1, kbdim         !  horizontal kbdim in the slab
        DO jj = zpbl(ii), klev !  vertical grid
           pnuc_rate(ii,jj) = zKhet*pc_sa(ii,jj)*pc_org(ii,jj)*1.E6 ![#/m3/s]

           !-- Both Organic compounds and H2SO4 are involved when
           !-- heteromolecular nucleation is assumed.
           pn_crit_sa(ii,jj) = 1.
           pn_crit_ocnv(ii,jj) = 1.
           pk_sa(ii,jj) = 1.
           pk_ocnv(ii,jj) = 1.
        END DO
     END DO

     pd_crit = 1.5E-9 ![m]

   END SUBROUTINE hetnucl

   !********************************************************************
   !********************************************************************

   ! ---------------------------------------------------------------
   ! Nucleation by Paasonen et al. 2010
   !
   ! SANUCL scheme takes into account the homomolecular nucleation of
   ! sulphuric acid H2SO4 with both of the available vapours
   ! ---------------------------------------------------------------

   SUBROUTINE SAnucl(kbdim,        klev,          &
                     pc_sa,      pc_org,                      &
                     pnuc_rate,  pd_crit,      ppbl,          &
                     pn_crit_sa, pn_crit_ocnv, pk_sa, pk_ocnv)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) :: &
          kbdim,   &     ! dimension for arrays
          klev           ! number of vertical klev

     REAL, INTENT(IN) ::    &
          pc_org(kbdim,klev),   & ! organic vapour concentration [#/cm3]
          pc_sa(kbdim,klev)       ! sulphuric acid concentration [#/cm3]

     !-- Output variables ------------------
     REAL, INTENT(OUT) ::      &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! diameter of critical cluster [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

     INTEGER :: ppbl(kbdim)             !boundary layer top

     !-- Local variables --------------------
     INTEGER :: ii, jj     ! loop indices
     INTEGER :: zpbl(kbdim)
     REAL    :: zKsa1 = 0.82e-14, zKsa2 = 7.0e-14 ![cm3/s] (Paasonen et al. Table 3.)

     !-------------------------------------
     zpbl(:) = ppbl(:)
     DO ii = 1, kbdim !  horizontal kbdim in the slab
        DO jj = zpbl(ii), klev !  vertical grid
           pnuc_rate(ii,jj) = (zKsa1*pc_sa(ii,jj)**2 + zKsa2*pc_sa(ii,jj)*pc_org(ii,jj))*1.E6 ![#/m3/s]

           !-- Both Organic compounds and H2SO4 are involved when SAnucleation is assumed.
           pn_crit_sa(ii,jj) = 3.
           pn_crit_ocnv(ii,jj) = 1.
           pk_sa(ii,jj) = 1.
           pk_ocnv(ii,jj) = 1.
        END DO
     END DO
     pd_crit = 1.5E-9 ![m]

   END SUBROUTINE SAnucl

   !********************************************************************
   !********************************************************************

   ! ---------------------------------------------------------------
   ! Nucleation by Paasonen et al. 2010
   !
   ! SAORGNUCL scheme takes into account the homomolecular nucleation of
   ! both sulphuric acid and organic with heteromolecular nucleation
   ! ---------------------------------------------------------------

   SUBROUTINE SAORGnucl(kbdim,        klev,          &
                        pc_sa,      pc_org,                      &
                        pnuc_rate,  pd_crit,      ppbl,          &
                        pn_crit_sa, pn_crit_ocnv, pk_sa, pk_ocnv)

     IMPLICIT NONE

     !-- Input variables -------------------
     INTEGER, INTENT(IN) :: &
          kbdim,   &     ! dimension for arrays
          klev           ! number of vertical klev

     REAL, INTENT(IN) ::    &
          pc_org(kbdim,klev),   & ! organic vapour concentration [#/cm3]
          pc_sa(kbdim,klev)       ! sulphuric acid concentration [#/cm3]

     !-- Output variables ------------------
     REAL, INTENT(OUT) ::      &
          pnuc_rate(kbdim,klev),   & ! nucleation rate [#/(m3 s)]
          pd_crit(kbdim,klev),     & ! diameter of critical cluster [m]
          pn_crit_sa(kbdim,klev),  & ! number of sulphuric acid molecules in cluster [1]
          pn_crit_ocnv(kbdim,klev),& ! number of organic molecules in cluster [1]
          pk_sa(kbdim,klev),       & ! Lever: If pk_sa = 1, h2so4 is involved in nucleation.
          pk_ocnv(kbdim,klev)        ! Lever: If pk_ocnv = 1, organic compounds are involved in nucleation.

     INTEGER :: ppbl(kbdim)             !boundary layer top

     !-- Local variables --------------------
     INTEGER :: ii, jj     ! loop indices
     INTEGER :: zpbl(kbdim)
     REAL    :: zKs1 = 0.92e-14, zKs2 = 7.0e-14, zKs3 = 0.098e-14 ![cm3/s] (Paasonen et al. Table 3.)
     !-------------------------------------
     zpbl(:) = ppbl(:)
     DO ii = 1, kbdim !  horizontal points in the slab
        DO jj = zpbl(ii), klev !  vertical grid
           pnuc_rate(ii,jj) = (zKs1*pc_sa(ii,jj)**2 + &
                               zKs2*pc_sa(ii,jj)*pc_org(ii,jj)+ zKs3*pc_org(ii,jj)**2)*1.E6  ![#/m3/s]

           !-- Organic compounds not involved when kinetic nucleation is assumed.
           pn_crit_sa(ii,jj) = 3.
           pn_crit_ocnv(ii,jj) = 3.
           pk_sa(ii,jj) = 1.
           pk_ocnv(ii,jj) = 1.

        END DO
     END DO

     pd_crit = 1.5E-9 ![m]

   END SUBROUTINE SAORGnucl

   !********************************************************************
   !********************************************************************

   !----------------------------------------------------------------
   ! Function zNnuc_tayl is connected to the calculation of
   ! self-coagualtion of small particles. It calculates number of
   ! the particles in the size range [zdcrit,dx] using Taylor-expansion
   ! (please note that the expansion is not valid for certain
   ! rational numbers, e.g. -4/3 and -3/2)
   !----------------------------------------------------------------

   FUNCTION zNnuc_tayl(d1,dx,zm_para,zjnuc_t,zeta,zGRtot) RESULT(zNnuc_taylor)
     IMPLICIT NONE
     INTEGER :: i
     REAL :: &
          d1, dx,zjnuc_t, zeta, &
          term1, term2, term3, term4, term5, &
          zNnuc_taylor, zGRtot, zm_para

     zNnuc_taylor= 0.

     DO i = 0, 29

        IF (i == 0 .OR. i == 1) THEN
           term1 = 1.
        ELSE
           term1 = term1*REAL(i)
        END IF
        term2 = (REAL(i)*(zm_para+1.)+1.)*term1
        term3 = zeta**i
        term4 = term3/term2
        term5 = REAL(i)*(zm_para+1.)+1.
        zNnuc_taylor = zNnuc_taylor + term4*(dx**term5 - d1**term5)

     END DO

     zNnuc_taylor = zNnuc_taylor*zjnuc_t*exp(-zeta*(d1**(zm_para+1)))/zGRtot

   END FUNCTION zNnuc_tayl

   !--------------------------------------------------------

END MODULE mo_salsa_nucleation
