MODULE mo_salsa_cloud
   IMPLICIT NONE

  !*********************************************************
  !  MO_SALSA_CLOUD
  !*********************************************************
  !
  ! Purpose: Calculates the number of activated cloud droplets
  !
  !
  ! Interface:
  ! ----------
  ! Called from the main aerosol model
  !
  !
  ! Coded by:
  ! ---------
  ! T. Anttila (FMI)     2007
  ! H. Kokkola (FMI)     2007
  ! A.-I. Partanen (FMI) 2007
  !
  !*********************************************************

CONTAINS

   SUBROUTINE cloud_activation(kbdim, klev,   &
                               temp,   pres,  rv,     &
                               rs,     w,     paero,  &
                               pcloud, pactd          )

      USE mo_submctl, ONLY : t_section, fn2b, ncld, &
                             lsactintst, lsactbase

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::              &
         kbdim,                     & ! dimension for arrays
         klev                       ! number of vertical levels

      REAL, INTENT(in) ::             &
         pres(kbdim,klev),          &
         temp(kbdim,klev),          &
         w(kbdim,klev)

      REAL, INTENT(inout) :: rv(kbdim,klev) ! Water vapor mixing ratio
      REAL, INTENT(in)    :: rs(kbdim,klev) ! Saturation vapor mixing ratio

      TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld),  &
                                        paero(kbdim,klev,fn2b)

      ! Properties of newly activate particles
      TYPE(t_section), INTENT(out) :: pactd(kbdim,klev,ncld)

      INTEGER :: ii, jj, kk

      ! This is needed for cloud base activation, but must be set to zero for interstitial activation
      DO jj = 1, klev    ! vertical grid
         DO ii = 1, kbdim ! horizontal grid
            ! Reset activated
            DO kk = 1, ncld
               pactd(ii,jj,kk)%volc(:) = 0.
               pactd(ii,jj,kk)%numc = 0.
            END DO
         END DO
      END DO

      ! -------------------------------------
      ! Interstitial activation
      ! -------------------------------------
      IF ( lsactintst ) THEN

         CALL actInterst(kbdim,klev,paero,pcloud,rv,rs,temp)

      END IF

      ! -----------------------------------
      ! Activation at cloud base
      ! -----------------------------------
      IF ( lsactbase ) THEN

         CALL ActCloudBase(kbdim,klev,paero,pres,temp,w,pactd)

      END IF

   END SUBROUTINE cloud_activation


   ! -----------------------------------------------------------------
   ! Calculates the number of moles of dissolved solutes in one particle
   !
   SUBROUTINE getSolute(kbdim,klev,paero,pns)
      USE mo_submctl, ONLY : t_section,nlim,       &
                             fn2b,           &
                             rhosu, rhooc,   &
                             rhonh, rhono,   &
                             rhoss,          &
                             msu, moc,       &
                             mnh, mno,       &
                             mss
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: kbdim,klev
      TYPE(t_section), INTENT(IN) :: paero(kbdim,klev,fn2b)
      REAL, INTENT(OUT) :: pns(kbdim,klev,fn2b)

      INTEGER :: ii,jj,kk

      pns = 0.

      DO jj = 1, klev

         DO ii = 1, kbdim

            !-- subranges 1a, 2a and 2b

            DO kk = 1, fn2b
               IF (paero(ii,jj,kk)%numc > nlim) THEN

                  !-- number of moles of solute in one particle [mol]
                  !   BC and dust are insoluble - ignored
                  !   SS or NaCl produces 2 ions
                  !   SO or H2SO4 produces 3 ions
                  pns(ii,jj,kk) = (3.*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &
                                  paero(ii,jj,kk)%volc(2) * rhooc/moc +      &
                                  paero(ii,jj,kk)%volc(6) * rhono/mno +      &
                                  paero(ii,jj,kk)%volc(7) * rhonh/mnh +      &
                                  2.*paero(ii,jj,kk)%volc(5) * rhoss/mss)/   &
                                  paero(ii,jj,kk)%numc

               END IF
            END DO

         END DO

      END DO


   END SUBROUTINE getSolute


   ! -----------------------------------------------------------------

   SUBROUTINE ActCloudBase(kbdim,klev,paero,pres,temp,w,pactd)
      ! Cloud base activation following:
      !
      ! Abdul-Razzak et al: "A parameterization of aerosol activation -
      !                      3. Sectional representation"
      !                      J. Geophys. Res. 107, 10.1029/2001JD000483, 2002.
      !                      [Part 3]
      !
      ! Abdul Razzak et al: "A parameterization of aerosol activation -
      !                      1. Single aerosol type"
      !                      J. Geophys. Res. 103, 6123-6130, 1998.
      !                      [Part 1]
      !
      !
      ! Note: updates pactd, but does not change pcloud?
      ! Note: insoluble species are not properly accounted for
      !

      USE mo_submctl, ONLY :     &
         rg,                             & ! molar gas constant [J/(mol K)]
         surfw0,                       & ! surface tension of water [J/m2]
         nlim,                           & ! lowest possible particle conc. in a bin [#/m3]
         rhowa, mwa,                & ! Density and molar mass of water
         pi,                             &
         cpa, mair,                    & ! Air properties
         in1a,in2b,fn2a, fn2b,      & ! size regime bin indices
         t_section,                    & ! Data type for cloud/rain drops
         grav,                      & ! Standard acceleration due to gravity
         ncld                               ! Total number of cloud bins

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::              &
         kbdim,                     & ! dimension for arrays
         klev                       ! number of vertical levels

      REAL, INTENT(in) ::             &
         pres(kbdim,klev),          &
         temp(kbdim,klev),          &
         w(kbdim,klev)

      TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,fn2b)

      ! Properties of newly activate particles
      TYPE(t_section), INTENT(out) :: pactd(kbdim,klev,ncld)

      !-- local variables --------------
      INTEGER :: ii, jj, kk             ! loop indices

      INTEGER :: bcrita(kbdim,klev),bcritb(kbdim,klev) ! Index of the critical aerosol bin for regimes a and b

      REAL ::                         &
         sil,                       & ! critical supersaturation
                                      !     at the upper bound of the bin
         siu,                       & !  "  at the lower bound of the bin
         scrit(fn2b),               & !  "  at the center of the bin
         aa,                        & ! curvature (Kelvin) effect [m]
         bb,                        & ! solute (Raoult) effect [m3]
         ns(kbdim,klev,fn2b),                  & ! number of moles of solute
         nshi,                      & !             " at the upper bound of the bin
         nslo,                      & !             " at the lower bound of the bin
         s_max,                     & ! maximum supersaturation
         s_eff,                     & ! effective supersaturation
         x, x1, a1, sum1,       & ! technical variables
         ka1,                       & ! thermal conductivity
         dv1,                       & ! diffusion coefficient
         Gc,                        & ! growth coefficient
         alpha,                     & ! see Abdul-Razzak and Ghan, part 3
         gamma,                     & ! see Abdul-Razzak and Ghan, part 3
         L,                         & ! latent heat of evaporation
         ps,                        & ! saturation vapor pressure of water [Pa]
         khi,                       & ! see Abdul-Razzak and Ghan, part 3
         theta,                     & ! see Abdul-Razzak and Ghan, part 3
         frac(kbdim,klev,fn2b),                & ! fraction of activated droplets in a bin
         ntot,                      & ! total number conc of particles [#/m3]
         zdcrit(kbdim,klev,fn2b),   & ! critical diameter [m]
         zdcstar(kbdim,klev),                   & ! Critical diameter corresponding to Smax
         zdcrhi(kbdim,klev,fn2b),   & ! Critical diameter at the high end of a bin
         zdcrlo(kbdim,klev,fn2b),   & ! Critical diameter at the low end of a bin
         V,                         & ! updraft velocity [m/s]
         rref                       ! reference radius [m]

      ! ------------------------------------------------------------------
      ! Initialization
      ! ------------------------------------------------------------------

      bb = 6.*mwa/(pi*rhowa)             ! Raoult effect [m3]
                                            ! NOTE!
                                            ! bb must be multiplied
                                            ! by the number of moles of
                                            ! solute
      zdcrit(:,:,:) = 0.
      zdcrlo(:,:,:) = 0.
      zdcrhi(:,:,:) = 0.
      frac(:,:,:) = 0.
      bcrita(:,:) = fn2a
      bcritb(:,:) = fn2b

      ! Get moles of solute at the middle of the bin
      CALL getSolute(kbdim,klev,paero,ns)

      ! ----------------------------------------------------------------

      DO jj = 1, klev    ! vertical grid
         DO ii = 1, kbdim ! horizontal grid
            ! Reset activated - done already
            !   pactd(ii,jj,kk)%volc(:) = 0.
            !   pactd(ii,jj,kk)%numc = 0.

            ! Positive updraft velocity required
            IF (w(ii,jj) <= 0.) CYCLE

            aa = 4.*mwa*surfw0/(rg*rhowa*temp(ii,jj)) ! Kelvin effect [m]

            x  = 4.*aa**3/(27.*bb)

            ! Get the critical supersaturation for aerosol bins for necessary summation terms (sum1 & ntot)

            scrit(in1a:fn2b) = exp(sqrt(x/max(epsilon(1.0),ns(ii,jj,in1a:fn2b)))) - 1.

            !-- sums in equation (8), part 3
            ntot = SUM(paero(ii,jj,in1a:fn2b)%numc)
            sum1 = SUM(paero(ii,jj,in1a:fn2b)%numc/scrit(in1a:fn2b)**(2./3.))

            IF(ntot < nlim) CYCLE
            V = w(ii,jj)

            !-- latent heat of evaporation [J/kg]
            L = 2.501e6-2370.*(temp(ii,jj)-273.15)

            !-- saturation vapor pressure of water [Pa]
            a1 = 1.-(373.15/temp(ii,jj))
            ps = 101325.*                                                 &
                 exp(13.3185*a1-1.976*a1**2-0.6445*a1**3-0.1299*a1**4)

            !-- part 1, eq (11)
            alpha = grav*mwa*L/(cpa*rg*temp(ii,jj)**2)-                            &
                    grav*mair/(rg*temp(ii,jj))

            !-- part 1, eq (12)
            gamma = rg*temp(ii,jj)/(ps*mwa) &
                    + mwa*L**2/(cpa*pres(ii,jj)*mair*temp(ii,jj))

            !-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
            x1  = pres(ii,jj) / 101325.
            dv1 = 1.e-4 * (0.211/x1) * ((temp(ii,jj)/273.)**1.94)

            rref = 10.e-9

            !-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
            ka1 = 1.e-3 * (4.39 + 0.071 * temp(ii,jj))

            !-- growth coefficient, part 1, eq (16)
            !-- (note: here uncorrected diffusivities and conductivities are used
            !    based on personal communication with H. Abdul-Razzak, 2007)
            Gc = 1./(rhowa*rg*temp(ii,jj)/(ps*dv1*mwa) +                      &
                 L*rhowa/(ka1*temp(ii,jj)) * (L*mwa/(temp(ii,jj)*rg)-1.))

            !-- effective critical supersaturation: part 3, eq (8)

            s_eff = sqrt( (ntot/sum1) **3)

            !-- part 3, equation (5)

            theta = ( sqrt( (alpha*V/Gc)**3 ) ) /(2.*pi*rhowa*gamma*ntot)

            !-- part 3, equation (6)
            khi = (2./3.)*aa*sqrt(alpha*V/Gc)

            !-- maximum supersaturation of the air parcel: part 3, equation (9)

            s_max = s_eff / sqrt(  0.5* sqrt( (khi/theta)**3 ) &
                    +  ((s_eff**2)/(theta+3.*khi))**(0.75)  )

            !-- Juha: Get the critical diameter corresponding to the maximum supersaturation
            zdcstar = 2.*aa/(3.*s_max)

            DO kk = in1a, fn2b

               IF (paero(ii,jj,kk)%numc < nlim) CYCLE

               !-- moles of solute in particle at the upper bound of the bin
               nshi = ns(ii,jj,kk)*paero(ii,jj,kk)%vratiohi

               !-- critical supersaturation
               sil = exp(sqrt(x/nshi)) - 1.

               IF(s_max < sil) CYCLE

               !-- moles of solute at the lower bound of the bin:
               nslo = ns(ii,jj,kk)*paero(ii,jj,kk)%vratiolo

               !-- critical supersaturation
               siu = exp(sqrt(x/nslo)) - 1.

               !-- fraction of activated in a bin, eq (13), part 3
               frac(ii,jj,kk) = min(1.,log(s_max/sil)/log(siu/sil))

               !-- Critical diameters for each bin and bin edges
               zdcrlo(ii,jj,kk) = 2.*sqrt(3.*nslo*bb/aa)
               zdcrhi(ii,jj,kk) = 2.*sqrt(3.*nshi*bb/aa)
               zdcrit(ii,jj,kk) = 2.*sqrt(3.*ns(ii,jj,kk)*bb/aa)

            END DO ! kk

            ! Find critical bin
            DO kk = in1a, fn2a
               IF (frac(ii,jj,kk) < 1. .AND. frac(ii,jj,kk) > 0.) THEN
                  bcrita(ii,jj) = kk
                  EXIT
               END IF
            END DO
            DO kk = in2b, fn2b
               IF (frac(ii,jj,kk) < 1. .AND. frac(ii,jj,kk) > 0.) THEN
                  bcritb(ii,jj) = kk
                  EXIT
               END IF
            END DO

         END DO ! ii

      END DO ! jj

      CALL activate3(kbdim, klev, paero, bcrita,  &
                     zdcrit, zdcrlo, zdcrhi, zdcstar, pactd  )

   END SUBROUTINE ActCloudBase


   SUBROUTINE actInterst(kbdim,klev,paero,pcloud,prv,prs,temp)
      !
      ! Activate interstitial aerosols if they've reached their critical size
      !
      ! 1. Formulate the profiles of Dwet within bins
      !      - Get the slopes between adjacent bin mids (known)
      !      - Interpolate Dwet to bin edges
      ! 2. Based on Dcrit and Dwet slopes, estimate the Ddry where Dwet becomes > Dcrit
      ! 3. Formulate the slopes for number concentration
      ! 4. Use the Dry limits from (2) as the integration limits if they are defined
      !
      ! Note: insoluble species are not properly accounted for
      !
      USE mo_submctl, ONLY :     &
         rg,                             & ! molar gas constant [J/(mol K)]
         surfw0,                       & ! surface tension of water [J/m2]
         t_section,nlim,pi6,         &
         ica,fca,icb,fcb,             &
         mwa, rhowa,                &
         in1a,fn2a,in2b,fn2b,       &
         nbins,ncld
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: kbdim, klev
      TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins),  &
                                        pcloud(kbdim,klev,ncld)
      REAL, INTENT(IN) :: prv(kbdim,klev), prs(kbdim,klev)  ! Water vapour and saturation mixin ratios
      REAL, INTENT(in) :: temp(kbdim,klev)  ! Absolute temperature

      TYPE(t_section) :: pactd(ncld) ! Local variable

      REAL :: paa        ! Coefficient for Kelvin effect

      REAL :: zdcstar, zvcstar  ! Critical diameter/volume corresponding to S_LES
      REAL :: zactvol           ! Total volume of the activated particles

      REAL :: Nact, Vact(8)        ! Helper variables for transferring the activated particles

      REAL :: Nmid, Nim1, Nip1     ! Bin number concentrations in current and adjacent bins
      REAL :: dNmid, dNim1, dNip1  ! Density function value of the number distribution for current and adjacent bins

      REAL :: Vmid, Vim1, Vip1     ! Dry particle volume in the middle of the bin
      REAL :: Vlo, Vhi             ! Dry particle volume scaled to bin edges
      REAL :: Vlom1, Vhim1         ! - '' - For adjacent bins
      REAL :: Vlop1, Vhip1         !

      REAL :: Vwmid, Vwim1, Vwip1  ! Wet particle volume in the middle of the bin
      REAL :: Vwlo, Vwhi           ! Wet particle volume at bin edges

      REAL :: zs1, zs2          ! Slopes for number concetration distributions within bins

      REAL :: N01, N02          ! Origin values for number distribution slopes
      REAL :: V01, V02          ! Origin values for wet particle volume slopes
      REAL :: Nnorm, Vnorm      ! Normalization factors for number and volume integrals

      REAL    :: vcut, vint1, vint2  ! cut volume, integration limit volumes
      LOGICAL :: intrange(4)         ! Logical table for integration ranges depending on the shape of the wet size profile:
                                     ! [Vlo -- vint1][vint1 -- Vmid][Vmid -- vint2][vint1 -- Vhi]
      INTEGER :: cb, ab, ii, jj, ss


      DO jj = 1, klev
         DO ii = 1, kbdim
            IF ( prv(ii,jj)/prs(ii,jj) <= 1.000 ) CYCLE

            paa = 4.*mwa*surfw0/(rg*rhowa*temp(ii,jj)) ! Kelvin effect [m]

            ! Determine Dstar == critical diameter corresponding to the host model S
            zdcstar = 2.*paa/( 3.*( (prv(ii,jj)/prs(ii,jj))-1. ) )
            zvcstar = pi6*zdcstar**3

            ! For insolubles zdstar = paa/(prv(ii,jj)/prs(ii,jj)-1.)

            ! Loop over cloud droplet (and aerosol) bins
            DO cb = ica%cur, fcb%cur
               IF (cb <= fca%cur) THEN
                  ! a-bins
                  ab = ica%par + (cb-ica%cur)
               ELSE
                  ! b-bins
                  ab = icb%par + (cb-icb%cur)
               END IF
               pactd(cb)%numc = 0.d0
               pactd(cb)%volc(:) =0.d0
               IF ( paero(ii,jj,ab)%numc < nlim) CYCLE
               intrange = .FALSE.

               ! Define some parameters
               Nmid = paero(ii,jj,ab)%numc     ! Number concentration at the current bin center
               Vwmid = SUM(paero(ii,jj,ab)%volc(:))/Nmid  ! Wet volume at the current bin center
               Vmid = SUM(paero(ii,jj,ab)%volc(1:7))/Nmid ! Dry volume at the current bin center
               Vlo = Vmid*paero(ii,jj,ab)%vratiolo        ! Dry vol at low limit
               Vhi = Vmid*paero(ii,jj,ab)%vratiohi        ! Dry vol at high limit

               ! Number concentrations and volumes at adjacent bins (check for sizedistribution boundaries)
               IF (ab == in1a .OR. ab == in2b) THEN
                  Nim1 = nlim
                  Vim1 = Vlo/2.
                  Vlom1 = 0.
                  Vhim1 = Vlo
                  Vwim1 = Vwmid/3.
               ELSE
                  Nim1 = paero(ii,jj,ab-1)%numc
                  IF (Nim1 > nlim) THEN
                     Vim1 = SUM(paero(ii,jj,ab-1)%volc(1:7))/Nim1
                     Vwim1 = SUM(paero(ii,jj,ab-1)%volc(:))/Nim1
                  ELSE
                     Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                     Vwim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                  END IF
                  Vlom1 = Vim1*paero(ii,jj,ab-1)%vratiolo
                  Vhim1 = Vim1*paero(ii,jj,ab-1)%vratiohi
               END IF
               IF (ab == fn2a .OR. ab == fn2b ) THEN
                  Nip1 = nlim
                  Vip1 = Vhi + 0.5*(Vhi-Vlo)
                  Vlop1 = Vhi
                  Vhip1 = Vhi + (Vhi-Vlo)
                  Vwip1 = Vhip1
               ELSE
                  Nip1 = paero(ii,jj,ab+1)%numc
                  IF (Nip1 > nlim) THEN
                     Vip1 = SUM(paero(ii,jj,ab+1)%volc(1:7))/Nip1
                     Vwip1 = SUM(paero(ii,jj,ab+1)%volc(:))/Nip1
                  ELSE
                     Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                     Vwip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                  END IF
                  Vlop1 = Vip1*paero(ii,jj,ab+1)%vratiolo
                  Vhip1 = Vip1*paero(ii,jj,ab+1)%vratiohi
               END IF

               ! Keeping things smooth...
               Vip1 = MAX(Vhi,Vip1)
               Vim1 = MIN(Vlo,Vim1)

               ! First, make profiles of particle wet radius in
               ! order to determine the integration boundaries
               zs1 = (Vwmid - MAX(Vwim1,0.))/(Vmid - Vim1)
               zs2 = (MAX(Vwip1,0.) - Vwmid)/(Vip1 - Vmid)

               ! Get the origin values for slope equations
               V01 = Vwmid - zs1*Vmid
               V02 = Vwmid - zs2*Vmid

               ! Get the wet sizes at bins edges
               Vwlo = MAX(V01 + zs1*Vlo, 0.)
               Vwhi = MAX(V02 + zs2*Vhi, 0.)

               ! Find out dry vol integration boundaries based on *zvcstar*:
               IF ( zvcstar < Vwlo .AND. zvcstar < Vwmid .AND. zvcstar < Vwhi ) THEN
                  ! Whole bin activates
                  vint1 = Vlo
                  vint2 = Vhi

                  intrange(1:4) = .TRUE.

               ELSE IF ( zvcstar > Vwlo .AND. zvcstar > Vwmid .AND. zvcstar > Vwhi) THEN
                  ! None activates
                  vint1 = 999.
                  vint2 = 999.

                  intrange(1:4) = .FALSE.

               ELSE
                  ! Partial activation:
                   ! Slope1
                  vcut = (zvcstar - V01)/zs1  ! Where the wet size profile intersects the critical size (slope 1)
                  IF (vcut < Vlo .OR. vcut > Vmid) THEN
                     ! intersection volume outside the current size range -> set as the lower limit
                     vint1 = Vlo
                  ELSE
                     vint1 = vcut
                  END IF

                  ! Slope2
                  vcut = (zvcstar - V02)/zs2  ! Where the wet size profile intersects the critical size (slope 2)
                  IF (vcut < Vmid .OR. vcut > Vhi) THEN
                     ! Intersection volume outside the current size range -> set as the lower limit
                     vint2 = Vmid
                  ELSE
                     vint2 = vcut
                  END IF

                  ! Determine which size ranges have wet volume larger than the critical
                  intrange(1) = ( Vwlo > zvcstar )
                  intrange(2) = ( Vwmid > zvcstar )
                  intrange(3) = ( Vwmid > zvcstar )
                  intrange(4) = ( Vwhi > zvcstar )

               END IF

               ! Number concentration profiles within bins and integration for number of activated:
               ! -----------------------------------------------------------------------------------
               ! get density distribution values for number concentration
               dNim1 = Nim1/(Vhim1-Vlom1)
               dNip1 = Nip1/(Vhip1-Vlop1)
               dNmid = Nmid/(Vhi-Vlo)

               ! Get slopes
               zs1 = ( dNmid - dNim1 )/( Vmid - Vim1 )
               zs2 = ( dNip1 - dNmid )/( Vip1 - Vmid )

               N01 = dNmid - zs1*Vmid  ! Origins
               N02 = dNmid - zs2*Vmid  !

               ! Define normalization factors
               Nnorm = intgN(zs1,N01,Vlo,Vmid) + intgN(zs2,N02,Vmid,Vhi)
               Vnorm = intgV(zs1,N01,Vlo,Vmid) + intgV(zs2,N02,Vmid,Vhi)

               ! Accumulated variables
               zactvol = 0.
               Nact    = 0.
               Vact(:) = 0.

               ! Integration over each size range within a bin
               IF ( intrange(1) ) THEN
                  Nact = Nact + (Nmid/Nnorm)*intgN(zs1,N01,Vlo,vint1)
                  zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs1,N01,Vlo,vint1)
               END IF

               IF ( intrange(2) ) THEN
                  Nact = Nact + (Nmid/Nnorm)*intgN(zs1,N01,vint1,Vmid)
                  zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs1,N01,vint1,Vmid)
               END IF

               IF ( intrange(3) ) THEN
                  Nact = Nact + (Nmid/Nnorm)*intgN(zs2,N02,Vmid,vint2)
                  zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs2,N02,Vmid,vint2)
               END IF

               IF ( intrange(4) ) THEN
                  Nact = Nact + (Nmid/Nnorm)*intgN(zs2,N02,vint2,Vhi)
                  zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs2,N02,vint2,Vhi)
               END IF

               DO ss = 1, 8
                  Vact(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
               END DO

               ! Store the number concentration and mass of activated particles for current bins
               pactd(cb)%numc = MIN(Nact,Nmid)
               pactd(cb)%volc(:) = MIN(Vact(:),paero(ii,jj,ab)%volc(:))

            END DO ! cb

            ! Apply the number and mass activated to aerosol and cloud bins
            paero(ii,jj,ica%par:fcb%par)%numc =   &
               MAX(0., paero(ii,jj,ica%par:fcb%par)%numc - pactd(ica%cur:fcb%cur)%numc)
            pcloud(ii,jj,ica%cur:fcb%cur)%numc = pcloud(ii,jj,ica%cur:fcb%cur)%numc + pactd(ica%cur:fcb%cur)%numc
            DO ss = 1, 8
               paero(ii,jj,ica%par:fcb%par)%volc(ss) =  &
                  MAX(0., paero(ii,jj,ica%par:fcb%par)%volc(ss) - pactd(ica%cur:fcb%cur)%volc(ss))
               pcloud(ii,jj,ica%cur:fcb%cur)%volc(ss) = pcloud(ii,jj,ica%cur:fcb%cur)%volc(ss) + pactd(ica%cur:fcb%cur)%volc(ss)
            END DO

         END DO ! ii

      END DO ! jj

   END SUBROUTINE actInterst


   ! ----------------------------------------------

   SUBROUTINE activate3(kbdim,klev,paero,pbcrita, &
                        pdcrit, pdcrlo, pdcrhi, pdcstar, pactd   )
      !
      ! Gets the number and mass activated in the critical aerosol size bin
      !
      USE mo_submctl, ONLY : t_section, pi6, nlim, fn2b, ncld,  &
                             in1a, fn2a, ica, fca, icb, fcb, in2b, fn2b
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: kbdim,klev
      TYPE(t_section), INTENT(IN) :: paero(kbdim,klev,fn2b)
      INTEGER, INTENT(IN) :: pbcrita(kbdim,klev)          ! Index of the critical aerosol bin in regime a

      REAL, INTENT(IN) :: pdcrit(kbdim,klev,fn2b),    & ! Bin middle critical diameter
                          pdcrlo(kbdim,klev,fn2b),    & ! Critical diameter at low limit
                          pdcrhi(kbdim,klev,fn2b)       ! Critical diameter at high limit
      REAL, INTENT(IN) :: pdcstar(kbdim,klev)           ! Critical diameter corresponding to Smax
      TYPE(t_section), INTENT(OUT) :: pactd(kbdim,klev,ncld) ! Properties of the maximum amount of newly activated droplets


      REAL :: zvcstar, zvcint
      REAL :: zs1,zs2             ! Slopes
      REAL :: Nmid, Nim1,Nip1
      REAL :: dNmid, dNim1, dNip1
      REAL :: Vmid,Vlo,Vhi
      REAL :: Vim1,Vlom1,Vhim1
      REAL :: Vip1,Vlop1,Vhip1
      REAL :: Nnorm,Vnorm,N01,N02
      REAL :: zactvol

      REAL :: zhlp

      INTEGER :: ii,jj,ss,cb,ab

      DO jj = 1, klev
         DO ii = 1, kbdim

            ! This means in practice that vertical velocity is <= 0 or Ntot == 0
            IF ( ALL(pdcrit(ii,jj,:) < epsilon(1.0)) ) CYCLE

            zvcstar = 0.

            IF ( paero(ii,jj,pbcrita(ii,jj))%numc < nlim ) THEN
               Vmid = pi6*paero(ii,jj,pbcrita(ii,jj))%dmid**3
            ELSE
               Vmid = SUM( paero(ii,jj,pbcrita(ii,jj))%volc(1:7) )/MAX(nlim,paero(ii,jj,pbcrita(ii,jj))%numc)
            END IF
            Vhi = paero(ii,jj,pbcrita(ii,jj))%vhilim
            Vlo = paero(ii,jj,pbcrita(ii,jj))%vlolim


            IF ( pdcstar(ii,jj) >= pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN
               zhlp = ( pi6*pdcstar(ii,jj)**3 - pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 ) / &
                      MAX( epsilon(1.0), pi6*pdcrhi(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 )

               zvcstar = Vmid + zhlp*(Vhi-Vmid)
            ELSE IF (pdcstar(ii,jj) < pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN

               zhlp = ( pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcstar(ii,jj)**3 ) / &
                      MAX( epsilon(1.0), pi6*pdcrit(ii,jj,pbcrita(ii,jj))**3 - pi6*pdcrlo(ii,jj,pbcrita(ii,jj))**3 )

               zvcstar = Vmid - zhlp*(Vmid-Vlo)
            END IF

            zvcstar = MAX( zvcstar, paero(ii,jj,pbcrita(ii,jj))%vlolim )
            zvcstar = MIN( zvcstar, paero(ii,jj,pbcrita(ii,jj))%vhilim )

            ! Loop over cloud droplet (and aerosol) bins
            DO cb = ica%cur, fcb%cur
               IF (cb <= fca%cur) THEN
                  ab = ica%par + (cb-ica%cur)
               ELSE
                  ab = icb%par + (cb-icb%cur)
               END IF

               IF ( paero(ii,jj,ab)%numc < nlim) CYCLE

               ! Formulate a slope for Wet particle size within bins and integrate over
               ! the particles larger than zvcstar

               Nmid = MAX(paero(ii,jj,ab)%numc, nlim)
               Vmid = SUM(paero(ii,jj,ab)%volc(1:7))/Nmid ! Dry bin mid volume
               Vlo = paero(ii,jj,ab)%vlolim      ! Mid dry volume scaled to bin low limit (this is mostly an educated guess... )
               Vhi = paero(ii,jj,ab)%vhilim      ! Same for high limit

               IF (ab == in1a .OR. ab == in2b) THEN
                  Nim1 = nlim
                  Vim1 = Vlo/2.
                  Vlom1 = 0.
                  Vhim1 = Vlo
               ELSE
                  Nim1 = MAX(paero(ii,jj,ab-1)%numc, nlim)
                  IF (Nim1 > nlim) THEN
                     Vim1 = SUM(paero(ii,jj,ab-1)%volc(1:7))/Nim1
                  ELSE
                     Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                  END IF
                  Vlom1 = paero(ii,jj,ab-1)%vlolim
                  Vhim1 = paero(ii,jj,ab-1)%vhilim
               END IF
               IF (ab == fn2a .OR. ab == fn2b) THEN
                  Nip1 = nlim
                  Vip1 = Vhi + 0.5*(Vhi-Vlo)
                  Vlop1 = Vhi
                  Vhip1 = Vhi + (Vhi-Vlo)
               ELSE
                  Nip1 = MAX(paero(ii,jj,ab+1)%numc, nlim)
                  IF (Nip1 > nlim) THEN
                     Vip1 = SUM(paero(ii,jj,ab+1)%volc(1:7))/Nip1
                  ELSE
                     Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                  END IF
                  Vlop1 = paero(ii,jj,ab+1)%vlolim
                  Vhip1 = paero(ii,jj,ab+1)%vhilim
               END IF

               Vip1 = MAX(Vlop1,MIN(Vip1,Vhip1))
               Vim1 = MAX(Vlom1,MIN(Vim1,Vhim1))
             

               ! get density distribution values for
               dNim1 = Nim1/(Vhim1-Vlom1)
               dNip1 = Nip1/(Vhip1-Vlop1)
               dNmid = Nmid/(Vhi-Vlo)

               ! Get slopes
               zs1 = ( dNmid - dNim1 )/( Vmid - Vim1 )
               zs2 = ( dNip1 - dNmid )/( Vip1 - Vmid )

               N01 = dNmid - zs1*Vmid  ! Origins
               N02 = dNip1 - zs2*Vip1  !

               ! Define normalization factors
               Nnorm = intgN(zs1,N01,Vlo,Vmid) + intgN(zs2,N02,Vmid,Vhi)
               Vnorm = intgV(zs1,N01,Vlo,Vmid) + intgV(zs2,N02,Vmid,Vhi)

               IF (zvcstar < Vmid) THEN

                  ! Use actual critical volume only in the critical bin, otherwise current bin limits
                  zvcint = MAX(zvcstar, Vlo)

                  pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs1,N01,zvcint,Vmid) + intgN(zs2,N02,Vmid,Vhi) )
                  ! For different species, assume the mass distribution identical in particles within the bin
                  zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs1,N01,zvcint,Vmid) + intgV(zs2,N02,Vmid,Vhi) )
                  DO ss = 1, 8
                     pactd(ii,jj,cb)%volc(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                  END DO

               ELSE IF (zvcstar >= Vmid) THEN

                  ! Use actual critical volume only in the critical bin, otherwise current bin limits
                  zvcint = MIN(zvcstar,Vhi)

                  pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs2,N02,zvcint,Vhi) )
                  zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs2,N02,zvcint,Vhi) )
                  DO ss = 1, 8
                     pactd(ii,jj,cb)%volc(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                  END DO

               END IF

               pactd(ii,jj,cb)%numc = MAX(0., pactd(ii,jj,cb)%numc)
               DO ss = 1, 8
                  pactd(ii,jj,cb)%volc(ss) = MAX(0., pactd(ii,jj,cb)%volc(ss))
               END DO

               ! "Artificially" adjust the wet size of newly activated a little bit to prevent them from being
               ! evaporated immediately
               pactd(ii,jj,cb)%volc(8) = pactd(ii,jj,cb)%numc*pi6*(pdcrit(ii,jj,ab)**3) *  &
                                         MIN(2.,(3.e-6/max(epsilon(1.0),pdcrit(ii,jj,ab)))**2)

            END DO ! cb

         END DO ! ii
      END DO ! jj

   END SUBROUTINE activate3

   ! ------------------------------------------------
   REAL FUNCTION intgN(ikk,icc,ilow,ihigh)
      ! Gets the integral over a (linear) number concentration distribution
    
      IMPLICIT NONE
      REAL, INTENT(in) :: ikk,icc,ilow,ihigh
      intgN = 0.5*ikk*MAX(ihigh**2 - ilow**2,0.) + icc*MAX(ihigh - ilow,0.)

   END FUNCTION intgN

   ! ------------------------------------------------
   REAL FUNCTION intgV(ikk,icc,ilow,ihigh)
      ! Gets the integral over a volume volume distribution based on a linear
      ! number concentration distribution
    
      IMPLICIT NONE
      REAL, INTENT(in) :: ikk,icc,ilow,ihigh
      intgV = (1./3.)*ikk*MAX(ihigh**3 - ilow**3,0.) + 0.5*icc*MAX(ihigh**2 - ilow**2,0.)

   END FUNCTION intgV


   !-----------------------------------------
   SUBROUTINE autoconv2(kbdim,klev,   &
                        pcloud,pprecp, ptstep  )
      !
      ! Uses a more straightforward method for converting cloud droplets to drizzle.
      ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
      ! parameter and is set to 1.2 by default
      !
    
      USE mo_submctl, ONLY : t_section,   &
                             ncld,        &
                             nprc,        &
                             pi6,         &
                             nlim, prlim
      IMPLICIT NONE

      INTEGER, INTENT(in) :: kbdim,klev
      REAL, INTENT(in)    :: ptstep
      TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld)
      TYPE(t_section), INTENT(inout) :: pprecp(kbdim,klev,nprc)

      REAL :: Vrem, Nrem, Vtot, Ntot
      REAL :: dvg,dg

      REAL, PARAMETER :: zd0 = 50.e-6
      REAL, PARAMETER :: sigmag = 1.2

      INTEGER :: ii,jj,cc,ss

      DO jj = 1, klev
         DO ii = 1, kbdim
            DO cc = 1, ncld ! Start from the largest drops
               ! Autoconversion rate can be limited

               Ntot = pcloud(ii,jj,cc)%numc
               Vtot = SUM(pcloud(ii,jj,cc)%volc(:))

               IF ( Ntot > nlim .AND. Vtot > 0. ) THEN

                  ! Volume geometric mean diameter
                  dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(sigmag)**2)/2. )
                  dg = dvg*EXP( -3.*LOG(sigmag)**2 )

                  Vrem = Vtot*( 1. - cumlognorm(dvg,sigmag,zd0) )
                  Nrem = Ntot*( 1. - cumlognorm(dg,sigmag,zd0) )

                  IF ( Vrem > 0. .AND. Nrem > prlim) THEN

                     ! Put the mass and number to the first precipitation bin and remove from
                     ! cloud droplets
                     DO ss = 1, 7
                        pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + pcloud(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                        pcloud(ii,jj,cc)%volc(ss) = pcloud(ii,jj,cc)%volc(ss)*(1. - (Nrem/Ntot))
                     END DO
                   
                     pprecp(ii,jj,1)%volc(8) = pprecp(ii,jj,1)%volc(8) + pcloud(ii,jj,cc)%volc(8)*(Vrem/Vtot)
                     pcloud(ii,jj,cc)%volc(8) = pcloud(ii,jj,cc)%volc(8)*(1. - (Vrem/Vtot))

                     pprecp(ii,jj,1)%numc = pprecp(ii,jj,1)%numc + Nrem
                     pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc - Nrem

                  END IF ! Nrem Vrem

               END IF ! Ntot Vtot

            END DO ! cc
         END DO ! ii
      END DO ! jj

   END SUBROUTINE autoconv2


   !***********************************************
   !***********************************************
   ! Ice nucleation
   ! Sources
   !   Mor05   Morrison et al. (Journal of Aerosol Science, 62, 1665-1677, 2005)
   !   KC00    Khvorostyanov and Curry (Geophysical Research letters, 27, 4081-4084, 2000)
   !   KS98  Khvorostyanov and Sassen (Geophysical Research letters, 25, 3155-3158, 1998)
   !   PK97  Pruppacher and Klett, Microphysics of Clouds and Precipitation, 1997

   SUBROUTINE ice_nucl_driver(kbdim,klev,   &
                              paero,pcloud,pprecp,pice,psnow, &
                              ptemp,ppres,prv,prsi,ptstep )

    
      USE mo_submctl, ONLY : t_section,   &
                             fn2b,        &
                             ncld,        &
                             nice,        &
                             rhowa,       &
                             rhoic,       &
                             planck,      &
                             pi,          &
                             nlim,        &
                             prlim,       &
                             rd, avog,    &
                             calc_Sw_eq, nprc, nsnw

      IMPLICIT NONE

      INTEGER, INTENT(in) :: kbdim,klev
      REAL, INTENT(in)    :: ptstep
      REAL, INTENT(in)    :: ptemp(kbdim,klev),  &
                             prv(kbdim,klev),    &
                             prsi(kbdim,klev),   &
                             ppres(kbdim,klev)

      TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,fn2b), &
                                        pcloud(kbdim,klev,ncld), &
                                        pprecp(kbdim,klev,nprc), &
                                        pice(kbdim,klev,nice),  &
                                        psnow(kbdim,klev,nsnw)

      INTEGER :: ii, jj, kk, ss, ba
      REAL :: pf, jf, Sw, Si, rn, rw, frac, pw, zvol, zvola, zvolnew
      LOGICAL :: isdry
      ! Flags (adjust with the nucleation parameters):
      !   Which species can nucleate ice
      LOGICAL :: aerosol_ice = .FALSE., cloud_ice = .TRUE., precip_ice = .TRUE.
      !   Which nucleation modes are used (contact freezing not implemented)
      LOGICAL :: ice_hom = .TRUE., ice_imm=.TRUE., ice_dep=.FALSE.
      !   Other options
      !REAL :: tmax_homog=243.15 ! Maximum temperature for homogenous nucleation

      DO ii = 1, kbdim
         DO jj = 1, klev
            IF (ptemp(ii,jj) > 273.15 ) CYCLE

            ! Precipitation
            DO kk = 1, nprc

               IF (pprecp(ii,jj,kk)%numc<prlim .OR. .NOT.precip_ice) CYCLE

               ! Radius of the insoluble portion of the droplet
               rn = MAX(0., (3.*sum(pprecp(ii,jj,kk)%volc(3:4))/pprecp(ii,jj,kk)%numc/4./pi)**(1./3.) )
               ! Droplet radius
               rw = (3.*sum(pprecp(ii,jj,kk)%volc(:))/pprecp(ii,jj,kk)%numc/4./pi)**(1./3.)
               ! Equilibrium saturation ratio
               Sw = calc_Sw_eq(pprecp(ii,jj,kk),ptemp(ii,jj))

               IF (rn > 1.e-10 .AND. ice_imm) THEN
                  ! Immersion and condensation freezing
                  jf = calc_Jcf(rn,ptemp(ii,jj),Sw)
                  pf = 1. - exp( -jf*ptstep )
               ELSE IF (rn < 1.e-10 .AND. ice_hom) THEN
                  ! Homogeneous freezing
                  jf = calc_Jhf(ptemp(ii,jj),Sw)
                  pf = 1. - exp( -jf*4./3.*pi*(rw**3-rn**3)*ptstep )
               ELSE
                  pf = 0.
               END IF

               frac = MIN(1.,pf)
               IF (pprecp(ii,jj,kk)%numc*frac <prlim) CYCLE

               ! Move to ice bin with closest matching dry volume. Ain't perfect but
               ! the bin update subroutine in SALSA will take care of the rest.
               zvol = SUM( pprecp(ii,jj,kk)%volc(1:7) ) ! Dry volume

               ba = 0
               zvola = -1.
               DO ss = 1, nice
                  IF (pice(kk,ii,ss)%numc>prlim) THEN
                     zvolnew = SUM( pice(ii,jj,ss)%volc(1:7) ) ! Dry volume
                     IF (abs(zvolnew-zvol)<abs(zvola-zvol)) THEN
                        ! New closest match
                        ba=ss
                        zvola=zvolnew
                     END IF
                  END IF
               END DO

               IF (ba==0) STOP 'FAIL: no target ice for freezing rain drops'

               DO ss = 1, 7
                  pice(ii,jj,ba)%volc(ss) = max(0.,pice(ii,jj,ba)%volc(ss) + pprecp(ii,jj,kk)%volc(ss)*frac)
                  pprecp(ii,jj,kk)%volc(ss) = max(0.,pprecp(ii,jj,kk)%volc(ss) - pprecp(ii,jj,kk)%volc(ss)*frac)
               END DO

               ss = 8
               pice(ii,jj,ba)%volc(ss) = max(0.,pice(ii,jj,ba)%volc(ss) + pprecp(ii,jj,kk)%volc(ss)*frac*rhowa/rhoic)
               pprecp(ii,jj,kk)%volc(ss) = max(0.,pprecp(ii,jj,kk)%volc(ss) - pprecp(ii,jj,kk)%volc(ss)*frac)

               pice(ii,jj,ba)%numc = max(0.,pice(ii,jj,ba)%numc + pprecp(ii,jj,kk)%numc*frac)
               pprecp(ii,jj,kk)%numc = max(0.,pprecp(ii,jj,kk)%numc-pprecp(ii,jj,kk)%numc*frac)
            END DO !kk


            ! Cloud droplets
            DO kk = 1, ncld
               IF (pcloud(ii,jj,kk)%numc < nlim .OR. .NOT. cloud_ice) CYCLE

               ! Radius of the insoluble portion of the droplet
               rn = MAX(0., (3.*sum(pcloud(ii,jj,kk)%volc(3:4))/pcloud(ii,jj,kk)%numc/4./pi)**(1./3.) )
               ! Droplet radius
               rw = (3.*sum(pcloud(ii,jj,kk)%volc(:))/pcloud(ii,jj,kk)%numc/4./pi)**(1./3.)
               ! Equilibrium saturation ratio
               Sw = calc_Sw_eq(pcloud(ii,jj,kk),ptemp(ii,jj))
               ! Is it dry?
               isdry = (pcloud(ii,jj,kk)%volc(8)<1e-20)

               IF (isdry .AND. ice_dep) THEN
                  ! Deposition freezing
                  Si = prv(ii,jj)/prsi(ii,jj)
                  pw = prv(ii,jj)*ppres(ii,jj)/(0.622-prv(ii,jj))
                  jf = calc_Jdep(rn,ptemp(ii,jj),Si,pw)
                  pf = 1. - exp( -jf*ptstep )
               ELSE IF (rn > 1.e-10 .AND. ice_imm) THEN
                  ! Immersion and condensation freezing
                  jf = calc_Jcf(rn,ptemp(ii,jj),Sw)
                  pf = 1. - exp( -jf*ptstep )
               ELSE IF (rn<1.e-10 .AND. ice_hom) THEN
                  ! Homogeneous freezing
                  jf = calc_Jhf(ptemp(ii,jj),Sw)
                  pf = 1. - exp( -jf*4./3.*pi*(rw**3-rn**3)*ptstep )
               ELSE
                  pf = 0.
               END IF

               frac = MIN(1.,pf)
               IF (pcloud(ii,jj,kk)%numc*frac <prlim) CYCLE
  
               DO ss = 1, 7
                  pice(ii,jj,kk)%volc(ss) = max(0.,pice(ii,jj,kk)%volc(ss) + pcloud(ii,jj,kk)%volc(ss)*frac)
                  pcloud(ii,jj,kk)%volc(ss) = max(0.,pcloud(ii,jj,kk)%volc(ss) - pcloud(ii,jj,kk)%volc(ss)*frac)
               END DO
               ss = 8
               pice(ii,jj,kk)%volc(ss) = max(0.,pice(ii,jj,kk)%volc(ss) + pcloud(ii,jj,kk)%volc(ss)*frac*rhowa/rhoic)
               pcloud(ii,jj,kk)%volc(ss) = max(0.,pcloud(ii,jj,kk)%volc(ss) - pcloud(ii,jj,kk)%volc(ss)*frac)

               pice(ii,jj,kk)%numc = max(0.,pice(ii,jj,kk)%numc + pcloud(ii,jj,kk)%numc*frac)
               pcloud(ii,jj,kk)%numc = max(0.,pcloud(ii,jj,kk)%numc-pcloud(ii,jj,kk)%numc*frac)
            END DO !kk

            ! Aerosol
            DO kk = 1, fn2b
               IF (paero(ii,jj,kk)%numc<nlim .OR. .NOT.aerosol_ice) CYCLE

               ! Radius of the insoluble portion of the droplet
               rn = MAX(0., (3.*sum(paero(ii,jj,kk)%volc(3:4))/paero(ii,jj,kk)%numc/4./pi)**(1./3.) )
               ! Droplet radius
               rw = (3.*sum(paero(ii,jj,kk)%volc(:))/paero(ii,jj,kk)%numc/4./pi)**(1./3.)
               ! Equilibrium saturation ratio
               Sw = calc_Sw_eq(paero(ii,jj,kk),ptemp(ii,jj))
               ! Is it dry?
               isdry = (paero(ii,jj,kk)%volc(8)<1e-20)

               IF (isdry .AND. ice_dep) THEN
                  ! Deposition freezing
                  Si = prv(ii,jj)/prsi(ii,jj)
                  pw = prv(ii,jj)*ppres(ii,jj)/(0.622-prv(ii,jj))
                  jf = calc_Jdep(rn,ptemp(ii,jj),Si,pw)
                  pf = 1. - exp( -jf*ptstep )
               ELSE IF (rn>1.e-10 .AND. ice_imm) THEN
                  ! Immersion and condensation freezing
                  jf = calc_Jcf(rn,ptemp(ii,jj),Sw)
                  pf = 1. - exp( -jf*ptstep )
               ELSE IF (rn<1.e-10 .AND. ice_hom) THEN
                  ! Homogeneous freezing
                  jf = calc_Jhf(ptemp(ii,jj),Sw)
                  pf = 1. - exp( -jf*4./3.*pi*(rw**3-rn**3)*ptstep )
               ELSE
                  pf = 0.
               END IF

               frac = MIN(1.,pf)
               IF (paero(ii,jj,kk)%numc*frac <prlim) CYCLE
  
               DO ss = 1, 7
                  pice(ii,jj,kk)%volc(ss) = max(0.,pice(ii,jj,kk)%volc(ss) + paero(ii,jj,kk)%volc(ss)*frac)
                  paero(ii,jj,kk)%volc(ss) = max(0.,paero(ii,jj,kk)%volc(ss) - paero(ii,jj,kk)%volc(ss)*frac)
               END DO
               ss = 8
               pice(ii,jj,kk)%volc(ss) = max(0.,pice(ii,jj,kk)%volc(ss) + paero(ii,jj,kk)%volc(ss)*frac*rhowa/rhoic)
               paero(ii,jj,kk)%volc(ss) = max(0.,paero(ii,jj,kk)%volc(ss) - paero(ii,jj,kk)%volc(ss)*frac)

               pice(ii,jj,kk)%numc = max(0.,pice(ii,jj,kk)%numc + paero(ii,jj,kk)%numc*frac)
               paero(ii,jj,kk)%numc = max(0.,paero(ii,jj,kk)%numc-paero(ii,jj,kk)%numc*frac)
            END DO !kk
         END DO !jj
      END DO !ii

   END SUBROUTINE ice_nucl_driver
   !***********************************************
   !
   ! homogenous nucleation according to Morrison et al. 2005 (JAS 62:1665-1677) eq. (27)
   ! as of referenced as [Mor05]
   !
   !***********************************************
   REAL FUNCTION calc_Jcf(rn,temp,Sw)
      ! The rate of germ formation (#/s) through condensation (or immersion) freezing
      !   Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
      ! - "Freezing" mode includes condensation, contact and immersion freezing modes
      USE mo_submctl, ONLY : boltz, planck, pi, rg, mwa, surfw0, avog
      IMPLICIT NONE
      REAL, INTENT(in) :: rn,temp,Sw

      REAL :: Tc, sigma_is, Lefm, act_energy, GG, r_g, x, fii, fpsi, sf, crit_energy
      REAL, PARAMETER ::   & ! More or less constant parameters
               epsi = 0.,   & ! Elastic strain produced in ice embryo by the insoluble substrate (set to zero)
               C = 1.7e10,  & ! Constant (1.7e11 dyn cm^-2 = 1.7e11*1e-5/1e-4 N m^-2 = 1.7e10 N m^-2)
               mis = 0.5,   & ! Cosine of the contact angle (wettability parameter) at the solution-ice interface (adjustable)
               alpha = 0.0, & ! Relative area of active sites (adjustable)
               rho_ice = 900., & ! Density of ice (kg/m^3)
               c_1s = 1e19, & ! The concentration of water molecules adsorbed on 1 m^2 a surface (#/m^2)
               T0 = 273.15    ! 0 C in Kelvins

      Tc = temp-T0 ! Temperature in Celsius

      IF (rn < 1e-10) THEN
         ! Homogenous freezing not considered here
         calc_Jcf = 0.
         RETURN
      END IF

      ! Activation energy at the solution-ice interface
      IF (Tc <= -30.0) THEN
         !  KS98, when T < -30 C, homogenous nucleation (erg = 1e-7 J)
         act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))
      ELSE
         ! Li et al. (JGR, 118, 11213-11227, doi:0.1002/jgrd.50861, 2013)
         act_energy = 5.55*exp(-8.423e-3*Tc+6.384e-4*Tc**2+7.891e-6*Tc**3)/avog*4.1868e3
      END IF

      ! Critical energy of germ formation (eq. 2.10 in KC00)
      ! a) Ice germ radius (eq. 2.6 in KC00)
      Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6 in KS98; cal/g=4.1868e3 J/kg)
      GG = rg*temp/mwa/Lefm ! Eq 2.7 in KC00
      sigma_is = 28e-3+0.25e-3*Tc ! Surface tension of ice at the ice-solution interface (from KS98)
      r_g = 2.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG)-C*epsi**2)
      IF (r_g<=0. .OR. (T0/temp)*Sw**GG<1.0001) THEN
         ! Must be positive
         calc_JCF = 0.
         RETURN
      END IF

      ! b) Shape factor (eq 2.9 in KC00)
      x = rn/r_g
      fii = sqrt(1.-2.*mis*x+x**2)
      fpsi = (x-mis)/fii
      sf = 0.5* (1. + ( ( 1.-mis*x)/fii)**3 + (2.-3.*fpsi+fpsi**3)*x**3 + 3.*mis*(fpsi-1.)*x**2)
      ! Critical energy (eq. 2.10)
      crit_energy = 4.*pi/3.*sigma_is*r_g**2*sf-alpha*(1-mis)*rn**2
      ! Eq 2.1 in KC00
      calc_JCF= boltz*temp/planck*c_1s*4.*pi*rn**2*exp((-act_energy-crit_energy)/(boltz*temp))

   END FUNCTION calc_JCF

   REAL FUNCTION calc_Jdep(rn,temp,Si,pw)
      ! The rate of germ formation (#/s) through deposition freezing
      !   Hoose et al., J. Atmos. Sci., 67, 2483-2503, 2010
      ! See also:
      !   Chen et al., Atmos. Chem. Phys., 8, 7431-7449, 2008
      !   Khvorostyanov and Curry, Geophys. Res. Lett., 27, 4081-4084, 2000 [KC00]
      ! - Requires solid surface and supersaturation
      !
      USE mo_submctl, ONLY : boltz, planck, pi, rg, mwa, rhowa, surfw0, avog
      IMPLICIT NONE
      REAL, INTENT(in) :: rn,temp,Si,pw
      REAL :: Tc, act_energy, sigma_iv, vw, mw, r_g , x, fii, fpsi, sf, crit_energy
      REAL, PARAMETER :: & ! More or less constant parameters
               mis = 0.5, & ! Cosine of the contact angle (wettability parameter) at the solution-ice interface (adjustable)
               T0 = 273.15, & ! 0 C in Kelvins
               vs = 1e13 ! Frequency of vibration of water vapor molecule adsorbed on solid substrate (1/s)

      Tc = temp-T0 ! Temperature in Celsius

      IF (rn<=0. .OR. Si<=1.0001) THEN
         ! A core is required and supersaturaturation over ice
         calc_Jdep = 0.
         RETURN
      END IF

      ! Activation energy for deposition nucleation
      act_energy = 0.

      ! Critical energy of germ formation (eq. 2.12 in KC00)
      ! a) Ice germ radius (eq. 2.12 in KC00)
      sigma_iv = ( (76.1-0.155*Tc) + (28.5+0.25*Tc) )*1e-3 ! Hoose et al. (2010)
      vw = mwa/rhowa ! Volume of mole of water (m^3/mol)
      mw = mwa/avog ! Mass of water molecule (kg)
      r_g = 2.*vw*sigma_iv/( rg*temp*log(Si) ) ! Eq. 8 (Hoose et al., 2010)

      IF (r_g <= 0.) THEN
         ! Must be positive
         calc_Jdep = 0.
         RETURN
      END IF

      ! b) Shape factor (eq 2.9 in KC00 or Eq 2 in Hoose et al.)
      x = rn/r_g
      fii = sqrt(1.-2.*mis*x+x**2)
      fpsi = (x-mis)/fii
      sf = 0.5* (1. + ( ( 1.-mis*x)/fii)**3 + (2.-3.*fpsi+fpsi**3)*x**3 + 3.*mis*(fpsi-1.)*x**2)
      ! Critical energy
      crit_energy = 4.*pi/3.*sigma_iv*r_g**2*sf

      ! Eqs 1 and 10
      calc_Jdep = (pw**2*vw/(mw*boltz*temp*vs)*sqrt(sigma_iv/boltz/temp))*rn**2*sqrt(sf)*exp((-act_energy-crit_energy)/(boltz*temp))

   END FUNCTION calc_Jdep


   REAL FUNCTION calc_Jhf(temp,Sw)
      ! Homogeneous freezing based on Khovosrotyanov and Sassen, Geophys. Res. Lett., 25, 3155-3158, 1998 [KS98]
      ! - Must have an aqueous phase
      ! - Important for low temperatures (<243 K)
      !
      USE mo_submctl, ONLY : boltz, planck, pi, mwa, rhowa, rg
      IMPLICIT NONE
      REAL, intent(in) :: temp, Sw ! Temperature (K) and water vapor saturation ratio
      REAL :: Tc, act_energy, Lefm, GG, sigma_is, r_g, crit_energy
      REAL, PARAMETER :: &
         Nc = 5.85e16, & ! The number of water molecules contacting unit area of ice germ (#/m^2) [from KC00]
         T0 = 273.15, & ! 0 C in Kelvins
         rho_ice = 900. ! Density of ice (kg/m^3)

      Tc = temp-T0 ! Temperature in Celsius

      ! Activation energy at the solution-ice interface (erg = 1e-7 J)
      act_energy = 0.694e-19 * (1.+ 0.027*(Tc+30.))

      ! Critical energy of germ formation (eq. 9b in KS98)
      ! a) Ice germ radius (eq. 9a in KS98)
      Lefm = (79.7+0.708*Tc-2.5e-3*Tc**2)*4.1868e3 ! Effective latent heat of melting (eq. 6 in KS98; cal/g=4.1868e3 J/kg)
      GG = rg*temp/Lefm/mwa
      sigma_is = 28e-3+0.25e-3*Tc ! Surface tension of ice at the ice-solution interface (from KS98)
      r_g = 2.*sigma_is/( rho_ice*Lefm*log((T0/temp)*Sw**GG) )

      IF (r_g<=0. .OR. (T0/temp)*Sw**GG<1.0001) THEN
         ! Critical size must be positive
         calc_Jhf = 0.
         RETURN
      END IF

      ! Critical energy
      crit_energy = 4.*pi/3.*sigma_is*r_g**2

      ! Eq. 1 in KS98
      calc_Jhf = 2.0*Nc*(rhowa/rho_ice)*boltz*temp/planck*sqrt(sigma_is/boltz/temp)*exp((-crit_energy-act_energy)/(boltz*temp))

      ! Nucleation rate [m^-3 s^-1] is for droplet liquid volume, so freezing propability is 1-exp( -J*V*tstep), where V=4/3*pi*(rwet**3-rdry**3)

   END FUNCTION calc_Jhf
   !***********************************************
   !
   ! heterogenous immersion nucleation according to Dieh & Wurzler 2006 JAS 61:2063-2072
   ! as of referenced as [Mor05]
   !
   !***********************************************
   SUBROUTINE ice_fixed_NC(  kbdim,  klev,   &
                             pcloud,  pice,   &
                             ptemp,  ppres,  prv,  prsi)


      USE mo_submctl, ONLY : t_section,   &
                             ncld,        &
                             nice,        &
                             rhowa,       &
                             rhoic,       &
                             rda,         &
                             nlim,        &
                             fixinc   ! Target ice number concentration (#/kg)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: kbdim,klev

      REAL, INTENT(in) :: &
                             ptemp(kbdim,klev),    &
                             ppres(kbdim,klev),    &
                             prv(kbdim,klev),    &
                             prsi(kbdim,klev)

      TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                        pice(kbdim,klev,nice)

      INTEGER :: ii,jj,kk,ss

      REAL :: pdn, iceSupSat, rc_tot, Ni0,  &
              sumICE, iceTendecyNumber, liqToIceFrac


      DO ii = 1, kbdim
         DO jj = 1, klev
            pdn = ppres(ii,jj)/(rda*ptemp(ii,jj)) ! Air density (kg/m^3)

            iceSupSat = prv(ii,jj) / prsi(ii,jj)  - 1.0 ! ice supersaturation
            rc_tot = sum( pcloud(ii,jj,:)%volc(8) )*rhowa/pdn ! cloud water mixing ratio (kg/kg)

            ! conditions for ice nucleation
            IF ( icesupsat < 0.05 .OR. rc_tot < 0.001e-3  ) CYCLE

            ! target number concentration of ice, converted to #/m^3
            Ni0 = fixinc * pdn

            ! current ice number concentration (#/m^3)
            sumICE = sum( pice(ii,jj,:)%numc )

            IF ( sumICE > Ni0 ) CYCLE

            DO kk = nice,1,-1 ! Assuming nice=ncld
               IF( sumICE < Ni0 .AND. pcloud(ii,jj,kk)%numc > nlim) THEN

                  iceTendecyNumber = max( 0.0, min( Ni0 - pice(ii,jj,kk)%numc , pcloud(ii,jj,kk)%numc )  )

                  pice(ii,jj,kk)%numc = pice(ii,jj,kk)%numc + iceTendecyNumber
                  sumICE = sumICE + iceTendecyNumber

                  liqToIceFrac   = MAX( 0.0, MIN( 1.0, iceTendecyNumber/pcloud(ii,jj,kk)%numc ) )
                  pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc - iceTendecyNumber

                  DO ss = 1,7
                     pice(ii,jj,kk)%volc(ss) = pice(ii,jj,kk)%volc(ss) + max(0., pcloud(ii,jj,kk)%volc(ss)*liqToIceFrac )
                     pcloud(ii,jj,kk)%volc(ss) = pcloud(ii,jj,kk)%volc(ss) - max(0., pcloud(ii,jj,kk)%volc(ss)*liqToIceFrac )
                  END DO
                  ss = 8 ! Water
                  pice(ii,jj,kk)%volc(ss) = pice(ii,jj,kk)%volc(ss) + max(0., pcloud(ii,jj,kk)%volc(ss)*liqToIceFrac*rhowa/rhoic )
                  pcloud(ii,jj,kk)%volc(ss) = pcloud(ii,jj,kk)%volc(ss) - max(0., pcloud(ii,jj,kk)%volc(ss)*liqToIceFrac )
               END IF
            END DO !kk
         END DO !jj
      END DO !ii

   END SUBROUTINE ice_fixed_NC

   ! ------------------------------------------------------------

   SUBROUTINE ice_melt(kbdim,klev,   &
                       pcloud,pice,pprecp,psnow, &
                       ptemp )

      USE mo_submctl, ONLY : t_section,   &
                             ncld,        &
                             nice,        &
                             nsnw,        &
                             nprc,        &
                             rhowa, rhoic, rhosn,      &
                             pi6,         &
                             prlim

      IMPLICIT NONE

      INTEGER, INTENT(in) :: kbdim,klev
      REAL, INTENT(in)    :: ptemp(kbdim,klev)

      TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), &
                                        pice(kbdim,klev,nice),   &
                                        psnow(kbdim,klev,nsnw),  &
                                        pprecp(kbdim,klev,nprc)

      INTEGER :: ii,jj,kk,ss

      DO ii = 1, kbdim
         DO jj = 1, klev
            ! Ice and snow melt when temperature above 273.15 K
            ! => should add the effect of freezing point depression
            IF (ptemp(ii,jj) <= 273.15 ) CYCLE

            DO kk = 1, nice
               ! Ice => cloud water (parallel bin)
               IF (pice(ii,jj,kk)%numc < prlim) CYCLE
               DO ss = 1, 7
                  pcloud(ii,jj,kk)%volc(ss) = pcloud(ii,jj,kk)%volc(ss) + pice(ii,jj,kk)%volc(ss)
                  pice(ii,jj,kk)%volc(ss) = 0.
               END DO
               ss = 8 ! Water
               pcloud(ii,jj,kk)%volc(ss) = pcloud(ii,jj,kk)%volc(ss) + pice(ii,jj,kk)%volc(ss)*rhoic/rhowa
               pice(ii,jj,kk)%volc(ss) = 0.

               pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc + pice(ii,jj,kk)%numc
               pice(ii,jj,kk)%numc = 0.
            END DO

            DO kk = 1, nsnw
               ! Snow => precipitation (bin 1)
               IF (psnow(ii,jj,kk)%numc < prlim) CYCLE
               DO ss = 1, 7
                  pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + psnow(ii,jj,kk)%volc(ss)
                  psnow(ii,jj,kk)%volc(ss) = 0.
               END DO
               ss = 8 ! Water
               pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + psnow(ii,jj,kk)%volc(ss)*rhosn/rhowa
               psnow(ii,jj,kk)%volc(ss) = 0.

               pprecp(ii,jj,1)%numc = pprecp(ii,jj,1)%numc + psnow(ii,jj,kk)%numc
               psnow(ii,jj,kk)%numc = 0.
            END DO
         END DO
      END DO

   END SUBROUTINE ice_melt


   SUBROUTINE autosnow(kbdim,klev,   &
                       pice,psnow ,ptstep        )
      !
      ! Uses a more straightforward method for converting ice to snow.
      ! Assume a lognormal ice distribution for each bin. Sigma_g is an adjustable
      ! parameter and is set to 1.2 by default
      !
    
      USE mo_submctl, ONLY : t_section,   &
                             nice,        &
                             nsnw,        &
                             pi6,         &
                             rhosn, rhoic,       &
                             prlim,          &
                             rd
      IMPLICIT NONE

      INTEGER, INTENT(in) :: kbdim,klev
      REAL, INTENT(in) :: ptstep
      TYPE(t_section), INTENT(inout) :: pice(kbdim,klev,nice)
      TYPE(t_section), INTENT(inout) :: psnow(kbdim,klev,nsnw)

      REAL :: Vrem, Nrem, Vtot, Ntot
      REAL :: dvg,dg

      REAL, PARAMETER :: zd0 = 250.e-6
      REAL, PARAMETER :: sigmag = 1.2

      INTEGER :: ii,jj,cc,ss

      DO jj = 1, klev
         DO ii = 1, kbdim
            DO cc = 1, nice

               Ntot = pice(ii,jj,cc)%numc
               Vtot = SUM(pice(ii,jj,cc)%volc(:))

               IF ( Ntot > prlim .AND. Vtot > 0. ) THEN
                  ! Volume geometric mean diameter
                  dvg = ((Vtot/Ntot/pi6)**(1./3.))*EXP( (3.*LOG(sigmag)**2)/2. )
                  dg  = dvg*EXP( -3.*LOG(sigmag)**2 )

                  Vrem = Vtot*( 1. - cumlognorm(dvg,sigmag,zd0) )
                  Nrem = Ntot*( 1. - cumlognorm(dg,sigmag,zd0) )

                  IF ( Vrem > 0. .AND. Nrem > prlim) THEN
                     ! Put the mass and number to the first snow bin and remover from cloud droplets

                     DO ss = 1, 7
                        psnow(ii,jj,1)%volc(ss) = psnow(ii,jj,1)%volc(ss) + pice(ii,jj,1)%volc(ss)*(Nrem/Ntot)
                        pice(ii,jj,cc)%volc(ss) = pice(ii,jj,cc)%volc(ss)*(1. - Nrem/Ntot)
                     END DO
                     ! From ice to snow volume
                     psnow(ii,jj,1)%volc(8) = psnow(ii,jj,1)%volc(8) + pice(ii,jj,cc)%volc(8)*(Vrem/Vtot)*rhoic/rhosn
                     pice(ii,jj,cc)%volc(8) = pice(ii,jj,cc)%volc(8)*(1. - (Vrem/Vtot))

                     psnow(ii,jj,1)%numc = psnow(ii,jj,1)%numc + Nrem
                     pice(ii,jj,cc)%numc = pice(ii,jj,cc)%numc - Nrem

                  END IF ! Nrem Vrem

               END IF ! Ntot Vtot

            END DO ! cc
         END DO ! ii
      END DO ! jj

   END SUBROUTINE autosnow

   !
   ! -----------------------------------------------------------------
   !
   REAL FUNCTION cumlognorm(dg,sigmag,dpart)
    
      IMPLICIT NONE
      ! Cumulative lognormal function
      REAL, INTENT(in) :: dg
      REAL, INTENT(in) :: sigmag
      REAL, INTENT(in) :: dpart

      REAL :: hlp1,hlp2

      hlp1 = ( LOG(dpart) - LOG(dg) )
      hlp2 = SQRT(2.)*LOG(sigmag)
      cumlognorm = 0.5 + 0.5*ERF( hlp1/hlp2 )

   END FUNCTION cumlognorm



END MODULE mo_salsa_cloud
