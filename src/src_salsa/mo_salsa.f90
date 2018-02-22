MODULE mo_salsa
   IMPLICIT NONE

   ! --------------------------------------------------------------------------
   ! The SALSA subroutine
   !
   ! Modified for the new aerosol datatype,
   ! Juha Tonttila, FMI, 2014.
   ! --------------------------------------------------------------------------


   PRIVATE

   ! -- subroutines
   PUBLIC :: salsa

CONTAINS

   SUBROUTINE salsa(kbdim,    klev,           &
                    ppres,    prv, prs, prsi,    ptemp, ptstep,     &
                    pc_h2so4, pc_ocnv,  pc_ocsv, pc_hno3,    &
                    pc_nh3,   paero,    pcloud,  pprecp,     &
                    pice,     psnow,                         &
                    pactd,    pw,       prtcl,   level )

      USE mo_salsa_dynamics, ONLY : coagulation, condensation
      USE mo_salsa_update,   ONLY : distr_update
      USE mo_salsa_cloud,    ONLY : cloud_activation, autoconv2, &
                                    autosnow, ice_fixed_NC, ice_nucl_driver, ice_melt

      USE mo_submctl, ONLY :      &
         fn2b,                      & ! size section and composition indices
         t_section,                 & ! For cloud bins
         ncld,                      &
         nprc,                      &
         nice,                      &
         nsnw,                      &
         lscoag,                    &
         lscnd,                     &
         lsauto,                    &
         lsautosnow,                &
         lsactiv,                   &
         lsicenucl,                 &
         lsfixinc,                  &
         lsicmelt,                  &
         lsdistupdate

      USE class_componentIndex, ONLY : ComponentIndex

      IMPLICIT NONE

      !-- Input parameters and variables --------------------------------------------------
      INTEGER, INTENT(in) ::      &
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev                         ! number of vertical levels (='klev')


      REAL, INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         ptstep                          ! time step [s]

      TYPE(ComponentIndex), INTENT(in) :: prtcl

      REAL, INTENT(in) ::            & ! Vertical velocity
         pw(kbdim,klev)

      !-- Input variables that are changed within --------------------------------------
      REAL, INTENT(inout) ::      & ! gas phase concentrations at each grid point [#/m3]
         pc_h2so4(kbdim,klev),      & ! sulphuric acid
         pc_hno3 (kbdim,klev),      & ! nitric acid
         pc_nh3  (kbdim,klev),      & ! ammonia
         pc_ocnv (kbdim,klev),      & ! nonvolatile organic compound
         pc_ocsv (kbdim,klev),      & ! semivolatile organic compound
         prv(kbdim,klev),           & ! Water vapour mixing ratio  [kg/kg]
         prs(kbdim,klev),           & ! Saturation mixing ratio    [kg/kg]
         prsi(kbdim,klev)              ! Saturation mixing ratio over ice   [kg/kg]

      TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     &
         paero(kbdim,klev,fn2b),      &
         pprecp(kbdim,klev,nprc),     &
         pice(kbdim,klev,nice),       &
         psnow(kbdim,klev,nsnw)

      TYPE(t_section), INTENT(out) :: &
         pactd(kbdim,klev,ncld)

      INTEGER, INTENT(in) :: level                         ! thermodynamical level

      !-- Local variables ------------------------------------------------------------------

      INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level

      zpbl(:) = 1

      ! Coagulation
      IF (lscoag) &
         CALL coagulation( kbdim,  klev,                   &
                           paero,  pcloud, pprecp, pice, psnow,    &
                           ptstep, ptemp,  ppres   )

      ! Condensation
      IF (lscnd) &
         CALL condensation(kbdim,    klev,               &
                           paero,   pcloud,   pprecp,             &
                           pice,    psnow,                        &
                           pc_h2so4, pc_ocnv, pc_ocsv,  pc_hno3,  &
                           pc_nh3, prv, prs, prsi, ptemp, ppres,  &
                           ptstep, zpbl, prtcl                    )

      ! Autoconversion (liquid)
      IF (lsauto) &
         CALL autoconv2(kbdim,klev, &
                        pcloud, pprecp, ptstep     )

      ! Cloud activation
      IF (lsactiv )  &
         CALL cloud_activation(kbdim, klev,   &
                               ptemp,  ppres, prv,    &
                               prs,    pw,    paero,  &
                               pcloud, pactd          )
    ! Fixed ice number concentration

    IF (lsfixinc) &
          CALL  ice_fixed_NC(kbdim, klev,   &
                             pcloud,   pice,   &
                             ptemp,  ppres,  prv,  prsi)

    ! Ice nucleation

      IF (lsicenucl) &
        CALL ice_nucl_driver(kbdim,klev,   &
                             paero,pcloud,pprecp,pice,psnow, &
                             ptemp,ppres,prv,prsi,ptstep)

      ! Melting of ice and snow

      IF (lsicmelt) &
         CALL ice_melt(kbdim,klev,              &
                       pcloud,pice,pprecp,psnow,ptemp)

      ! ! Snow formation ~ autoconversion from ice
      IF (lsautosnow) &
         CALL autosnow(kbdim,klev, &
                        pice, psnow, ptstep )

      ! Size distribution bin update
      IF (lsdistupdate ) &
         CALL distr_update(kbdim, klev,     &
                           paero,  pcloud, pprecp,  &
                           pice, psnow, level       )

   END SUBROUTINE salsa


END MODULE mo_salsa
