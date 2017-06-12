
MODULE mo_submctl


   IMPLICIT NONE

   PRIVATE

   PUBLIC :: lsdistupdate,nsnucl

   PUBLIC :: debug

   !SALSA:
   PUBLIC :: act_coeff,nj3

   PUBLIC :: in1a,in2a,in2b,fn1a,fn2a,fn2b,nbins
   PUBLIC :: nbin, nbin2, nbin3,reglim,nlim,prlim,nreg
   PUBLIC :: vhilim,vlolim,vratiohi,vratiolo,dpmid
   PUBLIC :: pi, pi6, rg, avog, planck, boltz, cpa, mair, grav, eps
   PUBLIC :: rhosu,rhooc, rhobc,rhoss, rhodu, rhowa, rhonh, rhono, rhoic,rhosn
   PUBLIC :: msu,mdu,mno,mnh,n3,massacc,d_sa,pstand,mss,mbc,moc,epsoc,mwa,ions,&
             mvsu,mvoc,mvss,surfw0,surfi0,mvwa,mvno,mvnh

   PUBLIC :: t_section,t_parallelbin
   PUBLIC :: ncldbin,ica,fca,icb,fcb,ira,fra,ncld,nprc,dmincld

   PUBLIC :: nicebin,iia,fia,iib,fib,isa,fsa,nice,nsnw,dminice

   PUBLIC :: aerobins, cloudbins, precpbins, icebins, snowbins

   PUBLIC :: liqFracA, iceFracA, liqFracB, iceFracB

   PUBLIC :: rd,rv,alv,als,alf


   PUBLIC :: nlcoag,                &
             nlcgaa,nlcgcc,nlcgpp,  &
             nlcgca,nlcgpa,nlcgpc,  &
             nlcgia,nlcgic,nlcgii,  &
             nlcgip,nlcgsa,nlcgsc,  &
             nlcgsi,nlcgsp,nlcgss,  &
             nlcnd,                 &
             nlcndgas,              &
             nlcndh2oae,nlcndh2ocl, &
             nlcndh2oic,            & ! ice'n'snow
             nlauto,                &
             nlautosnow,            &
             nlactiv,nlactbase,     &
             nlactintst,            &
             nlichom,               &
             nlichet,               &
             nlicimmers,            &
             nlicmelt,              &
             nldebug,               &
             lscoag,                &
             lscgaa,lscgcc,lscgpp,  &
             lscgca,lscgpa,lscgpc,  &
             lscgia,lscgic,lscgii,  &
             lscgip,lscgsa,lscgsc,  &
             lscgsi,lscgsp,lscgss,  &
             lscnd,                 &
             lscndgas,              &
             lscndh2oae,lscndh2ocl, &
             lscndh2oic,            & ! ice'n'snow
             lsauto,                &
             lsautosnow,            &
             lsactiv,lsactbase,     &
             lsactintst,            &
             lsichom,               &
             lsichet,               &
             lsicimmers,            &
             lsicmelt

   PUBLIC :: nspec, listspec, maxspec, nmod
   PUBLIC :: sigmag, dpg, n, volDistA, volDistB, nf2a
   PUBLIC :: rhlim
   PUBLIC :: isdtyp
   PUBLIC :: initliqice


   ! Datatype used to store information about the binned size distributions of aerosols,cloud,drizzle and ice
   ! ---------------------------------------------------------------------------------------------------------
   TYPE t_section
      REAL :: vhilim,     & ! bin volume at the high limit
              vlolim,     & ! - '' - at the low limit
              vratiohi,   & ! volume ratio between the center and high limit
              vratiolo,   & ! - '' - and the low limit
              dmid,       & ! bin middle diameter
         !******************************************************
         ! ^ Do NOT change the stuff above after initialization !
         !******************************************************
              dwet,       & ! Wet diameter or mean droplet diameter

              volc(8),    & ! Volume concentrations of aerosol species + water
                       ! Since most of the stuff in SALSA is hard coded, these
                       ! *have to be* in the order: 1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O

              veqh2o,     & ! Equilibrium h2o concentration for each particle
              numc,       & ! Number concentration of particles/droplets
              core          ! Volume of dry particle
   END TYPE t_section
   ! ---------------------------------------------------------------------------------------------------

   TYPE t_parallelbin
      ! Map bin indices between parallel size distributions
      INTEGER :: cur  ! Index for current distribution
      INTEGER :: par  ! Index for corresponding parallel distribution
   END TYPE t_parallelbin


   !Switches for SALSA aerosol microphysical processes
   LOGICAL :: nldebug = .FALSE., & ! debuggin output
              debug

   ! Process switches: nl* is read from the NAMELIST and NOT changed during runtime.
   !                   ls* is the switch actually used and will get the value of nl*
   !                   except for special circumstances such as spinup period etc.

   LOGICAL :: nlcoag  = .TRUE., & ! Coagulation master switch
              lscoag

   LOGICAL :: nlcgaa  = .TRUE., & ! Coagulation between aerosols
              lscgaa

   LOGICAL :: nlcgcc  = .TRUE., & ! Collision-coalescence between cloud droplets
              lscgcc

   LOGICAL :: nlcgca  = .TRUE., & ! Cloud collection of aerosols
              lscgca

   LOGICAL :: nlcgpc  = .TRUE., & ! Collection of cloud droplets by rain
              lscgpc

   LOGICAL :: nlcgpa  = .TRUE., & ! Collection of aerosols by rain
              lscgpa

   LOGICAL :: nlcgpp  = .TRUE., & ! Collision between rain drops
              lscgpp

   LOGICAL :: nlcgia  = .TRUE., & ! Ice collection of aerosols
              lscgia

   LOGICAL :: nlcgic  = .TRUE., & ! Collection of cloud droplets by ice particles
              lscgic

   LOGICAL :: nlcgii  = .TRUE., & ! Collision-coalescence between ice particles
              lscgii

   LOGICAL :: nlcgip  = .TRUE., & ! Collection of precipitation by ice particles
              lscgip

   LOGICAL :: nlcgsa  = .TRUE., & ! Collection of aerosols by snow
              lscgsa

   LOGICAL :: nlcgsc  = .TRUE., & ! Collection of cloud droplets by snow
              lscgsc

   LOGICAL :: nlcgsi  = .TRUE., & ! Collection of ice by snow
              lscgsi

   LOGICAL :: nlcgsp  = .TRUE., & ! Collection of precipitation by snow
              lscgsp

   LOGICAL :: nlcgss  = .TRUE., & ! Collision-coalescence between snow particles
              lscgss


   LOGICAL :: nlcnd      = .TRUE., & ! Condensation
              lscnd

   LOGICAL :: nlcndgas   = .TRUE., & ! Condensation of precursor gases
              lscndgas

   LOGICAL :: nlcndh2ocl = .TRUE., & ! Condensation of water vapour on clouds (drizzle)
              lscndh2ocl

   LOGICAL :: nlcndh2oae = .TRUE., & ! Condensation of water vapour on aerosol particles (FALSE -> equilibrium calc.)
              lscndh2oae

   LOGICAL :: nlcndh2oic = .TRUE., & ! Condensation of water vapour on ice particles
              lscndh2oic


   LOGICAL :: nlauto      = .TRUE.,   & ! Autoconversion of cloud droplets (needs activation)
              lsauto

   LOGICAL :: nlautosnow  = .TRUE.,   & ! Autoconversion of ice particles to snow (needs activation)
              lsautosnow

   LOGICAL :: nlactiv     = .TRUE.,   & ! Cloud droplet activation
              lsactiv

   LOGICAL :: nlactintst  = .TRUE.,   & ! Switch for interstitial activation: use particle wet size determined by
              lsactintst                ! codensation equations and supersaturation directly from the host model

   LOGICAL :: nlactbase   = .TRUE.,   & ! Switch for cloud base activation: use the regular parameterized method
              lsactbase                 ! for maximum supersaturation and cloud activation.


   LOGICAL :: nlichom     = .TRUE., & ! homogenous ice nucleation
              lsichom

   LOGICAL :: nlichet     = .TRUE., & ! heterogenous ice nucleation
              lsichet

   LOGICAL :: nlicimmers  = .TRUE., & ! ice nucleation by immersion
              lsicimmers

   LOGICAL :: nlicmelt    = .TRUE., & ! ice melting
              lsicmelt

   LOGICAL :: lsdistupdate = .TRUE.  ! Perform the size distribution update

   ! 1) Switches for aerosol microphysical processes ------------------------
   INTEGER, PARAMETER :: nmod = 7

   INTEGER :: nsnucl = 0             ! Choice of the nucleation scheme:
                                     ! 0 = off
                                     ! 1 = binary nucleation
                                     ! 2 = activation type nucleation
                                     ! 3 = kinetic nucleation
                                     ! 4 = ternary nucleation
                                     ! 5 = nucleation with ORGANICs
                                     ! 6 = activation type of nucleation with H2SO4+ORG
                                     ! 7 = heteromolecular nucleation with H2SO4*ORG
                                     ! 8 = homomolecular nucleation of  H2SO4 +
                                     !           heteromolecular nucleation with H2SO4*ORG
                                     ! 9 = homomolecular nucleation of  H2SO4 and ORG +
                                     !           heteromolecular nucleation with H2SO4*ORG
  
   ! 1) Switches for SALSA aerosol microphysical processes ------------------------

   INTEGER ::                    & ! J3 parametrization
      nj3 = 1              ! 1 = condensational sink (Kerminen&Kulmala, 2002)
                                   ! 2 = coagulational sink (Lehtinen et al. 2007)
                                   ! 3 = coagS+self-coagulation (Anttila et al. 2010)
   REAL :: act_coeff = 1.e-7  ! activation coefficient

   ! RH Limit: used for initialization and spinup within SALSA to limit the water vapour mixing ratio.
   ! Prevents unrealistically high RH in cloud activation and condensation procedures that is often assigned
   ! in the LES input files to immediately generate cloud. Given in %/100.
   REAL :: rhlim = 1.20
  
   ! Define which aerosol species used and initial size distributions
   INTEGER :: nspec = 1
   INTEGER, PARAMETER :: maxspec = 7
   CHARACTER(len=3) :: listspec(maxspec) = (/'SO4','   ','   ','   ','   ','   ','   '/)

   ! Volume fractions between aerosol species for A and B-bins
   REAL :: volDistA(maxspec) = (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
   REAL :: volDistB(maxspec) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
   ! Number fraction allocated to a-bins in regime 2 (b-bins will get 1-nf2a)
   REAL :: nf2a = 1.0

   ! Should not be necessary!
   LOGICAL :: initliqice = .FALSE. ! initialize ice and liquid cloud particles from aerosol bins
   REAL :: liqFracA = 0.0
   REAL :: iceFracA = 0.0
   REAL :: liqFracB = 0.0
   REAL :: iceFracB = 0.0

   INTEGER :: isdtyp = 0  ! Type of input aerosol size distribution: 0 - Uniform
                          !                                          1 - Read vertical profile of the mode
                          !                                              number concentrations from an input file

   REAL :: sigmag(nmod) = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0/),               & ! Stdev
           dpg(nmod) = (/0.03, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2/),    & ! Mode diam in um
           n(nmod) = (/1600.,640.,0.,0.,0.,0.,0./)                   ! #/mg ~ #/cm3

   INTEGER, PARAMETER ::            &
      nreg = 2                          ! number of main size regimes

   REAL ::                       &
      reglim(nreg+2) =                            & ! low/high diameter limits
      (/ 3.e-9, 5.e-8, 7.e-7, 1.e-5 /) ! of main size regimes [m]

   INTEGER :: &
      nbin(nreg) = (/ 3, 7 /)   ! number of bins in each main regime

   INTEGER ::      &
      nbin2,         & ! number of bins in former 2-region
      nbin3            ! number of bins in former 3-region

   INTEGER ::      & ! start index for bin regimes
      in1a,          & ! regime 1a
      in2a,          & ! regime 2a
      in2b,          & ! regime 2b
                       ! last index for bin regimes
      fn1a,          & ! regime 1a
      fn2a,          & ! regime 2a
      fn2b,          & ! regime 2b
      nbins            ! total number of size bins

    ! Juha: Cloud and rain bins:
   INTEGER :: ncldbin(2) = (/7,7/)        ! Number of bins for cloud bins in regime a and b
                                         
   TYPE(t_parallelbin) ::   ica, & ! cloud droplets (first, regime a)
                            fca, & ! cloud droplets (last, regime a)
                            icb, & ! cloud droplets (first, regime b)
                            fcb    ! cloud droplets (last, regime b)

   INTEGER             ::   ira,fra! Rain/drizzle bin indices
   INTEGER             ::   ncld   ! Total number of cloud bins
   INTEGER             ::   nprc   ! Total number of precipitation bins

   REAL                ::   dmincld = 5.e-8   ! Minimum diameter for the cloud droplet regime in terms of the
                                              ! ccn dry radius. The first cloud droplet bin is taken to coincide
                                              ! with the smallest full aerosol bin that conforms with this diameter.

    ! Jaakko: ice bins:
   INTEGER :: nicebin(2) = (/7,7/)        ! Number of bins for ice bins in regime a and b

   TYPE(t_parallelbin) ::   iia, & ! ice particles (first, regime a)
                            fia, & ! ice particles (last, regime a)
                            iib, & ! ice particles (first, regime b)
                            fib    ! ice particles (last, regime b)

   INTEGER             ::   isa,fsa! snow bin indices
   INTEGER             ::   nice   ! Total number of ice bins
   INTEGER             ::   nsnw   ! Total number of snow bins

   REAL                :: dminice = 5.e-8    ! Minimum diameter for the ice particle regime in terms of the
                                             ! ccn dry radius. The first cloud droplet bin is taken to coincide
                                             ! with the smallest full aerosol bin that conforms with this diameter

   REAL, ALLOCATABLE :: aerobins(:),  &  ! These are just to deliver information about the bin diameters if the
                        cloudbins(:), &  ! host model needs it (lower limits).
                        precpbins(:), &
                        icebins(:), &
                        snowbins(:)

   REAL, ALLOCATABLE :: vhilim(:),         &
                        vlolim(:),         &
                        vratiohi(:),       &
                        vratiolo(:),       &
                        dpmid(:)

   REAL, PARAMETER ::     &
      avog   = 6.0221e+23,   & ! Avogadro number (#/mol)
      boltz  = 1.3807e-23,   & ! Boltzmann constant (J/K)
      planck = 6.626070040e-34, & ! Planck constant (J*s)
      grav   = 9.81,         & ! gravitational acceleration (m/s^2)
      pstand = 1.01325e+5,   & ! standard pressure (Pa)
      rg     = 8.314,        & ! molar gas constant (J/(mol K))
      pi     = 3.1415927,    & ! self explanatory
      pi6    = 0.5235988,    & ! pi/6
      cpa    = 1010.,        & ! specific heat of dry air, constant P (J/kg/K)
      mair   = 28.97e-3,     & ! molar mass of air (mol/kg)
      deltav = 1.096e-7,     & ! vapor jump length (m)
      deltaT = 2.16e-7,      & ! thermal jump length (m)
      alphaT = 0.96,         & ! thermal accomodation coefficient
      alphac = 1.0,          & ! condensation coefficient
      eps    = epsilon(1.0)       ! epsilon


   REAL, PARAMETER ::     & ! molar mass [kg/mol]
      msu = 98.08e-3,        & ! sulphate
      !msu = 132.14e-3,        & ! ammonium sulphate (for ASCOS simulations) TR
      mno = 62.01e-3,        & ! HNO3
      mnh = 18.04e-3,        & ! NH3
      moc = 150.e-3,         & ! organic carbon
      mbc = 12.e-3,          & ! black carbon
      mss = 58.44e-3,        & ! sea salt (NaCl)
      mdu = 100.e-3,         & ! mineral dust
      mwa = 18.016e-3,       & ! water
      mas = 132.14e-3,       & ! ammoniums sulphate ((NH4)2SO4)
                                  !
                                  ! densities [kg/m3]
      rhosu = 1830.,         & ! sulphate
      !rhosu = 1770.,         & ! ammoniun sulphate (for ASCOS simulations) TR
      rhono = 1479.,         & ! HNO3
      rhonh = 1530.,         & ! NH3
      rhooc = 2000.,         & ! organic carbon
      rhobc = 2000.,         & ! black carbon
      rhoss = 2165.,         & ! sea salt (NaCl)
      rhodu = 2650.,         & ! mineral dust
      rhowa = 1000.,         & ! water
      rhoic = 917.,          & ! ice
      rhosn = 300.,          & ! snow !!new snow density
                                  !
                                  ! volume of molecule [kg/#]
      mvsu = msu/avog/rhosu,    & ! sulphate
      mvno = mno/avog/rhono,    & ! HNO3
      mvnh = mnh/avog/rhonh,    & ! NH3
      mvoc = moc/avog/rhooc,    & ! organic carbon
      mvss = mss/avog/rhoss,    & ! sea salt
      mvwa = mwa/avog/rhowa,    &
                                 !
      volratio =                & ! ratio of molecular volumes for
      (msu*rhoss)/(rhosu*mss), & ! sulphate and sea salt
                                  !
      n3 = 158.79               ! number of H2SO4 molecules in 3 nm cluster
                                !  assuming d_sa = 5.54 ???

   ! Dry air and water vapour thermodynamic constants
   REAL, PARAMETER ::  &
      rd  = 287.05,    &        ! gas constant for dry air in J/K/kg
      rv  = 461.51,    &        ! gas constant for water vapour in J/K/kg
      alv = 2.5008e6,  &        ! latent heat for vaporisation in J/kg
      als = 2.8345e6,  &        ! latent heat for sublimation in J/kg
      alf = als-alv             ! latent heat for fusion in J/kg

   !-- 4.3) Properties of condensing vapours

   REAL, PARAMETER :: & ! diameter of condensing molecule [m]
      d_sa  = 5.539376964394570e-10,  &
      d_oc  = 6.195906936656752e-10,  &
      d_h2o = 3.851565216195334e-10

   REAL, PARAMETER :: &
      ions   = 3.0,       & ! van't Hoff factor (ions produced upon dissociation)
      surfw0 = 0.073,     & ! surface tension of pure water @ ~ 293 K [J/m2]
      surfi0 = 0.105,     & ! surface tension of ice
      epsoc  = 0.15         ! water uptake of organic material


   !-- 7) Parameters for cloud activation

   REAL, ALLOCATABLE :: massacc(:)


   REAL, PARAMETER :: &
      nlim  = 1.,  & ! Number conc. limit (#/kg) for aerosol and cloud droplets
      prlim = 1.e-6  ! The same for precipitation and ice species for which concentrations are normally much lower [#/m3]
  

  !--- 12) Service routines for initialization and auxiliary computations ----------

END MODULE mo_submctl
