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
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE grid

   USE ncio, ONLY : open_nc, define_nc
   USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
   USE classFieldArray
   USE class_componentIndex, ONLY : componentIndex

   IMPLICIT NONE

   INTEGER :: nxp = 132           ! number of x points
   INTEGER :: nyp = 132           ! number of y points
   INTEGER :: nzp = 105           ! number of z points


   LOGICAL :: nxpart = .TRUE.     ! number of processors in x

   REAL    :: deltax = 35.        ! dx for basic grid
   REAL    :: deltay = 35.        ! dy for basic grid
   REAL    :: deltaz = 17.5       ! dz for basic grid
   REAL    :: dzrat  = 1.02       ! grid stretching ratio
   REAL    :: dzmax  = 1200.      ! height to start grid-stretching
   REAL    :: dtlong = 10.0       ! long timestep
   REAL    :: th00   = 288.       ! basic state temperature

   REAL    :: CCN = 150.e6

   LOGICAL :: lbinanl = .FALSE.   ! Whether to write binned data to analysis files (takes a lot of space + mainly used for debugging)
   LOGICAL :: lbinprof = .TRUE.   ! Same for profile statistics
   INTEGER :: igrdtyp = 1         ! vertical grid type
   INTEGER :: isgstyp = 1         ! sgs model type
   INTEGER :: iradtyp = 0         ! radiation model type
   INTEGER :: level   = 0         ! thermodynamic level
   INTEGER :: naddsc  = 0         ! number of additional scalars;
   INTEGER :: nsalsa  = 0         ! Number of tracers for SALSA
   INTEGER :: nfpt = 10           ! number of rayleigh friction points
   REAL    :: distim = 300.0      ! dissipation timescale

   REAL    :: sst = 283.   ! Surface temperature      added by Zubair Maalick
   REAL    :: W1  = 0.9    ! Water content
   REAL    :: W2  = 0.9
   REAL    :: W3  = 0.9


   !CHARACTER (len=7), ALLOCATABLE, SAVE :: sanal(:)
   CHARACTER (len=80) :: expnme = 'Default' ! Experiment name
   CHARACTER (len=80) :: filprf = 'x'       ! File Prefix
   CHARACTER (len=7)  :: runtype = 'INITIAL'! Run Type Selection

   REAL               :: Tspinup = 7200.    ! Spinup period in seconds (added by Juha)


   CHARACTER (len=7),  PRIVATE :: v_snm = 'sxx    '
   CHARACTER (len=80), PRIVATE :: fname

   INTEGER, PRIVATE, SAVE  ::  nrec0, nvar0, varnum, nbase=15

   INTEGER           :: nz, nxyzp, nxyp
   REAL              :: dxi, dyi, dtl, dtlv, dtlt, umean, vmean, psrf
   REAL, ALLOCATABLE :: xt(:), xm(:), yt(:), ym(:), zt(:), zm(:), dzt(:), dzm(:)
   REAL, ALLOCATABLE :: u0(:), v0(:), pi0(:), pi1(:), th0(:), dn0(:), rt0(:)
   REAL, ALLOCATABLE :: spng_wfct(:), spng_tfct(:)
   REAL, ALLOCATABLE, target :: tmp_icep(:,:,:,:), tmp_icet(:,:,:,:)
   !
   ! velocity variables (past, current and tendency)
   !
   REAL, ALLOCATABLE, TARGET :: a_up(:,:,:), a_uc(:,:,:), a_ut(:,:,:)
   REAL, ALLOCATABLE, TARGET :: a_vp(:,:,:), a_vc(:,:,:), a_vt(:,:,:)
   REAL, ALLOCATABLE, TARGET :: a_wp(:,:,:), a_wc(:,:,:), a_wt(:,:,:)
   !
   ! wsave variables used in fft in x and y directons
   !
   REAL, ALLOCATABLE :: wsavex(:), wsavey(:)
   !
   ! prognostic scalar variables
   !
   TYPE(FloatArray3D), TARGET :: a_tp, a_tt

   TYPE(FloatArray3D), TARGET :: a_rp, a_rt  !Juha: In standard version this is the TOTAL water content.
                                             !      With SALSA this is taken as just the water VAPOUR content,
                                             !      in order not to over-specify the problem.
   TYPE(FloatArray3D), TARGET :: a_rpp, a_rpt
   TYPE(FloatArray3D), TARGET :: a_npp, a_npt

   TYPE(FieldArray) :: fields

   TYPE(FloatArray3D), TARGET :: a_qp, a_qt
   REAL, POINTER :: a_sp(:,:,:), a_st(:,:,:)  !dont touch yet AZ

   ! Juha: SALSA tracers
   !---------------------------------------------------------------------------
   ! -- Masses given in kg/kg, number concentrations in #/kg
   ! -- Each size bin/species will be treated as a separate tracer.
   ! -- NOTE: Volume concentration arrays are reduced to 4 dims.
   !          The 4th dim contains all the size bins sequentially for
   !          each aerosol species  + water
   !
   !          Gas tracers are contained sequentially in dimension
   !          4 as: 1. SO4, 2. HNO3, 3. NH3, 4. OCNV, 5. OCSV

   ! Prognostic tracers
   ! -- Number concentrations
   TYPE(FloatArray4D), TARGET :: a_naerop,  a_naerot,  &
                                 a_ncloudp, a_ncloudt, &
                                 a_nprecpp, a_nprecpt, &
                                 a_nicep,   a_nicet,   &
                                 a_nsnowp,  a_nsnowt
   ! -- Volume concentrations
   TYPE(FloatArray4D), TARGET :: a_maerop,  a_maerot,  &
                                 a_mcloudp, a_mcloudt, &
                                 a_mprecpp, a_mprecpt, &
                                 a_micep,   a_micet,   &
                                 a_msnowp,  a_msnowt
   ! -- Gas compound tracers
   TYPE(FloatArray4D), TARGET :: a_gaerop, a_gaerot

   ! Some stuff for tendency formulation
   REAL, ALLOCATABLE :: a_vactd(:,:,:,:), a_nactd(:,:,:,:)

   ! Particle component index tables
   TYPE(componentIndex) :: prtcl ! Contains "getIndex" which gives the index for a given
                                ! aerosol component name, i.e. 1:SO4, 2:OC, 3:BC, 4:DU,
                                ! 5:SS, 6:NO, 7:NH, 8:H2O

   ! Physical properties for the selected aerosol components
   REAL, ALLOCATABLE :: dens(:), mws(:), diss(:), dens_ice(:), dens_snow(:)

   !---------------------------------------------------------------------------

   REAL, ALLOCATABLE, TARGET :: a_sclrp(:,:,:,:),a_sclrt(:,:,:,:)
    !
   ! 3d diagnostic quantities
   !
   TYPE(FloatArray3D), TARGET :: a_theta  ! dry potential temp (k)
   TYPE(FloatArray3D), TARGET :: a_pexnr  ! perturbation exner func
   TYPE(FloatArray3D), TARGET :: a_press  ! pressure (hpa)
   TYPE(FloatArray3D), TARGET :: a_rc     ! Total cloud water +rain (level<=3) or aerosol+cloud (level>=4) water mixing ratio
   TYPE(FloatArray3D), TARGET :: a_ri     ! Total ice water mixing ratio
   TYPE(FloatArray3D), TARGET :: a_rv     ! Water vapor mixing ratio  (used only for levels < 4!)
   TYPE(FloatArray3D), TARGET :: a_srp    ! Total rain water mixing ratio (levels >= 4)(Diagnostic scalar!!)
   TYPE(FloatArray3D), TARGET :: a_snrp   ! Total number drop number mixing ratio (levels >=4) (Diagnostic scalar!!)
   TYPE(FloatArray3D), TARGET :: a_srs    ! Total snow water mixing ratio (level 5)
   TYPE(FloatArray3D), TARGET :: a_snrs   ! Total snow number mixing ratio (level 5) (Diagnostic scalar!!)
   TYPE(FloatArray3D), TARGET :: a_rh     ! Relative humidity
   TYPE(FloatArray3D), TARGET :: a_rsl    ! Water vapor saturation mixing ratio
   TYPE(FloatArray3D), TARGET :: a_rhi    ! Relative humidity over ice
   TYPE(FloatArray3D), TARGET :: a_rsi    ! Water vapor saturation mixing ratio over ice
   TYPE(FloatArray3D), TARGET :: a_dn     ! Air density (for normalizing concentrations according to mass, levels < 4!)
   !
   ! scratch arrays
   !
   REAL, ALLOCATABLE, DIMENSION (:,:,:) :: a_rflx, a_sflx, &
                                           a_fus, a_fds, a_fuir, a_fdir, &
                                           a_temp
   !
   !
   REAL, ALLOCATABLE :: a_ustar(:,:)
   REAL, ALLOCATABLE :: a_tstar(:,:)
   REAL, ALLOCATABLE :: a_rstar(:,:)
   REAL, ALLOCATABLE :: uw_sfc(:,:)
   REAL, ALLOCATABLE :: vw_sfc(:,:)
   REAL, ALLOCATABLE :: ww_sfc(:,:)
   REAL, ALLOCATABLE :: wt_sfc(:,:)
   REAL, ALLOCATABLE :: wq_sfc(:,:)
   REAL, ALLOCATABLE :: precip(:,:,:), snowin(:,:,:), albedo(:,:)

   ! Juha:
   ! Diagnostic variables needed to track mass conservation (of water).
   ! These are reset at every statistical output timestep (better use quite long statistical periods...)
   REAL :: mc_Mtot           ! Initial mass of water normalized by domain volume (== domain mean concentration)
   REAL :: mc_Matm           ! Atmospheric water content (instantaneous, domain mean concentration)
   REAL :: mc_Mevap          ! Evaporated water content, accumulated, Normalized by *domain surface area*/*domain volume*
   REAL :: mc_Mprec          ! Precipitated water content, - '' -
   REAL :: mc_Vdom           ! Domain volume
   REAL :: mc_Adom           ! Domain surface area
   REAL :: mc_ApVdom         ! Volume/Area

   CLASS(*), POINTER :: dummy, dummyt, dummyp, dummy2, dummy3 !dummy pointers for field arrays

   TYPE(FloatArray3d), POINTER :: testi_diag

   TYPE(FieldArray) :: ProgF
   TYPE(FieldArray) :: DiagF

   ! Targets for output variables
   TYPE(FloatArray3D), TARGET :: a_u, a_v, a_w, a_the, a_p, a_q, a_l, a_r, a_f, a_i, a_s, a_n, a_stke, a_rfl

   TYPE(FloatArray3D), TARGET :: a_Nac, a_Na,   a_Rwaa, a_Nb, a_Rwab, &
                                 a_Nc,  a_Rwca, a_Rwcb, a_Np, a_Rwpa, &
                                 a_Ni,  a_Rwia, a_Rwib, a_Ns, a_Rwsa, &
                                 a_aSO4a, a_aNH3a, a_aNO3a, a_aOCa, a_aBCa, a_aDUa, a_aSSa, &
                                 a_aSO4b, a_aNH3b, a_aNO3b, a_aOCb, a_aBCb, a_aDUb, a_aSSb, &
                                 a_cSO4a, a_cNH3a, a_cNO3a, a_cOCa, a_cBCa, a_cDUa, a_cSSa, &
                                 a_cSO4b, a_cNH3b, a_cNO3b, a_cOCb, a_cBCb, a_cDUb, a_cSSb, &
                                 a_iSO4a, a_iNH3a, a_iNO3a, a_iOCa, a_iBCa, a_iDUa, a_iSSa, &
                                 a_iSO4b, a_iNH3b, a_iNO3b, a_iOCb, a_iBCb, a_iDUb, a_iSSb

   TYPE(FloatArray4D), TARGET :: a_Naba, a_Rwaba, a_Nabb,  a_Rwabb, a_Ncba, a_Ncbb,  a_Rwcba, a_Rwcbb, &
                                 a_Npba, a_Rwpba, a_Rwiba, a_Niba,  a_Nibb, a_Rwibb, a_Nsba,  a_Rwsba


   !
   INTEGER :: nscl = 1
   INTEGER, SAVE :: ncid0,ncid_s
   !
CONTAINS
   !
   !----------------------------------------------------------------------
   ! Subroutine define_vars
   !
   ! Modified for level 4
   ! Juha Tonttila, FMI, 2014.
   !
   SUBROUTINE define_vars

      USE mpi_interface, ONLY : myid
      USE mo_submctl, ONLY : nbins, ncld, nprc,  & ! Number of aerosol and hydrometeor size bins for SALSA
                             nice,  nsnw,        & ! number of ice and snow size bins for SALSA
                             nspec, maxspec, listspec, &
                             rhosu, rhooc, rhobc, rhono, rhonh, rhoss, rhodu, rhowa, rhoic, rhosn,  &
                             msu,   moc,   mbc,   mno,   mnh,   mss,   mdu,   mwa

      USE class_ComponentIndex, ONLY : ComponentIndexConstructor,  &
                                       GetNcomp, IsUsed

      INTEGER :: memsize
      INTEGER :: zz
      INTEGER :: nc
      REAL :: zeros3d(nzp,nxp,nyp)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros4d(nzp,nxp,nyp,nbins+ncld+nprc+nice+nsnw)

      zeros3d(:,:,:) = 0.
      zeros4d(:,:,:,:) = 0.

      ProgF = FieldArray()
      DiagF = FieldArray()

      ! Juha: Number of prognostic tracers for SALSA
      !            Aerosol bins + Cloud bins + gas compound tracers

      IF (level >= 4) THEN
         ! Create index tables for different aerosol components (can be fetched by name using getIndex)
         CALL ComponentIndexConstructor(prtcl, nspec, maxspec, listspec)
         nc = GetNcomp(prtcl)

         nsalsa = (nc+2)*nbins + (nc+2)*ncld + (nc+2)*nprc + 5
         IF (level>=5) nsalsa = nsalsa + (nc+2)*nice + (nc+2)*nsnw

      END IF

      ! Juha: Stuff that's allocated for all configurations
      !----------------------------------------------------------
      ALLOCATE (u0(nzp),v0(nzp),pi0(nzp),pi1(nzp),th0(nzp),dn0(nzp),rt0(nzp))

      memsize = 2*nxyzp ! complexarray in pressure solver

      ALLOCATE (a_up(nzp,nxp,nyp),a_vp(nzp,nxp,nyp),a_wp(nzp,nxp,nyp))
      a_up(:,:,:) = 0.
      a_vp(:,:,:) = 0.
      a_wp(:,:,:) = 0.

      ALLOCATE (a_uc(nzp,nxp,nyp),a_vc(nzp,nxp,nyp),a_wc(nzp,nxp,nyp))
      a_uc(:,:,:) = 0.
      a_vc(:,:,:) = 0.
      a_wc(:,:,:) = 0.

      ALLOCATE (a_ut(nzp,nxp,nyp),a_vt(nzp,nxp,nyp),a_wt(nzp,nxp,nyp))
      a_ut(:,:,:) = 0.
      a_vt(:,:,:) = 0.
      a_wt(:,:,:) = 0.

      a_theta = FloatArray3D(zeros3d,store=.TRUE.)
      dummy => a_theta
      CALL DiagF%NewField("theta","Potential temperature","K","tttt",.TRUE.,dummy)

      a_pexnr = FloatArray3D(zeros3d,store=.TRUE.)
      dummy => a_pexnr
      CALL DiagF%NewField("pexnr"," "," "," ",.FALSE.,dummy)

      a_press = FloatArray3D(zeros3d,store=.TRUE.)
      dummy => a_press
      CALL DiagF%NewField("press","Pressure","Pa","tttt",.TRUE.,dummy)

      memsize = memsize + nxyzp*13 !

      IF (iradtyp > 0 ) THEN
         ALLOCATE (a_rflx(nzp,nxp,nyp))
         a_rflx(:,:,:) = 0.
         memsize = memsize + nxyzp
      END IF
      IF (iradtyp >= 3) THEN
         ALLOCATE (a_sflx(nzp,nxp,nyp),albedo(nxp,nyp))
         a_sflx(:,:,:) = 0.
         albedo(:,:) = 0.
         ALLOCATE (a_fus(nzp,nxp,nyp),a_fds(nzp,nxp,nyp),a_fuir(nzp,nxp,nyp),a_fdir(nzp,nxp,nyp))
         a_fus(:,:,:) = 0.
         a_fds(:,:,:) = 0.
         a_fuir(:,:,:) = 0.
         a_fdir(:,:,:) = 0.
         memsize = memsize + nxyzp + nxyp + 4*nxyp
      END IF


      ALLOCATE (a_temp(nzp,nxp,nyp))
      a_temp(:,:,:) = 0.

      a_rsl = FloatArray3D(zeros3d,store=.TRUE.)
      dummy => a_rsl
      CALL DiagF%NewField("rsl"," "," "," ",.FALSE.,dummy)

      memsize = memsize + nxyzp*3

      ! Juha: Stuff that's allocated if SALSA is not used
      !-----------------------------------------------------
      IF (level < 4) THEN

         IF (level >= 0) THEN
            a_rv = FloatArray3D(zeros3d,store=.TRUE.)
            dummy => a_rv
            CALL DiagF%NewField("rv"," "," "," ",.FALSE.,dummy)

            memsize = memsize + nxyzp
            IF (level > 1) THEN
               a_rc = FloatArray3D(zeros3d,store=.TRUE.)
               dummy => a_rc
               CALL DiagF%NewField("rc"," "," "," ",.FALSE.,dummy)

               memsize = memsize + nxyzp
            END IF
         END IF

         nscl = nscl+naddsc
         IF (level   > 0) nscl = nscl+1
         IF (level   > 2) nscl = nscl+2
         IF (isgstyp > 1) nscl = nscl+1

         ALLOCATE (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
         a_sclrp(:,:,:,:) = 0.
         a_sclrt(:,:,:,:) = 0.

         a_tp = FloatArray3D(a_sclrp(:,:,:,1),store=.FALSE.)
         a_tt = FloatArray3D(a_sclrt(:,:,:,1),store=.FALSE.)
         dummyp => a_tp
         dummyt => a_tt
         CALL ProgF%NewField("t","Liquid water potential temperature","K","tttt",.FALSE.,dummyp,dummyt)

         IF (level >= 0) THEN
            a_rp = FloatArray3D(a_sclrp(:,:,:,2),store=.FALSE.)
            a_rt = FloatArray3D(a_sclrt(:,:,:,2),store=.FALSE.)
            dummyp => a_rp
            dummyt => a_rt
            CALL ProgF%NewField("r","Rain-water mixing ratio","kg/kg","tttt",.FALSE.,dummyp,dummyt)
         END IF

         IF (level >= 3) THEN
            a_rpp = FloatArray3D(a_sclrp(:,:,:,3),store=.FALSE.)
            a_rpt = FloatArray3D(a_sclrt(:,:,:,3),store=.FALSE.)
            dummyp => a_rpp
            dummyt => a_rpt
            CALL ProgF%NewField("rp"," "," "," ",.FALSE.,dummyp,dummyt)

            a_npp = FloatArray3D(a_sclrp(:,:,:,4),store=.FALSE.)
            a_npt = FloatArray3D(a_sclrt(:,:,:,4),store=.FALSE.)
            dummyp => a_npp
            dummyt => a_npt
            CALL ProgF%NewField("np"," "," "," ",.FALSE.,dummyp,dummyt)
         END IF
         IF (isgstyp > 1) THEN
            a_qp = FloatArray3D(a_sclrp(:,:,:,nscl - naddsc),store=.FALSE.)
            a_qt = FloatArray3D(a_sclrt(:,:,:,nscl - naddsc),store=.FALSE.)
            dummyp => a_qp
            dummyt => a_qt
            CALL ProgF%NewField("q"," "," "," ",.FALSE.,dummyp,dummyt)
         END IF


      ! Allocate float array targets of the field arrays
         a_u = FloatArray3D(zeros3d,store=.TRUE.)
         a_v = FloatArray3D(zeros3d,store=.TRUE.)
         a_w = FloatArray3D(zeros3d,store=.TRUE.)
         a_the = FloatArray3D(zeros3d,store=.TRUE.)
         a_p = FloatArray3D(zeros3d,store=.TRUE.)
         a_q = FloatArray3D(zeros3d,store=.TRUE.)
         a_l = FloatArray3D(zeros3d,store=.TRUE.)
         a_r = FloatArray3D(zeros3d,store=.TRUE.)
         a_f = FloatArray3D(zeros3d,store=.TRUE.)
         a_i = FloatArray3D(zeros3d,store=.TRUE.)
         a_s = FloatArray3D(zeros3d,store=.TRUE.)
         a_n = FloatArray3D(zeros3d,store=.TRUE.)
         a_stke = FloatArray3D(zeros3d,store=.TRUE.)
         a_rfl = FloatArray3D(zeros3d,store=.TRUE.)

      !Juha: Stuff that's allocated when SALSA is used
      !---------------------------------------------------
      ELSE IF (level >= 4) THEN

         nc = GetNcomp(prtcl) ! number of aerosol components used. For allocations + 1 for water.

         ALLOCATE (a_nactd(nzp,nxp,nyp,ncld), a_vactd(nzp,nxp,nyp,(nc+1)*ncld)  )

         a_rc = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_rc
         CALL DiagF%NewField("rc"," "," "," ",.FALSE.,dummy)

         a_srp = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_srp
         CALL DiagF%NewField("srp"," "," "," ",.FALSE.,dummy)

         a_snrp = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_snrp
         CALL DiagF%NewField("snrp"," "," "," ",.FALSE.,dummy)

         a_rh = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_rh
         CALL DiagF%NewField("rh","Relative humidity","%","tttt",.TRUE.,dummy)

         a_dn = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_dn
         CALL DiagF%NewField("dn"," "," "," ",.FALSE.,dummy)

         a_nactd(:,:,:,:) = 0.
         a_vactd(:,:,:,:) = 0.
         memsize = memsize + 4*nxyzp + 3*nbins*nxyzp + 3*ncld*nxyzp + nxyzp*(nc+1)*ncld + 2*nprc*nxyzp

         ! ice'n'snow
         a_ri = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_ri
         CALL DiagF%NewField("ri"," "," "," ",.FALSE.,dummy)

         a_rsi = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_rsi
         CALL DiagF%NewField("rsi"," "," "," ",.FALSE.,dummy)

         a_rhi = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_rhi
         CALL DiagF%NewField("rhi"," "," "," ",.FALSE.,dummy)

         a_srs = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_srs
         CALL DiagF%NewField("srs"," "," "," ",.FALSE.,dummy)

         a_snrs = FloatArray3D(zeros3d,store=.TRUE.)
         dummy => a_snrs
         CALL DiagF%NewField("snrs"," "," "," ",.FALSE.,dummy)

         memsize = memsize + 5*nxyzp + 2*nice*nxyzp + 2*nsnw*nxyzp

         ! Total number of prognostic scalars: temp + total water + SALSA + tke(?)
         nscl = 2 + nsalsa
         IF (isgstyp > 1) nscl = nscl+1

         ALLOCATE (a_sclrp(nzp,nxp,nyp,nscl), a_sclrt(nzp,nxp,nyp,nscl))
         a_sclrp(:,:,:,:) = 0.
         a_sclrt(:,:,:,:) = 0.

         a_tp = FloatArray3D(a_sclrp(:,:,:,1),store=.FALSE.)
         a_tt = FloatArray3D(a_sclrt(:,:,:,1),store=.FALSE.)
         dummyp => a_tp
         dummyt => a_tt
         CALL ProgF%NewField("t","Liquid water potential temperature","K","tttt",.TRUE.,dummyp,dummyt)
         a_rp = FloatArray3D(a_sclrp(:,:,:,2),store=.FALSE.)
         a_rt = FloatArray3D(a_sclrt(:,:,:,2),store=.FALSE.)
         dummyp => a_rp
         dummyt => a_rt
         CALL ProgF%NewField("r","Rain-water mixing ratio","kg/kg","tttt",.TRUE.,dummyp,dummyt)

         IF (isgstyp > 1) THEN
            a_qp = FloatArray3D(a_sclrp(:,:,:,nscl - nsalsa),store=.FALSE.)
            a_qt = FloatArray3D(a_sclrt(:,:,:,nscl - nsalsa),store=.FALSE.)
            dummyp => a_qp
            dummyt => a_qt
            CALL ProgF%NewField("q","Total water mixing ratio","kg/kg","tttt",.TRUE.,dummyp,dummyt)
         END IF

         !JT: Set the pointers for prognostic SALSA variables (levels 4 & 5)
         zz = nscl-nsalsa
         a_naerop = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+nbins),store=.FALSE.)
         a_naerot = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+nbins),store=.FALSE.)
         dummyp => a_naerop
         dummyt => a_naerot
         CALL ProgF%NewField("naero"," "," "," ",.FALSE.,dummyp,dummyt)


         zz = zz+nbins
         a_maerop = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+(nc+1)*nbins),store=.FALSE.)
         a_maerot = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+(nc+1)*nbins),store=.FALSE.)
         dummyp => a_maerop
         dummyt => a_maerot
         CALL ProgF%NewField("maero"," "," "," ",.FALSE.,dummyp,dummyt)

         zz = zz+(nc+1)*nbins
         a_ncloudp = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+ncld),store=.FALSE.)
         a_ncloudt = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+ncld),store=.FALSE.)
         dummyp => a_ncloudp
         dummyt => a_ncloudt
         CALL ProgF%NewField("ncloud"," "," "," ",.FALSE.,dummyp,dummyt)

         zz = zz+ncld
         a_mcloudp = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+(nc+1)*ncld),store=.FALSE.)
         a_mcloudt = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+(nc+1)*ncld),store=.FALSE.)
         dummyp => a_mcloudp
         dummyt => a_mcloudt
         CALL ProgF%NewField("mcloud"," "," "," ",.FALSE.,dummyp,dummyt)

         zz = zz+(nc+1)*ncld
         a_nprecpp = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+nprc),store=.FALSE.)
         a_nprecpt = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+nprc),store=.FALSE.)
         dummyp => a_nprecpp
         dummyt => a_nprecpt
         CALL ProgF%NewField("nprecp"," "," "," ",.FALSE.,dummyp,dummyt)

         zz = zz+nprc
         a_mprecpp = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+(nc+1)*nprc),store=.FALSE.)
         a_mprecpt = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+(nc+1)*nprc),store=.FALSE.)
         dummyp => a_mprecpp
         dummyt => a_mprecpt
         CALL ProgF%NewField("mprecp"," "," "," ",.FALSE.,dummyp,dummyt)

         zz = zz+(nc+1)*nprc
         a_gaerop = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+5),store=.FALSE.)
         a_gaerot = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+5),store=.FALSE.)
         dummyp => a_gaerop
         dummyt => a_gaerot
         CALL ProgF%NewField("gaero"," "," "," ",.FALSE.,dummyp,dummyt)

         IF (level >= 5) THEN
            ! Level 5
            zz = zz+5
            a_nicep = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+nice),store=.FALSE.)
            a_nicet = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+nice),store=.FALSE.)
            dummyp => a_nicep
            dummyt => a_nicet
            CALL ProgF%NewField("nice"," "," "," ",.FALSE.,dummyp,dummyt)

            zz = zz+nice

            a_micep = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+(nc+1)*nice),store=.FALSE.)
            a_micet = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+(nc+1)*nice),store=.FALSE.)
            dummyp => a_micep
            dummyt => a_micet
            CALL ProgF%NewField("mice"," "," "," ",.FALSE.,dummyp,dummyt)

            zz = zz+(nc+1)*nice

            a_nsnowp = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+nsnw),store=.FALSE.)
            a_nsnowt = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+nsnw),store=.FALSE.)
            dummyp => a_nsnowp
            dummyt => a_nsnowt
            CALL ProgF%NewField("nsnow"," "," "," ",.FALSE.,dummyp,dummyt)

            zz = zz+nsnw

            a_msnowp = FloatArray4D(a_sclrp(:,:,:,zz+1:zz+(nc+1)*nsnw),store=.FALSE.)
            a_msnowt = FloatArray4D(a_sclrt(:,:,:,zz+1:zz+(nc+1)*nsnw),store=.FALSE.)
            dummyp => a_msnowp
            dummyt => a_msnowt
            CALL ProgF%NewField("msnow"," "," "," ",.FALSE.,dummyp,dummyt)
         ELSE
            ! Ice not included so allocate zero arrays for ice pointers
            ALLOCATE (tmp_icep(nzp,nxp,nyp,(nc+1)*MAX(nice,nsnw)), &
                          tmp_icet(nzp,nxp,nyp,(nc+1)*MAX(nice,nsnw)))
            tmp_icep = 0.
            tmp_icet = 0.
            a_nicep = FloatArray4D(tmp_icep(:,:,:,1:nice),store=.FALSE.)
            a_nicet = FloatArray4D(tmp_icet(:,:,:,1:nice),store=.FALSE.)
            dummyp => a_nicep
            dummyt => a_nicet
            CALL ProgF%NewField("nice"," "," "," ",.FALSE.,dummyp,dummyt)

            a_micep = FloatArray4D(tmp_icep(:,:,:,1:(nc+1)*nice),store=.FALSE.)
            a_micet = FloatArray4D(tmp_icet(:,:,:,1:(nc+1)*nice),store=.FALSE.)
            dummyp => a_micep
            dummyt => a_micet
            CALL ProgF%NewField("mice"," "," "," ",.FALSE.,dummyp,dummyt)

            a_nsnowp = FloatArray4D(tmp_icep(:,:,:,1:nsnw),store=.FALSE.)
            a_nsnowt = FloatArray4D(tmp_icet(:,:,:,1:nsnw),store=.FALSE.)
            dummyp => a_nsnowp
            dummyt => a_nsnowt
            CALL ProgF%NewField("nsnow"," "," "," ",.FALSE.,dummyp,dummyt)

            a_msnowp = FloatArray4D(tmp_icep(:,:,:,1:(nc+1)*nsnw),store=.FALSE.)
            a_msnowt = FloatArray4D(tmp_icet(:,:,:,1:(nc+1)*nsnw),store=.FALSE.)
            dummyp => a_msnowp
            dummyt => a_msnowt
            CALL ProgF%NewField("msnow"," "," "," ",.FALSE.,dummyp,dummyt)

         END IF


         ! Density, molecular weight, dissociation factor and molar volume arrays for the used species
         !   1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O
         ALLOCATE ( dens(nc+1), mws(nc+1), diss(nc+1), dens_ice(nc+1), dens_snow(nc+1) )
         zz = 0
         IF (IsUsed(prtcl,'SO4')) THEN
            zz = zz+1
            dens(zz) = rhosu
            diss(zz) = 3.  ! H2SO4
            mws(zz) = msu
         END IF
         IF (IsUsed(prtcl,'OC')) THEN
            zz = zz+1
            dens(zz) = rhooc
            diss(zz) = 1.
            mws(zz) = moc
         END IF
         IF (IsUsed(prtcl,'BC')) THEN
            zz = zz+1
            dens(zz) = rhobc
            diss(zz) = 0.  ! Insoluble
            mws(zz) = mbc
         END IF
         IF (IsUsed(prtcl,'DU')) THEN
            zz = zz+1
            dens(zz) = rhodu
            diss(zz) = 0.  ! Insoluble
            mws(zz) = mdu
         END IF
         IF (IsUsed(prtcl,'SS')) THEN
            zz = zz+1
            dens(zz) = rhoss
            diss(zz) = 2.  ! NaCl
            mws(zz) = mss
         END IF
         IF (IsUsed(prtcl,'NO')) THEN
            zz = zz+1
            dens(zz) = rhono
            diss(zz) = 1.  ! NO3- ??
            mws(zz) = mno
         END IF
         IF (IsUsed(prtcl,'NH')) THEN
            zz = zz+1
            dens(zz) = rhonh
            diss(zz) = 1.  ! NH4+ ??
            mws(zz) = mnh
         END IF

         ! Water is always there
         zz = zz+1
         dens(zz) = rhowa
         diss(zz) = 1.
         mws(zz) = mwa

         ! .. but sometimes it is frozen
         dens_ice=dens; dens_ice(zz)=rhoic
         dens_snow=dens; dens_snow(zz)=rhosn

         IF (zz /= nc+1) THEN
            WRITE(*,*) 'Physical properties not found for all species!'
            STOP
         END IF

         !Initialize float arrays for SALSA output variables

         a_u = FloatArray3D(zeros3d,store=.TRUE.)
         a_v = FloatArray3D(zeros3d,store=.TRUE.)
         a_w = FloatArray3D(zeros3d,store=.TRUE.)
         a_the = FloatArray3D(zeros3d,store=.TRUE.)
         a_p = FloatArray3D(zeros3d,store=.TRUE.)
         a_q = FloatArray3D(zeros3d,store=.TRUE.)
         a_l = FloatArray3D(zeros3d,store=.TRUE.)
         a_r = FloatArray3D(zeros3d,store=.TRUE.)
         a_f = FloatArray3D(zeros3d,store=.TRUE.)
         a_i = FloatArray3D(zeros3d,store=.TRUE.)
         a_s = FloatArray3D(zeros3d,store=.TRUE.)

         a_Nac  = FloatArray3D(zeros3d,store=.TRUE.)
         a_Na = FloatArray3D(zeros3d,store=.TRUE.)
         a_Naba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwaa = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwaba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Nb = FloatArray3D(zeros3d,store=.TRUE.)
         a_Nabb = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwab = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwabb = FloatArray4D(zeros4d,store=.TRUE.)

         a_Nc = FloatArray3D(zeros3d,store=.TRUE.)
         a_Ncba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Ncbb = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwca = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwcb = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwcba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwcbb = FloatArray4D(zeros4d,store=.TRUE.)

         a_Np = FloatArray3D(zeros3d,store=.TRUE.)
         a_Npba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwpa = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwpba = FloatArray4D(zeros4d,store=.TRUE.)

         a_Ni = FloatArray3D(zeros3d,store=.TRUE.)
         a_Niba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Nibb = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwia = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwib = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwiba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwibb = FloatArray4D(zeros4d,store=.TRUE.)

         a_Ns = FloatArray3D(zeros3d,store=.TRUE.)
         a_Nsba = FloatArray4D(zeros4d,store=.TRUE.)
         a_Rwsa  = FloatArray3D(zeros3d,store=.TRUE.)
         a_Rwsba = FloatArray4D(zeros4d,store=.TRUE.)

         a_aSO4a = FloatArray3D(zeros3d,store=.TRUE.)
         a_aNH3a = FloatArray3D(zeros3d,store=.TRUE.)
         a_aNO3a = FloatArray3D(zeros3d,store=.TRUE.)
         a_aOCa = FloatArray3D(zeros3d,store=.TRUE.)
         a_aBCa = FloatArray3D(zeros3d,store=.TRUE.)
         a_aDUa = FloatArray3D(zeros3d,store=.TRUE.)
         a_aSSa = FloatArray3D(zeros3d,store=.TRUE.)

         a_aSO4b = FloatArray3D(zeros3d,store=.TRUE.)
         a_aNH3b = FloatArray3D(zeros3d,store=.TRUE.)
         a_aNO3b = FloatArray3D(zeros3d,store=.TRUE.)
         a_aOCb = FloatArray3D(zeros3d,store=.TRUE.)
         a_aBCb = FloatArray3D(zeros3d,store=.TRUE.)
         a_aDUb = FloatArray3D(zeros3d,store=.TRUE.)
         a_aSSb = FloatArray3D(zeros3d,store=.TRUE.)

         a_cSO4a = FloatArray3D(zeros3d,store=.TRUE.)
         a_cNH3a = FloatArray3D(zeros3d,store=.TRUE.)
         a_cNO3a = FloatArray3D(zeros3d,store=.TRUE.)
         a_cOCa = FloatArray3D(zeros3d,store=.TRUE.)
         a_cBCa = FloatArray3D(zeros3d,store=.TRUE.)
         a_cDUa = FloatArray3D(zeros3d,store=.TRUE.)
         a_cSSa = FloatArray3D(zeros3d,store=.TRUE.)

         a_cSO4b = FloatArray3D(zeros3d,store=.TRUE.)
         a_cNH3b = FloatArray3D(zeros3d,store=.TRUE.)
         a_cNO3b = FloatArray3D(zeros3d,store=.TRUE.)
         a_cOCb = FloatArray3D(zeros3d,store=.TRUE.)
         a_cBCb = FloatArray3D(zeros3d,store=.TRUE.)
         a_cDUb  = FloatArray3D(zeros3d,store=.TRUE.)
         a_cSSb = FloatArray3D(zeros3d,store=.TRUE.)

         a_iSO4a = FloatArray3D(zeros3d,store=.TRUE.)
         a_iNH3a = FloatArray3D(zeros3d,store=.TRUE.)
         a_iNO3a = FloatArray3D(zeros3d,store=.TRUE.)
         a_iOCa = FloatArray3D(zeros3d,store=.TRUE.)
         a_iBCa = FloatArray3D(zeros3d,store=.TRUE.)
         a_iDUa = FloatArray3D(zeros3d,store=.TRUE.)
         a_iSSa = FloatArray3D(zeros3d,store=.TRUE.)

         a_iSO4b = FloatArray3D(zeros3d,store=.TRUE.)
         a_iNH3b = FloatArray3D(zeros3d,store=.TRUE.)
         a_iNO3b = FloatArray3D(zeros3d,store=.TRUE.)
         a_iOCb = FloatArray3D(zeros3d,store=.TRUE.)
         a_iBCb = FloatArray3D(zeros3d,store=.TRUE.)
         a_iDUb = FloatArray3D(zeros3d,store=.TRUE.)
         a_iSSb = FloatArray3D(zeros3d,store=.TRUE.)

      END IF ! level

      !----------------------------------------------------

      ALLOCATE (a_ustar(nxp,nyp),a_tstar(nxp,nyp),a_rstar(nxp,nyp))
      ALLOCATE (uw_sfc(nxp,nyp),vw_sfc(nxp,nyp),ww_sfc(nxp,nyp))
      ALLOCATE (wt_sfc(nxp,nyp),wq_sfc(nxp,nyp))
      IF (level >= 3) THEN
         ALLOCATE(precip(nzp,nxp,nyp))
         precip = 0.
         memsize = memsize + nxyzp
      END IF

      ALLOCATE(snowin(nzp,nxp,nyp))
      memsize = memsize + nxyzp

      a_ustar(:,:) = 0.
      a_tstar(:,:) = 0.
      a_rstar(:,:) = 0.
      uw_sfc(:,:)  = 0.
      vw_sfc(:,:)  = 0.
      ww_sfc(:,:)  = 0.
      wt_sfc(:,:)  = 0.
      wq_sfc(:,:)  = 0.
      umean = 0.
      vmean = 0.

      memsize = memsize +  nxyzp*nscl*2 + 3*nxyp + nxyp*10

      IF(myid == 0) THEN
         PRINT "(//' ',49('-')/,' ',/3x,i3.3,' prognostic scalars')", nscl
         PRINT "('   memory to be allocated  -  ',f8.3,' mbytes')", &
               memsize*1.e-6*kind(0.0)
      END IF

   END SUBROUTINE define_vars
   !
   !----------------------------------------------------------------------
   !
   SUBROUTINE define_grid

      USE mpi_interface, ONLY : xoffset, yoffset, wrxid, wryid, nxpg, nypg,   &
                                appl_abort, myid

      INTEGER :: i, j, k, kmax, nchby
      REAL    :: dzrfm, dz, zb, dzmin
      REAL    :: zmnvc(-1:nzp+1)
      CHARACTER (len=51) :: &
         fm1 = '(//" ",49("-")/,"   grid DIMENSIONs:"/)            ',      &
         fm2 = '("   nxp-4 = ",i3,", dx, dx = ",f8.1,",",f8.1," m")',      &
         fm3 = '("   nyp-4 = ",i3,", dy, dy = ",f8.1,",",f8.1," m")',      &
         fm4 = '("   nzp   = ",i3,", dz, dz = ",f8.1,",",f8.1," m")',      &
         fm5 = '("   timestep: ",f7.3,"s ")                        ',      &
         fm6 = '("   thermo level: ",i3)                        '

      nxyzp = nxp*nyp*nzp
      nxyp  = nxp*nyp

      nz = nzp-1

      dxi = 1./deltax
      dyi = 1./deltay
      ALLOCATE(wsavex(4*nxpg+100),wsavey(4*nypg+100))
      wsavex = 0.0
      wsavey = 0.0

      !
      ! define xm array for grid 1 from deltax
      !
      ALLOCATE (xm(nxp))

      xm(1) = -float(max(nxpg-2,1))*.5*deltax+xoffset(wrxid)*deltax
      DO i = 2, nxp-1
         xm(i) = xm(i-1)+deltax
      END DO
      xm(nxp) = 2*xm(nxp-1)-xm(nxp-2)

      !
      ! define ym array for grid 1 from deltay
      !
      ALLOCATE (ym(nyp))

      ym(1) = -float(max(nypg-2,1))*.5*deltay+yoffset(wryid)*deltay
      DO j = 2, nyp-1
         ym(j) = ym(j-1)+deltay
      END DO
      ym(nyp) = 2*ym(nyp-1)-ym(nyp-2)

        !
        !      define where the momentum points will lie in vertical
        !
      ALLOCATE (zm(nzp))

      SELECT CASE (abs(igrdtyp))
            !
            ! Read in grid spacings from a file
            !
         CASE(3)
            OPEN(1,file='zm_grid_in',status='old',form='formatted')

            DO k = 1, nzp
               READ(1,*) zm(k)
            END DO

            CLOSE(1)

            IF (zm(1) /= 0.) THEN
               IF (myid == 0) PRINT *, 'ABORTING:  Error in input grid'
               CALL appl_abort(0)
            END IF
            !
            ! Tschebyschev Grid with vertical size given by dzmax
            !
         CASE(2)

            zm(1) = 0.
            nchby = nzp-3

            DO k = 1, nzp-1
               zm(k+1) = cos( ((2.*nchby - 1. - 2.*(k-1))*2.*asin(1.))/(2.*nchby))
               zm(k+1) = (zm(k+1)+1.)*dzmax/2.
            END DO

            zm(nzp-1) = dzmax
            zm(nzp)   = dzmax + zm(2)*zm(2)/(zm(3)-zm(2))
            !
            ! define zm array for grid 1 from deltaz and dzrat, if dzrat is
            ! negative compress grid so that dzmin is the grid spacing in a 100m
            ! interval below dzmax.  In both CASEs stretcvh grid uniformly by the
            ! ration |dzrat| above dzmax
            !
         CASE(1)
            zm(1) = 0.
            zm(2) = deltaz
            zb = dzmax+100.

            IF (dzrat < 0.) THEN
               dzmin = -float(int(dzrat))
               dzrat =  dzrat+dzmin-1
               kmax = int(log(deltaz/dzmin)/log(abs(dzrat)))
               zb = dzmax-100.

               DO k = 1, kmax
                  zb = zb-dzmin*abs(dzrat)**k
               END DO

            END IF

            dz = deltaz

            DO k = 3, nzp

               IF(zm(k-1) > zb .AND. zm(k-1) < dzmax)then
                  dz = max(dzmin,dz/abs(dzrat))
               ELSE IF (zm(k-1) >= dzmax) THEN
                  dz = dz*abs(dzrat)
               END IF

               zm(k) = zm(k-1)+dz

            END DO

         CASE DEFAULT
            zm(1) = 0.

            DO k = 2, nzp ! Fixed: used to start from 1
               zm(k) = zm(k-1)+deltaz
            END DO

      END SELECT
      !
      ! Grid Points for Thermal Points (T-Grid):
      !
      ALLOCATE (xt(nxp))

      DO i = 2, nxp
         xt(i) = .5*(xm(i)+xm(i-1))
      END DO

      xt(1) = 1.5*xm(1)-.5*xm(2)
      !
      ALLOCATE (yt(nyp))

      DO j = 2, nyp
         yt(j) = .5*(ym(j)+ym(j-1))
      END DO

      yt(1) = 1.5*ym(1)-.5*ym(2)
      !
      ALLOCATE (zt(nzp))

      IF (igrdtyp < 0) THEN
         !
         ! Read in grid spacings from a file
         !
         OPEN(2,file='zt_grid_in',status='old',form='formatted')

         DO k = 1, nzp
            READ(2,*) zt(k)
         END DO

         CLOSE(2)

      ELSE
         !
         ! calculate where the thermo points will lie based on geometric
         ! interpolation from the momentum points
         !
         DO k = 1, nzp
            zmnvc(k) = zm(k)
         END DO

         zmnvc(0) = -(zmnvc(2)-zmnvc(1))**2 /(zmnvc(3)-zmnvc(2))
         zmnvc(-1) = zmnvc(0)-(zmnvc(1)-zmnvc(0))**2 /(zmnvc(2)-zmnvc(1))
         zmnvc(nzp+1) = zmnvc(nzp)+(zmnvc(nzp)-zmnvc(nzp-1))**2              &
                        /(zmnvc(nzp-1)-zmnvc(nzp-2))

         DO k = 1, nzp
            dzrfm = sqrt(sqrt((zmnvc(k+1)-zmnvc(k)) /(zmnvc(k-1)-zmnvc(k-2))))
            zt(k) = zmnvc(k-1)+(zmnvc(k)-zmnvc(k-1))/(1.+dzrfm)
         END DO

      END IF
      !
      ! compute other arrays based on the vertical grid.
      !   dzm: inverse of distance between thermal points k+1 and k
      !   dzt: inverse of distance between momentum points k and k-1
      !
      ALLOCATE (dzm(nzp))

      DO k = 1, nzp-1
         dzm(k) = 1./(zt(k+1)-zt(k))
      END DO

      dzm(nzp) = dzm(nzp-1)*dzm(nzp-1)/dzm(nzp-2)

      ALLOCATE (dzt(nzp))

      DO k = 2, nzp
         dzt(k) = 1./(zm(k)-zm(k-1))
      END DO

      dzt(1) = dzt(2)*dzt(2)/dzt(3)
      !
      ! set timesteps
      !
      dtl = dtlong
      dtlv = 2.*dtl
      dtlt = dtl
      !
      IF(myid == 0) THEN
         WRITE(6,fm1)
         WRITE(6,fm2) nxpg-4, deltax, 2.*xt(nxp-2)
         WRITE(6,fm3) nypg-4, deltay, 2.*yt(nyp-2)

         WRITE(6,fm4) nzp,zm(2)-zm(1),zm(nzp)
         WRITE(6,fm5) dtl
         WRITE(6,fm6) level
      END IF

   END SUBROUTINE define_grid
   !
   ! ----------------------------------------------------------------------
   ! Subroutine init_anal:  Defines the netcdf Analysis file
   !
   ! Modified for level 4.
   ! Juha Tonttila, FMI, 2014
   !
   !
   SUBROUTINE init_anal(time,salsa_b_bins)
      USE classFieldArray
      USE mo_structured_datatypes
      USE mpi_interface, ONLY : myid, ver, author, info
      USE mo_submctl, ONLY : fn2a, fn2b, fca, fcb, fra, &
                             fia, fib, fsa
      USE class_ComponentIndex, ONLY : IsUsed

      REAL, INTENT (in)    :: time

      LOGICAL, INTENT (in) :: salsa_b_bins
      INTEGER              :: nbeg, nend, i, nnum
      CHARACTER (len=10), ALLOCATABLE, SAVE :: sbase(:)
      LOGICAL, ALLOCATABLE, SAVE :: salsabool(:)

      ! dummy arrays to help allocate FloatArrays to zero -AZ
      REAL :: zeros1d(1)
      REAL :: zeros2d(2,2)
      REAL :: zeros3d(3,3,3)
      REAL :: zeros4d(4,4,4,4)

      LOGICAL :: l1, l2, l3, l4, l5           ! booleans for levels
      LOGICAL :: lso4, lnh3, lno3, lss, ldu, loc, lbc ! booleans for particles
      LOGICAL :: lsb !boolean for b bins

      ! dummy float arrays required to initialize field arrays -AZ
      TYPE(FloatArray1d), TARGET :: dumarr1d
      TYPE(FloatArray2d), TARGET :: dumarr2d
      TYPE(FloatArray3d), TARGET :: dumarr3d
      TYPE(FloatArray4d), TARGET :: dumarr4d
      TYPE(ArrayElement), POINTER :: ArrEl

      ! Set logical variables for levels -AZ
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      l4 = .FALSE.
      l5 = .FALSE.

      IF (level >= 1) l1 = .TRUE.
      IF (level >= 2) l2 = .TRUE.
      IF (level >= 3) l3 = .TRUE.
      IF (level >= 4) l4 = .TRUE.
      IF (level >= 5) l5 = .TRUE.
      lsb = salsa_b_bins

      ! Set logical variables for aerosol species if SALSA is used
      IF (level >= 4) THEN
         lso4 = IsUsed(prtcl,'SO4')
         loc = IsUsed(prtcl,'OC')
         lbc = IsUsed(prtcl,'BC')
         ldu = IsUsed(prtcl,'DU')
         lss = IsUsed(prtcl,'SS')
         lno3 = IsUsed(prtcl,'NO')
         lnh3 = IsUsed(prtcl,'NH')
      END IF

      ! Create a field array for output variables
      fields = FieldArray()

      ! Initialize dummy arrays and dummy pointer
      zeros1d(:) = 0.
      zeros2d(:,:) = 0.
      zeros3d(:,:,:) = 0.
      zeros4d(:,:,:,:) = 0.
      dumarr3d = FloatArray3D(zeros3d,store=.TRUE.)
      dumarr4d = FloatArray4D(zeros4d,store=.TRUE.)
      dummy => dumarr3d
      dummy2 => dumarr4d

      ! Create fields for output variables
      CALL fields%NewField("time","Time","s","time",.TRUE.,dummy)
      CALL fields%NewField("zt","Vertical displacement of cell centers","m","zt",.TRUE.,dummy)
      CALL fields%NewField("zm","Vertical displacement of cell edges","m","zm",.TRUE.,dummy)
      CALL fields%NewField("xt","East-west displacement of cell centers","m","xt",.TRUE.,dummy)
      CALL fields%NewField("xm","East-west displacement of cell edges","m","xm",.TRUE.,dummy)
      CALL fields%NewField("yt","North-south displacement of cell centers","m","yt",.TRUE.,dummy)
      CALL fields%NewField("ym","North-south displacement of cell edges","m","ym",.TRUE.,dummy)

      ! Size bins
      IF (level >= 4) THEN

         CALL fields%NewField("aea","Aerosol size bins, regime a","m","aea",lbinanl,dummy)
         CALL fields%NewField("aeb","Aerosol size bins, regime b","m","aeb",(lbinanl .AND. lsb),dummy)
         CALL fields%NewField("cla","Cloud droplet size bins, regime a","m","cla",lbinanl,dummy)
         CALL fields%NewField("clb","Cloud droplet size bins, regime b","m","clb",(lbinanl .AND. lsb),dummy)
         CALL fields%NewField("prc","Precipitation size bins","m","prc",lbinanl,dummy)

         CALL fields%NewField("ica","Ice cloud droplet size bins, regime a","m","ica",(lbinanl .AND. l5),dummy)
         CALL fields%NewField("icb","Ice cloud droplet size bins","m","icb",(lbinanl .AND. l5 .AND. lsb),dummy)
         CALL fields%NewField("snw","Snow size bins","m","snw",(lbinanl .AND. l5),dummy)

      END IF

      CALL fields%NewField("u0","Geostrophic zonal wind","m/s","zt",.TRUE.,dummy)
      CALL fields%NewField("v0","Geostrophic meridional wind","m/s","zt",.TRUE.,dummy)
      CALL fields%NewField("dn0","Base-state density","kg/m^3","zt",.TRUE.,dummy)

      !Point the dummy pointers to float arrays containing the numerical values, then create the fields
      dummy => a_u
      CALL fields%NewField("u","'Zonal wind'","m/s","mttt",.TRUE.,dummy)
      dummy => a_v
      CALL fields%NewField("v","Meridional wind","m/s","tmtt",.TRUE.,dummy)
      dummy => a_w
      CALL fields%NewField("w","Vertical velocity","m/s","ttmt",.TRUE.,dummy)
      dummy => a_theta
      CALL fields%NewField("theta","Potential temperature","K","tttt",.TRUE.,dummy)
      dummy => a_press
      CALL fields%NewField("p","Pressure","Pa","tttt",.TRUE.,dummy)

      dummy => a_l
      CALL fields%NewField("l","Liquid water mixing ratio","kg/kg","tttt",l2,dummy)
      dummy => a_q
      CALL fields%NewField("q","Total water mixing ratio","kg/kg","tttt",l1,dummy)
      dummy => a_r
      CALL fields%NewField("r","Rain-water mixing ratio","kg/kg","tttt",l3,dummy)


      IF (level >= 5) THEN

         dummy => a_f
         CALL fields%NewField("f","Total ice mixing ratio","kg/kg","tttt",l5,dummy)
         dummy => a_i
         CALL fields%NewField("i","Ice mixing ratio","kg/kg","tttt",l5,dummy)
         dummy => a_s
         CALL fields%NewField("s","Snow mixing ratio","kg/kg","tttt",l5,dummy)

      END IF


      IF (level < 4) THEN

         dummy => a_n
         IF (level == 3) THEN
            CALL fields%NewField("n","Rain-drop number mixing ratio","#/kg","tttt",.TRUE.,dummy)
         ELSE
            dummy => a_n
            CALL fields%NewField("n","Rain-drop number mixing ratio","#/kg","tttt",.FALSE.,dummy)
         END IF

         dummy => a_stke
         IF (isgstyp>1) THEN
            CALL fields%NewField("stke","Sub-filter scale TKE","J/kg","mttt",.TRUE.,dummy)
         ELSE
            CALL fields%NewField("stke","Sub-filter scale TKE","J/kg","mttt",.FALSE.,dummy)
         END IF

         dummy => a_rfl
         IF (iradtyp>2) THEN
            CALL fields%NewField("rflx","Total Radiative flux","W/m^2","ttmt",.TRUE.,dummy)
         ELSE
            CALL fields%NewField("rflx","Total Radiative flux","W/m^2","ttmt",.FALSE.,dummy)
         END IF

      END IF


      ! Create salsa exclusive output variables
      IF (level >= 4) THEN


         dummy => a_rh
         CALL fields%NewField("S_RH","SALSA Relative humidity","1","tttt",.TRUE.,dummy) !
         dummy => a_rhi
         CALL fields%NewField("S_RHI","SALSA Relative humidity over ice","1","tttt",l5,dummy) !
         dummy => a_nac
         CALL fields%NewField("S_Nact","SALSA Number of newly activated droplets","kg^-1","tttt",.TRUE.,dummy) !

         ! AEROSOL
         dummy => a_na
         CALL fields%NewField("S_Na","SALSA total number of soluble aerosols, (regime A)","kg^-1","tttt",.TRUE.,dummy) !
         dummy2 => a_naba
         CALL fields%NewField("S_Naba","Aerosol size distribution, regime A","kg^-1","ttttaea",lbinanl,dummy2) !
         dummy => a_rwaa
         CALL fields%NewField("S_Rwaa","SALSA number mean wet radius of aerosols, regime A","m","tttt",.TRUE.,dummy) !
         dummy2 => a_rwaba
         CALL fields%NewField("S_Rwaba","SALSA bin aerosol wet radius, regime A","m","ttttaea",lbinanl,dummy2) !
         dummy => a_nb
         CALL fields%NewField("S_Nb","SALSA total number of insoluble aerosols, (regime B)","kg^-1","tttt",lsb,dummy) !
         dummy2 => a_nabb
         CALL fields%NewField("S_Nabb","Aerosol size distribution, regime B","kg^-1","ttttaeb",(lbinanl .AND. lsb),dummy2) !
         dummy => a_rwab
         CALL fields%NewField("S_Rwab","SALSA number mean wet radius of aerosols, regime B","m","tttt",lsb,dummy) !
         dummy2 => a_rwabb
         CALL fields%NewField("S_Rwabb","SALSA bin aerosol wet radius, regime B","m","ttttaeb",(lbinanl .AND. lsb),dummy2) !

         ! CLOUD
         dummy => a_nc
         CALL fields%NewField("S_Nc","SALSA cdnc","kg^-1","tttt",.TRUE.,dummy) !
         dummy2 => a_ncba
         CALL fields%NewField("S_Ncba","SALSA cloud droplet size distribution, regime a","kg^-1","ttttcla",lbinanl,dummy2) !
         dummy2 => a_ncbb
         CALL fields%NewField("S_Ncbb","SALSA cloud droplet size distribution, regime b","kg^-1","ttttclb",&
                              (lbinanl .AND. lsb),dummy2) !
         dummy => a_rwca
         CALL fields%NewField("S_Rwca","SALSA number mean radius of cloud droplets, regime a","m","tttt",.TRUE.,dummy) !
         dummy => a_rwcb
         CALL fields%NewField("S_Rwcb","SALSA number mean radius of cloud droplets, regime b","m","tttt",lsb,dummy) !
         dummy2 => a_rwcba
         CALL fields%NewField("S_Rwcba","SALSA bin cloud droplet radius, regime a","m","ttttcla",lbinanl,dummy2) !
         dummy2 => a_rwcbb
         CALL fields%NewField("S_Rwcbb","SALSA bin cloud droplet radius, regime b","m","ttttclb",(lbinanl .AND. lsb),dummy2) !

         ! PRECIPITATION
         dummy => a_np
         CALL fields%NewField("S_Np","SALSA rdnc","kg^-1","tttt",.TRUE.,dummy) !
         dummy2 => a_npba
         CALL fields%NewField("S_Npba","SALSA precipitation size distribution","kg^-1","ttttprc",lbinanl,dummy2) !
         dummy => a_rwpa
         CALL fields%NewField("S_Rwpa","SALSA number mean radius of precipitation particles","m","tttt",.TRUE.,dummy) !
         dummy2 => a_rwpba
         CALL fields%NewField("S_Rwpba","SALSA bin precipitation particle radius","m","ttttprc",lbinanl,dummy2) !

         ! ICE
         dummy => a_ni
         CALL fields%NewField("S_Ni","SALSA ice nuclei","m^-3","tttt",l5,dummy) ! !
         dummy2 => a_niba
         CALL fields%NewField("S_Niba","SALSA ice particle size distribution, regime a","m^-3","ttttica",&
                              (lbinanl .AND. l5),dummy2) !
         dummy2 => a_nibb
         CALL fields%NewField("S_Nibb","SALSA ice particle size distribution, regime b","m^-3","tttticb",&
                              (lbinanl .AND. l5 .AND. lsb),dummy2)
         dummy => a_rwia
         CALL fields%NewField("S_Rwia","SALSA number mean radius of ice particles, regime a","m","tttt",l5,dummy) !
         dummy => a_rwib
         CALL fields%NewField("S_Rwib","SALSA number mean radius of ice particles, regime b","m","tttt",(l5 .AND. lsb),dummy) !
         dummy2 => a_rwiba
         CALL fields%NewField("S_Rwiba","SALSA bin ice particle radius, regime a","m","ttttica",(lbinanl .AND. l5),dummy2) !
         dummy2 => a_rwibb
         CALL fields%NewField("S_Rwibb","SALSA bin ice particle radius, regime b","m","tttticb",&
                              (lbinanl .AND. l5 .AND. lsb),dummy2) !

         ! SNOW
         dummy => a_ns
         CALL fields%NewField("S_Ns","SALSA sdnc","m^-3","tttt",l5,dummy) !
         dummy2 => a_nsba
         CALL fields%NewField("S_Nsba","SALSA snow size distribution","m^-3","ttttsnw",(lbinanl .AND. l5),dummy2) !
         dummy => a_rwsa
         CALL fields%NewField("S_Rwsa","SALSA number mean radius of snow particles","m","tttt",l5,dummy) !
         dummy2 => a_rwsba
         CALL fields%NewField("S_Rwsba","SALSA bin snow particle radius","m","ttttsnw",(lbinanl .AND. l5),dummy2) !


         ! AEROSOL SPECIES

         ! aerosols, regime A
         dummy => a_aso4a
         CALL fields%NewField("S_aSO4a","SALSA aerosol mass concentration of SO4, regime A","kg/kg","tttt",lso4,dummy) !
         dummy => a_anh3a
         CALL fields%NewField("S_aNH3a","SALSA aerosol mass concentration of NH3, regime A","kg/kg","tttt",lnh3,dummy) !
         dummy => a_ano3a
         CALL fields%NewField("S_aNO3a","SALSA aerosol mass concentration of NO3, regime A","kg/kg","tttt",lno3,dummy) !
         dummy => a_aoca
         CALL fields%NewField("S_aOCa","SALSA aerosol mass concentration of OC, regime A","kg/kg","tttt",loc,dummy) !
         dummy => a_abca
         CALL fields%NewField("S_aBCa","SALSA aerosol mass concentration of BC, regime A","kg/kg","tttt",lbc,dummy) !
         dummy => a_adua
         CALL fields%NewField("S_aDUa","SALSA aerosol mass concentration of DU, regime A","kg/kg","tttt",ldu,dummy) !
         dummy => a_assa
         CALL fields%NewField("S_aSSa","SALSA aerosol mass concentration of SS, regime A","kg/kg","tttt",lss,dummy) !

         ! aerosols, regime B
         dummy => a_aso4b
         CALL fields%NewField("S_aSO4b","SALSA aerosol mass concentration of SO4, regime B","kg/kg","tttt",(lso4 .AND. lsb),dummy) !
         dummy => a_anh3b
         CALL fields%NewField("S_aNH3b","SALSA aerosol mass concentration of NH3, regime B","kg/kg","tttt",(lnh3 .AND. lsb),dummy) !
         dummy => a_ano3b
         CALL fields%NewField("S_aNO3b","SALSA aerosol mass concentration of NO3, regime B","kg/kg","tttt",(lno3 .AND. lsb),dummy) !
         dummy => a_aocb
         CALL fields%NewField("S_aOCb","SALSA aerosol mass concentration of OC, regime B","kg/kg","tttt",(loc .AND. lsb),dummy) !
         dummy => a_abcb
         CALL fields%NewField("S_aBCb","SALSA aerosol mass concentration of BC, regime B","kg/kg","tttt",(lbc .AND. lsb),dummy) !
         dummy => a_adub
         CALL fields%NewField("S_aDUb","SALSA aerosol mass concentration of DU, regime B","kg/kg","tttt",(ldu .AND. lsb),dummy) !
         dummy => a_assb
         CALL fields%NewField("S_aSSb","SALSA aerosol mass concentration of SS, regime B","kg/kg","tttt",(lss .AND. lsb),dummy) !

         ! cloud aerosols, regime A
         dummy => a_cso4a
         CALL fields%NewField("S_cSO4a","SALSA CCN mass concentration of SO4, regime A","kg/kg","tttt",lso4,dummy)
         dummy => a_cnh3a
         CALL fields%NewField("S_cNH3a","SALSA CCN mass concentration of NH3, regime A","kg/kg","tttt",lnh3,dummy)
         dummy => a_cno3a
         CALL fields%NewField("S_cNO3a","SALSA CCN mass concentration of NO3, regime A","kg/kg","tttt",lno3,dummy)
         dummy => a_coca
         CALL fields%NewField("S_cOCa","SALSA CCN mass concentration of OC, regime A","kg/kg","tttt",loc,dummy)
         dummy => a_cbca
         CALL fields%NewField("S_cBCa","SALSA CCN mass concentration of BC, regime A","kg/kg","tttt",lbc,dummy)
         dummy => a_cdua
         CALL fields%NewField("S_cDUa","SALSA CCN mass concentration of DU, regime A","kg/kg","tttt",ldu,dummy)
         dummy => a_cssa
         CALL fields%NewField("S_cSSa","SALSA CCN mass concentration of SS, regime A","kg/kg","tttt",lss,dummy)

         ! cloud aerosols, regime B
         dummy => a_cso4b
         CALL fields%NewField("S_cSO4b","SALSA CCN mass concentration of SO4, regime B","kg/kg","tttt",(lso4 .AND. lsb),dummy)
         dummy => a_cnh3b
         CALL fields%NewField("S_cNH3b","SALSA CCN mass concentration of NH3, regime B","kg/kg","tttt",(lnh3 .AND. lsb),dummy)
         dummy => a_cno3b
         CALL fields%NewField("S_cNO3b","SALSA CCN mass concentration of NO3, regime B","kg/kg","tttt",(lno3 .AND. lsb),dummy)
         dummy => a_cocb
         CALL fields%NewField("S_cOCb","SALSA CCN mass concentration of OC, regime B","kg/kg","tttt",(loc .AND. lsb),dummy)
         dummy => a_cbcb
         CALL fields%NewField("S_cBCb","SALSA CCN mass concentration of BC, regime B","kg/kg","tttt",(lbc .AND. lsb),dummy)
         dummy => a_cdub
         CALL fields%NewField("S_cDUb","SALSA CCN mass concentration of DU, regime B","kg/kg","tttt",(ldu .AND. lsb),dummy)
         dummy => a_cssb
         CALL fields%NewField("S_cSSb","SALSA CCN mass concentration of SS, regime B","kg/kg","tttt",(lss .AND. lsb),dummy)

         ! ice aerosols, regime A
         dummy => a_iso4a
         CALL fields%NewField("S_iSO4a","SALSA IN mass concentration of SO4, regime A","kg/kg","tttt",(lso4 .AND. l5),dummy)
         dummy => a_inh3a
         CALL fields%NewField("S_iNH3a","SALSA IN mass concentration of NH3, regime A","kg/kg","tttt",(lnh3 .AND. l5),dummy)
         dummy => a_ino3a
         CALL fields%NewField("S_iNO3a","SALSA IN mass concentration of NO3, regime A","kg/kg","tttt",(lno3 .AND. l5),dummy)
         dummy => a_ioca
         CALL fields%NewField("S_iOCa","SALSA IN mass concentration of OC, regime A","kg/kg","tttt",(loc .AND. l5),dummy)
         dummy => a_ibca
         CALL fields%NewField("S_iBCa","SALSA IN mass concentration of BC, regime A","kg/kg","tttt",(lbc .AND. l5),dummy)
         dummy => a_idua
         CALL fields%NewField("S_iDUa","SALSA IN mass concentration of DU, regime A","kg/kg","tttt",(ldu .AND. l5),dummy)
         dummy => a_issa
         CALL fields%NewField("S_iSSa","SALSA IN mass concentration of SS, regime A","kg/kg","tttt",(lss .AND. l5),dummy)

         ! ice aerosols, regime B
         dummy => a_iso4b
         CALL fields%NewField("S_iSO4b","SALSA IN mass concentration of SO4, regime B","kg/kg","tttt",&
                              (lso4 .AND. l5 .AND. lsb),dummy)
         dummy => a_inh3b
         CALL fields%NewField("S_iNH3b","SALSA IN mass concentration of NH3, regime B","kg/kg","tttt",&
                              (lnh3 .AND. l5 .AND. lsb),dummy)
         dummy => a_ino3b
         CALL fields%NewField("S_iNO3b","SALSA IN mass concentration of NO3, regime B","kg/kg","tttt",&
                              (lno3 .AND. l5 .AND. lsb),dummy)
         dummy => a_iocb
         CALL fields%NewField("S_iOCb","SALSA IN mass concentration of OC, regime B","kg/kg","tttt",&
                              (loc .AND. l5 .AND. lsb),dummy)
         dummy => a_ibcb
         CALL fields%NewField("S_iBCb","SALSA IN mass concentration of BC, regime B","kg/kg","tttt",&
                              (lbc .AND. l5 .AND. lsb),dummy)
         dummy => a_idub
         CALL fields%NewField("S_iDUb","SALSA IN mass concentration of DU, regime B","kg/kg","tttt",&
                              (ldu .AND. l5 .AND. lsb),dummy)
         dummy => a_issb
         CALL fields%NewField("S_iSSb","SALSA IN mass concentration of SS, regime B","kg/kg","tttt",&
                              (lss .AND. l5 .AND. lsb),dummy)

      END IF

      !Number of total variables
      varnum = fields%count

      fname =  trim(filprf)
      IF(myid == 0) PRINT                                                  &
         "(//' ',49('-')/,' ',/,'   Initializing: ',A20)",trim(fname)

      CALL open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid0, nrec0, ver, author, info)

      IF (level < 4 .OR. .NOT. lbinanl) THEN

         CALL define_nc( ncid0, nrec0, varnum, fields, n1=nzp, n2=nxp-4, n3=nyp-4)

      ELSE IF (level == 4 .AND. lbinanl) THEN

         CALL define_nc( ncid0, nrec0, varnum, fields, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                         inae_a=fn2a, inae_b=fn2b-fn2a, incld_a=fca%cur,          &
                         incld_b=fcb%cur-fca%cur, inprc=fra )

      ELSE IF (level == 5 .AND. lbinanl) THEN

         CALL define_nc( ncid0, nrec0, varnum, fields, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                         inae_a=fn2a,  inae_b =fn2b-fn2a, incld_a=fca%cur,        &
                         incld_b=fcb%cur-fca%cur, inprc=fra, inice_a=fia%cur,     &
                         inice_b=fib%cur-fia%cur, insnw=fsa )

      END IF

      IF (myid == 0) PRINT *,'   ...starting record: ', nrec0

   END SUBROUTINE init_anal

   !Subroutine for calculating output values
   SUBROUTINE calc_anal()
      USE class_ComponentIndex, ONLY : IsUsed
      USE mo_submctl, ONLY : in1a, fn2a, in2b, fn2b,            &
                             ica, fca, icb, fcb, ira, fra,        &
                             iia, fia, iib, fib, isa, fsa,        &
                             aerobins, cloudbins, precpbins, &
                             icebins, snowbins, nlim, prlim, &
                             nspec, nbins, ncld, nice, nprc, nsnw

      INTEGER :: n, i
      INTEGER :: ibeg(4), icnt(4), i1, i2, j1, j2
      INTEGER :: ibegsd(5),  icntaea(5), icntaeb(5), icntcla(5), icntclb(5), icntpra(5), & ! Juha: For sizedistribution variables
                 icntica(5), icnticb(5), icntsna(5)
      REAL :: zsum(nzp,nxp,nyp) ! Juha: Helper for computing bulk output diagnostics
      REAL :: zvar(nzp,nxp,nyp)
      REAL :: a_Rawet(nzp,nxp,nyp,nbins), a_Rcwet(nzp,nxp,nyp,ncld), a_Rpwet(nzp,nxp,nyp,nprc), &
              a_Riwet(nzp,nxp,nyp,nice),  a_Rswet(nzp,nxp,nyp,nsnw)
      TYPE(ArrayElement), POINTER :: ArrEl

      icnt =    (/nzp, nxp-4, nyp-4, 1/)
      icntaea = (/nzp,nxp-4,nyp-4, fn2a, 1 /)
      icntaeb = (/nzp,nxp-4,nyp-4, fn2b-fn2a, 1/)
      icntcla = (/nzp,nxp-4,nyp-4, fca%cur, 1/)
      icntclb = (/nzp,nxp-4,nyp-4, fcb%cur-fca%cur, 1/)
      icntpra = (/nzp,nxp-4,nyp-4, fra, 1/)
      icntica = (/nzp,nxp-4,nyp-4, fia%cur, 1/)
      icnticb = (/nzp,nxp-4,nyp-4, fib%cur-fia%cur, 1/)
      icntsna = (/nzp,nxp-4,nyp-4, fsa, 1/)
      ibeg =   (/1,1,1,nrec0/)
      ibegsd = (/1,1,1,1,nrec0/)

      i1 = 3
      i2 = nxp-2
      j1 = 3
      j2 = nyp-2

      a_u%data(:,:,:) = a_up(:,:,:)
      a_v%data(:,:,:) = a_vp(:,:,:)
      a_w%data(:,:,:) = a_wp(:,:,:)
      a_p%data(:,:,:) = a_press%data(:,:,:)
      a_l%data(:,:,:) = a_rc%data(:,:,:)

      IF(level < 4) THEN

         n = 2
         CALL newsclr(n)
         a_q%data(:,:,:) = a_sp(:,:,:)
         n=n+1

         CALL newsclr(n)
         a_r%data = a_sp(:,:,:)
         n=n+1

         CALL newsclr(n)
         a_n%data(:,:,:) = a_sp(:,:,:)
         n=n+1

         CALL newsclr(n)
         a_stke%data(:,:,:) = a_sp(:,:,:)
         n=n+1

         a_rfl%data(:,:,:) = a_rflx(:,:,:)

      END IF


      IF(level >= 4) THEN
         a_q%data(:,:,:) = a_rp%data(:,:,:) + a_rc%data(:,:,:) + a_srp%data(:,:,:)+ &   ! Rain water
                           a_ri%data(:,:,:) + & ! Ice water (level 5)
                           a_srs%data(:,:,:)      ! Snow water (level 5)
         a_r%data(:,:,:) = a_srp%data(:,:,:)

         !Number of newly activated droplets
         a_nac%data(:,:,:) = SUM(a_nactd(:,:,:,:),DIM=4)

         !Calculate cloud droplet number concentration
         CALL bulkNumc('cloud','a',zvar(:,:,:))
         zsum(:,:,:) = zvar(:,:,:)
         CALL bulkNumc('cloud','b',zvar(:,:,:))
         zsum(:,:,:) = zsum(:,:,:) + zvar(:,:,:)
         a_nc%data(:,:,:) = zsum(:,:,:)

         !Cloud droplet radius, regime A
         CALL meanRadius('cloud','a',zvar(:,:,:))
         a_rwca%data(:,:,:) = zvar(:,:,:)

         !Cloud droplet radius, regime B
         CALL meanRadius('cloud','b',zvar(:,:,:))
         a_rwcb%data(:,:,:) = zvar(:,:,:)

         IF (lbinanl) THEN

            !SALSA bin cloud droplet radius for both regimes
            CALL getBinRadius(ncld,nspec+1,a_ncloudp%data(:,:,:,:),a_mcloudp%data(:,:,:,:),nlim,a_Rcwet,2)
            a_rwcba%data(:,:,:,:) = a_Rcwet(:,:,:,:)
            a_rwcbb%data(:,:,:,:) = a_Rcwet(:,:,:,:)

            !Cloud droplet size distribution for both regimes
            a_ncba%data(:,:,:,:) = a_ncloudp%data(:,:,:,:)
            a_ncbb%data(:,:,:,:) = a_ncloudp%data(:,:,:,:)

            !Precipitation size distribution
            a_npba%data(:,:,:,:) = a_nprecpp%data(:,:,:,:)

            !Bin precipitation particle radius
            CALL getBinRadius(nprc,nspec+1,a_nprecpp%data(:,:,:,:),a_mprecpp%data(:,:,:,:),prlim,a_Rpwet,3)
            a_rwpba%data(:,:,:,:) = a_Rpwet(:,:,:,:)

            !Aerosol size distribution, both regimes
            a_naba%data(:,:,:,:) = a_naerop%data(:,:,:,:)
            a_nabb%data(:,:,:,:) = a_naerop%data(:,:,:,:)

            !Bin aerosol wet radius, both regimes
            CALL getBinRadius(nbins,nspec+1,a_naerop%data(:,:,:,:),a_maerop%data(:,:,:,:),nlim,a_Rawet,1)
            a_rwaba%data(:,:,:,:) = a_Rawet(:,:,:,:)
            a_rwabb%data(:,:,:,:) = a_Rawet(:,:,:,:)

         END IF

         !Rain droplet number concentration
         CALL bulkNumc('precp','a',zvar(:,:,:))
         a_np%data(:,:,:) = zvar(:,:,:)

         !Number mean radius of precipitation particles
         CALL meanRadius('precp','a',zvar(:,:,:))
         a_rwpa%data(:,:,:) = zvar(:,:,:)

         !Total number of soluble aerosols, regime A
         CALL bulkNumc('aerosol','a',zvar(:,:,:))
         a_na%data(:,:,:) = zvar(:,:,:)

         !Total number of soluble aerosols, regime B
         CALL bulkNumc('aerosol','b',zvar(:,:,:))
         a_nb%data(:,:,:) = zvar(:,:,:)

         !Number mean wet radius of aerosols, regime A
         CALL meanRadius('aerosol','a',zvar(:,:,:))
         a_rwaa%data(:,:,:) = zvar(:,:,:)

         !Number mean wet radius of aerosols, regime A
         CALL meanRadius('aerosol','b',zvar(:,:,:))
         a_rwab%data(:,:,:) = zvar(:,:,:)

         IF (level == 5) THEN

            !a_f%data(:,:,:) = a_ri%data(:,:,:) + a_srs%data(:,:,:)
            a_i%data(:,:,:) = a_ri%data(:,:,:)
            a_s%data(:,:,:) = a_srs%data(:,:,:)

            !Number of ice nuclei
            CALL bulkNumc('ice','a',zvar(:,:,:))
            zsum(:,:,:) = zvar(:,:,:)
            CALL bulkNumc('ice','b',zvar(:,:,:))
            zsum(:,:,:) = zsum(:,:,:) + zvar(:,:,:)
            a_ni%data(:,:,:) = zsum(:,:,:)

            !Snow density number concentration
            CALL bulkNumc('snow','a',zvar(:,:,:))
            a_ns%data(:,:,:) = zvar(:,:,:)

            !Number mean radius of ice particles, regime A
            CALL meanRadius('ice','a',zvar(:,:,:))
            a_rwia%data(:,:,:) = zvar(:,:,:)

            !Number mean radius of ice particles, regime B
            CALL meanRadius('ice','b',zvar(:,:,:))
            a_rwib%data(:,:,:) = zvar(:,:,:)

            !Number mean radius of snow particles
            CALL meanRadius('snow','a',zvar(:,:,:))
            a_rwsa%data(:,:,:) = zvar(:,:,:)

            IF(lbinanl) THEN

               !Ice particle size distribution, regime A
               a_niba%data(:,:,:,:) = a_nicep%data(:,:,:,:)
               !Ice particle size distribution, regime B
               a_nibb%data(:,:,:,:) = a_nicep%data(:,:,:,:)

               CALL getBinRadius(nice,nspec+1,a_nicep%data,a_micep%data,prlim,a_Riwet,4)
               !Bin ice particle radius, regime A
               a_rwiba%data(:,:,:,:) = a_Riwet(:,:,:,:)
               !Bin ice particle radius, regime B
               a_rwibb%data(:,:,:,:) = a_Riwet(:,:,:,:)

               !Bin snow particle size distribution
               a_nsba%data(:,:,:,:) = a_nsnowp%data(:,:,:,:)

               !Bin snow particle radius
               CALL getBinRadius(nsnw,nspec+1,a_nsnowp%data,a_msnowp%data,prlim,a_Rswet,5)
               a_rwsba%data(:,:,:,:) = a_Rswet(:,:,:,:)

            END IF
         END IF

         ! AEROSOL OUTPUT

         IF (IsUsed(prtcl,'SO4')) THEN

            ! --Sulphate (aerosol, regime A)
            CALL bulkMixrat('SO4','aerosol','a',zvar(:,:,:))
            a_aSO4a%data(:,:,:) = zvar(:,:,:)

            ! --Sulphate (aerosol, regime B)
            CALL bulkMixrat('SO4','aerosol','b',zvar(:,:,:))
            a_aSO4b%data(:,:,:) = zvar(:,:,:)

            ! --Sulphate (clouds, regime A)
            CALL bulkMixrat('SO4','cloud','a',zvar(:,:,:))
            a_cSO4a%data(:,:,:) = zvar(:,:,:)

            ! --Sulphate (clouds, regime B)
            CALL bulkMixrat('SO4','cloud','b',zvar(:,:,:))
            a_cSO4b%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN

               ! --Sulphate (ice, regime A)
               CALL bulkMixrat('SO4','ice','a',zvar(:,:,:))
               a_iSO4a%data(:,:,:) = zvar(:,:,:)

               ! --Sulphate (ice, regime B)
               CALL bulkMixrat('SO4','ice','b',zvar(:,:,:))
               a_iSO4b%data(:,:,:) = zvar(:,:,:)
            END IF ! level 5

         END IF

         IF (IsUSed(prtcl,'NH')) THEN

            !-- Ammonium (aerosol, regime A)
            CALL bulkMixrat('NH','aerosol','a',zvar(:,:,:))
            a_aNH3a%data(:,:,:) = zvar(:,:,:)

            !-- Ammonium (aerosol, regime B)
            CALL bulkMixrat('NH','aerosol','b',zvar(:,:,:))
            a_aNH3b%data(:,:,:) = zvar(:,:,:)

            !-- Ammonium (clouds, regime A)
            CALL bulkMixrat('NH','cloud','a',zvar(:,:,:))
            a_cNH3a%data(:,:,:) = zvar(:,:,:)

            !-- Ammonium (clouds, regime B)
            CALL bulkMixrat('NH','cloud','b',zvar(:,:,:))
            a_cNH3b%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN
               ! --Ammonium (ice, regime A)
               CALL bulkMixrat('NH','ice','a',zvar(:,:,:))
               a_iNH3a%data(:,:,:) = zvar(:,:,:)

               ! --Ammonium (ice, regime B)
               CALL bulkMixrat('NH','ice','b',zvar(:,:,:))
               a_iNH3b%data(:,:,:) = zvar(:,:,:)
            END IF ! level 5

         END IF

         IF (IsUsed(prtcl,'NO')) THEN

            !-- Nitrate (aerosol, regime A)
            CALL bulkMixrat('NO','aerosol','a',zvar(:,:,:))
            a_aNO3a%data(:,:,:) = zvar(:,:,:)

            !-- Nitrate (aerosol, regime B)
            CALL bulkMixrat('NO','aerosol','b',zvar(:,:,:))
            a_aNO3b%data(:,:,:) = zvar(:,:,:)

            !-- Nitrate (clouds, regime A)
            CALL bulkMixrat('NO','cloud','a',zvar(:,:,:))
            a_cNO3a%data(:,:,:) = zvar(:,:,:)

            !-- Nitrate (clouds, regime B)
            CALL bulkMixrat('NO','cloud','b',zvar(:,:,:))
            a_cNO3b%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN
               ! --Nitrate (ice, regime A)
               CALL bulkMixrat('NO','ice','a',zvar(:,:,:))
               a_iNO3a%data(:,:,:) = zvar(:,:,:)

               ! --Nitrate (ice, regime B)
               CALL bulkMixrat('NO','ice','b',zvar(:,:,:))
               a_iNO3b%data(:,:,:) = zvar(:,:,:)
            END IF ! level 5

         END IF

         IF (IsUsed(prtcl,'OC')) THEN

            !-- Organic Carbon (aerosol, regime A)
            CALL bulkMixrat('OC','aerosol','a',zvar(:,:,:))
            a_aOCa%data(:,:,:) = zvar(:,:,:)

            !-- Organic Carbon (aerosol, regime B)
            CALL bulkMixrat('OC','aerosol','b',zvar(:,:,:))
            a_aOCb%data(:,:,:) = zvar(:,:,:)

            !-- Organic Carbon (clouds, regime A)
            CALL bulkMixrat('OC','cloud','a',zvar(:,:,:))
            a_cOCa%data(:,:,:) = zvar(:,:,:)

            !-- Organic Carbon (clouds, regime B)
            CALL bulkMixrat('OC','cloud','b',zvar(:,:,:))
            a_cOCb%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN
               ! --Organic Carbon (ice, regime A)
               CALL bulkMixrat('OC','ice','a',zvar(:,:,:))
               a_iOCa%data(:,:,:) = zvar(:,:,:)

               ! --Organic Carbon (ice, regime B)
               CALL bulkMixrat('OC','ice','b',zvar(:,:,:))
               a_iOCb%data(:,:,:) = zvar(:,:,:)

            END IF ! level 5

         END IF

         IF (IsUsed(prtcl,'BC')) THEN

            !-- Black Carbon (aerosol, regime A)
            CALL bulkMixrat('BC','aerosol','a',zvar(:,:,:))
            a_aBCa%data(:,:,:) = zvar(:,:,:)

            !-- Black Carbon (aerosol, regime B)
            CALL bulkMixrat('BC','aerosol','b',zvar(:,:,:))
            a_aBCb%data(:,:,:) = zvar(:,:,:)

            !-- Black Carbon (clouds, regime A)
            CALL bulkMixrat('BC','cloud','a',zvar(:,:,:))
            a_cBCa%data(:,:,:) = zvar(:,:,:)

            !-- Black Carbon (clouds, regime B)
            CALL bulkMixrat('BC','cloud','b',zvar(:,:,:))
            a_cBCb%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN
               ! --Black Carbon (ice, regime A)
               CALL bulkMixrat('BC','ice','a',zvar(:,:,:))
               a_iBCa%data(:,:,:) = zvar(:,:,:)

               ! --Black Carbon (ice, regime B)
               CALL bulkMixrat('BC','ice','b',zvar(:,:,:))
               a_iBCb%data(:,:,:) = zvar(:,:,:)
            END IF ! level 5

         END IF

         IF (IsUsed(prtcl,'DU')) THEN

            !-- Dust (aerosol, regime A)
            CALL bulkMixrat('DU','aerosol','a',zvar(:,:,:))
            a_aDUa%data(:,:,:) = zvar(:,:,:)

            !-- Dust (aerosol, regime B)
            CALL bulkMixrat('DU','aerosol','b',zvar(:,:,:))
            a_aDUb%data(:,:,:) = zvar(:,:,:)

            !-- Dust (clouds, regime A)
            CALL bulkMixrat('DU','cloud','a',zvar(:,:,:))
            a_cDUa%data(:,:,:) = zvar(:,:,:)

            !-- Dust (clouds, regime B)
            CALL bulkMixrat('DU','cloud','b',zvar(:,:,:))
            a_cDUb%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN
               ! --Dust (ice, regime A)
               CALL bulkMixrat('DU','ice','a',zvar(:,:,:))
               a_iDUa%data(:,:,:) = zvar(:,:,:)

               ! --Dust (ice, regime B)
               CALL bulkMixrat('DU','ice','b',zvar(:,:,:))
               a_iDUb%data(:,:,:) = zvar(:,:,:)
            END IF ! level 5

         END IF

         IF (IsUsed(prtcl,'SS')) THEN

            !-- Sea Salt (aerosol, regime A)
            CALL bulkMixrat('SS','aerosol','a',zvar(:,:,:))
            a_aSSa%data(:,:,:) = zvar(:,:,:)

            !-- Sea Salt (aerosol, regime B)
            CALL bulkMixrat('SS','aerosol','b',zvar(:,:,:))
            a_aSSb%data(:,:,:) = zvar(:,:,:)

            !-- Sea Salt (clouds, regime A)
            CALL bulkMixrat('SS','cloud','a',zvar(:,:,:))
            a_cSSa%data(:,:,:) = zvar(:,:,:)

            !-- Sea Salt (clouds, regime B)
            CALL bulkMixrat('SS','cloud','b',zvar(:,:,:))
            a_cSSb%data(:,:,:) = zvar(:,:,:)

            IF (level == 5) THEN
               ! -- Sea Salt (ice, regime A)
               CALL bulkMixrat('SS','ice','a',zvar(:,:,:))
               a_iSSa%data(:,:,:) = zvar(:,:,:)

               ! -- Sea Salt (ice, regime B)
               CALL bulkMixrat('SS','ice','b',zvar(:,:,:))
               a_iSSb%data(:,:,:) = zvar(:,:,:)
            END IF ! level 5

         END IF

      END IF ! level 4


   END SUBROUTINE calc_anal

   !
   ! ----------------------------------------------------------------------
   ! Subroutine close_anal:  Closes netcdf anal file
   !
   INTEGER FUNCTION close_anal()

      USE netcdf

      close_anal = nf90_close(ncid0)

   END FUNCTION close_anal
   !
   ! ----------------------------------------------------------------------
   ! Subroutine Write_anal:  Writes the netcdf Analysis file
   !
   ! Modified for levels 4 and 5
   ! Juha Tonttila, FMI, 2014
   !
   !
   SUBROUTINE write_anal(time)
      USE netcdf
      USE mpi_interface, ONLY : myid, appl_abort
      USE class_ComponentIndex, ONLY : IsUsed
      USE mo_submctl, ONLY : in1a,fn2a,in2b,fn2b,            &
                             ica,fca,icb,fcb,ira,fra,        &
                             iia,fia,iib,fib,isa,fsa,        &
                             aerobins, cloudbins, precpbins, &
                             icebins, snowbins, nlim, prlim, &
                             nspec, nbins, ncld, nice, nprc, nsnw

      REAL, INTENT (in) :: time

      INTEGER :: iret, VarID, n, i
      INTEGER :: ibeg(4), icnt(4), i1, i2, j1, j2
      INTEGER :: ibegsd(5), icntaea(5), icntaeb(5), icntcla(5), icntclb(5), icntpra(5), & ! Juha: For sizedistribution variables
                 icntica(5), icnticb(5), icntsna(5)
      REAL :: zsum(nzp,nxp,nyp) ! Juha: Helper for computing bulk output diagnostics
      REAL :: zvar(nzp,nxp,nyp)
      REAL :: zvar2(nzp,nxp,nyp,fn2a-in1a)

      REAL :: a_Rawet(nzp,nxp,nyp,nbins), a_Rcwet(nzp,nxp,nyp,ncld), a_Rpwet(nzp,nxp,nyp,nprc), &
              a_Riwet(nzp,nxp,nyp,nice), a_Rswet(nzp,nxp,nyp,nsnw)
      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray3d), POINTER :: testi
      TYPE(FloatArray4d), POINTER :: testi2

      icnt =    (/nzp, nxp-4, nyp-4, 1/)
      icntaea = (/nzp, nxp-4, nyp-4, fn2a, 1 /)
      icntaeb = (/nzp, nxp-4, nyp-4, fn2b-fn2a, 1/)
      icntcla = (/nzp, nxp-4, nyp-4, fca%cur, 1/)
      icntclb = (/nzp, nxp-4, nyp-4, fcb%cur-fca%cur, 1/)
      icntpra = (/nzp, nxp-4, nyp-4, fra, 1/)
      icntica = (/nzp, nxp-4, nyp-4, fia%cur, 1/)
      icnticb = (/nzp, nxp-4, nyp-4, fib%cur-fia%cur, 1/)
      icntsna = (/nzp, nxp-4, nyp-4, fsa, 1/)
      ibeg =   (/1,1,1,nrec0/)
      ibegsd = (/1,1,1,1,nrec0/)

      i1 = 3
      i2 = nxp-2
      j1 = 3
      j2 = nyp-2

      iret = nf90_inq_varid(ncid0, "time", VarID) !time
      iret = nf90_put_var(ncid0, VarID, time, start=(/nrec0/))

      ! Calculate output variables
      CALL calc_anal()

      IF (nrec0 == 1) THEN

         iret = nf90_inq_varid(ncid0, "zt", VarID) !zt
         iret = nf90_put_var(ncid0, VarID, zt, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, "zm", VarID) !zm
         iret = nf90_put_var(ncid0, VarID, zm, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, "xt", VarID) !xt
         iret = nf90_put_var(ncid0, VarID, xt(i1:i2), start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, "xm", VarID) !xm
         iret = nf90_put_var(ncid0, VarID, xm(i1:i2), start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, "yt", VarID) !yt
         iret = nf90_put_var(ncid0, VarID, yt(j1:j2), start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, "ym", VarID) !ym
         iret = nf90_put_var(ncid0, VarID, ym(j1:j2), start = (/nrec0/))

         IF (level >= 4 .AND. lbinanl) THEN

            iret = nf90_inq_varid(ncid0,"aea", VarID) !aea
            iret = nf90_put_var(ncid0, VarID, aerobins(in1a:fn2a), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,"aeb", VarID) !aeb
            iret = nf90_put_var(ncid0,VarID, aerobins(in2b:fn2b), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,"cla", VarID) !cla
            iret = nf90_put_var(ncid0,VarID, cloudbins(ica%cur:fca%cur), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,"clb", VarID) !clb
            iret = nf90_put_var(ncid0,VarID, cloudbins(icb%cur:fcb%cur), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,"prc", VarID) !prc
            iret = nf90_put_var(ncid0,VarID, precpbins(ira:fra), start = (/nrec0/))

            IF (level == 5) THEN

               iret = nf90_inq_varid(ncid0,"ica", VarID) !ica
               iret = nf90_put_var(ncid0,VarID, icebins(iia%cur:fia%cur), start = (/nrec0/))
             
               iret = nf90_inq_varid(ncid0,"icb", VarID) !icb
               iret = nf90_put_var(ncid0,VarID, icebins(iib%cur:fib%cur), start = (/nrec0/))

               iret = nf90_inq_varid(ncid0,"snw", VarID) !snw
               iret = nf90_put_var(ncid0,VarID, snowbins(isa:fsa), start = (/nrec0/))

            END IF
         END IF
      END IF

      iret = nf90_inq_varid(ncid0, "u0", VarID) !u0
      iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
      iret = nf90_inq_varid(ncid0, "v0", VarID) !v0
      iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
      iret = nf90_inq_varid(ncid0, "dn0", VarID) !dn0
      iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))

      IF (level < 4) THEN ! Operation without SALSA

         DO n = 11, varnum

            CALL fields%getField(n,ArrEl)

            ! If output status is true, print variable to output file
            IF (ArrEl%outputstatus) THEN
               IF ((ArrEl%dimension == "tttt") .OR. (ArrEl%dimension == "mttt") .OR. &
                   (ArrEl%dimension == "tmtt") .OR. (ArrEl%dimension == "ttmt") ) THEN

                  CALL fields%getData(1,testi,name=ArrEl%name)
                  iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                  iret = nf90_put_var(ncid0,VarID,testi%data(:,i1:i2,j1:j2),start=ibeg, &
                                         count=icnt)
                  CYCLE

              END IF
           END IF

         END DO

      END IF

      IF (level >= 4) THEN ! Operation with SALSA


         ! Loop through SALSA output variables
            DO n = 19, varnum
               CALL fields%getField(n,ArrEl)

         ! If output status is true, print variable to output file
               IF (ArrEl%outputstatus) THEN

                  ! Output differently based on dimension
                  IF ((ArrEl%dimension == "tttt") .OR. (ArrEl%dimension == "mttt") .OR. &
                      (ArrEl%dimension == "tmtt") .OR. (ArrEl%dimension == "ttmt") ) THEN
                     CALL fields%getData(1,testi,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi%data(:,i1:i2,j1:j2),start=ibeg, &
                                         count=icnt)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttaea") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,in1a:fn2a),start=ibegsd, &
                                         count=icntaea)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttaeb") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,in2b:fn2b),start=ibegsd, &
                                         count=icntaeb)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttcla") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,ica%cur:fca%cur),start=ibegsd, &
                                         count=icntcla)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttclb") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,icb%cur:fcb%cur),start=ibegsd, &
                                         count=icntclb)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttprc") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,ira:fra),start=ibegsd, &
                                         count=icntpra)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttica") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,iia%cur:fia%cur),start=ibegsd, &
                                         count=icntica)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "tttticb") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,iib%cur:fib%cur),start=ibegsd, &
                                         count=icnticb)
                     CYCLE
                  END IF

                  IF (ArrEl%dimension == "ttttsnw") THEN
                     CALL fields%getData(1,testi2,name=ArrEl%name)
                     iret = nf90_inq_varid(ncid0,ArrEl%name,VarID)
                     iret = nf90_put_var(ncid0,VarID,testi2%data(:,i1:i2,j1:j2,isa:fsa),start=ibegsd, &
                                         count=icntsna)
                     CYCLE
                  END IF

               END IF

            END DO
         END IF ! level

         IF (myid == 0) PRINT "(//' ',12('-'),'   Record ',I3,' to: ',A60)",    &
            nrec0,fname

         iret  = nf90_sync(ncid0)
         nrec0 = nrec0+1

   END SUBROUTINE write_anal

   !
   ! ----------------------------------------------------------------------
   ! Subroutine write_hist:  This subroutine writes a binary history file
   !
   SUBROUTINE write_hist(htype, time)

      USE mpi_interface, ONLY : appl_abort, myid, wrxid, wryid
      INTEGER :: errcode = -17

      INTEGER, INTENT (in) :: htype
      REAL, INTENT (in)    :: time

      CHARACTER (len=80) :: hname

      INTEGER :: n, iblank
      !
      ! create and open a new output file.
      !
      WRITE(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      hname = trim(hname)//'.'//trim(filprf)

      SELECT CASE(htype)
         CASE DEFAULT
            hname = trim(hname)//'.iflg'
         CASE(0)
            hname = trim(hname)//'.R'
         CASE(1)
            hname = trim(hname)//'.rst'
         CASE(2)
            iblank=index(hname,' ')
            WRITE(hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(time), 's'
      END SELECT
      !
      ! Write fields
      !
      IF (myid == 0) PRINT "(//' ',49('-')/,' ',/,'   History write to: ',A30)" &
                            ,hname
      OPEN(10,file=trim(hname), form='unformatted')

      WRITE(10) time, th00, umean, vmean, dtl, level, isgstyp, iradtyp, nzp, nxp, nyp, nscl
      WRITE(10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf, sst, W1, W2, W3 ! added by Zubair

      WRITE(10) a_ustar, a_tstar, a_rstar

      WRITE(10) a_pexnr%data
      WRITE(10) a_press%data
      WRITE(10) a_theta%data

      WRITE(10) a_up
      WRITE(10) a_vp
      WRITE(10) a_wp
      WRITE(10) a_uc
      WRITE(10) a_vc
      WRITE(10) a_wc

      DO n = 1, nscl
         CALL newsclr(n)
         WRITE(10) a_sp
      END DO

      IF ( allocated(a_rv%alloc)   ) WRITE(10) a_rv%data
      IF ( allocated(a_rc%alloc)   ) WRITE(10) a_rc%data
      IF ( allocated(a_rflx)       ) WRITE(10) a_rflx
      CLOSE(10)

      IF (myid == 0 .AND. htype < 0) THEN
         PRINT *, 'CFL Violation'
         CALL appl_abort(errcode)
      END IF

      RETURN
   END SUBROUTINE write_hist
   !
   ! ----------------------------------------------------------------------
   ! Subroutine read_hist:  This subroutine reads a binary history file
   !
   !                        Modified for level 4
   !                Juha Tonttila, FMI, 20140828
   !

   SUBROUTINE read_hist(time, hfilin)

      USE mpi_interface, ONLY : appl_abort, myid, wrxid, wryid

      CHARACTER(len=80), INTENT(in) :: hfilin
      REAL, INTENT(out)             :: time

      CHARACTER (len=80) :: hname
      INTEGER :: n, nxpx, nypx, nzpx, nsclx, iradx, isgsx, lvlx
      LOGICAL :: exans
      REAL    :: umx, vmx, thx
      !
      ! open input file.
      !

      WRITE(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      hname = trim(hname)//'.'//trim(hfilin)

      INQUIRE(file=trim(hname),exist=exans)
      IF (.NOT. exans) THEN
         PRINT *,'ABORTING: History file', trim(hname),' not found'
         CALL appl_abort(0)
      ELSE
         OPEN(10,file=trim(hname),status='old',form='unformatted')
         READ(10) time,thx,umx,vmx,dtl,lvlx,isgsx,iradx,nzpx,nxpx,nypx,nsclx

         IF (nxpx /= nxp .OR. nypx /= nyp .OR. nzpx /= nzp)  THEN
            IF (myid == 0) PRINT *, nxp, nyp, nzp, nxpx, nypx, nzpx
            CALL appl_abort(-1)
         END IF

         READ(10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf,sst,W1,W2,W3

         READ(10) a_ustar, a_tstar, a_rstar

         READ(10) a_pexnr%data
         READ(10) a_press%data
         READ(10) a_theta%data

         READ(10) a_up
         READ(10) a_vp
         READ(10) a_wp
         READ(10) a_uc
         READ(10) a_vc
         READ(10) a_wc

         DO n = 1, nscl
            CALL newsclr(n)
            IF (n <= nsclx) READ(10) a_sp
         END DO
         DO n = nscl+1, nsclx
            READ(10)
         END DO

         IF (lvlx > 0 .AND. lvlx < 4) THEN
            IF (level > 0 .AND. lvlx < 4) THEN
               READ(10) a_rv%data
            ELSE
               READ(10)
            END IF
         END IF
         IF (lvlx > 1) THEN
            IF (level > 1) THEN
               READ(10) a_rc%data
            ELSE
               READ(10)
            END IF
         END IF
         IF (iradx > 0) THEN
            IF (iradtyp > 0) THEN
               READ(10) a_rflx
            ELSE
               READ(10)
            END IF
         END IF

         CLOSE(10)
         !
         ! adjust namelist and basic state appropriately
         !
         IF (thx /= th00) THEN
            IF (myid == 0) PRINT "('  th00 changed  -  ',2f8.2)",th00,thx
            a_tp%data = a_tp%data + thx - th00
         END IF
         IF (umx /= umean) THEN
            IF (myid == 0) PRINT "('  umean changed  -  ',2f8.2)",umean,umx
            a_up = a_up + umx - umean
         END IF
         IF (vmx /= vmean) THEN
            IF (myid == 0) PRINT "('  vmean changed  -  ',2f8.2)",vmean,vmx
            a_vp = a_vp + vmx - vmean
         END IF
         dtlv = 2.*dtl
         dtlt = dtl

      END IF

   END SUBROUTINE read_hist
   !
   ! ----------------------------------------------------------------------
   ! Subroutine newsclr:  This routine updates the scalar pointer to the
   ! value corresponding to the next scalar in the scalar table
   !
   SUBROUTINE newsclr(iscnum)

      INTEGER, INTENT(in) :: iscnum

      a_sp => a_sclrp(:,:,:,iscnum)
      a_st => a_sclrt(:,:,:,iscnum)

      RETURN
   END SUBROUTINE newsclr
   !
   ! -----------------------------------
   ! Subroutine bulkMixrat: Find and calculate
   ! the total mixing ratio of a given compound
   ! in aerosol particles or hydrometeors
   !
   ! Juha Tonttila, FMI, 2015
   ! Jaakko Ahola, FMI, 2015
   SUBROUTINE bulkMixrat(icomp,ipart,itype,mixrat)
      USE mo_submctl, ONLY : ncld,nbins,nprc,   &
                             ica,fca,icb,fcb,   &
                             ira,fra,           &
                             nice,nsnw,         &
                             iia,fia,iib,fib,   &
                             isa,fsa,           &
                             in1a,in2b,         &
                             fn2a,fn2b

      USE class_ComponentIndex, ONLY : GetIndex, IsUsed

      CHARACTER(len=*), INTENT(in) :: icomp  ! This should be either:
                                             ! SO4,OC,NO,NH,BC,DU,SS,H2O.

      CHARACTER(len=*), INTENT(in) :: ipart  ! This should be either:
                                             ! aerosol,cloud,rain,ice,snow
      CHARACTER(len=*), INTENT(in) :: itype  ! Select bin regime: a or b

      REAL, INTENT(out) :: mixrat(nzp,nxp,nyp)

      INTEGER :: istr,iend,mm

      mixrat = 0.

      ! Determine multipliers

      mm = GetIndex(prtcl,icomp)

      ! Given in kg/kg
      SELECT CASE(ipart)
         CASE('aerosol')
            IF (itype == 'a') THEN
               istr = (mm-1)*nbins + in1a
               iend = (mm-1)*nbins + fn2a
            ELSE IF (itype == 'b') THEN
               istr = (mm-1)*nbins + in2b
               iend = (mm-1)*nbins + fn2b
            ELSE
               STOP 'bulkMixrat: Invalid aerosol bin regime SELECTion'
            END IF
            mixrat(:,:,:) = SUM(a_maerop%data(:,:,:,istr:iend),DIM=4)
         CASE('cloud')
            IF (itype == 'a') THEN
               istr = (mm-1)*ncld + ica%cur
               iend = (mm-1)*ncld + fca%cur
            ELSE IF (itype == 'b') THEN
               istr = (mm-1)*ncld + icb%cur
               iend = (mm-1)*ncld + fcb%cur
            ELSE
               STOP 'bulkMixrat: Invalid cloud bin regime SELECTion'
            END IF
            mixrat(:,:,:) = SUM(a_mcloudp%data(:,:,:,istr:iend),DIM=4)
         CASE('precp')
            istr = (mm-1)*nprc + ira
            iend = (mm-1)*nprc + fra
            mixrat(:,:,:) = SUM(a_mprecpp%data(:,:,:,istr:iend),DIM=4)
         CASE('ice')
            IF (itype == 'a') THEN
               istr = (mm-1)*nice + iia%cur
               iend = (mm-1)*nice + fia%cur
            ELSE IF (itype == 'b') THEN
               istr = (mm-1)*nice + iib%cur
               iend = (mm-1)*nice + fib%cur
            ELSE
               STOP 'bulkMixrat: Invalid ice bin regime SELECTion'
            END IF
            mixrat(:,:,:) = SUM(a_micep%data(:,:,:,istr:iend),DIM=4)
         CASE('snow')
            istr = (mm-1)*nsnw + isa
            iend = (mm-1)*nsnw + fsa
            mixrat(:,:,:) = SUM(a_msnowp%data(:,:,:,istr:iend),DIM=4)
      END SELECT

   END SUBROUTINE bulkMixrat
   !
   ! ----------------------------------------------
   ! Subroutine binSpecMixrat: Calculate the mixing
   ! ratio of selected aerosol species in individual
   ! bins.
   !
   ! Juha Tonttila, FMI, 2015
   SUBROUTINE binSpecMixrat(ipart,icomp,ibin,mixr)
      USE mo_submctl, ONLY : ncld, nbins, nprc, nice, nsnw
      USE class_componentIndex, ONLY : GetIndex

      CHARACTER(len=*), INTENT(in) :: icomp  ! This should be either:
                                             ! SO4,OC,NO,NH,BC,DU,SS,H2O.

      CHARACTER(len=*), INTENT(in) :: ipart  ! This should be either:
                                             ! aerosol,cloud,rain,ice,snow
      INTEGER, INTENT(in) :: ibin

      REAL, INTENT(out)   :: mixr(nzp,nxp,nyp)

      INTEGER :: mm

      ! Determine multipliers
      mm = GetIndex(prtcl,icomp)

      SELECT CASE(ipart)
         CASE('aerosol')
            mixr(:,:,:) = a_maerop%data(:,:,:,(mm-1)*nbins+ibin)
         CASE('cloud')
            mixr(:,:,:) = a_mcloudp%data(:,:,:,(mm-1)*ncld+ibin)
         CASE('precp')
            mixr(:,:,:) = a_mprecpp%data(:,:,:,(mm-1)*nprc+ibin)
         CASE('ice')
            mixr(:,:,:) = a_micep%data(:,:,:,(mm-1)*nice+ibin)
         CASE('snow')
            mixr(:,:,:) = a_msnowp%data(:,:,:,(mm-1)*nsnw+ibin)
      END SELECT

   END SUBROUTINE binSpecMixrat

   !
   ! ----------------------------------------------
   ! Subroutine bulkNumc: Calculate the total number
   ! concentration of particles of given type
   !
   ! Juha Tonttila, FMI, 2015
   !
   SUBROUTINE bulkNumc(ipart,itype,numc)
      USE mo_submctl, ONLY : ica,fca,icb,fcb,ira,fra, &
         iia,fia,iib,fib,isa,fsa, &
         in1a,in2b,fn2a,fn2b

      CHARACTER(len=*), INTENT(in) :: ipart
      CHARACTER(LEN=*), INTENT(in) :: itype
      REAL, INTENT(out) :: numc(nzp,nxp,nyp)
      INTEGER :: istr,iend

      istr = 0
      iend = 0

      ! Outputs #/kg
      ! No concentration limits (nlim or prlim) for number

      SELECT CASE(ipart)
         CASE('aerosol')
            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = in1a
               iend = fn2b
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = in1a
               iend = fn2a
            ELSE IF (itype == 'b') THEN ! 2b
               istr = in2b
               iend = fn2b
            END IF
            numc(:,:,:) = SUM(a_naerop%data(:,:,:,istr:iend),DIM=4)
         CASE('cloud')
            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = ica%cur
               iend = fcb%cur
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = ica%cur
               iend = fca%cur
            ELSE IF (itype == 'b') THEN
               istr = icb%cur
               iend = fcb%cur
            END IF
            numc(:,:,:) = SUM(a_ncloudp%data(:,:,:,istr:iend),DIM=4)
         CASE('precp')
            istr = ira
            iend = fra
            numc(:,:,:) = SUM(a_nprecpp%data(:,:,:,istr:iend),DIM=4)
         CASE('ice')
            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = iia%cur
               iend = fib%cur
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = iia%cur
               iend = fia%cur
            ELSE IF (itype == 'b') THEN
               istr = iib%cur
               iend = fib%cur
            END IF
            numc(:,:,:) = SUM(a_nicep%data(:,:,:,istr:iend),DIM=4)
         CASE('snow')
            istr = isa
            iend = fsa
            numc(:,:,:) = SUM(a_nsnowp%data(:,:,:,istr:iend),DIM=4)
      END SELECT

   END SUBROUTINE bulkNumc


   !
   ! -------------------------------------------------
   ! Subroutine meanRadius
   ! Gets the mean wet (water=nspec+1) radius for particles.
   !
   SUBROUTINE meanRadius(ipart,itype,rad)
      USE mo_submctl, ONLY : nbins,ncld,nprc,               &
                             nice,nsnw,                     &
                             ica,fca,icb,fcb,ira,fra,       &
                             iia,fia,iib,fib,isa,fsa,       &
                             in1a,fn2a,in2b,fn2b,           &
                             nlim,prlim,nspec
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(in) :: ipart
      CHARACTER(len=*), INTENT(in) :: itype
      REAL, INTENT(out) :: rad(nzp,nxp,nyp)

      INTEGER :: istr,iend

      rad = 0.

      SELECT CASE(ipart)
         CASE('aerosol')

            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = in1a
               iend = fn2b
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = in1a
               iend = fn2a
            ELSE IF (itype == 'b') THEN
               istr = in2b
               iend = fn2b
            ELSE
               STOP 'meanRadius: Invalid bin regime SELECTion (aerosol)'
            END IF

            CALL getRadius(istr,iend,nbins,nspec+1,a_naerop%data,a_maerop%data,nlim,rad,1)

         CASE('cloud')

            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = ica%cur
               iend = fcb%cur
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = ica%cur
               iend = fca%cur
            ELSE IF (itype == 'b') THEN
               istr = icb%cur
               iend = fcb%cur
            ELSE
               STOP 'meanRadius: Invalid bin regime SELECTion (cloud)'
            END IF

            CALL getRadius(istr,iend,ncld,nspec+1,a_ncloudp%data,a_mcloudp%data,nlim,rad,2)

         CASE('precp')

            istr = ira
            iend = fra

            CALL getRadius(istr,iend,nprc,nspec+1,a_nprecpp%data,a_mprecpp%data,prlim,rad,3)

         CASE('ice')

            IF (itype == 'ab') THEN ! Note: 1a and 2a, 2b combined
               istr = iia%cur
               iend = fib%cur
            ELSE IF (itype == 'a') THEN ! Note: 1a and 2a combined
               istr = iia%cur
               iend = fia%cur
            ELSE IF (itype == 'b') THEN
               istr = iib%cur
               iend = fib%cur
            ELSE
               STOP 'meanRadius: Invalid bin regime SELECTion (ice)'
            END IF

            CALL getRadius(istr,iend,nice,nspec+1,a_nicep%data,a_micep%data,prlim,rad,4)

         CASE('snow')

            istr = isa
            iend = fsa

            CALL getRadius(istr,iend,nsnw,nspec+1,a_nsnowp%data,a_msnowp%data,prlim,rad,5)

      END SELECT

 !  END SUBROUTINE meanRadius
   !
   ! ---------------------------------------------------
   ! Subroutine getRadius
   ! Calculates number mean wet radius (over selected bins) for the whole domain
  CONTAINS

   SUBROUTINE getRadius(zstr,zend,nn,n4,numc,mass,numlim,zrad,flag)
      USE mo_submctl, ONLY : pi6

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nn, n4 ! Number of bins (nn) and aerosol species (n4)
      INTEGER, INTENT(in) :: zstr,zEND  ! Start and end index for averaging
      REAL, INTENT(in) :: numc(nzp,nxp,nyp,nn)
      REAL, INTENT(in) :: mass(nzp,nxp,nyp,nn*n4)
      REAL, INTENT(in) :: numlim
      INTEGER, INTENT(in) :: flag

      REAL, INTENT(out) :: zrad(nzp,nxp,nyp)

      INTEGER :: k,i,j,bin
      REAL :: tot, rwet, tmp(n4)

      zrad(:,:,:) = 0.

      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 1, nzp

               tot = 0.
               rwet = 0.

               DO bin = zstr, zend
                  IF (numc(k,i,j,bin) > numlim) THEN
                     tot = tot+numc(k,i,j,bin)
                     tmp(:) = mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                     rwet = rwet+calc_eff_radius(n4,numc(k,i,j,bin),tmp,flag)*numc(k,i,j,bin)
                  END IF
               END DO

               IF (tot > numlim) THEN
                  zrad(k,i,j) = rwet/tot
               END IF

            END DO
         END DO
      END DO


   END SUBROUTINE getRadius
   END SUBROUTINE meanRadius

   ! SUBROUTINE getBinRadius
   ! Calculates wet radius for each bin in the whole domain

   SUBROUTINE getBinRadius(nn,n4,numc,mass,numlim,zrad,flag)
      USE mo_submctl, ONLY : pi6

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nn, n4 ! Number of bins (nn) and aerosol species (n4)
      REAL, INTENT(in)  :: numc(nzp,nxp,nyp,nn)
      REAL, INTENT(in)  :: mass(nzp,nxp,nyp,nn*n4)
      REAL, INTENT(in)  :: numlim
      REAL, INTENT(out) :: zrad(nzp,nxp,nyp,nn)
      INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)

      INTEGER :: k,i,j,bin
      REAL :: tmp(n4)

      zrad(:,:,:,:) = 0.
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 1, nzp

               DO bin = 1, nn
                  IF (numc(k,i,j,bin) > numlim) THEN
                     tmp(:) = mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                     zrad(k,i,j,bin) = calc_eff_radius(n4,numc(k,i,j,bin),tmp,flag)
                  END IF
               END DO

            END DO
         END DO
      END DO



   END SUBROUTINE getBinRadius


  !********************************************************************
  !
  ! Function for calculating effective (wet) radius for any particle type
  ! - Aerosol, cloud and rain are spherical
  ! - Snow and ice can be irregular and their densities can be size-dependent
  !
  ! Edit this function when needed (also update CalcDimension in mo_submctl.f90)
  !
  ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
  ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
  !
   REAL FUNCTION calc_eff_radius(n,numc,mass,flag)
      USE mo_submctl, ONLY : pi6
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n ! Number of species
      INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
      REAL, INTENT(IN) :: numc, mass(n)

      ! Don't calculate if very low number concentration
      IF (numc < 1e-15) RETURN

      calc_eff_radius=0.
      IF (flag == 4) THEN   ! Ice
         ! Spherical ice
         calc_eff_radius = 0.5*( SUM(mass(:)/dens_ice(:))/numc/pi6)**(1./3.)
      ELSE IF (flag == 5) THEN   ! Snow
         ! Spherical snow
         calc_eff_radius = 0.5*( SUM(mass(:)/dens_snow(:))/numc/pi6)**(1./3.)
      ELSE
         ! Radius from total volume of a spherical particle or aqueous droplet
         calc_eff_radius = 0.5*( SUM(mass(:)/dens(:))/numc/pi6)**(1./3.)
      END IF


   END FUNCTION calc_eff_radius
   !********************************************************************
END MODULE grid

