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
   !
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


   CHARACTER (len=7), ALLOCATABLE, SAVE :: sanal(:)
   CHARACTER (len=80) :: expnme = 'Default' ! Experiment name
   CHARACTER (len=80) :: filprf = 'x'       ! File Prefix
   CHARACTER (len=7)  :: runtype = 'INITIAL'! Run Type SELECTion

   REAL               :: Tspinup = 7200.    ! Spinup period in seconds (added by Juha)


   CHARACTER (len=7),  PRIVATE :: v_snm = 'sxx    '
   CHARACTER (len=80), PRIVATE :: fname

   INTEGER, PRIVATE, SAVE  ::  nrec0, nvar0, nbase=15

   INTEGER           :: nz, nxyzp, nxyp
   REAL              :: dxi, dyi, dtl, dtlv, dtlt, umean, vmean, psrf
   REAL, ALLOCATABLE :: xt(:), xm(:), yt(:), ym(:), zt(:), zm(:), dzt(:), dzm(:)
   REAL, ALLOCATABLE :: u0(:), v0(:), pi0(:), pi1(:), th0(:), dn0(:), rt0(:)
   REAL, ALLOCATABLE :: spng_wfct(:), spng_tfct(:)
   !
   ! velocity variables (past, current and tendency)
   !
   REAL, ALLOCATABLE, TARGET :: a_up(:,:,:),a_uc(:,:,:),a_ut(:,:,:)
   REAL, ALLOCATABLE, TARGET :: a_vp(:,:,:),a_vc(:,:,:),a_vt(:,:,:)
   REAL, ALLOCATABLE, TARGET :: a_wp(:,:,:),a_wc(:,:,:),a_wt(:,:,:)
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

   TYPE(FloatArray3D), TARGET :: a_qp, a_qt
   REAL, POINTER :: a_sp(:,:,:),a_st(:,:,:)  !dont touch yet AZ

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
   REAL, ALLOCATABLE :: dens(:), mws(:), mvol(:), diss(:)

   !---------------------------------------------------------------------------

   REAL, ALLOCATABLE, TARGET :: a_sclrp(:,:,:,:),a_sclrt(:,:,:,:)
    !
   ! 3d diagnostic quantities
   !
   TYPE(FloatArray3D), TARGET :: a_theta  ! dry potential temp (k)
   TYPE(FloatArray3D), TARGET :: a_pexnr  ! perturbation exner func
   TYPE(FloatArray3D), TARGET :: a_press  ! pressure (hpa)
   TYPE(FloatArray3D), TARGET :: a_rc     ! Total cloud water
   TYPE(FloatArray3D), TARGET :: a_ri     ! Total ice cloud content
   TYPE(FloatArray3D), TARGET :: a_rv     ! water vapor (used only for levels < 4!)
   TYPE(FloatArray3D), TARGET :: a_srp    ! Total rain water for use with LEVEL 4 (Diagnostic scalar!!)
   TYPE(FloatArray3D), TARGET :: a_snrp   ! Total number of rain drops for use with LEVEL 4 (Diagnostic scalar!!)
   TYPE(FloatArray3D), TARGET :: a_srs    ! Total snow for use with SALSA
   TYPE(FloatArray3D), TARGET :: a_snrs   ! Total number of snow particles for use with LEVEL 5 (Diagnostic scalar!!)
   TYPE(FloatArray3D), TARGET :: a_rh     ! Relative humidity
   TYPE(FloatArray3D), TARGET :: a_rsl    ! water saturation vapor mixing ratio
   TYPE(FloatArray3D), TARGET :: a_rhi    ! Relative humidity over ice
   TYPE(FloatArray3D), TARGET :: a_rsi    ! ice saturation vapor mixing ratio
   TYPE(FloatArray3D), TARGET :: a_dn     ! Air density (for normalizing concentrations according to mass, levels < 4!)
   !
   ! scratch arrays
   !
   REAL, ALLOCATABLE, DIMENSION (:,:,:) :: a_rflx, a_sflx, &
                                           a_fus, a_fds, a_fuir, a_fdir, &
                                           a_temp, a_temp0 ! store temperatures of previous timestep
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

   CLASS(*), POINTER :: dummy, dummyt, dummyp !dummy pointers for field arrays

   TYPE(FloatArray3d), POINTER :: testi_diag

   TYPE(FieldArray) :: ProgF
   TYPE(FieldArray) :: DiagF

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
      USE mo_submctl, ONLY : nbins,ncld,nprc,  & ! Number of aerosol and hydrometeor size bins for SALSA
                             nice,nsnw,        & ! number of ice and snow size bins for SALSA
                             nspec, maxspec, listspec, &
                             rhosu,rhooc,rhobc,rhono,rhonh,rhoss,rhodu,rhowa,  &
                             msu,moc,mbc,mno,mnh,mss,mdu,mwa

      USE class_ComponentIndex, ONLY : ComponentIndexConstructor,  &
                                       GetNcomp, IsUsed

      INTEGER :: memsize
      INTEGER :: zz
      INTEGER :: nc
      REAL :: zeros3d(nzp,nxp,nyp)  ! array to help allocate 3d FloatArrays to zero -AZ
      zeros3d(:,:,:) = 0.

      ProgF = FieldArray()
      DiagF = FieldArray()

      ! Juha: Number of prognostic tracers for SALSA
      !            Aerosol bins + Cloud bins + gas compound tracers

      IF (level >= 4) THEN
         ! Create index tables for different aerosol components (can be fetched by name using getIndex)
         CALL ComponentIndexConstructor(prtcl, nspec, maxspec, listspec)
         nc = GetNcomp(prtcl)

         nsalsa = (nc+2)*nbins + (nc+2)*ncld + (nc+2)*nprc + (nc+2)*nice + (nc+2)*nsnw + 5

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

      ALLOCATE (a_temp(nzp,nxp,nyp),a_temp0(nzp,nxp,nyp))
      a_temp(:,:,:) = 0.
      a_temp0(:,:,:) = 0.

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
         CALL ProgF%NewField("t","Liquid water potential temperature","K","tttt",.TRUE.,dummyp,dummyt)

         IF (level >= 0) THEN
            a_rp = FloatArray3D(a_sclrp(:,:,:,2),store=.FALSE.)
            a_rt = FloatArray3D(a_sclrt(:,:,:,2),store=.FALSE.)
            dummyp => a_rp
            dummyt => a_rt
            CALL ProgF%NewField("r","Rain-water mixing ratio","kg/kg","tttt",.TRUE.,dummyp,dummyt)
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

         zz = zz+(nc+1)*nsnw
         ! Density, molecular weight, dissociation factor and molar volume arrays for the used species
         !   1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O
         ALLOCATE ( dens(nc+1), mws(nc+1), diss(nc+1), mvol(nc+1) )
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

         IF (zz /= nc+1) THEN
            WRITE(*,*) 'Physical properties not found for all species!'
            STOP
         END IF

         ! Molar volume volume, (kg/mol)/(kg/m^3)=m^3/mol
         mvol(:) = mws(:)/dens(:)

      END IF ! level

      !----------------------------------------------------

      ALLOCATE (a_ustar(nxp,nyp),a_tstar(nxp,nyp),a_rstar(nxp,nyp))
      ALLOCATE (uw_sfc(nxp,nyp),vw_sfc(nxp,nyp),ww_sfc(nxp,nyp))
      ALLOCATE (wt_sfc(nxp,nyp),wq_sfc(nxp,nyp))
      !ALLOCATE (ra(nxp,nyp))
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

      INTEGER :: i,j,k,kmax,nchby
      REAL    :: dzrfm,dz,zb,dzmin
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

      USE mpi_interface, ONLY : myid, ver, author, info
      USE mo_submctl, ONLY : fn2a,fn2b,fca,fcb,fra, &
                             fia,fib,fsa
      USE class_ComponentIndex, ONLY : IsUsed
      INTEGER, PARAMETER :: nnames = 24
      INTEGER, PARAMETER :: salsa_nn = 104
      CHARACTER (len=7), SAVE :: sbase(nnames) =  (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     '   ,& ! 1
         'ym     ','u0     ','v0     ','dn0    ','u      ','v      '   ,& ! 7
         'w      ','theta  ','p      ','q      ','l      ','r      '   ,& ! 13
         'f      ','i      ','s      '                                 ,& ! 19 ice'n'snow
         'n      ','stke   ','rflx   '/)                                  ! 22 total 24
      ! Added for SALSA
      CHARACTER(len=7), SAVE :: salsa_sbase(salsa_nn) = (/ &
         'time   ','zt     ','zm     ','xt     ','xm     ','yt     ',  &  ! 1 
         'ym     ','aea    ','aeb    ','cla    ','clb    ','prc    ',  &  ! 7
         'ica    ','icb    ','snw    ','u0     ','v0     ','dn0    ',  &  ! 13
         'u      ','v      ','w      ','theta  ','p      ','q      ',  &  ! 19
         'l      ','r      ','f      ','i      ','s      ',            &  ! 25
         'S_RH   ','S_RHI  ','S_Nact ','S_Na   ','S_Naba ','S_Rwaa ',  &  ! 30
         'S_Rwaba','S_Nb   ','S_Nabb ','S_Rwab ','S_Rwabb','S_Nc   ',  &  ! 36
         'S_Ncba ','S_Ncbb ','S_Rwca ','S_Rwcb ','S_Rwcba','S_Rwcbb',  &  ! 42
         'S_Np   ','S_Npba ','S_Rwpa ','S_Rwpba','S_Nic  ','S_Niba ',  &  ! 48
         'S_Nibb ','S_Rwia ','S_Rwib ','S_Rwiba','S_Rwibb','S_Ns   ',  &  ! 54
         'S_Nsba ','S_Rwsa ','S_Rwsba','S_aSO4a','S_aNHa ','S_aNOa ',  &  ! 60
         'S_aOCa ','S_aBCa ','S_aDUa ','S_aSSa ','S_aSO4b','S_aNHb ',  &  ! 66
         'S_aNOb ','S_aOCb ','S_aBCb ','S_aDUb ','S_aSSb ','S_cSO4a',  &  ! 72
         'S_cNHa ','S_cNOa ','S_cOCa ','S_cBCa ','S_cDUa ','S_cSSa ',  &  ! 78
         'S_cSO4b','S_cNHb ','S_cNOb ','S_cOCb ','S_cBCb ','S_cDUb ',  &  ! 84
         'S_cSSb ','S_iSO4a','S_iNHa ','S_iNOa ','S_iOCa ','S_iBCa ',  &  ! 90
         'S_iDUa ','S_iSSa ','S_iSO4b','S_iNHb ','S_iNOb ','S_iOCb ',  &  ! 96
         'S_iBCb ','S_iDUb ','S_iSSb ' /)
           ! total 104

      LOGICAL, SAVE :: salsabool(salsa_nn)

      REAL, INTENT (in)    :: time
      LOGICAL, INTENT (in) :: salsa_b_bins
      INTEGER              :: nbeg, nend

      IF (level < 4) THEN  ! Standard operation for levels 1-3

         nvar0 = nbase + naddsc
         IF (level >= 1) nvar0 = nvar0+1
         IF (level >= 2) nvar0 = nvar0+1
         IF (level == 3) nvar0 = nvar0+2
         IF (isgstyp > 1) nvar0 = nvar0+1
         IF (iradtyp > 1) nvar0 = nvar0+1

         ALLOCATE (sanal(nvar0))
         sanal(1:nbase) = sbase(1:nbase)

         nvar0 = nbase
         !
         ! add liquid water, which is a diagnostic variable, first
         !
         IF (level >= 2) THEN
            nvar0 = nvar0+1
            sanal(nvar0) = sbase(nbase+2)
         END IF
         !
         ! add additional scalars, in the order in which they appear in scalar
         ! table
         !
         IF (level >= 1) THEN
            nvar0 = nvar0+1
            sanal(nvar0) = sbase(nbase+1)
         END IF

         IF (level == 3) THEN
            nvar0 = nvar0+1
            sanal(nvar0) = sbase(nbase+3)
            nvar0 = nvar0+1
            sanal(nvar0) = sbase(nbase+4+3)
         END IF

         IF (isgstyp > 1) THEN
            nvar0 = nvar0+1
            sanal(nvar0) = sbase(nbase+5+3)
         END IF

         IF (iradtyp > 2) THEN
            nvar0 = nvar0+1
            sanal(nvar0) = sbase(nbase+6+3)
         END IF

         nbeg = nvar0+1
         nend = nvar0+naddsc
         DO nvar0 = nbeg, nend
            WRITE(v_snm(2:3),'(i2.2)') nvar0-nbeg
            sanal(nvar0) = v_snm
         END DO
         nvar0 = nend

      ELSE IF (level >= 4) THEN ! Operation with SALSA

         ! Make a boolean array masking the output variables that are used.
         ! This mainly concerns unused aerosol species.
         salsabool(:) = .TRUE.
         IF (.NOT. lbinanl) THEN
            salsabool(8:15) = .FALSE.
            salsabool((/34,36,38,40,42,43,46,47,49,51,53,54,57,58,60,62/)) = .FALSE.
         END IF

         IF (level < 5 ) THEN
            salsabool(13:15) = .FALSE. ! ica, icb, snw
            salsabool(27:29) = .FALSE.
            salsabool(91:104) = .FALSE. ! aerosols in ice particles
            salsabool(31) = .FALSE.
            salsabool(52:62) = .FALSE.
         END IF

         IF (.NOT. IsUsed(prtcl,'SO4')) &
            salsabool((/ 63, 70, 77, 84, 91,  98 /)) = .FALSE.

         IF (.NOT. IsUsed(prtcl,'NH'))  &
            salsabool((/ 64, 71, 78, 85, 92,  99 /)) = .FALSE.

         IF (.NOT. IsUsed(prtcl,'NO'))  &
            salsabool((/ 65, 72, 79, 86, 93,  100 /)) = .FALSE.

         IF (.NOT. IsUsed(prtcl,'OC'))  &
            salsabool((/ 66, 73, 80, 87, 94,  101 /)) = .FALSE.

         IF (.NOT. IsUsed(prtcl,'BC'))  &
            salsabool((/ 67, 74, 81, 88, 95,  102 /)) = .FALSE.

         IF (.NOT. IsUsed(prtcl,'DU'))  &
            salsabool((/ 68, 75, 82, 89, 96,  103 /)) = .FALSE.

         IF (.NOT. IsUsed(prtcl,'SS'))  &
            salsabool((/ 69, 76, 83, 90, 97,  104 /)) = .FALSE.

            ! b-bins are not always saved
         IF (.NOT. salsa_b_bins) THEN
             salsabool((/ 37, 38, 39, 40, 43, 45, 47, 54, 56, 58 /)) = .FALSE.
             salsabool(70:76) = .FALSE.    ! Aerosol species
             salsabool(84:90) = .FALSE.    ! Cloud species
             salsabool(98:104) = .FALSE.   ! Ice species
         END IF

         nvar0 = COUNT(salsabool) + naddsc
         ALLOCATE(sanal(nvar0))

         sanal = PACK(salsa_sbase,salsabool)

      END IF

      fname =  trim(filprf)
      IF(myid == 0) PRINT                                                  &
         "(//' ',49('-')/,' ',/,'   Initializing: ',A20)",trim(fname)
      CALL open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid0, nrec0, ver, author, info)

      IF (level < 4 .OR. .NOT. lbinanl) THEN
         CALL define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4)

      ELSE IF (level == 4 .AND. lbinanl) THEN
         CALL define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                         inae_a=fn2a, inae_b=fn2b-fn2a, incld_a=fca%cur,          &
                         incld_b=fcb%cur-fca%cur, inprc=fra )
      ELSE IF (level == 5 .AND. lbinanl) THEN
         CALL define_nc( ncid0, nrec0, nvar0, sanal, n1=nzp, n2=nxp-4, n3=nyp-4,  &
                         inae_a=fn2a,  inae_b =fn2b-fn2a, incld_a=fca%cur,        &
                         incld_b=fcb%cur-fca%cur, inprc=fra, inice_a=fia%cur,     &
                         inice_b=fib%cur-fia%cur, insnw=fsa )
      END IF
      IF (myid == 0) PRINT *,'   ...starting record: ', nrec0


   END SUBROUTINE init_anal
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

      INTEGER :: iret, VarID, nn, n
      INTEGER :: ibeg(4), icnt(4), i1, i2, j1, j2
      INTEGER :: ibegsd(5), icntaea(5), icntaeb(5), icntcla(5), icntclb(5), icntpra(5), & ! Juha: For sizedistribution variables
                 icntica(5), icnticb(5), icntsna(5)
      REAL :: zsum(nzp,nxp,nyp) ! Juha: Helper for computing bulk output diagnostics
      REAL :: zvar(nzp,nxp,nyp)
      REAL :: a_Rawet(nzp,nxp,nyp,nbins), a_Rcwet(nzp,nxp,nyp,ncld),a_Rpwet(nzp,nxp,nyp,nprc), &
              a_Riwet(nzp,nxp,nyp,nice),a_Rswet(nzp,nxp,nyp,nsnw)

      icnt = (/nzp, nxp-4, nyp-4, 1/)
      icntaea = (/nzp,nxp-4,nyp-4, fn2a, 1 /)
      icntaeb = (/nzp,nxp-4,nyp-4, fn2b-fn2a, 1/)
      icntcla = (/nzp,nxp-4,nyp-4, fca%cur, 1/)
      icntclb = (/nzp,nxp-4,nyp-4, fcb%cur-fca%cur, 1/)
      icntpra = (/nzp,nxp-4,nyp-4, fra, 1/)
      icntica = (/nzp,nxp-4,nyp-4, iia%cur, 1/)
      icnticb = (/nzp,nxp-4,nyp-4, iib%cur-iia%cur, 1/)
      icntsna = (/nzp,nxp-4,nyp-4, fsa, 1/)
      ibeg = (/1  ,1  ,1  ,nrec0/)
      ibegsd = (/1,1,1,1,nrec0/)

      i1 = 3
      i2 = nxp-2
      j1 = 3
      j2 = nyp-2

      iret = nf90_inq_varid(ncid0, sanal(1), VarID)
      iret = nf90_put_var(ncid0, VarID, time, start=(/nrec0/))

      IF (nrec0 == 1) THEN
         iret = nf90_inq_varid(ncid0, sanal(2), VarID)
         iret = nf90_put_var(ncid0, VarID, zt, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(3), VarID)
         iret = nf90_put_var(ncid0, VarID, zm, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(4), VarID)
         iret = nf90_put_var(ncid0, VarID, xt(i1:i2), start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(5), VarID)
         iret = nf90_put_var(ncid0, VarID, xm(i1:i2), start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(6), VarID)
         iret = nf90_put_var(ncid0, VarID, yt(j1:j2), start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(7), VarID)
         iret = nf90_put_var(ncid0, VarID, ym(j1:j2), start = (/nrec0/))

         IF (level >= 4 .AND. lbinanl) THEN

            iret = nf90_inq_varid(ncid0,sanal(8), VarID)
            iret = nf90_put_var(ncid0, VarID, aerobins(in1a:fn2a), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,sanal(9), VarID)
            iret = nf90_put_var(ncid0,VarID, aerobins(in2b:fn2b), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,sanal(10), VarID)
            iret = nf90_put_var(ncid0,VarID, cloudbins(ica%cur:fca%cur), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,sanal(11), VarID)
            iret = nf90_put_var(ncid0,VarID, cloudbins(icb%cur:fcb%cur), start = (/nrec0/))

            iret = nf90_inq_varid(ncid0,sanal(12), VarID)
            iret = nf90_put_var(ncid0,VarID, precpbins(ira:fra), start = (/nrec0/))

            IF (level == 5) THEN
               iret = nf90_inq_varid(ncid0,sanal(13), VarID)
               iret = nf90_put_var(ncid0,VarID, icebins(iia%cur:fia%cur), start = (/nrec0/))
             
               iret = nf90_inq_varid(ncid0,sanal(14), VarID)
               iret = nf90_put_var(ncid0,VarID, icebins(iib%cur:fib%cur), start = (/nrec0/))

               iret = nf90_inq_varid(ncid0,sanal(15), VarID)
               iret = nf90_put_var(ncid0,VarID, snowbins(isa:fsa), start = (/nrec0/))
            END IF

         END IF
      END IF

      IF (level < 4) THEN
         iret = nf90_inq_varid(ncid0, sanal(8), VarID)
         iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(9), VarID)
         iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, sanal(10), VarID)
         iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))

         iret = nf90_inq_varid(ncid0, sanal(11), VarID)
         iret = nf90_put_var(ncid0, VarID, a_up(:,i1:i2,j1:j2), start=ibeg,    &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, sanal(12), VarID)
         iret = nf90_put_var(ncid0, VarID, a_vp(:,i1:i2,j1:j2), start=ibeg,    &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, sanal(13), VarID)
         iret = nf90_put_var(ncid0, VarID, a_wp(:,i1:i2,j1:j2), start=ibeg,    &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, sanal(14), VarID)
         iret = nf90_put_var(ncid0, VarID, a_theta%data(:,i1:i2,j1:j2), start=ibeg, &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, sanal(15), VarID)
         iret = nf90_put_var(ncid0, VarID, a_press%data(:,i1:i2,j1:j2), start=ibeg, &
                             count=icnt)

      ELSE IF (level >= 4) THEN
         iret = nf90_inq_varid(ncid0, 'u0', VarID)
         iret = nf90_put_var(ncid0, VarID, u0, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, 'v0', VarID)
         iret = nf90_put_var(ncid0, VarID, v0, start = (/nrec0/))
         iret = nf90_inq_varid(ncid0, 'dn0', VarID)
         iret = nf90_put_var(ncid0, VarID, dn0, start = (/nrec0/))

         iret = nf90_inq_varid(ncid0, 'u', VarID)
         iret = nf90_put_var(ncid0, VarID, a_up(:,i1:i2,j1:j2), start=ibeg,    &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, 'v', VarID)
         iret = nf90_put_var(ncid0, VarID, a_vp(:,i1:i2,j1:j2), start=ibeg,    &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, 'w', VarID)
         iret = nf90_put_var(ncid0, VarID, a_wp(:,i1:i2,j1:j2), start=ibeg,    &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, 'theta', VarID)
         iret = nf90_put_var(ncid0, VarID, a_theta%data(:,i1:i2,j1:j2), start=ibeg, &
                             count=icnt)
         iret = nf90_inq_varid(ncid0, 'p', VarID)
         iret = nf90_put_var(ncid0, VarID, a_press%data(:,i1:i2,j1:j2), start=ibeg, &
                             count=icnt)

      END IF


      IF (level < 4) THEN ! Normal operation for levels 1-3

         nn = nbase
         IF (level >= 2)  THEN
            nn = nn+1
            iret = nf90_inq_varid(ncid0, 'l', VarID)
            iret = nf90_put_var(ncid0, VarID, a_rc%data(:,i1:i2,j1:j2), start=ibeg, &
                                count=icnt)
         END IF

         DO n = 2, nscl
            nn = nn+1
            CALL newsclr(n)
            iret = nf90_inq_varid(ncid0, sanal(nn), VarID)
            iret = nf90_put_var(ncid0,VarID,a_sp(:,i1:i2,j1:j2), start=ibeg,   &
                                count=icnt)
         END DO

         IF (isgstyp > 1)  THEN
            nn = nn+1
            iret = nf90_inq_varid(ncid0, 'stke', VarID)
            iret = nf90_put_var(ncid0, VarID, a_qp%data(:,i1:i2,j1:j2), start=ibeg, &
                                count=icnt)
         END IF

         IF (iradtyp > 1)  THEN
            nn = nn+1
            iret = nf90_inq_varid(ncid0, 'rflx', VarID)
            iret = nf90_put_var(ncid0, VarID, a_rflx(:,i1:i2,j1:j2), start=ibeg, &
                                count=icnt)
         END IF

         IF (nn /= nvar0) THEN
            IF (myid == 0) PRINT *, 'ABORTING:  Anal write error'
            CALL appl_abort(0)
         END IF

      ELSE IF (level >= 4) THEN ! Operation with SALSA

         IF (.TRUE.) THEN

            ! Relative humidity
            iret = nf90_inq_varid(ncid0,'S_RH',VarID)
            iret = nf90_put_var(ncid0,VarID,a_rh%data(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            IF ( level == 5) THEN !
               ! Relative humidity
               iret = nf90_inq_varid(ncid0,'S_RHI',VarID)
               iret = nf90_put_var(ncid0,VarID,a_rhi%data(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
            END IF
       
            ! Total water mixing ratio
            zvar(:,:,:) = a_rp%data + &   ! Water vapor
                          a_rc%data + &   ! Liquid water
                          a_srp%data      ! Rain
       
            iret = nf90_inq_varid(ncid0,'q',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            ! Liquid water mixing ratio
            zvar(:,:,:) = a_rc%data
            iret = nf90_inq_varid(ncid0,'l',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            ! Rain water mixing ratio
            zvar(:,:,:) = a_srp%data
            iret = nf90_inq_varid(ncid0,'r',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
     
            IF ( level == 5) THEN !
          
               ! Total ice mixing ratio
               zvar(:,:,:) = a_ri%data + & ! Ice cloud content
                  a_srs%data      ! Snow
               iret = nf90_inq_varid(ncid0,'f',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! Ice mixing ratio
               zvar(:,:,:) = a_ri%data
               iret = nf90_inq_varid(ncid0,'i',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! Snow ice mixing ratio
               zvar(:,:,:) = a_srs%data
               iret = nf90_inq_varid(ncid0,'s',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
            END IF

            ! Total number of cloud droplets
            CALL bulkNumc('cloud','a',zvar(:,:,:))
            zsum = zvar
            CALL bulkNumc('cloud','b',zvar(:,:,:))
            zsum = zsum + zvar
            iret = nf90_inq_varid(ncid0,'S_Nc',VarID)
            iret = nf90_put_var(ncid0,VarID,zsum(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            IF (lbinanl) THEN
               ! Cloud droplet size distribution (regime A)
               iret = nf90_inq_varid(ncid0,'S_Ncba',VarID)
               iret = nf90_put_var(ncid0,VarID,a_ncloudp%data(:,i1:i2,j1:j2,ica%cur:fca%cur), &
                                   start=ibegsd,count=icntcla)
          
               ! Cloud droplet size distribution (regime B)
               iret = nf90_inq_varid(ncid0,'S_Ncbb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,a_ncloudp%data(:,i1:i2,j1:j2,icb%cur:fcb%cur), &
                                   start=ibegsd,count=icntclb)
            END IF
       
            ! Mean cloud droplet wet radius (regime A)
            CALL meanRadius('cloud','a',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Rwca',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            ! Mean cloud droplet wet radius (regime B)
            CALL meanRadius('cloud','b',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Rwcb',VarID)
            IF (iret==NF90_NOERR) &
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            IF (lbinanl) THEN
               ! Cloud droplet bin wet radius (regime A)
               CALL getBinRadius(ncld,nspec+1,a_ncloudp%data,a_mcloudp%data,nlim,a_Rcwet)
               iret = nf90_inq_varid(ncid0,'S_Rwcba',VarID)
               iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,ica%cur:fca%cur), &
                                   start=ibegsd,count=icntcla)
          
               ! Cloud droplet bin wet radius (regime B)
               iret = nf90_inq_varid(ncid0,'S_Rwcbb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,a_Rcwet(:,i1:i2,j1:j2,icb%cur:fcb%cur), &
                                   start=ibegsd,count=icntclb)
            END IF
       
            ! Number of rain droplets
            CALL bulkNumc('precp','a',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Np',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            IF (lbinanl) THEN
               ! Precipitation size distribution
               iret = nf90_inq_varid(ncid0,'S_Npba',VarID)
               iret = nf90_put_var(ncid0,VarID,a_nprecpp%data(:,i1:i2,j1:j2,ira:fra), &
                                   start=ibegsd,count=icntpra)
            END IF
       
            ! Mean rain drop wet radius
            CALL meanRadius('precp','a',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Rwpa',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
                                count=icnt)
       
            IF (lbinanl) THEN
               ! Rain drop bin wet radius
               CALL getBinRadius(nprc,nspec+1,a_nprecpp%data,a_mprecpp%data,prlim,a_Rpwet)
               iret = nf90_inq_varid(ncid0,'S_Rwpba',VarID)
               iret = nf90_put_var(ncid0,VarID,a_Rpwet(:,i1:i2,j1:j2,ira:fra),  &
                                   start=ibegsd,count=icntpra)
            END IF
       
            ! Number of soluble aerosols (regime A)
            CALL bulkNumc('aerosol','a',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Na',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            ! Number of insoluble aerosols (regime B)
            CALL bulkNumc('aerosol','b',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Nb',VarID)
            IF (iret==NF90_NOERR) &
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)
       
            IF (lbinanl) THEN
               ! Aerosol size distribution (regime A)
               iret = nf90_inq_varid(ncid0,'S_Naba',VarID)
               iret = nf90_put_var(ncid0,VarId,a_naerop%data(:,i1:i2,j1:j2,in1a:fn2a), &
                                   start=ibegsd,count=icntaea)
          
               !Aerosol size distribution (regime B)
               iret = nf90_inq_varid(ncid0,'S_Nabb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,a_naerop%data(:,i1:i2,j1:j2,in2b:fn2b), &
                                   start=ibegsd,count=icntaeb)
            END IF
       
            ! Mean aerosol wet radius (regime A)
            CALL meanRadius('aerosol','a',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Rwaa',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
                                count=icnt)
       
            ! Mean aerosol wet radius (regime B)
            CALL meanRadius('aerosol','b',zvar(:,:,:))
            iret = nf90_inq_varid(ncid0,'S_Rwab',VarID)
            IF (iret==NF90_NOERR) &
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
                                count=icnt)
       
            IF (lbinanl) THEN
               ! Aerosol bin wet radius (regime A)
               CALL getBinRadius(nbins,nspec+1,a_naerop%data,a_maerop%data,nlim,a_Rawet)
               iret = nf90_inq_varid(ncid0,'S_Rwaba',VarID)
               iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in1a:fn2a),  &
                                   start=ibegsd,count=icntaea)
          
               ! Aerosol bin wet radius (regime B)
               iret = nf90_inq_varid(ncid0,'S_Rwabb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,a_Rawet(:,i1:i2,j1:j2,in2b:fn2b),  &
                                   start=ibegsd,count=icntaeb)
            END IF
       
            IF (level == 5) THEN
               ! Number of ice particles
               CALL bulkNumc('ice','a',zvar(:,:,:))
               zsum = zvar
               CALL bulkNumc('ice','b',zvar(:,:,:))
               zsum = zsum + zvar
               iret = nf90_inq_varid(ncid0,'S_Nic',VarID)
               iret = nf90_put_var(ncid0,VarID,zsum(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! Number of snow droplets
               CALL bulkNumc('snow','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_Ns',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! Mean ice particle wet radius (a)
               CALL meanRadius('ice','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_Rwia',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! Mean ice particle wet radius (b)
               CALL meanRadius('ice','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_Rwib',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! Mean snow drop wet radius
               CALL meanRadius('snow','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_Rwsa',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg,  &
                                   count=icnt)
          
               IF (lbinanl) THEN
                  ! Ice particle size distribution reg. a
                  iret = nf90_inq_varid(ncid0,'S_Niba',VarID)
                  iret = nf90_put_var(ncid0,VarID,a_nicep%data(:,i1:i2,j1:j2,iia%cur:fia%cur), &
                                      start=ibegsd,count=icntica)
             
                  ! Ice particle size distribution reg. b
                  iret = nf90_inq_varid(ncid0,'S_Nibb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,a_nicep%data(:,i1:i2,j1:j2,iib%cur:fib%cur), &
                                      start=ibegsd,count=icntica)
             
                  ! Ice particle bin wet radius regime a
                  CALL getBinRadius(nice,nspec+1,a_nicep%data,a_micep%data,prlim,a_Riwet)
                  iret = nf90_inq_varid(ncid0,'S_Rwiba',VarID)
                  iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,iia%cur:fia%cur), &
                                      start=ibegsd,count=icntica)
             
                  ! Ice particle bin wet radius regime b
                  iret = nf90_inq_varid(ncid0,'S_Rwibb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,a_Riwet(:,i1:i2,j1:j2,iib%cur:fib%cur), &
                                      start=ibegsd,count=icntica)
             
                  ! Snow size distribution
                  iret = nf90_inq_varid(ncid0,'S_Nsba',VarID)
                  iret = nf90_put_var(ncid0,VarID,a_nsnowp%data(:,i1:i2,j1:j2,isa:fsa), &
                                      start=ibegsd,count=icntsna)
             
                  ! Snow drop bin wet radius
                  CALL getBinRadius(nsnw,nspec+1,a_nsnowp%data,a_msnowp%data,prlim,a_Rswet)
                  iret = nf90_inq_varid(ncid0,'S_Rwsba',VarID)
                  iret = nf90_put_var(ncid0,VarID,a_Rswet(:,i1:i2,j1:j2,isa:fsa),  &
                                      start=ibegsd,count=icntsna)
             
               END IF !(lbinanl)

            END IF !level 5

            ! Number of newly activated
            zvar(:,:,:) = SUM(a_nactd(:,:,:,:), DIM=4)
            iret = nf90_inq_varid(ncid0,'S_Nact',VarID)
            iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                count=icnt)

            ! Mass mixing ratios
            IF (IsUsed(prtcl,'SO4')) THEN

               ! --Sulphate (aerosol, regime A)
               CALL bulkMixrat('SO4','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aSO4a',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)

               ! --Sulphate (aerosol, regime B)
               CALL bulkMixrat('SO4','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aSO4b',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! --Sulphate (clouds, regime A)
               CALL bulkMixrat('SO4','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cSO4a',VarID)
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               ! --Sulphate (clouds, regime B)
               CALL bulkMixrat('SO4','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cSO4b',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)

               IF (level == 5) THEN
                  ! --Sulphate (ice, regime A)
                  CALL bulkMixrat('SO4','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iSO4a',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! --Sulphate (ice, regime B)
                  CALL bulkMixrat('SO4','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iSO4b',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5

            END IF

            IF (IsUSed(prtcl,'NH')) THEN

               !-- Ammonium (aerosol, regime A)
               CALL bulkMixrat('NH','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aNH3a',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Ammonium (aerosol, regime B)
               CALL bulkMixrat('NH','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aNH3b',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Ammonium (clouds, regime A)
               CALL bulkMixrat('NH','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cNH3a',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Ammonium (clouds, regime B)
               CALL bulkMixrat('NH','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cNH3b',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)

               IF (level == 5) THEN
                  ! --Ammonium (ice, regime A)
                  CALL bulkMixrat('NH','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iNHa',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! --Ammonium (ice, regime B)
                  CALL bulkMixrat('NH','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iNHb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5

            END IF

            IF (IsUsed(prtcl,'NO')) THEN

               !-- Nitrate (aerosol, regime A)
               CALL bulkMixrat('NO','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aNO3a',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Nitrate (aerosol, regime B)
               CALL bulkMixrat('NO','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aNO3b',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Nitrate (clouds, regime A)
               CALL bulkMixrat('NO','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cNO3a',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Nitrate (clouds, regime B)
               CALL bulkMixrat('NO','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cNO3b',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               IF (level == 5) THEN
                  ! --Nitrate (ice, regime A)
                  CALL bulkMixrat('NO','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iNOa',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! --Nitrate (ice, regime B)
                  CALL bulkMixrat('NO','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iNOb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5
          
            END IF

            IF (IsUsed(prtcl,'OC')) THEN
        
               !-- Organic Carbon (aerosol, regime A)
               CALL bulkMixrat('OC','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aOCa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Organic Carbon (aerosol, regime B)
               CALL bulkMixrat('OC','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aOCb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Organic Carbon (clouds, regime A)
               CALL bulkMixrat('OC','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cOCa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Organic Carbon (clouds, regime B)
               CALL bulkMixrat('OC','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cOCb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               IF (level == 5) THEN
                  ! --Organic Carbon (ice, regime A)
                  CALL bulkMixrat('OC','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iOCa',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! --Organic Carbon (ice, regime B)
                  CALL bulkMixrat('OC','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iOCb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5

            END IF
       
            IF (IsUsed(prtcl,'BC')) THEN

               !-- Black Carbon (aerosol, regime A)
               CALL bulkMixrat('BC','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aBCa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Black Carbon (aerosol, regime B)
               CALL bulkMixrat('BC','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aBCb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Black Carbon (clouds, regime A)
               CALL bulkMixrat('BC','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cBCa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Black Carbon (clouds, regime B)
               CALL bulkMixrat('BC','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cBCb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)

               IF (level == 5) THEN
                  ! --Black Carbon (ice, regime A)
                  CALL bulkMixrat('BC','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iBCa',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! --Black Carbon (ice, regime B)
                  CALL bulkMixrat('BC','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iBCb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5

            END IF

            IF (IsUsed(prtcl,'DU')) THEN

               !-- Dust (aerosol, regime A)
               CALL bulkMixrat('DU','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aDUa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Dust (aerosol, regime B)
               CALL bulkMixrat('DU','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aDUb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Dust (clouds, regime A)
               CALL bulkMixrat('DU','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cDUa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Dust (clouds, regime B)
               CALL bulkMixrat('DU','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cDUb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               IF (level == 5) THEN
                  ! --Dust (ice, regime A)
                  CALL bulkMixrat('DU','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iDUa',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! --Dust (ice, regime B)
                  CALL bulkMixrat('DU','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iDUb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5

            END IF
       
            IF (IsUsed(prtcl,'SS')) THEN

               !-- Sea Salt (aerosol, regime A)
               CALL bulkMixrat('SS','aerosol','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aSSa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Sea Salt (aerosol, regime B)
               CALL bulkMixrat('SS','aerosol','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_aSSb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Sea Salt (aerosol, regime A)
               CALL bulkMixrat('SS','cloud','a',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cSSa',VarID)
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               !-- Sea Salt (aerosol, regime B)
               CALL bulkMixrat('SS','cloud','b',zvar(:,:,:))
               iret = nf90_inq_varid(ncid0,'S_cSSb',VarID)
               IF (iret==NF90_NOERR) &
               iret = nf90_put_var(ncid0,Varid,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                   count=icnt)
          
               IF (level == 5) THEN
                  ! -- Sea Salt (ice, regime A)
                  CALL bulkMixrat('SS','ice','a',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iSSa',VarID)
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
             
                  ! -- Sea Salt (ice, regime B)
                  CALL bulkMixrat('SS','ice','b',zvar(:,:,:))
                  iret = nf90_inq_varid(ncid0,'S_iSSb',VarID)
                  IF (iret==NF90_NOERR) &
                  iret = nf90_put_var(ncid0,VarID,zvar(:,i1:i2,j1:j2),start=ibeg, &
                                      count=icnt)
               END IF ! level 5
          
            END IF

         END IF
       
      END IF

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

      WRITE(10) time,th00,umean,vmean,dtl,level,isgstyp,iradtyp,nzp,nxp,nyp,nscl
      WRITE(10) xt, xm, yt, ym, zt, zm, dn0, th0, u0, v0, pi0, pi1, rt0, psrf,sst,W1,W2,W3 ! added by Zubair

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
      IF ( allocated(a_rflx) ) WRITE(10) a_rflx
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
   ! Subroutine binMixrat: Calculate the total dry or wet
   ! Mass concentration for individual bins
   !
   ! Juha Tonttila, FMI, 2015
   ! Tomi Raatikainen, FMI, 2016
   SUBROUTINE binMixrat(ipart,itype,ibin,ii,jj,kk,sumc)
      USE mo_submctl, ONLY : ncld,nbins,nprc,nice,nsnw
      USE class_ComponentIndex, ONLY : GetNcomp
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(in) :: ipart
      CHARACTER(len=*), INTENT(in) :: itype
      INTEGER, INTENT(in) :: ibin,ii,jj,kk
      REAL, INTENT(out) :: sumc

      INTEGER :: iend

      ! Number of components-1
      IF (itype == 'dry') THEN
         iend = GetNcomp(prtcl)-1 ! dry CASE
      ELSE IF (itype == 'wet') THEN
         iend = GetNcomp(prtcl) ! wet CASE
      ELSE
         STOP 'Error in binMixrat!'
      END IF

      SELECT CASE(ipart)
         CASE('aerosol')
            sumc = SUM( a_maerop%data(kk,ii,jj,ibin:iend*nbins+ibin:nbins) )
         CASE('cloud')
            sumc = SUM( a_mcloudp%data(kk,ii,jj,ibin:iend*ncld+ibin:ncld) )
         CASE('precp')
            sumc = SUM( a_mprecpp%data(kk,ii,jj,ibin:iend*nprc+ibin:nprc) )
         CASE('ice')
            sumc = SUM( a_micep%data(kk,ii,jj,ibin:iend*nice+ibin:nice) )
         CASE('snow')
            sumc = SUM( a_msnowp%data(kk,ii,jj,ibin:iend*nsnw+ibin:nsnw) )
         CASE DEFAULT
            STOP 'bin mixrat error'
      END SELECT

   END SUBROUTINE binMixrat

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

            CALL getRadius(istr,iend,nbins,nspec+1,a_naerop%data,a_maerop%data,nlim,rad)

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

            CALL getRadius(istr,iend,ncld,nspec+1,a_ncloudp%data,a_mcloudp%data,nlim,rad)

         CASE('precp')

            istr = ira
            iend = fra

            CALL getRadius(istr,iend,nprc,nspec+1,a_nprecpp%data,a_mprecpp%data,prlim,rad)

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

            CALL getRadius(istr,iend,nice,nspec+1,a_nicep%data,a_micep%data,prlim,rad)

         CASE('snow')

            istr = isa
            iend = fsa

            CALL getRadius(istr,iend,nsnw,nspec+1,a_nsnowp%data,a_msnowp%data,prlim,rad)

      END SELECT

   END SUBROUTINE meanRadius
   !
   ! ---------------------------------------------------
   ! Subroutine getRadius
   ! Calculates number mean wet radius (over selected bins) for the whole domain
   SUBROUTINE getRadius(zstr,zend,nn,n4,numc,mass,numlim,zrad)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: nn, n4 ! Number of bins (nn) and aerosol species (n4)
      INTEGER, INTENT(in) :: zstr,zEND  ! Start and end index for averaging
      REAL, INTENT(in) :: numc(nzp,nxp,nyp,nn)
      REAL, INTENT(in) :: mass(nzp,nxp,nyp,nn*n4)
      REAL, INTENT(in) :: numlim

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
                     rwet = rwet+calc_wet_radius(n4,numc(k,i,j,bin),tmp)*numc(k,i,j,bin)
                  END IF
               END DO

               IF (tot > numlim) THEN
                  zrad(k,i,j) = rwet/tot
               END IF

            END DO
         END DO
      END DO


   END SUBROUTINE getRadius

   SUBROUTINE getBinRadius(nn,n4,numc,mass,numlim,zrad)
      IMPLICIT NONE

      INTEGER, INTENT(in) :: nn, n4 ! Number of bins (nn) and aerosol species (n4)
      REAL, INTENT(in)  :: numc(nzp,nxp,nyp,nn)
      REAL, INTENT(in)  :: mass(nzp,nxp,nyp,nn*n4)
      REAL, INTENT(in)  :: numlim
      REAL, INTENT(out) :: zrad(nzp,nxp,nyp,nn)

      INTEGER :: k,i,j,bin
      REAL :: tmp(n4)

      zrad(:,:,:,:) = 0.
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 1, nzp

               DO bin = 1, nn
                  IF (numc(k,i,j,bin) > numlim) THEN
                     tmp(:) = mass(k,i,j,bin:(n4-1)*nn+bin:nn)
                     zrad(k,i,j,bin) = calc_wet_radius(n4,numc(k,i,j,bin),tmp)
                  END IF
               END DO

            END DO
         END DO
      END DO



   END SUBROUTINE getBinRadius

   ! Aerosol and cloud droplet composition
   !     Needs local density array dens
   REAL FUNCTION calc_wet_radius(n,numc,mass)
      USE mo_submctl, ONLY : pi6
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL, INTENT(IN)    :: numc, mass(n)

      calc_wet_radius = 0.

      ! Don't calculate if very low number concentration
      IF (numc < 1e-15) RETURN

      ! Radius from total volume per particle
      calc_wet_radius = 0.5*( SUM(mass(:)/dens(:))/numc/pi6)**((1./3.))

   END FUNCTION calc_wet_radius

END MODULE grid

