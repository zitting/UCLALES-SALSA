!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributsed in the hope that it will be useful,
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
MODULE stat

   USE classFieldArray
   USE ncio, ONLY : open_nc, define_nc, define_nc_cs
   USE grid, ONLY : level, lbinprof, lbinanl
   USE util, ONLY : get_avg3, get_cor3, get_var3, get_avg_ts, &
                    get_avg2dh, get_3rd3

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: nvar1 = 29,               &
                         nv1sbulk = 48,            &
                         nv1MB = 4,                &
                         nv1_lvl5 = 32,            &
                         nvar2 = 96,               &
                         nv2sbulk = 49,            &
                         nv2_lvl5 = 29,            &
                         nv2saa = 8, nv2sab = 8,   &
                         nv2sca = 8, nv2scb = 8,   &
                         nv2sp = 8,                 &
                         nv2sia = 8, nv2sib = 8,   &
                         nv2ss = 8

   ! All SALSA species
   CHARACTER(len=3), PARAMETER :: zspec(8) = (/'SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O'/)
   ! Active SALSA species
   CHARACTER (len=3), SAVE :: actspec(8)
   INTEGER, SAVE      :: nspec=0

   INTEGER, SAVE      :: nrec1, nrec2, nrec3, ncid1, ncid2, ncid3, nv1=nvar1, nv2 = nvar2
   REAL, SAVE         :: fsttm, lsttm, nsmp = 0
   REAL, SAVE         :: avgtime = 0

   ! Field arrays containing output variables
   TYPE(FieldArray) :: fields1,fields2
   TYPE(FieldArray) :: f1,f1sbulk,f1lvl5,f2,f2sbulk,f2lvl5,f2aea,f2aeb,f2cla,f2clb,f2prc,f2ica,f2icb,f2snw

   !Float arrays containing output variable data
   TYPE(FloatArray1D), TARGET  :: &
      a_time, a_cfl, a_maxdiv, a_zi1_bar, a_zi2_bar, a_zi3_bar, &
      a_vtke, a_sfcbflx, a_wmax, a_tsrf, a_ustar, a_shf_bar,    &
      a_lhf_bar, a_zi_bar, a_lwp_bar, a_lwp_var, a_zc, a_zb,    &
      a_cfrac, a_lmax, a_albedo, a_rwp_bar, a_prcp, a_pfrac,    &
      a_CCN, a_nrain, a_nrcnt, a_nccnt, a_prcp_bc

   TYPE(FloatArray1D), TARGET  :: &
      a_Nc_ic,   a_Na_int,  a_Na_oc,                                  &
      a_SO4_ic,  a_SO4_int, a_SO4_oc,                                 &
      a_OC_ic,   a_OC_int,  a_OC_oc,                                  &
      a_BC_ic,   a_BC_int,  a_BC_oc,                                  &
      a_DU_ic,   a_DU_int,  a_DU_oc,                                  &
      a_SS_ic,   a_SS_int,  a_SS_oc,                                  &
      a_NH_ic,   a_NH_int,  a_NH_oc,                                  &
      a_NO_ic,   a_NO_int,  a_NO_oc,                                  &
      a_rmH2Oae, a_rmH2Ocl, a_rmH2Opr,                                &
      a_rmSO4dr, a_rmSO4cl, a_rmSO4pr, a_rmSO4wt, a_rmSO4tt,          &
      a_rmOCdr,  a_rmOCcl,  a_rmOCpr,  a_rmOCwt,  a_rmOCtt,           &
      a_rmBCdr,  a_rmBCcl,  a_rmBCpr,  a_rmBCwt,  a_rmBCtt,           &
      a_rmDUdr,  a_rmDUcl,  a_rmDUpr,  a_rmDUwt,  a_rmDUtt,           &
      a_rmSSdr,  a_rmSScl,  a_rmSSpr,  a_rmSSwt,  a_rmSStt,           &
      a_rmNOdr,  a_rmNOcl,  a_rmNOpr,  a_rmNOwt,  a_rmNOtt,           &
      a_rmNHdr,  a_rmNHcl,  a_rmNHpr,  a_rmNHwt,  a_rmNHtt

   TYPE(FloatArray1D), TARGET :: &
      a_Ni_ic, a_Ni_ii,   a_Ni_is, a_Ns_ic, a_Ns_ii, a_Ns_is, &
      a_Ri_ii, a_iwp_bar, a_imax, a_nicnt, &
      a_Rs_is, a_swp_bar, a_smax, a_nscnt, &
      a_rmSO4ic, a_rmSO4sn, a_rmOCic,  a_rmOCsn, &
      a_rmBCic,  a_rmBCsn,  a_rmDUic,  a_rmDUsn, &
      a_rmNOic,  a_rmNOsn,  a_rmNHic,  a_rmNHsn, &
      a_rmSSic,  a_rmSSsn,  a_rmH2Oic, a_rmH2Osn, &
      a_sfrac, a_sprcp

   TYPE(FloatArray1D), TARGET  :: &
      a_time2,   a_zt,      a_zm,       a_dn0,     a_u0,      a_v0,      &
      a_fsttm,   a_lsttm,   a_nsmp,     a_u,       a_v,       a_theta_1, &
      a_p,       a_u_2,     a_v_2,      a_w_2,     a_theta_2, a_w_3,     &
      a_theta_3, a_tot_tw,  a_sfs_tw,   a_tot_uw,  a_sfs_uw,  a_tot_vw,  &
      a_sfs_vw,  a_tot_ww,  a_sfs_ww,   a_km,      a_kh,      a_lmbd,    &
      a_lmbde,   a_sfs_tke, a_sfs_boy,  a_sfs_shr, a_boy_prd, a_shr_prd, &
      a_trans,   a_diss,    a_dff_u,    a_dff_v,   a_dff_w,   a_adv_u,   &
      a_adv_v,   a_adv_w,   a_prs_u,    a_prs_v,   a_prs_w,   a_prd_uw,  &
      a_storage, a_q,       a_q_2,      a_q_3,     a_tot_qw,  a_sfs_qw,  &
      a_rflx1,   a_rflx2,   a_sflx1,    a_sflx2,   a_l,       a_l_2,     &
      a_l_3,     a_tot_lw,  a_sed_lw,   a_cs1,     a_cnt_cs1, a_w_cs1,   &
      a_t_cs1 ,  a_tv_cs1,  a_rt_cs1,   a_rc_cs1,  a_wt_cs1,  a_wtv_cs1, &
      a_wrt_cs1, a_cs2,     a_cnt_cs2,  a_w_cs2,   a_t_cs2,   a_tv_cs2,  &
      a_rt_cs2,  a_rc_cs2,  a_wt_cs2,   a_wtv_cs2, a_wrt_cs2, a_Nc,      &
      a_Nr,      a_rr,      a_rrate,    a_evap,    a_frc_prc, a_prc_prc, &
      a_frc_ran, a_hst_srf, a_sw_up,    a_sw_down, a_lw_up,   a_lw_down


   TYPE(FloatArray1D), TARGET  :: &
      a_aea,   a_aeb,   a_cla,  a_clb,  a_prc,          &
      a_Naa,   a_Nab,   a_Nca,  a_Ncb,  a_Np,           &
      a_Rwaa,  a_Rwab,  a_Rwca, a_Rwcb, a_Rwp,          &
      a_cSO4a, a_cSO4c, a_cSO4p,                        &
      a_cOCa,  a_cOCc,  a_cOCp,                         &
      a_cBCa,  a_cBCc,  a_cBCp,                         &
      a_cDUa,  a_cDUc,  a_cDUp,                         &
      a_cSSa,  a_cSSc,  a_cSSp,                         &
      a_cNOa,  a_cNOc,  a_cNOp,                         &
      a_cNHa,  a_cNHc,  a_cNHp,                         &
      a_cH2Oa, a_cH2Oc, a_cH2Op,                        &
      a_prl, a_Prr, a_prv, a_pRH,                       &
      a_Na_c, a_Nc_c, a_Np_c, a_Pcfrac, a_clw_c, a_thl_c

   TYPE(FloatArray1D), TARGET :: &
      a_ica,   a_icb,   a_snw, &
      a_Nia,   a_Nib,   a_Ns,    a_Rwia, a_Rwib, a_Rws, &
      a_cSO4i, a_cSO4s, a_cOCi,  a_cOCs,  &
      a_cBCi,  a_cBCs,  a_cDUi,  a_cDUs,  &
      a_cSSi,  a_cSSs,  a_cNOi,  a_cNOs,  &
      a_cNHi,  a_cNHs,  a_cH2Oi, a_cH2Os, &
      a_pri,   a_prs,   a_pRHi,  a_srate

   TYPE(FloatArray2D), TARGET  :: &
      a_Naba, a_SO4aa, a_OCaa, a_BCaa,      &
      a_DUaa, a_SSaa,  a_NOaa, a_NHaa,      &
      a_Nabb, a_SO4ab, a_OCab, a_BCab,      &
      a_DUab, a_SSab,  a_NOab, a_NHab,      &
      a_Ncba, a_SO4ca, a_OCca, a_BCca,      &
      a_DUca, a_SSca,  a_NOca, a_NHca,      &
      a_Ncbb, a_SO4cb, a_OCcb, a_BCcb,      &
      a_DUcb, a_SScb,  a_NOcb, a_NHcb,      &
      a_Npb,  a_SO4pb, a_OCpb, a_BCpb,      &
      a_DUpb, a_SSpb,  a_NOpb, a_NHpb

   TYPE(FloatArray2D), TARGET :: &
      a_Niba, a_SO4ia, a_OCia, a_BCia,      &
      a_DUia, a_SSia,  a_NOia, a_NHia,      &
      a_Nibb, a_SO4ib, a_OCib, a_BCib,      &
      a_DUib, a_SSib,  a_NOib, a_NHib,      &
      a_Nsb,  a_SO4sb, a_OCsb, a_BCsb,      &
      a_DUsb, a_SSsb,  a_NOsb, a_NHsb


   LOGICAL            :: sflg = .FALSE.
   LOGICAL            :: mcflg = .FALSE.
   LOGICAL            :: csflg = .FALSE.
   LOGICAL            :: salsa_b_bins = .FALSE.
   LOGICAL            :: cloudy_col_stats = .FALSE.
   REAL               :: ssam_intvl = 30.   ! statistical sampling interval
   REAL               :: savg_intvl = 1800. ! statistical averaging interval

   CLASS(*), POINTER  :: dummy,dummy2,dummy3
   ! New output variables:
   ! 1. Add them to one of the arrays below
   ! 2. Update the dimensions of the output arrays accordingly
   ! 3. Make a subroutine for accumulating the new data and to plane them
   !    in the output arrays
   ! 4. Add the new variables in ncio.f90 list of variables
   ! 5. Make sure stat_init is up to date (boolean arrays etc).

   CHARACTER (len=3), SAVE :: spec3(7)

   REAL, SAVE, ALLOCATABLE :: tke_sgs(:), tke_res(:), tke0(:), wtv_sgs(:),            &
                              wtv_res(:), wrl_sgs(:), thvar(:), svctr(:,:), ssclr(:), &
      ! Additional ssclr and svctr for BULK SALSA output
                              svctr_b(:,:), ssclr_b(:),                               &
      ! Additional ssclr and svctr for BINNED SALSA output.
                              svctr_aa(:,:,:), svctr_ca(:,:,:), svctr_p(:,:,:),       &
                              svctr_ab(:,:,:), svctr_cb(:,:,:),                       &
                              ! The same for SALSA level 5
       svctr_lvl5(:,:), ssclr_lvl5(:), &
       svctr_ia(:,:,:), svctr_ib(:,:,:), svctr_s(:,:,:), &
      ! Mass budget arrays
                              massbdg(:), scs_rm(:,:,:)

   PUBLIC :: sflg, ssam_intvl, savg_intvl, statistics, init_stat, write_ps,        &
             acc_tend, updtst, sfc_stat, close_stat, fill_scalar, tke_sgs, sgsflxs,&
             sgs_vel, comp_tke, get_zi, acc_removal, cs_rem_set, acc_massbudged, write_massbudged, &
             mcflg, csflg, salsa_b_bins, cloudy_col_stats

CONTAINS
   !
   ! ---------------------------------------------------------------------
   ! INIT_STAT:  This routine initializes the statistical arrays which
   ! are user/problem defined.  Note that svctr is given 100 elements, and
   ! elements 90 and above are used for computing the TKE budget. Hence
   ! if (nvar2 >= 90 the program stops
   !
   SUBROUTINE init_stat(time, filprf, expnme, nzp)

      USE classFieldArray
      USE mo_structured_datatypes
      USE grid,          ONLY : nxp, nyp, iradtyp, prtcl
      USE mpi_interface, ONLY : myid, ver, author, info
      USE mo_submctl,    ONLY : nprc, fn2a,fn2b,fca,fcb,fra,fia,fib,fsa
      USE class_ComponentIndex, ONLY : IsUsed

      CHARACTER (len=80), INTENT (in) :: filprf, expnme
      INTEGER, INTENT (in)            :: nzp
      REAL, INTENT (in)               :: time

      INTEGER :: i, e, n, nnum
      CHARACTER (len=80) :: fname
      CHARACTER (len=10) :: strn

      TYPE(ArrayElement), POINTER :: ArrEl

      ALLOCATE (wtv_sgs(nzp),wtv_res(nzp),wrl_sgs(nzp))
      ALLOCATE (tke_res(nzp),tke_sgs(nzp),tke0(nzp),thvar(nzp))

      wtv_sgs(:) = 0.
      wtv_res(:) = 0.
      wrl_sgs(:) = 0.
      tke_res(:) = 0.
      tke_sgs(:) = 0.
      tke0(:)    = 0.

      SELECT CASE(level)
         CASE (0)
            nv1 = 13
            nv2 = 58
         CASE (1)
            nv1 = 14
            nv2 = 58
         CASE (2)
            nv1 = 20
            nv2 = 83
            IF (iradtyp == 3) nv1 = 21
         ! For SALSA
         CASE (4,5)
            nv1 = nvar1
            nv2 = nvar2
         CASE DEFAULT

            nv1 = nvar1
            nv2 = nvar2
      END SELECT

      CALL initilize_arrays(nzp)
      CALL create_fields()

      nspec = nspec+1
      actspec(nspec)=zspec(e)

      IF ( level < 4 ) THEN
         IF (mcflg) &
            ALLOCATE ( massbdg(nv1MB) ) ! Mass budged array; ALL CALCULATIONS ARE NOT IMPLEMENTED FOR LEVEL < 4
         IF (mcflg) massbdg(:) = 0.

      ELSE IF ( level >= 4 ) THEN

         IF (mcflg) THEN
            ALLOCATE ( massbdg(nv1MB) ) ! Mass budged array
            massbdg(:) = 0.
         END IF

         IF (csflg .AND. level > 3) THEN
            ! Allocate array for level 4 removal rate column statistics
            ! Total number of outputs is 3 for warm (aerosol, cloud and precipitation)
            ! and 5 (add ice and snow) for each species including water
            IF (level == 4) THEN
               ALLOCATE( scs_rm(3*nspec,nxp,nyp) )
            ELSE
               ALLOCATE( scs_rm(5*nspec,nxp,nyp) )
            END IF
            scs_rm = 0. ! Set to zero (during spinup)
         END IF
      END IF ! If level >=4

      nnum = f1%count+f1sbulk%count+f1lvl5%count

      DO n = 1, nnum
      CALL fields1%NewField("","","","",.TRUE.,dummy)
      END DO

      DO n = 1, nnum

         IF (n <= f1%count) THEN
         fields1%list(n) = f1%list(n)
         ELSE IF (n <= f1%count+f1sbulk%count) THEN
         fields1%list(n) = f1sbulk%list(n-f1%count)
         ELSE IF (n > f1%count+f1sbulk%count) THEN
         fields1%list(n) = f1lvl5%list(n-f1%count-f1sbulk%count)
         END IF

         CALL fields1%getField(n,ArrEl)
         write(*,*) ArrEl%name,ArrEl%dimension,Arrel%outputstatus,n

      END DO


      nnum = f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count+f2prc%count+&
             f2lvl5%count+f2ica%count+f2icb%count+f2snw%count
      CALL fields2%ExN(nnum)

      DO n = 1, nnum

         IF (n <= f2%count) THEN
         fields2%list(n) = f2%list(n)

         ELSE IF (n <= f2%count+f2sbulk%count) THEN
         fields2%list(n) = f2sbulk%list(n-f2%count)

         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count) THEN
         fields2%list(n) = f2aea%list(n-f2%count-f2sbulk%count)

         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count) THEN
         fields2%list(n) = f2aeb%list(n-f2%count-f2sbulk%count-f2aea%count)

         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count) THEN
         fields2%list(n) = f2cla%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count)

         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count) THEN
         fields2%list(n) = f2clb%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count-f2cla%count)

         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count+f2prc%count) THEN
         fields2%list(n) = f2prc%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count-f2cla%count-f2clb%count)

         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count+f2prc%count&
                  +f2lvl5%count) THEN
         fields2%list(n) = f2lvl5%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count-f2cla%count-f2clb%count-f2prc%count)
         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count+f2prc%count&
                  +f2lvl5%count+f2ica%count) THEN
         fields2%list(n) = f2ica%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count-f2cla%count-f2clb%count-f2prc%count&
                           -f2lvl5%count)
         ELSE IF (n <= f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count+f2prc%count&
                  +f2lvl5%count+f2ica%count+f2icb%count) THEN
         fields2%list(n) = f2icb%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count-f2cla%count-f2clb%count-f2prc%count&
                           -f2lvl5%count-f2ica%count)
         ELSE IF (n > f2%count+f2sbulk%count+f2aea%count+f2aeb%count+f2cla%count+f2clb%count+f2prc%count&
                  +f2lvl5%count+f2ica%count+f2icb%count) THEN
         fields2%list(n) = f2snw%list(n-f2%count-f2sbulk%count-f2aea%count-f2aeb%count-f2cla%count-f2clb%count-f2prc%count&
                           -f2lvl5%count-f2ica%count-f2icb%count)
         END IF

         CALL fields2%getField(n,ArrEl)
         write(*,*) ArrEl%name,ArrEl%dimension,Arrel%outputstatus,n

      END DO

      nnum = fields1%count
      fname = trim(filprf)//'.ts'
      IF(myid == 0) PRINT                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
      CALL open_nc( fname, expnme, time, (nxp-4)*(nyp-4), ncid1, nrec1, ver, author, info)
      ! Juha: Modified for SALSA output
      CALL define_nc( ncid1, nrec1, nnum, fields1)
      IF (myid == 0) PRINT *, '   ...starting record: ', nrec1

      nnum = fields2%count
      fname = trim(filprf)//'.ps'
      IF(myid == 0) PRINT                                                  &
         "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
      CALL open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid2, nrec2, ver, author, info)
      ! Juha: Modified due to SALSA output
      CALL define_nc( ncid2, nrec2, nnum, fields2, n1=nzp, inae_a=fn2a, inae_b=fn2b-fn2a, &
                      incld_a=fca%cur, incld_b=fcb%cur-fca%cur, inice_a=fia%cur, inice_b=fib%cur-fia%cur, inprc=fra, insnw=fsa)

      IF (myid == 0) PRINT *, '   ...starting record: ', nrec2

      ! Optional column statistics
      IF (csflg) THEN
         fname = trim(filprf)//'.cs'
         IF(myid == 0) PRINT "(//' ',49('-')/,' ',/,'  Initializing: ',A20)",trim(fname)
         CALL open_nc( fname, expnme, time,(nxp-4)*(nyp-4), ncid3, nrec3, ver, author, info)
         IF (ncid3 >= 0) CALL define_nc_cs(ncid3, nrec3, nxp-4, nyp-4, level, iradtyp, actspec(1:nspec), nspec)
         IF (myid == 0) PRINT *, '   ...starting record: ', nrec3
      END IF

   END SUBROUTINE init_stat

   !
   ! ---------------------------------------------------------------------
   ! Subroutine Statistics:  This subroutine is the statistics driver
   ! it calls various other subroutines to compute and accumulate
   ! statistical quantities.  These are stored in two arrays:  SVCTR,
   ! and SSCLR (which accumulate scalar and vector statistics respectively
   !
   ! Modified for level 4: rxt contains a_rp if level < 4 and a_rp+a_rc
   ! if level == 4
   ! Juha Tonttila, FMI, 2014
   !
   ! Modified for level 5
   ! Jaakko Ahola, FMI, 2016
   SUBROUTINE statistics(time)

      USE grid, ONLY : a_up, a_vp, a_wp, a_rc, a_theta, a_rv,                                  &
                       a_rp, a_tp, a_press, nxp, nyp, nzp, dzm, dzt, zm, zt, th00, umean,      &
                       vmean, dn0, precip, a_rpp, a_npp, CCN, iradtyp, a_rflx, a_sflx,         &
                       a_fus, a_fds, a_fuir, a_fdir, albedo, a_srp, a_snrp, a_ncloudp, xt, yt, &
                       a_ri, a_nicep, a_srs, a_snrs, snowin


      REAL, INTENT (in) :: time

      REAL :: rxt(nzp,nxp,nyp), rnt(nzp,nxp,nyp), rxl(nzp,nxp,nyp), rxv(nzp,nxp,nyp)
      REAL :: xrpp(nzp,nxp,nyp), xnpp(nzp,nxp,nyp)
      INTEGER :: i
      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt

      SELECT CASE(level)
         CASE(1,2,3)
           rxt = a_rp%data ! Total water (vapor + condensed water and ice) = q
           rxl = a_rc%data-a_rpp%data ! Cloud water (+aerosol), but no precipitation or ice
           rxv = a_rv%data ! Water vapor
           xrpp = a_rpp%data ! Rain water
           xnpp = a_npp%data ! Rain number
         CASE(4,5)
           rxt = a_rp%data + a_rc%data + a_srp%data + a_ri%data + a_srs%data
           rxl = a_rc%data
           rxv = a_rp%data
           xrpp = a_srp%data
           xnpp = a_snrp%data
      END SELECT

      IF (nsmp == 0.) fsttm = time
      nsmp = nsmp+1.

      DO i=14,nvar1
         CALL fields1%getField(i,ArrEl)
         CALL fields1%getData(1,trgt,name=ArrEl%name)
         trgt%data = -999.
      END DO

      !
      ! profile statistics
      !
      CALL accum_stat(nzp, nxp, nyp, a_up, a_vp, a_wp, a_theta%data, a_press%data, umean, &
                      vmean)
      IF (iradtyp == 3) THEN
         CALL accum_rad(nzp, nxp, nyp, a_rflx, sflx=a_sflx, sup=a_fus, sdwn=a_fds, &
                        irup=a_fuir, irdwn=a_fdir, alb=albedo)
      ELSE IF (iradtyp > 0) THEN
         CALL accum_rad(nzp, nxp, nyp, a_rflx)
      END IF
      IF (level >= 1) CALL accum_lvl1(nzp, nxp, nyp, rxt)
      IF (level >= 2) CALL accum_lvl2(nzp, nxp, nyp, th00, dn0, zm, a_wp,        &
                                      a_theta%data, a_tp%data, rxv, rxl, rxt   )
      IF (level >= 3) CALL accum_lvl3(nzp, nxp, nyp, dn0, zm, rxl, xrpp,  &
                                      xnpp, precip, CCN                    )
      IF (level >= 4) CALL accum_lvl4(nzp, nxp, nyp)
      IF (level >= 5) CALL accum_lvl5(nzp, nxp, nyp, snowin)
       !for Salsa output in ps files .. by Zubair Maalick
      !
      ! scalar statistics
      !

      CALL set_ts(nzp, nxp, nyp, a_wp, a_theta%data, dn0, zt,zm,dzt,dzm,th00,time)
      IF ( level >= 1 ) CALL ts_lvl1(nzp, nxp, nyp, dn0, zt, dzm, rxt)
      IF ( level >= 2 ) CALL ts_lvl2(nzp, nxp, nyp, rxl, zt)
      IF ( level >= 4 ) CALL ts_lvl4(nzp, nxp, nyp, a_rc%data)
      IF ( level >= 5 ) CALL ts_lvl5(nzp, nxp, nyp, dn0, zt, a_rc%data, a_ri%data, a_srs%data, snowin)

      CALL write_ts

      !
      ! Column statistics
      !
      IF (csflg) THEN
         ! Radiation
         IF (iradtyp == 3) CALL set_cs_any(nxp,nyp,albedo,'albedo')

         ! Deposition statistics
         IF (level == 3) THEN
            CALL set_cs_any(nxp,nyp,precip(2,:,:),'prcp')
         ELSE IF (level > 3) THEN
            ! Surface deposition fluxes
            CALL cs_rem_save(nxp,nyp)
         END IF

         ! Warm cloud statistics
         rxt = 0.    ! Condensate
         rnt = 0.
         xrpp = 0.   ! Precipitate
         xnpp = 0.
         IF (level == 1) THEN
             ! No clouds or precipitation
         ELSE IF (level == 2) THEN
            ! Clouds available
            rxt = a_rc%data
            rnt = CCN
         ELSE IF (level == 3) THEN
            ! Clouds and precipitation available
            rxt = a_rc%data
            rnt = CCN
            xrpp = a_rpp%data
            xnpp = a_npp%data
         ELSE IF (level == 4 .OR. level == 5) THEN
            ! Levels 4 and 5
            rxt = a_rc%data
            rnt = SUM(a_ncloudp%data,DIM=4)
            xrpp = a_srp%data
            xnpp = a_snrp%data
         END IF

         CALL set_cs_warm(nzp,nxp,nyp,rxt,rnt,xrpp,xnpp,a_theta%data,dn0,zm,zt,dzm,xt,yt,time)

         ! Ice cloud statistics
         IF (level==5) THEN
             rxt = a_ri%data
             rnt = SUM(a_nicep%data,DIM=4)
             xrpp = a_srs%data
             xnpp = a_snrs%data
             CALL set_cs_cold(nzp,nxp,nyp,rxt,rnt,xrpp,xnpp,dn0,zm,xt,yt,time)
         END IF

         nrec3 = nrec3 + 1

      END IF

   END SUBROUTINE statistics
   !
   ! -----------------------------------------------------------------------
   ! Subroutines set_cs_warm, cs_rem_set, cs_rem_save and set_cs_any:
   ! write (and compute) column average statistics
   !
   ! Save named data (already available)
   SUBROUTINE set_cs_any(n2,n3,r,nam)
      USE netcdf
      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(IN)    :: r(n2,n3)
      CHARACTER (LEN=*) :: nam
      INTEGER :: iret, VarID

      IF (csflg) THEN
         iret = nf90_inq_varid(ncid3,trim(nam),VarID)
         IF (iret == NF90_NOERR) THEN
            iret = nf90_put_var(ncid3, VarID, r(3:n2-2,3:n3-2), start=(/1,1,nrec3/))
            iret = nf90_sync(ncid3)
         END IF
      END IF

   END SUBROUTINE set_cs_any
   !
   ! Removal statistics (level>3): calculate values for further use
   SUBROUTINE cs_rem_set(n2,n3,n4,raer,rcld,rprc,rice,rsnw)

      USE mo_submctl, ONLY : nbins, ncld, nprc, nice,  nsnw
      IMPLICIT NONE

      INTEGER, INTENT(in) :: n2,n3,n4   ! Grid dimensions
      REAL, INTENT(in) :: raer(n2,n3,n4*nbins), & ! Removal arrays
                          rcld(n2,n3,n4*ncld),  &
                          rprc(n2,n3,n4*nprc),  &
                          rice(n2,n3,n4*nice),  &
                          rsnw(n2,n3,n4*nsnw)

      INTEGER :: si, i, end,str

      IF (.NOT. csflg) RETURN

      ! Calculate all removal fluxes and save those to scs_rm for later use
      i = 1
      DO si = 1, nspec ! +1 for water
         ! Removal by sedimentation of aerosol
         str = (si-1)*nbins+1
         END = si*nbins
         scs_rm(i,:,:) = SUM(raer(:,:,str:end),DIM=3)
         i = i+1

         ! Removal by sedimentation of cloud droplets
         str = (si-1)*ncld+1
         END = si*ncld
         scs_rm(i,:,:) = SUM(rcld(:,:,str:end),DIM=3)
         i = i+1

         ! Removal by precipitation
         str = (si-1)*nprc+1
         END = si*nprc
         scs_rm(i,:,:) = SUM(rprc(:,:,str:end),DIM=3)
         i = i+1

         IF (level > 4) THEN
            ! Removal by sedimentation of ice particles
            str = (si-1)*nice+1
            END = si*nice
            scs_rm(i,:,:) = SUM(rice(:,:,str:end),DIM=3)
            i = i+1

            ! Removal by snow
            str = (si-1)*nsnw+1
            END = si*nsnw
            scs_rm(i,:,:) = SUM(rsnw(:,:,str:end),DIM=3)
            i = i+1
         END IF
      END DO

   END SUBROUTINE cs_rem_set
   !
   ! Removal statistics (level>3): save values
   SUBROUTINE cs_rem_save(n2,n3)

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n2, n3

      INTEGER :: si, i
      CHARACTER(LEN=3) :: nam

      IF (.NOT. csflg) RETURN

      ! Save all previously calculated removal fluxes
      !   Note: fluxes not calculated during spinup, so saving zeros
      i = 1
      DO si = 1, nspec

         nam=actspec(si)

         ! Removal by sedimentation of aerosol
         CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'dr') ! 'dr' should be for aerosol and 'ae' for water
         i = i+1

         ! Removal by sedimentation of cloud droplets
         CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'cl')
         i = i+1

         ! Removal by precipitation
         CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'pr')
         i = i+1

         IF (level > 4) THEN
            ! Removal by sedimentation of ice particles
            CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'ic')
            i = i+1

            ! Removal by snow
            CALL set_cs_any(n2,n3,scs_rm(i,:,:),'rm'//trim(nam)//'sn')
            i = i+1
         END IF
      END DO

   END SUBROUTINE cs_rem_save
   !
   ! Calculate warm cloud statistics
   SUBROUTINE set_cs_warm(n1,n2,n3,rc,nc,rp,np,th,dn0,zm,zt,dzm,xt,yt,time)

      USE netcdf

      INTEGER, INTENT(in) :: n1, n2, n3
      REAL, INTENT(in)    :: rc(n1,n2,n3), nc(n1,n2,n3), rp(n1,n2,n3), np(n1,n2,n3), th(n1,n2,n3)
      REAL, INTENT(in)    :: dn0(n1), zm(n1), zt(n1), dzm(n1), xt(n2), yt(n3), time
      REAL :: lwp(n2,n3), ncld(n2,n3), rwp(n2,n3), nrain(n2,n3), zb(n2,n3), zc(n2,n3), &
              th1(n2,n3), lmax(n2,n3)
      INTEGER :: ncloudy(n2,n3), nrainy(n2,n3)
      INTEGER :: i, j, k, iret, VarID
      REAL    :: cld, rn, sval, dmy

      ! No outputs for level 1
      IF (level < 2) RETURN

      ! Calculate stats
      lwp = 0.      ! LWP (kg/m^2)
      ncld = 0.     ! Average CDNC (#/kg)
      ncloudy = 0 ! Number of cloudy grid cells
      rwp = 0.      ! RWP (kg/m^2)
      nrain = 0.    ! Average RDNC (#/kg)
      nrainy = 0    ! Number of cloudy grid cells
      zb = zm(n1)+100.  ! Cloud base (m)
      zc = 0.           ! Cloud top (m)
      lmax = 0.     ! Liquid water mixing ratio (kg/kg)
      th1 = 0.      ! Height of the maximum theta gradient

      DO j = 3, n3-2
         DO i = 3, n2-2
            cld = 0.
            rn = 0.
            sval = 0.
            DO k = 2, n1

               IF (rc(k,i,j) > 0.01e-3) THEN
                  ! Cloudy grid
                  lwp(i,j) = lwp(i,j)+rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  ! Volume weighted average of the CDNC
                  ncld(i,j) = ncld(i,j)+nc(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  cld = cld+dn0(k)*(zm(k)-zm(k-1))
                  ! Number of cloudy pixels
                  ncloudy(i,j) = ncloudy(i,j)+1
                  ! Cloud base and top
                  zb(i,j) = min(zt(k),zb(i,j))
                  zc(i,j) = max(zt(k),zc(i,j))
               END IF

               IF (rp(k,i,j) > 0.001e-3) THEN
                  ! Rainy grid cell
                  rwp(i,j) = rwp(i,j)+rp(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  ! Volume weighted average of the RDNC
                  nrain(i,j) = nrain(i,j)+np(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  rn = rn+dn0(k)*(zm(k)-zm(k-1))
                  ! Number of rainy pixels
                  nrainy(i,j) = nrainy(i,j)+1
               END IF

               ! Maximum liquid water mixing ratio
               lmax(i,j) = max(lmax(i,j),rc(k,i,j))

               ! Height of the maximum theta gradient
               IF (k <= n1-5) THEN
                  dmy = (th(k+1,i,j)-th(k,i,j))*dzm(k)
                  IF (dmy > sval ) THEN
                     sval = dmy
                     th1(i,j) = zt(k)
                  END IF
               END IF

            END DO
            IF (cld > 0.) THEN
               ncld(i,j) = ncld(i,j)/cld
            ELSE
               zb(i,j) = -999.
               zc(i,j) = -999.
            END IF
            IF (rn > 0.) THEN
               nrain(i,j) = nrain(i,j)/rn
            END IF
         END DO
      END DO

      ! Save the data
      iret = nf90_inq_varid(ncid3,'time',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, time, start=(/nrec3/))

      IF (nrec3 == 1) THEN
         iret = nf90_inq_varid(ncid3, 'xt', VarID)
         iret = nf90_put_var(ncid3, VarID, xt(3:n2-2), start = (/nrec3/))
         iret = nf90_inq_varid(ncid3, 'yt', VarID)
         iret = nf90_put_var(ncid3, VarID, yt(3:n3-2), start = (/nrec3/))
      END IF

      iret = nf90_inq_varid(ncid3,'lwp',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, lwp(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'rwp',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, rwp(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'Nc',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ncld(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'Nr',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nrain(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'nccnt',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, ncloudy(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'nrcnt',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nrainy(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'zb',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zb(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'zc',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, zc(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'zi1',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, th1(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_inq_varid(ncid3,'lmax',VarID)
      IF (iret == NF90_NOERR) iret = nf90_put_var(ncid3, VarID, lmax(3:n2-2,3:n3-2), start = (/1,1,nrec3/))

      iret = nf90_sync(ncid3)
      nrec3 = nrec3 + 1

   END SUBROUTINE set_cs_warm


  ! Calculate cold cloud statistics
   SUBROUTINE set_cs_cold(n1,n2,n3,ri,ni,rs,ns,dn0,zm,xt,yt,time)

      USE netcdf
      USE mo_submctl, ONLY : prlim

      INTEGER, INTENT(in) :: n1, n2, n3
      REAL, INTENT(in)    :: ri(n1,n2,n3), ni(n1,n2,n3), rs(n1,n2,n3), ns(n1,n2,n3)
      REAL, INTENT(in)    :: dn0(n1), zm(n1), xt(n2), yt(n3), time
      REAL :: iwp(n2,n3), nice(n2,n3), swp(n2,n3), nsnow(n2,n3), imax(n2,n3)
      INTEGER :: nicy(n2,n3), nsnowy(n2,n3)
      INTEGER :: i, j, k, iret, VarID
      REAL    :: ice, sn

    ! No outputs for levels less than 5
      IF (level < 5) RETURN

    ! Calculate stats
      iwp = 0.  ! IWP (kg/m^2)
      nice = 0. ! Average ice number concentration (#/kg)
      nicy = 0  ! Number of icy grid cells
      swp = 0.  ! SWP (kg/m^2)
      nsnow = 0.    ! Average snow number concentration (#/kg)
      nsnowy = 0    ! Number of snowy grid cells
      imax = 0. ! Maximum ice water mixing ratio (kg/kg)

      DO j = 3, n3-2
         DO i = 3, n2-2
            ice = 0.
            sn = 0.
            DO k = 2, n1
               IF (ni(k,i,j) > prlim .AND. ri(k,i,j) > 1.e-15) THEN
                  ! Icy grid cell
                  iwp(i,j) = iwp(i,j)+ri(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  ! Volume weighted average of the ice number concentration
                  nice(i,j) = nice(i,j)+ni(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  ice = ice+dn0(k)*(zm(k)-zm(k-1))
                  ! Number of icy pixels
                  nicy(i,j) = nicy(i,j)+1
               END IF
               IF (ns(k,i,j) > prlim .AND. rs(k,i,j) > 1.e-20) THEN
                  ! Snowy grid cell
                  swp(i,j) = swp(i,j)+rs(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  ! Volume weighted average of the snow number concentration
                  nsnow(i,j) = nsnow(i,j)+ns(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
                  sn = sn+dn0(k)*(zm(k)-zm(k-1))
                  ! Number of snowy pixels
                 nsnowy(i,j) = nsnowy(i,j)+1
               END IF
               ! Maximum ice water mixing ratio
               imax(i,j) = max(imax(i,j),ri(k,i,j))
            END DO

            IF (ice > 0.) THEN
               nice(i,j) = nice(i,j)/ice
            END IF

            IF (sn > 0.) THEN
               nsnow(i,j) = nsnow(i,j)/sn
            END IF
       END DO
    END DO

    ! Save the data
    iret = nf90_inq_varid(ncid3,'time',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, time, start=(/nrec3/))

    IF (nrec3 == 1) THEN
       iret = nf90_inq_varid(ncid3, 'xt', VarID)
       iret = nf90_put_var(ncid3, VarID, xt(3:n2-2), start = (/nrec3/))
       iret = nf90_inq_varid(ncid3, 'yt', VarID)
       iret = nf90_put_var(ncid3, VarID, yt(3:n3-2), start = (/nrec3/))
    END IF

    iret = nf90_inq_varid(ncid3,'iwp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, iwp(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'swp',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, swp(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'Ni',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nice(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'Ns',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nsnow(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'nicnt',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nicy(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'nscnt',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, nsnowy(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_inq_varid(ncid3,'imax',VarID)
    IF (iret==NF90_NOERR) iret = nf90_put_var(ncid3, VarID, imax(3:n2-2,3:n3-2), start=(/1,1,nrec3/))

    iret = nf90_sync(ncid3)

   END SUBROUTINE set_cs_cold

   !
   ! -----------------------------------------------------------------------
   ! Subroutine set_ts: computes and writes time sequence stats
   !
   SUBROUTINE set_ts(n1,n2,n3,w,th,dn0,zt,zm,dzt,dzm,th00,time)

      USE defs, ONLY : cp

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: w(n1,n2,n3),th(n1,n2,n3)
      REAL, INTENT(in)    :: dn0(n1),zt(n1),zm(n1),dzt(n1),dzm(n1),th00,time

      INTEGER :: k
      REAL    :: bf(n1)

      a_time%data = time
      a_zi1_bar%data = get_zi(n1, n2, n3, 2, th, dzm, zt, 1.)   ! maximum gradient
      a_zi2_bar%data = get_zi(n1, n2, n3, 3, th, thvar, zt, 1.) ! maximum variance

      !
      ! buoyancy flux statistics
      !

      a_vtke%data = 0
      DO k = 2, n1-2

         bf(k) = wtv_res(k) + wtv_sgs(k)

         a_vtke%data = a_vtke%data + (tke_res(k)+tke_sgs(k))*dn0(k)/dzt(k)
         a_sfs_boy%data(k) = a_sfs_boy%data(k) + wtv_sgs(k)*9.8/th00

      END DO

      a_zi3_bar%data = get_zi(n1, n2, n3, 4, th, bf, zm, 1.) ! minimum buoyancy flux
      a_sfcbflx%data = bf(2)
      a_wmax%data = maxval(w)
      a_shf_bar%data = a_shf_bar%data*cp*(dn0(1)+dn0(2))*0.5

   END SUBROUTINE set_ts
   !
   ! -----------------------------------------------------------------------
   ! Subroutine ts_lvl1: computes and writes time sequence stats; for the
   ! zi calculation setting itype=1 selects a concentration threshold
   !
   SUBROUTINE ts_lvl1(n1,n2,n3,dn0,zt,dzm,q)

      USE defs, ONLY : alvl

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: q(n1,n2,n3)
      REAL, INTENT(in)    :: dn0(n1),zt(n1),dzm(n1)

      a_lhf_bar%data = a_lhf_bar%data*alvl*(dn0(1)+dn0(2))*0.5
      a_zi_bar%data = get_zi(n1, n2, n3, 2, q, dzm, zt, 0.5e-3)

   END SUBROUTINE ts_lvl1
   !
   ! -----------------------------------------------------------------------
   ! Subroutine ts_lvl2: computes and writes time sequence stats
   !
   SUBROUTINE ts_lvl2(n1,n2,n3,rc,zt)

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: rc(n1,n2,n3), zt(n1)

      INTEGER :: k,i,j
      REAL    :: cpnt, unit

      a_zb%data = zt(n1)
      a_cfrac%data = 0.
      a_lmax%data = 0.
      a_nccnt%data = 0.

      unit = 1./REAL((n2-4)*(n3-4))
      DO j = 3, n3-2
         DO i = 3, n2-2
            cpnt = 0.
            DO k = 2, n1-2
               IF (rc(k,i,j) > 1.e-5) THEN

                  a_zc%data = max(a_zc%data,zt(k))
                  a_zb%data = min(a_zb%data,zt(k))

                  cpnt = unit

                  a_lmax%data = max(a_lmax%data, rc(k,i,j))
                  a_nccnt%data = a_nccnt%data + 1.

               END IF
            END DO
            a_cfrac%data = a_cfrac%data + cpnt
         END DO
      END DO

      IF (a_zb%data(1) == zt(n1)) a_zb%data = -999.

   END SUBROUTINE ts_lvl2
   !
   ! -----------------------------------------------------------------------
   ! Subroutine ts_lvl4: computes and writes time sequence stats of Salsa variables --
   !  Implemented by Zubair Maalick 20/07/2015
   !  Some rewriting and adjusting by Juha Tonttila
   !
   SUBROUTINE ts_lvl4(n1,n2,n3,rc)
      USE mo_submctl, ONLY : nlim
      USE grid,       ONLY : prtcl, bulkNumc, bulkMixrat,dzt
      USE class_componentIndex, ONLY : IsUsed

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: rc(n1,n2,n3)

      REAL    :: a0(n1,n2,n3), a1(n1,n2,n3)
      INTEGER :: ii,ss
      LOGICAL :: cond_ic(n1,n2,n3), cond_oc(n1,n2,n3)

      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt

      CALL bulkNumc('cloud','a',a0)
      CALL bulkNumc('cloud','b',a1)
      cond_ic(:,:,:) = ( a0(:,:,:) + a1(:,:,:) > nlim .AND. rc(:,:,:) > 1.e-5 )
      cond_oc = .NOT. cond_ic

      a_Nc_ic%data = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
      CALL bulkNumc('aerosol','a',a0)
      CALL bulkNumc('aerosol','b',a1)

      a_Na_int%data = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
      a_Na_oc%data = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_oc)

      ii = 4
      DO ss = 1, 7 !Not including water

         IF (IsUsed(prtcl,zspec(ss))) THEN

            CALL bulkMixrat(zspec(ss),'cloud','a',a0)
            CALL bulkMixrat(zspec(ss),'cloud','b',a1)

            CALL f1sbulk%getField(ii,ArrEl)
            CALL f1sbulk%getData(1,trgt,index=ii)
            trgt%data = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
            ii = ii + 1

            CALL bulkMixrat(zspec(ss),'aerosol','a',a0)
            CALL bulkMixrat(zspec(ss),'aerosol','b',a1)

            CALL f1sbulk%getField(ii,ArrEl)
            CALL f1sbulk%getData(1,trgt,index=ii)
            trgt%data = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_ic)
            ii = ii + 1

            CALL f1sbulk%getField(ii,ArrEl)
            CALL f1sbulk%getData(1,trgt,index=ii)
            trgt%data = get_avg_ts(n1,n2,n3,a0+a1,dzt,cond_oc)
            ii = ii + 1

         ELSE
            ii = ii + 3
         END IF
      END DO

   END SUBROUTINE ts_lvl4

     ! -----------------------------------------------------------------------
  ! subroutine ts_lvl5: computes and writes time sequence stats of Salsa variables --
  !  Implemented by Jaakko Ahola 15/12/2016
 !
 SUBROUTINE ts_lvl5(n1,n2,n3,dn0,zm,rc,ri,rs,srate)
    USE mo_submctl, ONLY : nlim,prlim
    USE grid, ONLY : bulkNumc, bulkMixrat, meanRadius, dzt
    USE class_componentIndex, ONLY : IsUsed

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3
    REAL, INTENT(in) :: rc(n1,n2,n3), ri(n1,n2,n3), rs(n1,n2,n3), srate(n1,n2,n3), &
                        zm(n1) , dn0(n1)

    REAL :: a0(n1,n2,n3), a1(n1,n2,n3), totc(n1,n2,n3), toti(n1,n2,n3), tots(n1,n2,n3)
    REAL :: scr(n2,n3), scr2(n2,n3)
    REAL :: sscnt
    INTEGER :: i, j, k
    LOGICAL :: cond_ic(n1,n2,n3), cond_ii(n1,n2,n3), cond_is(n1,n2,n3)

    CALL bulkNumc('cloud','a',a0)
    CALL bulkNumc('cloud','b',a1)
    totc(:,:,:) = a0(:,:,:)+a1(:,:,:)
    ! In-cloud mask
    cond_ic(:,:,:) = ( totc(:,:,:) > nlim .AND. rc(:,:,:) > 1.e-5 )

    CALL bulkNumc('ice','a',a0)
    CALL bulkNumc('ice','b',a1)
    toti(:,:,:) = a0(:,:,:)+a1(:,:,:)
    ! In-ice mask (grid cells with ice)
    cond_ii(:,:,:) = ( toti(:,:,:) > prlim .AND. ri(:,:,:) > 1.e-15) ! Loose limits for ri

    CALL bulkNumc('snow','a',tots)
    ! In-snow mask (grid cells with snow)
    cond_is(:,:,:) = ( tots(:,:,:) > prlim .AND. rs(:,:,:) > 1.e-20 ) ! Loose limits for rs

    ! Outputs
    a_Ni_ic%data = get_avg_ts(n1,n2,n3,toti,dzt,cond_ic) ! Ice particles in liquid clouds
    a_Ni_ii%data = get_avg_ts(n1,n2,n3,toti,dzt,cond_ii) ! Ice particles in icy clouds
    a_Ni_is%data = get_avg_ts(n1,n2,n3,toti,dzt,cond_is) ! Ice particles in snowy clouds

    a_Ns_ic%data = get_avg_ts(n1,n2,n3,tots,dzt,cond_ic) ! The same for snow ...
    a_Ns_ii%data = get_avg_ts(n1,n2,n3,tots,dzt,cond_ii)
    a_Ns_is%data = get_avg_ts(n1,n2,n3,tots,dzt,cond_is)

    CALL meanRadius('ice','ab',a0)
    a_Ri_ii%data = get_avg_ts(n1,n2,n3,a0,dzt,cond_ii)
    CALL meanRadius('snow','ab',a0)
    a_Rs_is%data =  get_avg_ts(n1,n2,n3,a0,dzt,cond_is)

    ! Could include aerosol and cloud droplets in these regions, but maybe too much data?

    ! IWP, SWP, max(ice) and max(snow), and the number of ice/snow cells
    a_imax%data = 0.
    a_nicnt%data = 0.
    a_smax%data = 0.
    a_nscnt%data = 0.

    scr = 0.   ! IWP
    scr2 = 0.   ! SWP
    sscnt = 0
    DO j = 3, n3-2
       DO i = 3, n2-2
          DO k = 2, n1-2
             IF (cond_ii(k,i,j)) THEN
                a_imax%data = max(a_imax%data, ri(k,i,j))
                a_nicnt%data = a_nicnt%data + 1.
             END IF
             IF (cond_is(k,i,j)) THEN
                a_smax%data = max(a_smax%data, rs(k,i,j))
                a_nscnt%data = a_nscnt%data + 1.
             END IF
             !
             ! Ice-water path (without snow)
             scr(i,j) = scr(i,j)+ri(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
             !
             ! Snow-water path
             scr2(i,j) = scr2(i,j)+rs(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
          END DO
          !
          ! Surface snow rate
          IF (cond_is(2,i,j)) sscnt = sscnt+1
          !
       END DO
    END DO

    a_iwp_bar%data = get_avg2dh(n2,n3,scr(:,:))
    a_swp_bar%data = get_avg2dh(n2,n3,scr2(:,:))
    a_sfrac%data = REAL(sscnt)/REAL( (n3-4)*(n2-4) )

    scr2(:,:) = srate(2,:,:)
    a_sprcp%data = get_avg2dh(n2,n3,scr2)

  END SUBROUTINE ts_lvl5
  !

   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_STAT: Accumulates various statistics over an
   ! averaging period for base (level 0) version of model
   !
   SUBROUTINE accum_stat(n1,n2,n3,u,v,w,t,p,um,vm)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, DIMENSION (n1,n2,n3), INTENT (in) :: u, v, w, t, p
      REAL, INTENT (in) :: um, vm

      INTEGER :: k
      REAL    :: a1(n1), b1(n1), c1(n1), d1(n1), a3(n1), b3(n1), tmp(n1)

      CALL get_avg3(n1,n2,n3, u,a1)
      CALL get_avg3(n1,n2,n3, v,b1)
      CALL get_avg3(n1,n2,n3, t,c1)
      CALL get_avg3(n1,n2,n3, p,d1)
      CALL get_var3(n1,n2,n3, t, c1, thvar)
      CALL get_3rd3(n1,n2,n3, t, c1, b3) ! Used to be (t-a1)**3
      tmp(:) = 0.
      CALL get_3rd3(n1,n2,n3, w, tmp, a3) ! Now just w**3

      DO k = 1, n1
         a_u%data(k) = a_u%data(k) + a1(k) + um
         a_v%data(k) = a_v%data(k) + b1(k) + um
         a_theta_1%data(k) = a_theta_1%data(k) + c1(k)
         a_p%data(k) = a_p%data(k) + d1(k)
         a_theta_2%data(k) = a_theta_2%data(k) + thvar(k)
         a_w_3%data(k) = a_w_3%data(k) + a3(k)
         a_theta_3%data(k) = a_theta_3%data(k) + b3(k)
      END DO

   END SUBROUTINE accum_stat
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_STAT: Accumulates various statistics over an
   ! averaging period for radiation variables
   !
   SUBROUTINE accum_rad(n1,n2,n3,rflx,sflx,sup,sdwn,irup,irdwn,alb)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: rflx(n1,n2,n3)
      REAL, OPTIONAL, INTENT (in) :: sflx(n1,n2,n3), alb(n2,n3)
      REAL, OPTIONAL, INTENT (in) :: sup(n1,n2,n3), sdwn(n1,n2,n3), irup(n1,n2,n3), irdwn(n1,n2,n3)

      INTEGER :: k
      REAL    :: a1(n1),a2(n1)

      CALL get_avg3(n1,n2,n3,rflx,a1)
      CALL get_var3(n1,n2,n3,rflx,a1,a2)

      DO k = 1, n1
         a_rflx1%data(k) = a_rflx1%data(k) + a1(k)
         a_rflx2%data(k) = a_rflx2%data(k) + a2(k)
      END DO

      IF (present(sflx)) THEN
         CALL get_avg3(n1,n2,n3,sflx,a1)
         CALL get_var3(n1,n2,n3,sflx,a1,a2)

         DO k = 1, n1
            a_sflx1%data(k) = a_sflx1%data(k) + a1(k)
            a_sflx2%data(k) = a_sflx2%data(k) + a2(k)
         END DO

         a_albedo%data = get_avg2dh(n2,n3,alb)
      END IF

      IF (present(sup)) THEN
         CALL get_avg3(n1,n2,n3,sup,a1)
         a_sw_up%data(:) = a_sw_up%data(:) + a1(:)

         CALL get_avg3(n1,n2,n3,sdwn,a1)
         a_sw_down%data(:) = a_sw_down%data(:) + a1(:)

         CALL get_avg3(n1,n2,n3,irup,a1)
         a_lw_up%data(:) = a_lw_up%data(:) + a1(:)

         CALL get_avg3(n1,n2,n3,irdwn,a1)
         a_lw_down%data(:) = a_lw_down%data(:) + a1(:)
      END IF

   END SUBROUTINE accum_rad
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL1: Accumulates various statistics over an
   ! averaging period for moisture variable (smoke or total water)
   !
   SUBROUTINE accum_lvl1(n1,n2,n3,q)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: q(n1,n2,n3)

      INTEGER :: k
      REAL    :: a1(n1),a2(n1),a3(n1)

      CALL get_avg3(n1,n2,n3,q,a1)
      CALL get_var3(n1,n2,n3,q,a1,a2)
      CALL get_3rd3(n1,n2,n3,q,a1,a3)

      DO k = 1, n1
         a_q%data(k) = a_q%data(k) + a1(k)
         a_q_2%data(k) = a_q_2%data(k) + a2(k)
         a_q_3%data(k) = a_q_3%data(k) + a3(k)
      END DO

   END SUBROUTINE accum_lvl1
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL2: Accumulates specialized statistics that depend
   ! on level 2 variables.
   !
   SUBROUTINE accum_lvl2(n1, n2, n3, th00, dn0, zm, w, th, t, &
      rv, rl, rt)

      USE defs, ONLY : ep2

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)                      :: th00
      REAL, INTENT (in), DIMENSION(n1)       :: zm, dn0
      REAL, INTENT (in), DIMENSION(n1,n2,n3) :: w, th, t, rv, rl, rt

      REAL, DIMENSION(n1,n2,n3) :: tv    ! Local variable
      INTEGER                   :: k, i, j, kp1
      REAL, DIMENSION(n1)       :: a1, a2, a3, tvbar
      REAL, DIMENSION(n1,n2,n3) :: scr, xy1, xy2, tw, tvw, rtw
      LOGICAL :: cond(n1,n2,n3)

      !
      ! liquid water statistics
      !
      CALL get_avg3(n1,n2,n3,rl,a1)
      CALL get_var3(n1,n2,n3,rl,a1,a2)
      CALL get_3rd3(n1,n2,n3,rl,a1,a3)

      a_l%data(:) = a_l%data(:) + a1(:)
      a_l_2%data(:) = a_l_2%data(:) + a2(:)
      a_l_3%data(:) = a_l_3%data(:) + a3(:)

      !
      ! Do some conditional sampling statistics: cloud, cloud-core
      !
      tv(:,:,:) = th(:,:,:)*(1.+ep2*rv(:,:,:) - (rt(:,:,:)-rv(:,:,:))) ! Virtual potential temperature (K)
      CALL get_avg3(n1,n2,n3,tv,tvbar)
      !
      xy1 = 0.
      xy2 = 0.
      tvw = 0.
      DO k = 1, n1

         kp1 = k+1
         IF (k == n1) kp1 = k

         DO j = 3, n3-2
            DO i = 3, n2-2
               IF (rl(k,i,j) > 1e-5) THEN
                  xy1(k,i,j) = 1.
                  IF (tv(k,i,j) > tvbar(k)) THEN
                     xy2(k,i,j) = 1.
                  END IF
                  !
                  tw(k,i,j) = (.5*(t(k,i,j)+t(kp1,i,j))+th00)*w(k,i,j)
                  tvw(k,i,j) = (.5*(tv(k,i,j)+tv(kp1,i,j)))*w(k,i,j)
                  rtw(k,i,j) = (.5*(rt(k,i,j)+rt(kp1,i,j)))*w(k,i,j)
               END IF
            END DO
         END DO

      END DO

      CALL get_avg3(n1,n2,n3,xy1,a1)
      a_cs1%data(:) = a_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,xy1,a1,normalize=.FALSE.)
      a_cnt_cs1%data(:) = a_cnt_cs1%data(:) + a1(:)

      cond(:,:,:) = xy1(:,:,:)>0.5
      CALL get_avg3(n1,n2,n3,w,a1,cond=cond)
      a_w_cs1%data(:) = a_w_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,t+th00,a1,cond=cond)
      a_t_cs1%data(:) = a_t_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
      a_tv_cs1%data(:) = a_tv_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
      a_rt_cs1%data(:) = a_rt_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
      a_rc_cs1%data(:) = a_rc_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,tw,a1,cond=cond)
      a_wt_cs1%data(:) = a_wt_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
      a_wtv_cs1%data(:) = a_wtv_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
      a_wrt_cs1%data(:) = a_wrt_cs1%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,xy2,a1)
      a_cs2%data(:) = a_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,xy2,a1,normalize=.FALSE.)
      a_cnt_cs2%data(:) = a_cnt_cs2%data(:) + a1(:) ! Counts

      cond(:,:,:) = xy2(:,:,:)>0.5
      CALL get_avg3(n1,n2,n3,w,a1,cond=cond)
      a_w_cs2%data(:) = a_w_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,t+th00,a1,cond=cond)
      a_t_cs2%data(:) = a_t_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,tv,a1,cond=cond)
      a_tv_cs2%data(:) = a_tv_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,rt,a1,cond=cond)
      a_rt_cs2%data(:) = a_rt_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,rl,a1,cond=cond)
      a_rc_cs2%data(:) = a_rc_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,tw,a1,cond=cond)
      a_wt_cs2%data(:) = a_wt_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,tvw,a1,cond=cond)
      a_wtv_cs2%data(:) = a_wtv_cs2%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,rtw,a1,cond=cond)
      a_wrt_cs2%data(:) = a_wrt_cs2%data(:) + a1(:)

      !
      ! liquid water path (without precipitation)
      !
      DO j = 3, n3-2
         DO i = 3, n2-2
            scr(1,i,j) = 0.
            DO k = 2, n1
               scr(1,i,j) = scr(1,i,j)+rl(k,i,j)*dn0(k)*(zm(k)-zm(k-1))
            END DO
         END DO
      END DO

      a_lwp_bar%data = get_avg2dh(n2,n3,scr(1,:,:))

      scr(1,:,:) = (scr(1,:,:)-a_lwp_bar%data(1))**2
      a_lwp_var%data = get_avg2dh(n2,n3,scr(1,:,:))

   END SUBROUTINE accum_lvl2
   !
   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL3: Accumulates specialized statistics that depend
   ! on level 3 variables.
   !
   SUBROUTINE accum_lvl3(n1, n2, n3, dn0, zm, rc, rr, nr, rrate, CCN)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)                      :: CCN
      REAL, INTENT (in), DIMENSION(n1)       :: zm, dn0
      REAL, INTENT (in), DIMENSION(n1,n2,n3) :: rc, rr, nr, rrate

      INTEGER                :: k, i, j
      REAL                   :: nrsum, nrcnt, rrcnt, rrcb
      REAL                   :: rmax, rmin
      REAL, DIMENSION(n1)    :: a1
      REAL, DIMENSION(n2,n3) :: scr2
      REAL                   :: mask(n1,n2,n3), tmp(n1,n2,n3)
      LOGICAL :: below

      !
      ! Average rain water mixing ratio
      !
      CALL get_avg3(n1,n2,n3,rr,a1)
      a_rr%data(:) = a_rr%data(:) + a1(:)
      a_Nc%data(:) = a_Nc%data(:) + CCN ! nc (#/kg)

      !
      ! conditionally average rain droplet concentrations
      !
      WHERE (rr > 0.001e-3)
         mask = 1.
      ELSE WHERE
         mask = 0.
      END WHERE

      CALL get_avg3(n1,n2,n3,nr,a1,cond=(mask>0.5))
      a_Nr%data(:) = a_Nr%data(:)+a1(:)

      CALL get_avg3(n1,n2,n3,mask,a1)
      a_frc_ran%data(:) = a_frc_ran%data(:) +a1(:)

      !
      ! precipitation flux
      !
      CALL get_avg3(n1,n2,n3,rrate,a1)
      a_rrate%data(:) = a_rrate%data(:)+a1(:)

      !
      ! conditionally average precip fluxes
      !
      WHERE (rrate > 3.65e-5)
         mask = 1.
      ELSE WHERE
         mask = 0.
      END WHERE

      tmp(:,:,:) = mask(:,:,:)*rrate(:,:,:)
      CALL get_avg3(n1,n2,n3,tmp,a1)
      a_prc_prc%data(:) = a_prc_prc%data(:) + a1(:)

      CALL get_avg3(n1,n2,n3,mask,a1)
      a_frc_prc%data(:) = a_frc_prc%data(:) + a1(:)

      !
      ! Histogram of surface rain rates
      !
      mask = 0.
      DO k = 1, n1
         rmin = max(6.2e-8,(k-1)*3.421e-5)
         rmax =  k * 3.421e-5
         DO j = 3, n3-2
            DO i = 3, n2-2
               IF (rrate(2,i,j) > rmin .AND. rrate(2,i,j) <= rmax) mask(k,i,j) = 1.
            END DO
         END DO
      END DO

      CALL get_avg3(n1,n2,n3,mask,a1,normalize = .FALSE.)
      a_hst_srf%data(:) = a_hst_srf%data(:) + a1(:)

      !
      ! Temporal statistics
      !
      scr2(:,:) = 0.
      nrsum = 0.
      nrcnt = 0.
      rrcnt = 0.
      rrcb = 0.

      DO j = 3, n3-2
         DO i = 3, n2-2
            below = .TRUE.

            DO k = 2, n1

               ! RWP
               scr2(i,j) = scr2(i,j)+rr(k,i,j)*dn0(k)*(zm(k)-zm(k-1))

               ! Rainy grid cell
               IF (rr(k,i,j) > 0.001e-3) THEN
                  nrsum = nrsum + nr(k,i,j)
                  nrcnt = nrcnt + 1.
               END IF

               ! Surface precipitation for this column
               IF (k == 2 .AND. rrate(k,i,j) > 3.65e-5) rrcnt = rrcnt + 1.

               ! Precpitation at cloud base (no cloud = no precip.)
               IF (rc(k,i,j) > 1.e-5 .AND. below) THEN
                  ! Take precpitation from level k-1 (>=2), which is just below cloud base
                  rrcb = rrcb + rrate(max(2,k-1),i,j)
                  below = .FALSE.
               END IF
            END DO
         END DO
      END DO

      a_pfrac%data(:) = rrcnt/REAL( (n3-4)*(n2-4) )
      a_rwp_bar%data(:) = get_avg2dh(n2,n3,scr2)

      scr2(:,:) = rrate(2,:,:)
      a_prcp%data(:) = get_avg2dh(n2,n3,scr2)
      a_ccn%data(:) = CCN

      IF (nrcnt > 0.) a_nrain%data = nrsum/nrcnt

      a_nrcnt%data(:) = nrcnt
      a_prcp_bc%data(:) = rrcb/REAL( (n3-4)*(n2-4) )

   END SUBROUTINE accum_lvl3

   !---------------------------------------------------------------------
   ! Subroutine ACCUM_LVL4: Accumulates specialized statistics that depend
   ! on level 4 variables.
   !
   SUBROUTINE accum_lvl4(n1,n2,n3)

      USE mo_submctl, ONLY : in1a, in2b, fn2a, fn2b,     &
                             ica,  fca,  icb,  fcb,  ira,  fra, &
                             nprc, nlim, prlim
      USE grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, &
                       a_rc, a_srp, a_rp, a_rh, prtcl,    &
                       a_naerop, a_ncloudp, a_nprecpp, a_tp
      USE class_ComponentIndex, ONLY : IsUsed
      USE classFieldArray

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3
      INTEGER :: ii,ss,bb

      LOGICAL :: cloudmask(n1,n2,n3)
      LOGICAL :: drizzmask(n1,n2,n3)


      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt
      TYPE(FloatArray2d), POINTER :: trgt2

      REAL, DIMENSION(n1,n2,n3)           :: a1,a12
      REAL, DIMENSION(n1,5)               :: a2
      REAL, DIMENSION(n1,fn2a)            :: a3_a
      REAL, DIMENSION(n1,fn2b-fn2a)       :: a3_b
      REAL, DIMENSION(n1,fca%cur)         :: a4_a
      REAL, DIMENSION(n1,fcb%cur-fca%cur) :: a4_b
      REAL, DIMENSION(n1,nprc)            :: a5


      ! *************************
      ! Bulk output for SALSA
      ! *************************
      CALL bulkNumc('aerosol','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,1))

      CALL bulkNumc('aerosol','b',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,2))

      CALL bulkNumc('cloud','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,3))

      CALL bulkNumc('cloud','b',a12)
      CALL get_avg3(n1,n2,n3,a12,a2(:,4))

      ! In cloud mask
      WHERE (a1+a12 > nlim)
         cloudmask = .TRUE.
      ELSE WHERE
         cloudmask = .FALSE.
      END WHERE

      CALL bulkNumc('precp','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,5))

      ! In drizzle mask
      WHERE (a1 > prlim)
         drizzmask = .TRUE.
      ELSE WHERE
         drizzmask = .FALSE.
      END WHERE

      a_Naa%data(:) = a_Naa%data(:) + a2(:,1)
      a_Nab%data(:) = a_Nab%data(:) + a2(:,2)
      a_Nca%data(:) = a_Nca%data(:) + a2(:,3)
      a_Ncb%data(:) = a_Ncb%data(:) + a2(:,4)
      a_Np%data(:)  = a_Np%data(:) + a2(:,5)


      ! Particle radius
      CALL meanRadius('aerosol','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,1))

      CALL meanRadius('aerosol','b',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,2))

      ! In-cloud
      CALL meanRadius('cloud','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=cloudmask)

      CALL meanRadius('cloud','b',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,4),cond=cloudmask)

      ! In-drizzle
      CALL meanRadius('precp','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,5),cond=drizzmask)

      a_Rwaa%data(:) = a_Rwaa%data(:) + a2(:,1)
      a_Rwab%data(:) = a_Rwab%data(:) + a2(:,2)
      a_Rwca%data(:) = a_Rwca%data(:) + a2(:,3)
      a_Rwcb%data(:) = a_Rwcb%data(:) + a2(:,4)
      a_Rwp%data(:)  = a_Rwp%data(:) + a2(:,5)

      ! Bin number concentrations
      ! -------------------------------------------
      IF (lbinprof) THEN

         DO bb = in1a, fn2a
            CALL get_avg3(n1,n2,n3,a_naerop%data(:,:,:,bb),a3_a(:,bb))
         END DO

         DO bb = in2b, fn2b
            CALL get_avg3(n1,n2,n3,a_naerop%data(:,:,:,bb),a3_b(:,bb-fn2a))
         END DO

         DO bb = ica%cur, fca%cur
            CALL get_avg3(n1,n2,n3,a_ncloudp%data(:,:,:,bb),a4_a(:,bb))
         END DO

         DO bb = icb%cur, fcb%cur
            CALL get_avg3(n1,n2,n3,a_ncloudp%data(:,:,:,bb),a4_b(:,bb-fca%cur))
         END DO

         DO bb = ira, fra
            CALL get_avg3(n1,n2,n3,a_nprecpp%data(:,:,:,bb),a5(:,bb))
         END DO

         a_Naba%data(:,:) = a_Naba%data(:,:) + a3_a(:,:)
         a_Nabb%data(:,:) = a_Nabb%data(:,:) + a3_b(:,:)
         a_Ncba%data(:,:) = a_Ncba%data(:,:) + a4_a(:,:)
         a_Ncbb%data(:,:) = a_Ncbb%data(:,:) + a4_b(:,:)
         a_Npb%data(:,:)  = a_Npb%data(:,:) + a5(:,:)

      END IF


      ! Species mixing ratios
      ! -------------------------------------------
      ii = 16
      DO ss = 1, 8 !Include water ss = 8

         IF (ss == 8 .OR. IsUsed(prtcl,zspec(ss))) THEN
            ! Total mass mixing ratios
            CALL bulkMixrat(zspec(ss),'aerosol','a',a1)
            CALL bulkMixrat(zspec(ss),'aerosol','b',a12)
            CALL get_avg3(n1,n2,n3,a1+a12,a2(:,1))

            ! In-cloud
            CALL bulkMixrat(zspec(ss),'cloud','a',a1)
            CALL bulkMixrat(zspec(ss),'cloud','b',a12)
            CALL get_avg3(n1,n2,n3,a1+a12,a2(:,2),cond=cloudmask)

            ! In-drizzle
            CALL bulkMixrat(zspec(ss),'precp','a',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=drizzmask)

            CALL f2sbulk%getField(ii,Arrel)
            CALL f2sbulk%getData(1,trgt,name=Arrel%name)
            trgt%data(:) = trgt%data(:) + a2(:,1)

            CALL f2sbulk%getField(ii+1,Arrel)
            CALL f2sbulk%getData(1,trgt,name=Arrel%name)
            trgt%data(:) = trgt%data(:) + a2(:,2)

            CALL f2sbulk%getField(ii+2,Arrel)
            CALL f2sbulk%getData(1,trgt,name=Arrel%name)
            trgt%data(:) = trgt%data(:) + a2(:,3)

            ! Binned mixing ratios (not for water)
            IF (lbinprof .AND. ss<8) THEN
               DO bb = in1a, fn2a
                  CALL binSpecMixrat('aerosol',zspec(ss),bb,a1) ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a3_a(:,bb))         ! average profile for bin bb for species ss
               END DO

               DO bb = in2b, fn2b
                  CALL binSpecMixrat('aerosol',zspec(ss),bb,a1) ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a3_b(:,bb-fn2a))    ! average profile for bin bb for species ss
               END DO

               DO bb = ica%cur,fca%cur
                  CALL binSpecMixrat('cloud',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=cloudmask)        ! average profile for bin bb for species ss
               END DO

               DO bb = icb%cur,fcb%cur
                  CALL binSpecMixrat('cloud',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fca%cur),cond=cloudmask)! average profile for bin bb for species ss
               END DO

               ! Binned mixing ratios
               DO bb = 1,nprc
                  CALL binSpecMixrat('precp',zspec(ss),bb,a1)  ! z,x,y-field for bin bb
                  CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=drizzmask)          ! average profile for bin bb for species ss
               END DO

               CALL f2aea%getField(ss+1,Arrel)
               CALL f2aea%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a3_a(:,:)

               CALL f2aeb%getField(ss+1,Arrel)
               CALL f2aeb%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a3_b(:,:)

               CALL f2cla%getField(ss+1,Arrel)
               CALL f2cla%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a4_a(:,:)

               CALL f2clb%getField(ss+1,Arrel)
               CALL f2clb%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a4_b(:,:)

               CALL f2prc%getField(ss+1,Arrel)
               CALL f2prc%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a5(:,:)

            END IF
         END IF ! IsUsed

         ii = ii + 3

      END DO ! ss


      ! Liquid water mixing ratio
      CALL get_avg3(n1,n2,n3,a_rc%data,a2(:,1))

      ! Precipitation mixing ratio
      CALL get_avg3(n1,n2,n3,a_srp%data,a2(:,2))

      ! Water vapor mixing ratio
      CALL get_avg3(n1,n2,n3,a_rp%data,a2(:,3))

      ! Relative humidity
      CALL get_avg3(n1,n2,n3,a_rh%data,a2(:,4))
      a2(:,4) = a2(:,4)*100.0 ! RH in %

      a_prl%data(:) = a_prl%data + a2(:,1)
      a_prr%data(:) = a_prr%data(:) + a2(:,2)
      a_prv%data(:) = a_prv%data(:) + a2(:,3)
      a_prh%data(:) = a_prh%data(:) + a2(:,4)


      ! Stats for cloudy columns
      !   Cloudy column: LWC > 1e-5 kg/kg and CDNC>nlim anywhere in a column
      IF (cloudy_col_stats) THEN
         ! Total cloud droplets
         CALL bulkNumc('cloud','ab',a1)
         ! Which columns should be included
         cloudmask(1,:,:)=ANY( (a1>nlim .AND. a_rc%data>1.e-5), DIM=1)
         ! Fill array
         DO ii = 2, n1
            cloudmask(ii,:,:)=cloudmask(1,:,:)
         END DO

         ! Save the fraction of cloudy columns
         WHERE (cloudmask)
            a1 = 1.
         ELSE WHERE
            a1 = 0.
         END WHERE
         CALL get_avg3(n1,n2,n3,a1,a2(:,4),cond=cloudmask)

         ! Aerosol number concentration (a+b)
         CALL bulkNumc('aerosol','ab',a1)
         CALL get_avg3(n1,n2,n3,a1,a2(:,1),cond=cloudmask)

         ! Cloud droplet number concentration (a+b)
         CALL bulkNumc('cloud','ab',a1)
         CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=cloudmask)

         ! Rain drop number concentration
         CALL bulkNumc('precp','a',a1)
         CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=cloudmask)

         ! Save
         a_Na_c%data(:) = a_Na_c%data(:) + a2(:,1)
         a_Nc_c%data(:) = a_Nc_c%data(:) + a2(:,2)
         a_Np_c%data(:) = a_Np_c%data(:) + a2(:,3)
         a_cfrac%data(:) = a_cfrac%data(:) + a2(:,4)

         ! Cloud liquid water mixing ratio
         CALL get_avg3(n1,n2,n3,a_rc%data,a2(:,1),cond=cloudmask)

         ! Liquid water potential temperature
         CALL get_avg3(n1,n2,n3,a_tp%data,a2(:,2),cond=cloudmask)

         ! Save
         a_clw_c%data(:) = a_clw_c%data(:) + a2(:,1)
         a_thl_c%data(:) = a_thl_c%data(:) + a2(:,2)

      END IF

   END SUBROUTINE accum_lvl4

   !---------------------------------------------------------------------
   ! SUBROUTINE ACCUM_LVL5: Accumulates specialized statistics that depend
   ! on level 5 variables.
   !
   SUBROUTINE accum_lvl5(n1,n2,n3,srate)

      USE mo_submctl, ONLY : iia, fia, iib, fib, isa, fsa, nsnw, prlim
      USE grid, ONLY : bulkNumc, bulkMixrat, meanRadius, binSpecMixrat, &
                       a_ri, a_srs, a_rhi, prtcl, a_nicep, a_nsnowp
      USE class_ComponentIndex, ONLY : IsUsed

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT (in), DIMENSION(n1,n2,n3) :: srate
      INTEGER :: ii, ss, bb
      LOGICAL :: icemask(n1,n2,n3)
      LOGICAL :: snowmask(n1,n2,n3)

      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt
      TYPE(FloatArray2d), POINTER :: trgt2

      REAL, DIMENSION(n1,n2,n3) :: a1, a12
      REAL, DIMENSION(n1,3) :: a2
      REAL, DIMENSION(n1,fia%cur) :: a4_a
      REAL, DIMENSION(n1,fib%cur-fia%cur) :: a4_b
      REAL, DIMENSION(n1,nsnw) :: a5

      ! *************************
      ! Bulk output for SALSA
      ! *************************
      CALL bulkNumc('ice','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,1))

      CALL bulkNumc('ice','b',a12)
      CALL get_avg3(n1,n2,n3,a12,a2(:,2))

      ! In ice mask (grid cells with ice)
      icemask(:,:,:) = ( a1(:,:,:)+a12(:,:,:) > prlim .AND. a_ri%data(:,:,:) > 1.e-15) ! Loose limits for ri

      CALL bulkNumc('snow','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,3))

      ! In-snow mask (grid cells with snow)
      snowmask(:,:,:) = ( a1(:,:,:) > prlim .AND. a_srs%data(:,:,:) > 1.e-20 ) ! Loose limits for rs

      a_Nia%data(:) = a_Nia%data(:) + a2(:,1)
      a_Nib%data(:) = a_Nib%data(:) + a2(:,2)
      a_Ns%data(:)  = a_Ns%data(:) + a2(:,3)

      ! Particle radius

      ! Ice
      CALL meanRadius('ice','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,1),cond=icemask)

      CALL meanRadius('ice','b',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=icemask)

      ! Snow
      CALL meanRadius('snow','a',a1)
      CALL get_avg3(n1,n2,n3,a1,a2(:,3),cond=snowmask)

      a_Rwia%data(:) = a_Rwia%data(:) + a2(:,1)
      a_Rwib%data(:) = a_Rwib%data(:) + a2(:,2)
      a_Rws%data(:) = a_Rws%data(:) + a2(:,3)

      ! Bin number concentrations
      ! -------------------------------------------
      IF (lbinprof) THEN
         DO bb = iia%cur, fia%cur
            CALL get_avg3(n1,n2,n3,a_nicep%data(:,:,:,bb),a4_a(:,bb))
         END DO

         DO bb = iib%cur, fib%cur
            CALL get_avg3(n1,n2,n3,a_nicep%data(:,:,:,bb),a4_b(:,bb-fia%cur))
         END DO

         DO bb = isa, fsa
            CALL get_avg3(n1,n2,n3,a_nsnowp%data(:,:,:,bb),a5(:,bb))
         END DO

         a_Niba%data(:,:) = a_Niba%data(:,:) + a4_a(:,:)
         a_Nibb%data(:,:) = a_Nibb%data(:,:) + a4_b(:,:)
         a_Nsb%data(:,:) = a_Nsb%data(:,:) + a5(:,:)

      END IF

      ! Species mixing ratios
      ! -------------------------------------------
      ii = 10 ! 'P_cSO4i'
      DO ss = 1, 8  ! Including water
         IF (IsUsed(prtcl,zspec(ss)) .OR. (ss==8)) THEN
            ! Total mass mixing ratios

            ! In-ice
            CALL bulkMixrat(zspec(ss),'ice','a',a1)
            CALL bulkMixrat(zspec(ss),'ice','b',a12)
            CALL get_avg3(n1,n2,n3,a1+a12,a2(:,1),cond=icemask)

            ! In-snow
            CALL bulkMixrat(zspec(ss),'snow','a',a1)
            CALL get_avg3(n1,n2,n3,a1,a2(:,2),cond=snowmask)

            CALL f2lvl5%getField(ii,Arrel)
            CALL f2lvl5%getData(1,trgt,name=Arrel%name)
            trgt%data(:) = trgt%data(:) + a2(:,1)

            CALL f2lvl5%getField(ii+1,Arrel)
            CALL f2lvl5%getData(1,trgt,name=Arrel%name)
            trgt%data(:) = trgt%data(:) + a2(:,2)

            ! Binned mixing ratios
            IF (lbinprof .AND. ss<8) THEN
               DO bb = iia%cur, fia%cur
                  CALL binSpecMixrat('ice',zspec(ss),bb,a1)
                  CALL get_avg3(n1,n2,n3,a1,a4_a(:,bb),cond=icemask)
               END DO

               DO bb = iib%cur, fib%cur
                  CALL binSpecMixrat('ice',zspec(ss),bb,a1)
                  CALL get_avg3(n1,n2,n3,a1,a4_b(:,bb-fia%cur),cond=icemask)
               END DO

               ! Binned mixing ratios
               DO bb = 1, nsnw
                 CALL binSpecMixrat('snow',zspec(ss),bb,a1)
                 CALL get_avg3(n1,n2,n3,a1,a5(:,bb),cond=snowmask)
               END DO

               CALL f2ica%getField(ss+1,Arrel)
               CALL f2ica%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a4_a(:,:)

               CALL f2icb%getField(ss+1,Arrel)
               CALL f2icb%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a4_a(:,:)

               CALL f2snw%getField(ss+1,Arrel)
               CALL f2snw%getData(1,trgt2,name=Arrel%name)
               trgt2%data(:,:) = trgt2%data(:,:) + a5(:,:)
            END IF
         END IF ! IsUsed

         ii = ii + 2

      END DO ! ss

      ! Ice water mixing ratio
      CALL get_avg3(n1,n2,n3,a_ri%data,a2(:,1))

      ! Snow water mixing ratio
      CALL get_avg3(n1,n2,n3,a_srs%data,a2(:,2))

      ! Relative humidity ove ice
      CALL get_avg3(n1,n2,n3,a_rhi%data,a2(:,3))
      a2(:,3) = a2(:,3)*100.0 ! RH in %

      a_pri%data(:)  = a_pri%data(:) + a2(:,1)
      a_prs%data(:)  = a_prs%data(:) + a2(:,2)
      a_prhi%data(:) = a_prhi%data(:) + a2(:,3)

      ! Snow deposition flux
      CALL get_avg3(n1,n2,n3,srate,a2(:,1))
      a_srate%data(:) = a_srate%data(:) + a2(:,1)

      ! Could add statistics for ice-containing columns

   END SUBROUTINE accum_lvl5

   !
   !
   !
   SUBROUTINE comp_tke(n1,n2,n3,dzm,th00,u,v,w,s)

      INTEGER, INTENT (in) :: n1,n2,n3
      REAL, INTENT (in)    :: dzm(n1),th00,u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3)
      REAL, INTENT (inout) :: s(n1,n2,n3)

      INTEGER :: k,kp1
      REAL    :: x1(n1), x2(n1)

      !
      ! ------
      ! Calculates buoyancy forcing
      !
      CALL get_buoyancy(n1,n2,n3,s,w,th00)
      !
      ! ------
      ! Estimates shear component of TKE budget
      !
      CALL get_shear(n1,n2,n3,u,v,w,dzm)
      !
      ! ------
      ! Calculates horizontal variances and resolved TKE
      !
      CALL get_avg3(n1,n2,n3,u**2,x1)
      CALL get_avg3(n1,n2,n3,u,x2)
      a_u_2%data(:) = a_u_2%data(:) + (x1(:)-x2(:)**2)
      tke_res(:) = (x1(:)-x2(:)**2)

      CALL get_avg3(n1,n2,n3,v**2,x1)
      CALL get_avg3(n1,n2,n3,v,x2)
      a_v_2%data(:) = a_v_2%data(:) + (x1(:)-x2(:)**2)
      tke_res(:)  = tke_res(:) + (x1(:)-x2(:)**2)

      CALL get_avg3(n1,n2,n3,w**2,x1)
      a_w_2%data(:) = a_w_2%data(:) + x1(:)

      DO k = 1, n1
         kp1 = min(k+1,n1)
         tke_res(k)  = 0.5*(0.5*(tke_res(k)+tke_res(kp1)) + x1(k))
      END DO

      IF (nsmp == 0) tke0(:) = tke_res(:)

   END SUBROUTINE comp_tke
   !
   ! ---------------------------------------------------------------------
   ! get_buoyancy:  estimates buoyancy production term in tke budget
   !
   SUBROUTINE get_buoyancy(n1,n2,n3,b,w,th00)

      USE defs, ONLY : g

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: w(n1,n2,n3),th00
      REAL, INTENT(inout) :: b(n1,n2,n3)

      INTEGER :: i,j,k,kp1

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               kp1 = min(k+1,n1)
               b(k,i,j) = (b(k,i,j) + b(kp1,i,j))
            END DO
         END DO
      END DO

      CALL get_cor3(n1,n2,n3,b,w,wtv_res)

      DO k = 1, n1
         a_boy_prd%data(k) = a_boy_prd%data(k) + wtv_res(k)
         wtv_res(k) = wtv_res(k) * th00/g
      END DO

   END SUBROUTINE get_buoyancy
   !
   ! ---------------------------------------------------------------------
   ! get_shear:  estimates shear production term in tke budget
   !
   SUBROUTINE get_shear(n1,n2,n3,u,v,w,dzm)

      INTEGER, INTENT(in) :: n3,n2,n1
      REAL, INTENT(in)    :: w(n1,n2,n3),dzm(n1),u(n1,n2,n3),v(n1,n2,n3)

      REAL    :: ub(n1), vb(n1)
      INTEGER :: i,j,k
      REAL    :: fact, uw_shear, vw_shear

      fact = 0.25/float((n2-4)*(n3-4))

      CALL get_avg3(n1,n2,n3,u,ub)
      CALL get_avg3(n1,n2,n3,v,vb)

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1

               uw_shear = -(u(k,i,j)-ub(k))*fact*(                        &
                           (w(k,i,j)  +w(k,i+1,j)  )*(ub(k+1)-ub(k)  )*dzm(k) +    &
                           (w(k-1,i,j)+w(k-1,i+1,j))*(ub(k)  -ub(k-1))*dzm(k-1))

               IF (j > 1) vw_shear = -(v(k,i,j)-vb(k))*fact*(            &
                                      (w(k,i,j)  +w(k,i,j+1)  )*(vb(k+1)-vb(k)  )*dzm(k) +   &
                                      (w(k-1,i,j)+w(k-1,i,j+1))*(vb(k)  -vb(k-1))*dzm(k-1))

               a_prd_uw%data(k) = a_prd_uw%data(k) + uw_shear
               a_shr_prd%data(k) = a_shr_prd%data(k) + uw_shear + vw_shear

            END DO
         END DO
      END DO

   END SUBROUTINE get_shear
   !
   ! ----------------------------------------------------------------------
   ! Subroutine write_ts: writes the statistics file
   !
   SUBROUTINE write_ts()

      USE classFieldArray
      USE netcdf
      USE grid, ONLY : level

      INTEGER :: iret, n, VarID
      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt

      DO n = 1, f1%count
         CALL f1%getField(n,ArrEl)
         iret = nf90_inq_varid(ncid1, ArrEl%name, VarID)
         CALL f1%getData(1,trgt,name=ArrEl%name)
         iret = nf90_put_var(ncid1, VarID, trgt%data, start=(/nrec1/))
         trgt%data = 0.
      END DO

      IF (level >= 4) THEN
         DO n = 1, f1sbulk%count
            CALL f1sbulk%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid1, ArrEl%name, VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            CALL f1sbulk%getData(1,trgt,name=ArrEl%name)
            iret = nf90_put_var(ncid1, VarID, trgt%data, start=(/nrec1/))
            trgt%data = 0.
         END DO
      END IF

      IF (level >= 5) THEN
         DO n = 1, f1lvl5%count
            CALL f1lvl5%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid1, ArrEl%name, VarID)
            CALL f1lvl5%getData(1,trgt,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid1, VarID, trgt%data, start=(/nrec1/))
            trgt%data = 0.
         END DO
      END IF

      iret = nf90_sync(ncid1)
      nrec1 = nrec1 + 1

   END SUBROUTINE write_ts
   !
   ! ----------------------------------------------------------------------
   ! Subroutine write_ps: writes the time averaged elements of the
   ! statistics file
   !
   SUBROUTINE  write_ps(n1,dn0,u0,v0,zm,zt,time)

      USE netcdf
      USE defs, ONLY : alvl, cp
      USE mo_submctl, ONLY : in1a, in2b, fn2a, fn2b, fca, ica, fcb, icb, fra, ira, &
                             iia, fia, iib, fib, isa, fsa, &
                             aerobins, cloudbins, precpbins, icebins, snowbins

      INTEGER, INTENT (in) :: n1
      REAL, INTENT (in)    :: time
      REAL, INTENT (in)    :: dn0(n1), u0(n1), v0(n1), zm(n1), zt(n1)
      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt
      TYPE(FloatArray2d), POINTER :: trgt2

      INTEGER :: iret, VarID, k, n, kp1, i

      lsttm = time

      DO k = 1, n1
         kp1 = min(n1,k+1)
         a_tot_tw%data(k) = (a_tot_tw%data(k) + a_sfs_tw%data(k))*cp

         a_tot_uw%data(k) = a_tot_uw%data(k) + a_sfs_uw%data(k)

         a_tot_vw%data(k) = a_tot_vw%data(k) + a_sfs_vw%data(k)

         a_tot_ww%data(k) = a_tot_ww%data(k) + a_sfs_ww%data(k)

         a_tot_qw%data(k) = (a_tot_qw%data(k)+a_sfs_qw%data(k))*alvl

         a_sfs_tw%data(k) = a_sfs_tw%data(k)*cp

         a_sfs_qw%data(k) = a_sfs_qw%data(k)*alvl

         a_tot_lw%data(k) = a_tot_lw%data(k)*alvl

         a_trans%data(k) = a_adv_w%data(k) + a_prs_w%data(k) + ( &
                           a_prs_u%data(k) + a_prs_u%data(kp1) + a_prs_v%data(k) + a_prs_w%data(kp1) + &
                           a_adv_u%data(k) + a_adv_u%data(kp1) + a_adv_v%data(k) + a_adv_w%data(kp1) - &
                           a_shr_prd%data(k) - a_shr_prd%data(kp1))*0.5


         IF (lsttm > fsttm) THEN
            a_storage%data(k) = (tke_res(k) - tke0(k))/(lsttm-fsttm)
         ELSE
            a_storage%data(k) = 0.
         END IF

         DO i=10,f2%count
            CALL f2%getField(i,Arrel)
            CALL f2%getData(1,trgt,name=Arrel%name)
            trgt%data(k) = trgt%data(k)/nsmp
         END DO

         IF (level >= 4 ) THEN

            DO i = 1, f2sbulk%count
               CALL f2sbulk%getField(i,Arrel)
               CALL f2sbulk%getData(1,trgt,name=Arrel%name)
               trgt%data(k) = trgt%data(k)/nsmp
            END DO

            DO i = 1, f2aea%count
               CALL f2aea%getField(i,Arrel)
               CALL f2aea%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            DO i = 1, f2aeb%count
               CALL f2aeb%getField(i,Arrel)
               CALL f2aeb%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            DO i = 1, f2cla%count
               CALL f2cla%getField(i,Arrel)
               CALL f2cla%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            DO i = 1, f2clb%count
               CALL f2clb%getField(i,Arrel)
               CALL f2clb%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            DO i = 1, f2prc%count
               CALL f2prc%getField(i,Arrel)
               CALL f2prc%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            ! Replace level 3 CDNC=CCN with that from SALSA (a + b bins)
            a_Nc%data(k) = a_Nca%data(k) + a_Ncb%data(k)
         END IF

         IF (level >= 5) THEN
            DO i = 1, f2lvl5%count
               CALL f2lvl5%getField(i,Arrel)
               CALL f2lvl5%getData(1,trgt,name=Arrel%name)
               trgt%data(k) = trgt%data(k)/nsmp
            END DO

            DO i = 1, f2ica%count
               CALL f2ica%getField(i,Arrel)
               CALL f2ica%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            DO i = 1, f2icb%count
               CALL f2icb%getField(i,Arrel)
               CALL f2icb%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

            DO i = 1, f2snw%count
               CALL f2snw%getField(i,Arrel)
               CALL f2snw%getData(1,trgt2,name=Arrel%name)
               trgt2%data(k,:) = trgt2%data(k,:)/nsmp
            END DO

         END IF

      END DO

      iret = nf90_inq_VarID(ncid2, "time", VarID)
      iret = nf90_put_var(ncid2, VarID, time, start=(/nrec2/))

      IF (nrec2 == 1) THEN
         iret = nf90_inq_varid(ncid2, "zt", VarID)
         iret = nf90_put_var(ncid2, VarID, zt, start = (/nrec2/))
         iret = nf90_inq_varid(ncid2, "zm", VarID)
         iret = nf90_put_var(ncid2, VarID, zm, start = (/nrec2/))
         iret = nf90_inq_varid(ncid2, "dn0", VarID)
         iret = nf90_put_var(ncid2, VarID, dn0, start = (/nrec2/))
         ! Juha: For SALSA
         IF (level >= 4) THEN
            iret = nf90_inq_varid(ncid2,"aea",VarID)
            iret = nf90_put_var(ncid2,VarID,aerobins(in1a:fn2a),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,"aeb",VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,aerobins(in2b:fn2b),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,"cla",VarID)
            iret = nf90_put_var(ncid2,VarID,cloudbins(ica%cur:fca%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,"clb",VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,cloudbins(icb%cur:fcb%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,"prc",VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,precpbins(ira:fra),start=(/nrec2/))
         END IF
         IF (level >= 5) THEN
            iret = nf90_inq_varid(ncid2,"ica",VarID)
            iret = nf90_put_var(ncid2,VarID,cloudbins(iia%cur:fia%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,"icb",VarID)
            IF (iret == NF90_NOERR) iret = nf90_put_var(ncid2,VarID,icebins(iib%cur:fib%cur),start=(/nrec2/))
            iret = nf90_inq_varid(ncid2,"snw",VarID)
            iret = nf90_put_var(ncid2,VarID,snowbins(isa:fsa),start=(/nrec2/))
         END IF
         ! \\ SALSA
         iret = nf90_inq_varid(ncid2, "u0", VarID)
         iret = nf90_put_var(ncid2, VarID, u0, start = (/nrec2/))
         iret = nf90_inq_varid(ncid2, "v0", VarID)
         iret = nf90_put_var(ncid2, VarID, v0, start = (/nrec2/))
      END IF

      iret = nf90_inq_VarID(ncid2, "fsttm", VarID)
      iret = nf90_put_var(ncid2, VarID, fsttm, start=(/nrec2/))
      iret = nf90_inq_VarID(ncid2, "lsttm", VarID)
      iret = nf90_put_var(ncid2, VarID, lsttm, start=(/nrec2/))
      iret = nf90_inq_VarID(ncid2, "nsmp", VarID)
      iret = nf90_put_var(ncid2, VarID, nsmp,  start=(/nrec2/))

      DO n = 10, f2%count
         CALL f2%getField(n,ArrEl)
         iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
         IF (iret /= NF90_NOERR) CYCLE
         CALL f2%getData(1,trgt,name=ArrEl%name)
         iret = nf90_put_var(ncid2, VarID, trgt%data(:), start=(/1,nrec2/),count=(/n1,1/))
         trgt%data(:) = 0.
      END DO

      IF (level >= 4) THEN

         DO n = 6, f2sbulk%count
            CALL f2sbulk%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            CALL f2sbulk%getData(1,trgt,name=ArrEl%name)
            iret = nf90_put_var(ncid2, VarID, trgt%data(:), start=(/1,nrec2/),count=(/n1,1/))
            trgt%data(:) = 0.
         END DO

         DO n = 1, f2aea%count
            CALL f2aea%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2aea%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fn2a,1/))
            trgt2%data(:,:) = 0.
         END DO

         DO n = 1, f2aeb%count
            CALL f2aeb%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2aeb%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fn2b-fn2a,1/))
            trgt2%data(:,:) = 0.
         END DO

         DO n = 1, f2cla%count
            CALL f2cla%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2cla%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fca%cur,1/))
            trgt2%data(:,:) = 0.
         END DO

         DO n = 1, f2clb%count
            CALL f2clb%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2clb%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fcb%cur-fca%cur,1/))
            trgt2%data(:,:) = 0.
         END DO

         DO n = 1, f2prc%count
            CALL f2prc%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2prc%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fra,1/))
            trgt2%data(:,:) = 0.
         END DO
      END IF

      IF (level >= 5) THEN

         DO n = 4, f2lvl5%count
            CALL f2lvl5%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2lvl5%getData(1,trgt,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt%data(:), start=(/1,nrec2/),count=(/n1,1/))
            trgt%data(:) = 0.
         END DO

         DO n = 1, f2ica%count
            CALL f2ica%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2ica%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fia%cur,1/))
            trgt2%data(:,:) = 0.
         END DO

         DO n = 1, f2icb%count
            CALL f2icb%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2icb%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fib%cur-fia%cur,1/))
            trgt2%data(:,:) = 0.
         END DO

         DO n = 1, f2snw%count
            CALL f2snw%getField(n,ArrEl)
            iret = nf90_inq_varid(ncid2, ArrEl%name, VarID)
            CALL f2snw%getData(1,trgt2,name=ArrEl%name)
            IF (iret /= NF90_NOERR) CYCLE ! Probably due to aerosol species not being used
            iret = nf90_put_var(ncid2, VarID, trgt2%data(:,:), start=(/1,1,nrec2/), count=(/n1,fsa,1/))
            trgt2%data(:,:) = 0.
         END DO

      END IF

      iret  = nf90_sync(ncid2)
      nrec2 = nrec2+1
      nsmp  = 0.

   END SUBROUTINE write_ps
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: sfc_stat:  Updates statistical arrays with surface flux
   ! variables
   !
   SUBROUTINE sfc_stat(n2,n3,tflx,qflx,ustar,sst)

      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(in), DIMENSION(n2,n3) :: tflx, qflx, ustar
      REAL, INTENT(in)    :: sst

      a_tsrf%data = sst
      a_ustar%data =  get_avg2dh(n2,n3,ustar)
      a_shf_bar%data =  get_avg2dh(n2,n3,tflx)
      IF (level >= 1) a_lhf_bar%data = get_avg2dh(n2,n3,qflx)

   END SUBROUTINE sfc_stat
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: fills scalar array based on index
   ! 1: cfl; 2 max divergence
   !
   SUBROUTINE fill_scalar(index,xval)

      INTEGER, INTENT(in) :: index
      REAL, INTENT (in)   :: xval

      SELECT CASE(index)
         CASE(1)
            a_cfl%data(:) = xval
         CASE(2)
            a_maxdiv%data(:) = xval
      END SELECT

   END SUBROUTINE fill_scalar
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: calculates the dissipation for output diagnostics, if
   ! isgstyp equals 2 then le is passed in via diss
   !
   SUBROUTINE sgs_vel(n1,n2,n3,v1,v2,v3)

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: v1(n1),v2(n1),v3(n1)

      a_sfs_uw%data(:) = a_sfs_uw%data(:)+v1(:)/float((n2-2)*(n3-2))
      a_sfs_vw%data(:) = a_sfs_vw%data(:)+v2(:)/float((n2-2)*(n3-2))
      a_sfs_ww%data(:) = a_sfs_ww%data(:)+v3(:)/float((n2-2)*(n3-2))

   END SUBROUTINE sgs_vel
   !
   ! --------------------------------------------------------------------------
   ! SGSFLXS: estimates the sgs rl and tv flux from the sgs theta_l and sgs r_t
   ! fluxes
   !
   SUBROUTINE sgsflxs(n1,n2,n3,level,rl,rv,th,flx,type)

      USE defs, ONLY : alvl, cp, rm, ep2

      INTEGER, INTENT(in) :: n1,n2,n3,level
      REAL, INTENT(in)    :: rl(n1,n2,n3),rv(n1,n2,n3)
      REAL, INTENT(in)    :: th(n1,n2,n3),flx(n1,n2,n3)
      CHARACTER (len=2)   :: type

      INTEGER :: k,i,j
      REAL    :: rnpts      ! reciprical of number of points and
      REAL    :: fctl, fctt ! factors for liquid (l) and tv (t) fluxes

      IF (type == 'tl') THEN
         wrl_sgs(:) = 0.
         wtv_sgs(:) = 0.
      END IF
      rnpts = 1./REAL((n2-4)*(n3-4))
      !
      ! calculate fluxes assuming the possibility of liquid water.  if liquid
      ! water does not exist sgs_rl = 0.
      !
      IF ( level >= 2 ) THEN
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1-1
                  IF (rl(k+1,i,j) > 0.) THEN
                     fctt = rnpts*(1. + rv(k,i,j)*(1.+ep2 +ep2*rv(k,i,j)*alvl     &
                            /(rm*th(k,i,j))))                              &
                            /(1.+(rv(k,i,j)*(alvl/th(k,i,j))**2)/(rm*cp))
                     SELECT CASE (type)
                        CASE ('tl')
                           fctl = -rnpts/(rm*th(k,i,j)**2/(rv(k,i,j)*alvl)+alvl/cp)
                        CASE ('rt')
                           fctl = rnpts/(1.+(rv(k,i,j)*alvl**2)/(cp*rm*th(k,i,j)**2))
                           fctt = (alvl*fctt/cp - th(k,i,j)*rnpts)
                     END SELECT
                     wrl_sgs(k) = wrl_sgs(k) + fctl*flx(k,i,j)
                     wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
                  ELSE
                     SELECT CASE (type)
                        CASE ('tl')
                           fctt = rnpts*(1. + ep2*rv(k,i,j))
                        CASE ('rt')
                           fctt = rnpts*(ep2*th(k,i,j))
                     END SELECT
                     wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
                  END IF
               END DO
            END DO
         END DO
         !
         ! calculate fluxes for dry thermodynamics, i.e., wrl_sgs is by def
         ! zero
         !
      ELSE
         DO k = 1, n1
            wrl_sgs(k) = 0.
         END DO
         DO j = 3, n3-2
            DO i = 3, n2-2
               DO k = 1, n1-1
                  IF ( level >= 1) THEN
                     SELECT CASE (type)
                        CASE ('tl')
                           fctt = rnpts * (1. + ep2*rv(k,i,j))
                        CASE ('rt')
                           fctt = rnpts * ep2*th(k,i,j)
                     END SELECT
                  ELSE
                     fctt = rnpts
                  END IF
                  wtv_sgs(k) = wtv_sgs(k) + fctt*flx(k,i,j)
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE sgsflxs
   !
   ! ----------------------------------------------------------------------
   ! Subroutine fill_tend: fills arrays with current value of tendencies
   !
   SUBROUTINE acc_tend(n1,n2,n3,f1,f2,f3,t1,t2,t3,v1,v2,v3,ic,routine)

      INTEGER, INTENT(in) :: n1,n2,n3,ic
      REAL, INTENT(in)    :: f1(n1,n2,n3),f2(n1,n2,n3),f3(n1,n2,n3)
      REAL, INTENT(in)    :: t1(n1,n2,n3),t2(n1,n2,n3),t3(n1,n2,n3)
      REAL, INTENT(inout) :: v1(n1),v2(n1),v3(n1)
      CHARACTER (len=3)   :: routine

      INTEGER :: k,ii
      REAL    :: x1(n1),x2(n1),x3(n1)

      CALL get_cor3(n1,n2,n3,f1,t1,x1)
      CALL get_cor3(n1,n2,n3,f2,t2,x2)
      CALL get_cor3(n1,n2,n3,f3,t3,x3)

      SELECT CASE (routine)
         CASE ('sgs')
            ii = 39
         CASE ('adv')
            ii = 42
      END SELECT

      SELECT CASE (ic)
         CASE (1)
            DO k = 1, n1
               v1(k) = x1(k)
               v2(k) = x2(k)
               v3(k) = x3(k)
            END DO
         CASE (2)
            DO k = 1, n1

               IF (ii == 39) THEN
                  a_dff_u%data(k) = a_dff_u%data(k) + (x1(k)-v1(k))
                  a_dff_v%data(k)= a_dff_v%data(k) + (x2(k)-v2(k))
                  a_dff_w%data(k)=a_dff_w%data(k) + (x3(k)-v3(k))
               END IF

               IF (ii == 42) THEN
                  a_adv_u%data(k) = a_adv_u%data(k) + (x1(k)-v1(k))
                  a_adv_v%data(k)= a_adv_v%data(k) + (x2(k)-v2(k))
                  a_adv_w%data(k)=a_adv_w%data(k) + (x3(k)-v3(k))
               END IF
            END DO
      END SELECT

   END SUBROUTINE acc_tend
   !
   !---------------------------------------------------------------------
   ! Subroutine updtst: updates appropriate statistical arrays
   !
   SUBROUTINE updtst(n1,routine,nfld,values,ic)

      INTEGER, INTENT(in)            :: n1,nfld,ic
      REAL, INTENT (in)              :: values(n1)
      CHARACTER (len=3), INTENT (in) :: routine

      INTEGER :: nn,k
      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt
      SELECT CASE (routine)
         CASE("sgs")
            SELECT CASE (nfld)
               CASE (-6)
                  nn = 31 ! dissipation length-scale
               CASE (-5)
                  nn = 30 ! mixing length
               CASE (-4)
                  nn = 29 ! eddy diffusivity
               CASE (-3)
                  nn = 28 ! eddy viscosity
               CASE (-2)
                  nn = 38 ! dissipation
               CASE (-1)
                  nn = 32 ! estimated sgs energy
               CASE (1)
                  nn = 21 ! sgs tl flux
               CASE (2)
                  nn = 54 ! sgs rt flux
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE("adv")
            SELECT CASE (nfld)
               CASE (-3)
                  nn = 26 ! adv w flux
               CASE (-2)
                  nn = 24 ! adv v flux
               CASE (-1)
                  nn = 22 ! adv u flux
               CASE (0)
                  nn = 62 ! adv rl flux
               CASE (1)
                  nn = 20 ! adv tl flux
               CASE (2)
                  nn = 53 ! adv rt flux
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE("prs")
            SELECT CASE (nfld)
               CASE (1)
                  nn = 45 ! dpdx u corr
               CASE (2)
                  nn = 46 ! dpdy v corr
               CASE (3)
                  nn = 47 ! dpdz w corr
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE("prc")
            SELECT CASE (nfld)
               CASE (2)
                  nn = 88
               CASE (3)
                  nn = 63
               CASE DEFAULT
                  nn = 0
            END SELECT
         CASE DEFAULT
            nn = 0
      END SELECT

      IF (nn > 0) THEN
         CALL f2%getField(nn,Arrel)
         CALL f2%getData(1,trgt,name=Arrel%name)
         IF (ic == 0) trgt%data(:) = 0.
         DO k = 1, n1
            trgt%data(k) = trgt%data(k)+values(k)
         END DO
      END IF

   END SUBROUTINE updtst

   !
   ! -------------------------------------------------------------------------
   ! Similar to updtst but intended for making temporal statistics of the
   ! aerosol removal processes.
   ! Juha Tonttila, FMI, 2015

   ! Jaakko Ahola, FMI, 2016
   ! Modified for ice'n'snow
   !
   ! -------------------------------------------------------------------------
   !
   SUBROUTINE acc_removal(n2,n3,n4,raer,rcld,rprc,rice,rsnw)
      USE grid, ONLY : prtcl
      USE mo_submctl, ONLY : nbins, ncld, nprc, nice,  nsnw
      USE class_componentIndex, ONLY : IsUsed, GetIndex
      IMPLICIT NONE

      INTEGER, INTENT(in)        :: n2, n3, n4                     ! Grid dimensions
      REAL, INTENT(in)           :: raer(n2,n3,n4*nbins)        ! Array containing the binned 2d-field
      REAL, OPTIONAL, INTENT(in) :: rcld(n2,n3,n4*ncld), &
                                    rprc(n2,n3,n4*nprc), &     ! 2 optional arrays for calculating total removals
                                    rice(n2,n3,n4*nice), &
                                    rsnw(n2,n3,n4*nsnw)
      REAL :: zavg

      TYPE(ArrayElement), POINTER :: ArrEl
      TYPE(FloatArray1d), POINTER :: trgt

      INTEGER :: ss, si
      INTEGER :: tt
      INTEGER :: end,str

      DO ss = 1, 8
         IF ( .NOT. IsUsed(prtcl,zspec(ss)) .AND. (ss<8) ) CYCLE

         si = GetIndex(prtcl,zspec(ss))

         ! Index to ssclr_b and s1SalsaBulk
         tt = 25 +(ss-1)*3

         CALL f1sbulk%getField(tt,ArrEl)
         CALL f1sbulk%getData(1,trgt,index=tt)

         ! Removal by sedimentation of aerosol
         str = (si-1)*nbins+1
         end = si*nbins
         zavg = get_avg2dh( n2,n3,SUM(raer(:,:,str:end),DIM=3) )
         trgt%data = trgt%data + zavg

         tt = tt+1

         CALL f1sbulk%getField(tt,ArrEl)
         CALL f1sbulk%getData(1,trgt,index=tt)

         ! Removal by sedimentation of cloud droplets
         str = (si-1)*ncld+1
         end = si*ncld
         zavg = get_avg2dh( n2,n3,SUM(rcld(:,:,str:end),DIM=3) )
         trgt%data = trgt%data + zavg
         tt = tt+1

         CALL f1sbulk%getField(tt,ArrEl)
         CALL f1sbulk%getData(1,trgt,index=tt)

         ! Removal by precipitation
         str = (si-1)*nprc+1
         end = si*nprc
         zavg = get_avg2dh( n2,n3,SUM(rprc(:,:,str:end),DIM=3) )
         trgt%data = trgt%data + zavg
         tt = tt+1

         IF (level < 5) CYCLE

         ! Index to ssclr_lvl5 and s1_lvl5
         tt = 15 +(ss-1)*2

        ! Removal by sedimentation of ice particles
        str = (si-1)*nice+1
        end = si*nice
        zavg = get_avg2dh( n2,n3,SUM(rice(:,:,str:end),DIM=3) )
        CALL f1sbulk%getField(tt,ArrEl)
        CALL f1sbulk%getData(1,trgt,index=tt)
        trgt%data = trgt%data + zavg
        tt=tt+1

        ! Removal by snow
        str = (si-1)*nsnw+1
        end = si*nsnw
        zavg = get_avg2dh( n2,n3,SUM(rsnw(:,:,str:end),DIM=3) )
        CALL f1sbulk%getField(tt,ArrEl)
        CALL f1sbulk%getData(1,trgt,index=tt)
        trgt%data = trgt%data + zavg

      END DO

   END SUBROUTINE acc_removal
   !
   !--------------------------------------------------------------------------
   !
   ! Accumulate mass budget terms. Must be done for every timestep. Is called from step,srfc,mcrp
   ! Juha Tonttila, FMI, 2016
   !
   SUBROUTINE acc_massbudged(n1,n2,n3,type,tstep,dz,dn,    &
                             rv,rc,prc,revap,rdep,ApVdom   )

      IMPLICIT NONE

      INTEGER, INTENT(in) :: n1, n2, n3  !
      INTEGER, INTENT(in) :: type      ! 1: Atmospheric water, 2: Water evaporated to the atm,
                                       ! 3: Water removed by precip/deposition
      REAL, INTENT(in) :: tstep        ! Current timestep length

      REAL, INTENT(in) :: dz(n1)
      REAL, INTENT(in) :: dn(n1,n2,n3) ! Air density

      REAL, INTENT(in), OPTIONAL :: rv(n1,n2,n3)      ! Water vapor mixing ratio
      REAL, INTENT(in), OPTIONAL :: rc(n1,n2,n3), &   ! Cloud water mixing ratio
                                    prc(n1,n2,n3)     ! Precipitation mixing ratio

      REAL, INTENT(in), OPTIONAL :: revap(n2,n3), &   ! Gain of water through evaporation at the surface
                                    rdep(n2,n3)       ! Loss of water through deposition

      REAL, INTENT(in), OPTIONAL :: ApVdom            ! Domain surface area / Domain atm volume

      REAL :: a1, a2, a3

      SELECT CASE(type)
         CASE(0)
            IF ( .NOT. present(rv) .OR. .NOT. present(rc) .OR. .NOT. present(prc) ) &
               STOP 'acc_massbudget (stat): ERROR - for atm water q,rc and prc must be present'

            a1 = 0.; a2 = 0.; a3 = 0.
            a1 = get_avg_ts(n1,n2,n3,rv*dn,dz)
            a2 = get_avg_ts(n1,n2,n3,rc*dn,dz)
            a3 = get_avg_ts(n1,n2,n3,prc*dn,dz)
            massbdg(1) = (a1 + a2 + a3)

         CASE(1)
            IF ( .NOT. present(rv) .OR. .NOT. present(rc) .OR. .NOT. present(prc) ) &
               STOP 'acc_massbudget (stat): ERROR - for atm water q,rc and prc must be present'

            a1 = 0.; a2 = 0.; a3 = 0.
            a1 = get_avg_ts(n1,n2,n3,rv*dn,dz)
            a2 = get_avg_ts(n1,n2,n3,rc*dn,dz)
            a3 = get_avg_ts(n1,n2,n3,prc*dn,dz)
            massbdg(2) = (a1 + a2 + a3) ! Not accumulated

         CASE(2)
            IF ( .NOT. present(revap) .OR. .NOT. present(ApVdom) ) &
               STOP 'acc_massbudget (stat): ERROR - for evaporation stats revap must be present'

            a1 = 0.
            a1 = get_avg2dh(n2,n3,revap)
            massbdg(3) = massbdg(3) + a1*tstep*ApVdom

         CASE(3)
            IF ( .NOT. present(rdep) .OR. .NOT. present(ApVdom) ) &
               STOP 'acc_massbudget (stat): ERROR - for deposition stats rdep must be present'

            a1 = 0.
            a1 = get_avg2dh(n2,n3,rdep)
            massbdg(4) = massbdg(4) + a1*tstep*ApVdom

      END SELECT

      avgtime = avgtime + tstep ! Accumulate time for normalization

   END SUBROUTINE acc_massbudged
   !
   ! -------------------------------------------------------------------------
   !
   SUBROUTINE write_massbudged()
      IMPLICIT NONE

      OPEN(88,FILE='MASSBUDGED.TXT')
      WRITE(88,*) 'Initial mass (atm), Final mass (atm), Final-Initial, Total evpaoration, Total deposition, Evap-Dep  '
      WRITE(88,*) massbdg(1), massbdg(2), massbdg(2)-massbdg(1),  &
                  massbdg(3), massbdg(4), massbdg(3)-massbdg(4)
      CLOSE(88)

   END SUBROUTINE write_massbudged

   SUBROUTINE initilize_arrays(nzp)
      USE classFieldArray
      USE mo_structured_datatypes
      USE mo_submctl, ONLY : nprc, fn2a, fn2b, fca, fcb, fra, fia, fib, fsa

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nzp

      REAL :: zeros0d(1)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros1d(nzp)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_aa(nzp,fn2a)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_ab(nzp,fn2b-fn2a)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_ca(nzp,fca%cur)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_cb(nzp,fcb%cur-fca%cur)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_p(nzp,nprc)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_ia(nzp,fia%cur)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_ib(nzp,fib%cur-fia%cur)  ! array to help allocate 3d FloatArrays to zero -AZ
      REAL :: zeros2d_s(nzp,nprc)  ! array to help allocate 3d FloatArrays to zero -AZ

      zeros0d(:) = 0.
      zeros1d(:) = 0.
      zeros2d_aa(:,:) = 0.
      zeros2d_ab(:,:) = 0.
      zeros2d_ca(:,:) = 0.
      zeros2d_cb(:,:) = 0.
      zeros2d_p(:,:) = 0.
      zeros2d_ia(:,:) = 0.
      zeros2d_ib(:,:) = 0.
      zeros2d_s(:,:) = 0.


      a_time = FloatArray1D(zeros0d,store=.TRUE.)
      a_cfl = FloatArray1D(zeros0d,store=.TRUE.)
      a_maxdiv = FloatArray1D(zeros0d,store=.TRUE.)
      a_zi1_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_zi2_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_zi3_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_vtke = FloatArray1D(zeros0d,store=.TRUE.)
      a_sfcbflx = FloatArray1D(zeros0d,store=.TRUE.)
      a_wmax = FloatArray1D(zeros0d,store=.TRUE.)
      a_tsrf = FloatArray1D(zeros0d,store=.TRUE.)
      a_ustar = FloatArray1D(zeros0d,store=.TRUE.)
      a_shf_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_lhf_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_zi_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_lwp_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_lwp_var = FloatArray1D(zeros0d,store=.TRUE.)
      a_zc = FloatArray1D(zeros0d,store=.TRUE.)
      a_zb = FloatArray1D(zeros0d,store=.TRUE.)
      a_cfrac = FloatArray1D(zeros0d,store=.TRUE.)
      a_lmax = FloatArray1D(zeros0d,store=.TRUE.)
      a_albedo = FloatArray1D(zeros0d,store=.TRUE.)
      a_rwp_bar = FloatArray1D(zeros0d,store=.TRUE.)
      a_prcp = FloatArray1D(zeros0d,store=.TRUE.)
      a_pfrac = FloatArray1D(zeros0d,store=.TRUE.)
      a_CCN = FloatArray1D(zeros0d,store=.TRUE.)
      a_nrain = FloatArray1D(zeros0d,store=.TRUE.)
      a_nrcnt = FloatArray1D(zeros0d,store=.TRUE.)
      a_nccnt = FloatArray1D(zeros0d,store=.TRUE.)
      a_prcp_bc = FloatArray1D(zeros0d,store=.TRUE.)


      a_zt = FloatArray1D(zeros1d,store=.TRUE.)
      a_zm = FloatArray1D(zeros1d,store=.TRUE.)
      a_dn0 = FloatArray1D(zeros1d,store=.TRUE.)
      a_u0 = FloatArray1D(zeros1d,store=.TRUE.)
      a_v0 = FloatArray1D(zeros1d,store=.TRUE.)
      a_fsttm = FloatArray1D(zeros1d,store=.TRUE.)
      a_lsttm = FloatArray1D(zeros1d,store=.TRUE.)
      a_nsmp = FloatArray1D(zeros1d,store=.TRUE.)
      a_u = FloatArray1D(zeros1d,store=.TRUE.)
      a_v = FloatArray1D(zeros1d,store=.TRUE.)
      a_theta_1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_p = FloatArray1D(zeros1d,store=.TRUE.)
      a_u_2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_v_2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_w_2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_theta_2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_w_3 = FloatArray1D(zeros1d,store=.TRUE.)
      a_theta_3 = FloatArray1D(zeros1d,store=.TRUE.)
      a_tot_tw = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_tw = FloatArray1D(zeros1d,store=.TRUE.)
      a_tot_uw = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_uw = FloatArray1D(zeros1d,store=.TRUE.)
      a_tot_vw = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_vw = FloatArray1D(zeros1d,store=.TRUE.)
      a_tot_ww = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_ww = FloatArray1D(zeros1d,store=.TRUE.)
      a_km = FloatArray1D(zeros1d,store=.TRUE.)
      a_kh = FloatArray1D(zeros1d,store=.TRUE.)
      a_lmbd = FloatArray1D(zeros1d,store=.TRUE.)
      a_lmbde = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_tke = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_boy = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_shr = FloatArray1D(zeros1d,store=.TRUE.)
      a_boy_prd = FloatArray1D(zeros1d,store=.TRUE.)
      a_shr_prd = FloatArray1D(zeros1d,store=.TRUE.)
      a_trans = FloatArray1D(zeros1d,store=.TRUE.)
      a_diss = FloatArray1D(zeros1d,store=.TRUE.)
      a_dff_u = FloatArray1D(zeros1d,store=.TRUE.)
      a_dff_v = FloatArray1D(zeros1d,store=.TRUE.)
      a_dff_w = FloatArray1D(zeros1d,store=.TRUE.)
      a_adv_u = FloatArray1D(zeros1d,store=.TRUE.)
      a_adv_v = FloatArray1D(zeros1d,store=.TRUE.)
      a_adv_w = FloatArray1D(zeros1d,store=.TRUE.)
      a_prs_u = FloatArray1D(zeros1d,store=.TRUE.)
      a_prs_v = FloatArray1D(zeros1d,store=.TRUE.)
      a_prs_w = FloatArray1D(zeros1d,store=.TRUE.)
      a_prd_uw = FloatArray1D(zeros1d,store=.TRUE.)
      a_storage = FloatArray1D(zeros1d,store=.TRUE.)
      a_q = FloatArray1D(zeros1d,store=.TRUE.)
      a_q_2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_q_3 = FloatArray1D(zeros1d,store=.TRUE.)
      a_tot_qw = FloatArray1D(zeros1d,store=.TRUE.)
      a_sfs_qw = FloatArray1D(zeros1d,store=.TRUE.)
      a_rflx1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_rflx2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_sflx1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_sflx2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_l = FloatArray1D(zeros1d,store=.TRUE.)
      a_l_2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_l_3 = FloatArray1D(zeros1d,store=.TRUE.)
      a_tot_lw = FloatArray1D(zeros1d,store=.TRUE.)
      a_sed_lw = FloatArray1D(zeros1d,store=.TRUE.)
      a_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_cnt_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_w_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_t_cs1  = FloatArray1D(zeros1d,store=.TRUE.)
      a_tv_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_rt_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_rc_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_wt_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_wtv_cs1= FloatArray1D(zeros1d,store=.TRUE.)
      a_wrt_cs1 = FloatArray1D(zeros1d,store=.TRUE.)
      a_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_cnt_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_w_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_t_cs2  = FloatArray1D(zeros1d,store=.TRUE.)
      a_tv_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_rt_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_rc_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_wt_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_wtv_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_wrt_cs2 = FloatArray1D(zeros1d,store=.TRUE.)
      a_Nc = FloatArray1D(zeros1d,store=.TRUE.)
      a_Nr = FloatArray1D(zeros1d,store=.TRUE.)
      a_rr = FloatArray1D(zeros1d,store=.TRUE.)
      a_rrate = FloatArray1D(zeros1d,store=.TRUE.)
      a_evap = FloatArray1D(zeros1d,store=.TRUE.)
      a_frc_prc = FloatArray1D(zeros1d,store=.TRUE.)
      a_prc_prc = FloatArray1D(zeros1d,store=.TRUE.)
      a_frc_ran = FloatArray1D(zeros1d,store=.TRUE.)
      a_hst_srf = FloatArray1D(zeros1d,store=.TRUE.)
      a_sw_up = FloatArray1D(zeros1d,store=.TRUE.)
      a_sw_down = FloatArray1D(zeros1d,store=.TRUE.)
      a_lw_up = FloatArray1D(zeros1d,store=.TRUE.)
      a_lw_down = FloatArray1D(zeros1d,store=.TRUE.)

      IF(level >= 4) THEN

         a_Nc_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_Na_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_Na_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_SO4_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_SO4_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_SO4_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_OC_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_OC_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_OC_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_BC_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_BC_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_BC_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_DU_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_DU_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_DU_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_SS_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_SS_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_SS_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_NH_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_NH_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_NH_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_NO_ic  = FloatArray1D(zeros0d,store=.TRUE.)
         a_NO_int = FloatArray1D(zeros0d,store=.TRUE.)
         a_NO_oc  = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmH2Oae = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmH2Ocl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmH2Opr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmSO4dr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmSO4cl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmSO4pr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmOCdr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmOCcl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmOCpr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmBCdr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmBCcl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmBCpr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmDUdr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmDUcl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmDUpr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmSSdr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmSScl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmSSpr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmNHdr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmNHcl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmNHpr = FloatArray1D(zeros0d,store=.TRUE.)

         a_rmNOdr = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmNOcl = FloatArray1D(zeros0d,store=.TRUE.)
         a_rmNOpr = FloatArray1D(zeros0d,store=.TRUE.)

         a_aea     = FloatArray1D(zeros1d,store=.TRUE.)
         a_aeb     = FloatArray1D(zeros1d,store=.TRUE.)
         a_cla     = FloatArray1D(zeros1d,store=.TRUE.)
         a_clb     = FloatArray1D(zeros1d,store=.TRUE.)
         a_prc     = FloatArray1D(zeros1d,store=.TRUE.)

         a_Naa   = FloatArray1D(zeros1d,store=.TRUE.)
         a_Nab   = FloatArray1D(zeros1d,store=.TRUE.)
         a_Nca   = FloatArray1D(zeros1d,store=.TRUE.)
         a_Ncb   = FloatArray1D(zeros1d,store=.TRUE.)
         a_Np    = FloatArray1D(zeros1d,store=.TRUE.)

         a_Rwaa  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Rwab  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Rwca  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Rwcb  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Rwp   = FloatArray1D(zeros1d,store=.TRUE.)

         a_cSO4a = FloatArray1D(zeros1d,store=.TRUE.)
         a_cSO4c = FloatArray1D(zeros1d,store=.TRUE.)
         a_cSO4p = FloatArray1D(zeros1d,store=.TRUE.)

         a_cOCa  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cOCc  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cOCp  = FloatArray1D(zeros1d,store=.TRUE.)

         a_cBCa  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cBCc  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cBCp  = FloatArray1D(zeros1d,store=.TRUE.)

         a_cDUa  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cDUc  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cDUp  = FloatArray1D(zeros1d,store=.TRUE.)

         a_cSSa  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cSSc  = FloatArray1D(zeros1d,store=.TRUE.)
         a_cSSp  = FloatArray1D(zeros1d,store=.TRUE.)

         a_cNHa = FloatArray1D(zeros1d,store=.TRUE.)
         a_cNHc = FloatArray1D(zeros1d,store=.TRUE.)
         a_cNHp = FloatArray1D(zeros1d,store=.TRUE.)

         a_cNOa = FloatArray1D(zeros1d,store=.TRUE.)
         a_cNOc = FloatArray1D(zeros1d,store=.TRUE.)
         a_cNOp = FloatArray1D(zeros1d,store=.TRUE.)

         a_cH2Oa = FloatArray1D(zeros1d,store=.TRUE.)
         a_cH2Oc = FloatArray1D(zeros1d,store=.TRUE.)
         a_cH2Op = FloatArray1D(zeros1d,store=.TRUE.)

         a_prl    = FloatArray1D(zeros1d,store=.TRUE.)
         a_Prr   = FloatArray1D(zeros1d,store=.TRUE.)
         a_prv    = FloatArray1D(zeros1d,store=.TRUE.)
         a_pRH    = FloatArray1D(zeros1d,store=.TRUE.)

         a_Na_c  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Nc_c  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Np_c  = FloatArray1D(zeros1d,store=.TRUE.)
         a_Pcfrac= FloatArray1D(zeros1d,store=.TRUE.)
         a_clw_c = FloatArray1D(zeros1d,store=.TRUE.)
         a_thl_c = FloatArray1D(zeros1d,store=.TRUE.)


         a_Naba = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_SO4aa = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_OCaa = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_BCaa = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_DUaa = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_SSaa = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_NHaa = FloatArray2D(zeros2d_aa,store=.TRUE.)
         a_NOaa = FloatArray2D(zeros2d_aa,store=.TRUE.)

         a_Nabb = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_SO4ab = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_OCab = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_BCab = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_DUab = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_SSab = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_NHab = FloatArray2D(zeros2d_ab,store=.TRUE.)
         a_NOab = FloatArray2D(zeros2d_ab,store=.TRUE.)

         a_Ncba = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_SO4ca = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_OCca = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_BCca = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_DUca = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_SSca = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_NHca = FloatArray2D(zeros2d_ca,store=.TRUE.)
         a_NOca = FloatArray2D(zeros2d_ca,store=.TRUE.)

         a_Ncbb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_SO4cb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_OCcb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_BCcb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_DUcb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_SScb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_NHcb = FloatArray2D(zeros2d_cb,store=.TRUE.)
         a_NOcb = FloatArray2D(zeros2d_cb,store=.TRUE.)

         a_Npb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_SO4pb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_OCpb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_BCpb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_DUpb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_SSpb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_NHpb = FloatArray2D(zeros2d_p,store=.TRUE.)
         a_NOpb = FloatArray2D(zeros2d_p,store=.TRUE.)

         IF (level >= 5) THEN

            a_Ni_ic  = FloatArray1D(zeros0d,store=.TRUE.)
            a_Ni_ii  = FloatArray1D(zeros0d,store=.TRUE.)
            a_Ni_is = FloatArray1D(zeros0d,store=.TRUE.)
            a_Ns_ic  = FloatArray1D(zeros0d,store=.TRUE.)
            a_Ns_ii  = FloatArray1D(zeros0d,store=.TRUE.)
            a_Ns_is  = FloatArray1D(zeros0d,store=.TRUE.)

            a_Ri_ii  = FloatArray1D(zeros0d,store=.TRUE.)
            a_iwp_bar= FloatArray1D(zeros0d,store=.TRUE.)
            a_imax  = FloatArray1D(zeros0d,store=.TRUE.)
            a_nicnt  = FloatArray1D(zeros0d,store=.TRUE.)

            a_Rs_is  = FloatArray1D(zeros0d,store=.TRUE.)
            a_swp_bar= FloatArray1D(zeros0d,store=.TRUE.)
            a_smax  = FloatArray1D(zeros0d,store=.TRUE.)
            a_nscnt  = FloatArray1D(zeros0d,store=.TRUE.)

            a_rmSO4ic= FloatArray1D(zeros0d,store=.TRUE.)
            a_rmSO4sn= FloatArray1D(zeros0d,store=.TRUE.)
            a_rmOCic = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmOCsn = FloatArray1D(zeros0d,store=.TRUE.)

            a_rmBCic = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmBCsn = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmDUic = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmDUsn = FloatArray1D(zeros0d,store=.TRUE.)

            a_rmNOic = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmNOsn = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmNHic = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmNHsn = FloatArray1D(zeros0d,store=.TRUE.)

            a_rmSSic = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmSSsn = FloatArray1D(zeros0d,store=.TRUE.)
            a_rmH2Oic= FloatArray1D(zeros0d,store=.TRUE.)
            a_rmH2Osn= FloatArray1D(zeros0d,store=.TRUE.)

            a_sfrac  = FloatArray1D(zeros0d,store=.TRUE.)
            a_sprcp = FloatArray1D(zeros0d,store=.TRUE.)

            a_ica     = FloatArray1D(zeros1d,store=.TRUE.)
            a_icb     = FloatArray1D(zeros1d,store=.TRUE.)
            a_snw     = FloatArray1D(zeros1d,store=.TRUE.)

            a_Nia   = FloatArray1D(zeros1d,store=.TRUE.)
            a_Nib   = FloatArray1D(zeros1d,store=.TRUE.)
            a_Ns    = FloatArray1D(zeros1d,store=.TRUE.)
            a_Rwia  = FloatArray1D(zeros1d,store=.TRUE.)
            a_Rwib  = FloatArray1D(zeros1d,store=.TRUE.)
            a_Rws   = FloatArray1D(zeros1d,store=.TRUE.)

            a_cSO4i = FloatArray1D(zeros1d,store=.TRUE.)
            a_cSO4s = FloatArray1D(zeros1d,store=.TRUE.)
            a_cOCi  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cOCs  = FloatArray1D(zeros1d,store=.TRUE.)

            a_cBCi  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cBCs  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cDUi  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cDUs  = FloatArray1D(zeros1d,store=.TRUE.)

            a_cSSi  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cSSs  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cNOi  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cNOs  = FloatArray1D(zeros1d,store=.TRUE.)

            a_cNHi  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cNHs  = FloatArray1D(zeros1d,store=.TRUE.)
            a_cH2Oi = FloatArray1D(zeros1d,store=.TRUE.)
            a_cH2Os = FloatArray1D(zeros1d,store=.TRUE.)

            a_pri    = FloatArray1D(zeros1d,store=.TRUE.)
            a_prs    = FloatArray1D(zeros1d,store=.TRUE.)
            a_pRHi   = FloatArray1D(zeros1d,store=.TRUE.)
            a_srate  = FloatArray1D(zeros1d,store=.TRUE.)

            a_Niba = FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_SO4ia= FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_OCia = FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_BCia = FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_DUia = FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_SSia = FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_NOia = FloatArray2D(zeros2d_ia,store=.TRUE.)
            a_NHia = FloatArray2D(zeros2d_ia,store=.TRUE.)

            a_Nibb = FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_SO4ib= FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_OCib = FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_BCib = FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_DUib = FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_SSib = FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_NOib = FloatArray2D(zeros2d_ib,store=.TRUE.)
            a_NHib = FloatArray2D(zeros2d_ib,store=.TRUE.)

            a_Nsb  = FloatArray2D(zeros2d_s,store=.TRUE.)
            a_SO4sb= FloatArray2D(zeros2d_s,store=.TRUE.)
            a_OCsb = FloatArray2D(zeros2d_s,store=.TRUE.)
            a_BCsb = FloatArray2D(zeros2d_s,store=.TRUE.)

            a_DUsb = FloatArray2D(zeros2d_s,store=.TRUE.)
            a_SSsb = FloatArray2D(zeros2d_s,store=.TRUE.)
            a_NOsb = FloatArray2D(zeros2d_s,store=.TRUE.)
            a_NHsb = FloatArray2D(zeros2d_s,store=.TRUE.)

        END IF

      END IF

   END SUBROUTINE

   SUBROUTINE create_fields
    USE classFieldArray
      USE mo_structured_datatypes
      USE grid,          ONLY : nxp, nyp, iradtyp, prtcl
      USE mpi_interface, ONLY : myid, ver, author, info
      USE mo_submctl,    ONLY : nprc, fn2a,fn2b,fca,fcb,fra,fia,fib,fsa
      USE class_ComponentIndex, ONLY : IsUsed

      INTEGER :: i, e, n, nnum

      REAL :: zeros1d(1)  ! array to help allocate 1d FloatArrays to zero -AZ
      REAL :: zeros2d(2,2)  ! array to help allocate 2d FloatArrays to zero -AZ
      REAL :: zeros3d(3,3,3)  ! array to help allocate 3d FloatArrays to zero -AZ

      TYPE(FloatArray1d), TARGET :: dumarr1d
      TYPE(FloatArray2d), TARGET :: dumarr2d
      TYPE(FloatArray3d), TARGET :: dumarr3d
      TYPE(ArrayElement), POINTER :: ArrEl

      LOGICAL :: lso4, loc, lbc, lss, ldu, lnh3, lno3, lbb
      lbb = salsa_b_bins
      lso4 = .FALSE.
      ldu = .FALSE.
      loc = .FALSE.
      lbc = .FALSE.
      lss = .FALSE.
      lnh3 = .FALSE.
      lno3 = .FALSE.

      f1 = FieldArray()
      f2 = FieldArray()
      f1sbulk = FieldArray()
      f1lvl5 = FieldArray()
      f2sbulk = FieldArray()
      f2lvl5 = FieldArray()
      f2aea = FieldArray()
      f2aeb = FieldArray()
      f2cla = FieldArray()
      f2clb = FieldArray()
      f2prc = FieldArray()
      f2ica = FieldArray()
      f2icb = FieldArray()
      f2snw = FieldArray()

      !Initialize dummy arrays
      zeros1d(:) = 0.
      zeros2d(:,:) = 0.

      dumarr1d = FloatArray1D(zeros1d,store=.TRUE.)
      dumarr2d = FloatArray2D(zeros2d,store=.TRUE.)
      dumarr3d = FloatArray3D(zeros3d,store=.TRUE.)

      !Point dummy pointers to dummy arrays
      dummy => dumarr1d
      dummy2 => dumarr2d

      ! Create field arrays containing output values

      ! ** OUTPUT STATISTICS **
      dummy => a_time
      CALL f1%NewField("time","Time","s","time",.TRUE.,dummy)
      dummy => a_cfl
      CALL f1%NewField("cfl","Courant number","-","time",.TRUE.,dummy)
      dummy => a_maxdiv
      CALL f1%NewField("maxdiv","Maximum divergence","1/s","time",.TRUE.,dummy)
      dummy => a_zi1_bar
      CALL f1%NewField("zi1_bar","Height of maximum theta gradient","m","time",.TRUE.,dummy)
      dummy => a_zi2_bar
      CALL f1%NewField("zi2_bar","Height of maximum theta variance","m","time",.TRUE.,dummy)
      dummy => a_zi3_bar
      CALL f1%NewField("zi3_bar","Height of minimum buoyancy flux","m","time",.TRUE.,dummy)
      dummy => a_vtke
      CALL f1%NewField("vtke","Vertical integral of total TKE","kg/s","time",.TRUE.,dummy)
      dummy => a_sfcbflx
      CALL f1%NewField("sfcbflx","Surface Buoyancy Flux","m/s^2","time",.TRUE.,dummy)
      dummy => a_wmax
      CALL f1%NewField("wmax","Maximum vertical velocity","m/s","time",.TRUE.,dummy)
      dummy => a_tsrf
      CALL f1%NewField("tsrf","Surface temperature","K","time",.TRUE.,dummy)
      dummy => a_ustar
      CALL f1%NewField("ustar","Surface friction velocity","m/s","time",.TRUE.,dummy)
      dummy => a_shf_bar
      CALL f1%NewField("shf_bar","Sensible heat flux","W/m^2","time",.TRUE.,dummy)
      dummy => a_lhf_bar
      CALL f1%NewField("lhf_bar","Latent heat flux","W/m^2","time",.TRUE.,dummy)
      dummy => a_zi_bar
      CALL f1%NewField("zi_bar","Height of maximum scalar gradient","m","time",.TRUE.,dummy)
      dummy => a_lwp_bar
      CALL f1%NewField("lwp_bar","Liquid-water path","kg/m^2","time",.TRUE.,dummy)
      dummy => a_lwp_var
      CALL f1%NewField("lwp_var","Liquid-water path variance","kg^2/m^4","time",.TRUE.,dummy)
      dummy => a_zc
      CALL f1%NewField("zc","Cloud-top height","m","time",.TRUE.,dummy)
      dummy => a_zb
      CALL f1%NewField("zb","Cloud-base height","m","time",.TRUE.,dummy)
      dummy => a_cfrac
      CALL f1%NewField("cfrac","Cloud fraction","-","time",.TRUE.,dummy)
      dummy => a_lmax
      CALL f1%NewField("lmax","Maximum liquid water mixing ratio","kg/kg","time",.TRUE.,dummy)
      dummy => a_albedo
      CALL f1%NewField("albedo","Reflected (TOA) shortwave radiation","-","time",.TRUE.,dummy)
      dummy => a_rwp_bar
      CALL f1%NewField("rwp_bar","Rain-water path","kg/m^2","time",.TRUE.,dummy)
      dummy => a_prcp
      CALL f1%NewField("prcp","Surface precipitation rate","W/m^2","time",.TRUE.,dummy)
      dummy => a_pfrac
      CALL f1%NewField("pfrac","Surface precipitation fraction","-","time",.TRUE.,dummy)
      dummy => a_CCN
      CALL f1%NewField("CCN","Cloud condensation nuclei","kg^-1","time",.TRUE.,dummy)
      dummy => a_nrain
      CALL f1%NewField("nrain","Conditionally sampled rain number mixing ratio","kg^-1","time",.TRUE.,dummy)
      dummy => a_nrcnt
      CALL f1%NewField("nrcnt","Rain cell counts","#","time",.TRUE.,dummy)
      dummy => a_nccnt
      CALL f1%NewField("nccnt","Cloud cell counts","#","time",.TRUE.,dummy)
      dummy => a_prcp_bc
      CALL f1%NewField("prcp_bc","Below cloud precipitation rate","W/m^2","time",.TRUE.,dummy)

      ! ** PROFILE OUTPUT **
      dummy => a_time
      CALL f2%NewField("time","Time","s","time",.TRUE.,dummy)
      dummy => a_zt
      CALL f2%NewField("zt","Vertical displacement of cell centers","m","zt",.TRUE.,dummy)
      dummy => a_zm
      CALL f2%NewField("zm","Vertical displacement of cell edges","m","zm",.TRUE.,dummy)
      dummy => a_dn0
      CALL f2%NewField("dn0","Base-state density","kg/m^3","zt",.TRUE.,dummy)
      dummy => a_u0
      CALL f2%NewField("u0","Geostrophic zonal wind","m/s","zt",.TRUE.,dummy)
      dummy => a_v0
      CALL f2%NewField("v0","Geostrophic meridional wind","m/s","zt",.TRUE.,dummy)
      dummy => a_fsttm
      CALL f2%NewField("fsttm","First sample time","s","time",.TRUE.,dummy)
      dummy => a_lsttm
      CALL f2%NewField("lsttm","Last sample time","s","time",.TRUE.,dummy)
      dummy => a_nsmp
      CALL f2%NewField("nsmp","Number of samples","#","time",.TRUE.,dummy)
      dummy => a_u
      CALL f2%NewField("u","Zonal wind","m/s","mttt",.TRUE.,dummy)
      dummy => a_v
      CALL f2%NewField("v","Meridional wind","m/s","tmtt",.TRUE.,dummy)
      dummy => a_theta_1
      CALL f2%NewField("theta","Liquid water potential temperature","K","tttt",.TRUE.,dummy)
      dummy => a_p
      CALL f2%NewField("p","Pressure","Pa","tttt",.TRUE.,dummy)
      dummy => a_u_2
      CALL f2%NewField("u_2","Variance of u wind","m^2/s^2","tttt",.TRUE.,dummy)
      dummy => a_v_2
      CALL f2%NewField("v_2","Variance of v wind","m^2/s^2","tttt",.TRUE.,dummy)
      dummy => a_w_2
      CALL f2%NewField("w_2","Second raw moment of w wind","m^2/s^2","ttmt",.TRUE.,dummy)
      dummy => a_theta_2
      CALL f2%NewField("theta_2","Variance of theta_1","K^2","tttt",.TRUE.,dummy)
      dummy => a_w_3
      CALL f2%NewField("w_3","Third raw moment of w wind","m^3/s^3","ttmt",.TRUE.,dummy)
      dummy => a_theta_3
      CALL f2%NewField("theta_3","Third moment of theta_1","K^3","tttt",.TRUE.,dummy)
      dummy => a_tot_tw
      CALL f2%NewField("tot_tw","Total vertical flux of theta","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_sfs_tw
      CALL f2%NewField("sfs_tw","Sub-filter scale vertical flux of theta","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_tot_uw
      CALL f2%NewField("tot_uw","Total vertical flux of theta","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_sfs_uw
      CALL f2%NewField("sfs_uw","Sub-filter scale vertical flux of u-wind","m^2/s^2","ttmt",.TRUE.,dummy)
      dummy => a_tot_vw
      CALL f2%NewField("tot_vw","Total vertical flux of v-wind","m^2/s^2","ttmt",.TRUE.,dummy)
      dummy => a_sfs_vw
      CALL f2%NewField("sfs_vw","SGS vertical flux of v-wind","m^2/s^2","ttmt",.TRUE.,dummy)
      dummy => a_tot_ww
      CALL f2%NewField("tot_ww","Total vertical flux of w-wind","m^2/s^2","ttmt",.TRUE.,dummy)
      dummy => a_sfs_ww
      CALL f2%NewField("sfs_ww","SGS vertical flux of w-wind","m^2/s^2","ttmt",.TRUE.,dummy)

      dummy => a_km
      CALL f2%NewField("km","Eddy viscosity","m^2/s","ttmt",.TRUE.,dummy)
      dummy => a_kh
      CALL f2%NewField("kh","Eddy diffusivity","m^2/s","ttmt",.TRUE.,dummy)
      dummy => a_lmbd
      CALL f2%NewField("lmbd","Mixing lengthscale","m","ttmt",.TRUE.,dummy)
      dummy => a_lmbde
      CALL f2%NewField("lmbde","Dissipation lengthscale","m","ttmt",.TRUE.,dummy)
      dummy => a_sfs_tke
      CALL f2%NewField("sfs_tke","Subfilter scale TKE","m^2/s^2","ttmt",.TRUE.,dummy)
      dummy => a_sfs_boy
      CALL f2%NewField("sfs_boy","Subfilter buoyancy production of TKE","m^2/s^3","ttmt",.TRUE.,dummy)
      dummy => a_sfs_shr
      CALL f2%NewField("sfs_shr","Shear production of SGS TKE","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_boy_prd
      CALL f2%NewField("boy_prd","Buoyancy production of resolved TKE","m^2/s^3","ttmt",.TRUE.,dummy)
      dummy => a_shr_prd
      CALL f2%NewField("shr_prd","Shear production of resolved TKE","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_trans
      CALL f2%NewField("trans","Net transport of resolved TKE","m^2/s^3","ttmt",.TRUE.,dummy)
      dummy => a_diss
      CALL f2%NewField("diss","Dissipation rate of resolved TKE","m^2/s^3","ttmt",.TRUE.,dummy)
      dummy => a_dff_u
      CALL f2%NewField("dff_u","u(du/dt) from diffusion","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_dff_v
      CALL f2%NewField("dff_v","v(dv/dt) from diffusion","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_dff_w
      CALL f2%NewField("dff_w","w(dw/dt) from diffusion","m^2/s^3","ttmt",.TRUE.,dummy)
      dummy => a_adv_u
      CALL f2%NewField("adv_u","u(du/dt) from advection","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_adv_v
      CALL f2%NewField("adv_v","v(dv/dt) from advection","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_adv_w
      CALL f2%NewField("adv_w","w(dw/dt) from advection","m^2/s^3","ttmt",.TRUE.,dummy)
      dummy => a_prs_u
      CALL f2%NewField("prs_u","u(du/dt) from pressure","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_prs_v
      CALL f2%NewField("prs_v","v(dv/dt) from pressure","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_prs_w
      CALL f2%NewField("prs_w","w(dw/dt) from pressure","m^2/s^3","ttmt",.TRUE.,dummy)

      dummy => a_prd_uw
      CALL f2%NewField("prd_uw","uw shear production","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_storage
      CALL f2%NewField("storage","Rate of increase of resolved TKE","m^2/s^3","tttt",.TRUE.,dummy)
      dummy => a_q
      CALL f2%NewField("q","Total water mixing ratio","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_q_2
      CALL f2%NewField("q_2","Variance of total water","-","tttt",.TRUE.,dummy)
      dummy => a_q_3
      CALL f2%NewField("q_3","Third moment of total water","-","tttt",.TRUE.,dummy)
      dummy => a_tot_qw
      CALL f2%NewField("tot_qw","Total vertical flux of q","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_sfs_qw
      CALL f2%NewField("sfs_qw","Sub-filter scale vertical flux of q","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_rflx1
      CALL f2%NewField("rflx","Total radiative flux","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_rflx2
      CALL f2%NewField("rflx2","Variance of total radiative flux","W^2/m^4","ttmt",.TRUE.,dummy)
      dummy => a_sflx1
      CALL f2%NewField("sflx","Shortwave radiative flux","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_sflx2
      CALL f2%NewField("sflx2","Variance of shortwave radiative flux","W^2/m^4","ttmt",.TRUE.,dummy)
      dummy => a_l
      CALL f2%NewField("l","Liquid water mixing ratio","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_l_2
      CALL f2%NewField("l_2","Variance of liquid water mixing ratio","-","tttt",.TRUE.,dummy)
      dummy => a_l_3
      CALL f2%NewField("l_3","Third moment of liquid water mixing ratio","-","tttt",.TRUE.,dummy)
      dummy => a_tot_lw
      CALL f2%NewField("tot_lw","Resolved turbulent flux of liquid water mixing ratio","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_sed_lw
      CALL f2%NewField("sed_lw","Sedimentation flux of r_l","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_cs1
      CALL f2%NewField("cs1","Fraction of cloudy columns (cs1)","-","tttt",.TRUE.,dummy)
      dummy => a_cnt_cs1
      CALL f2%NewField("cnt_cs1","Number of cloudy columns (cs1)","#","tttt",.TRUE.,dummy)
      dummy => a_w_cs1
      CALL f2%NewField("w_cs1","Conditional average of w over cs1","m/s","ttmt",.TRUE.,dummy)
      dummy => a_t_cs1
      CALL f2%NewField("t_cs1 ","Conditional average of theta_l over cs1","K","ttmt",.TRUE.,dummy)
      dummy => a_tv_cs1
      CALL f2%NewField("tv_cs1","Conditional average of theta_v over cs1","K","tttt",.TRUE.,dummy)
      dummy => a_rt_cs1
      CALL f2%NewField("rt_cs1","Conditional average of rt over cs1","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_rc_cs1
      CALL f2%NewField("rc_cs1","Conditional average of rl over cs1","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_wt_cs1
      CALL f2%NewField("wt_cs1","Covariance of wtheta_l flux and cs1","K*m/s","ttmt",.TRUE.,dummy)
      dummy => a_wtv_cs1
      CALL f2%NewField("wtv_cs1","Covariance of wtheta_v flux and cs1","K*m/s","ttmt",.TRUE.,dummy)
      dummy => a_wrt_cs1
      CALL f2%NewField("wrt_cs1","Covariance of wr_t flux and cs1","kg/kg*m/s","ttmt",.TRUE.,dummy)

      dummy => a_cs2
      CALL f2%NewField("cs2","Fraction of cloud core columns (cs2)","-","tttt",.TRUE.,dummy)
      dummy => a_cnt_cs2
      CALL f2%NewField("cnt_cs2","Number of cloud core columns (cs2)","#","tttt",.TRUE.,dummy)
      dummy => a_w_cs2
      CALL f2%NewField("w_cs2","Conditional average of w over cs2","m/s","ttmt",.TRUE.,dummy)
      dummy => a_t_cs2
      CALL f2%NewField("t_cs2 ","Conditional average of theta_l over cs2","K","tttt",.TRUE.,dummy)
      dummy => a_tv_cs2
      CALL f2%NewField("tv_cs2","Conditional average of theta_v over cs2","K","tttt",.TRUE.,dummy)
      dummy => a_rt_cs2
      CALL f2%NewField("rt_cs2","Conditional average of rt over cs2","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_rc_cs2
      CALL f2%NewField("rc_cs2","Conditional average of rl over cs2","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_wt_cs2
      CALL f2%NewField("wt_cs2","Covariance of wtheta_l flux and cs2","K*m/s","ttmt",.TRUE.,dummy)
      dummy => a_wtv_cs2
      CALL f2%NewField("wtv_cs2","Covariance of wtheta_v flux and cs2","K*m/s","ttmt",.TRUE.,dummy)
      dummy => a_wrt_cs2
      CALL f2%NewField("wrt_cs2","Covariance of wr_t flux and cs2","kg/kg*m/s","ttmt",.TRUE.,dummy)


      dummy => a_Nc
      CALL f2%NewField("Nc","Cloud droplet number concentration","kg^-1","tttt",.TRUE.,dummy)
      dummy => a_Nr
      CALL f2%NewField("Nr","Rain drop number concentration","kg^-1","tttt",.TRUE.,dummy)
      dummy => a_rr
      CALL f2%NewField("rr","Rain water mixing ratio","kg/kg","tttt",.TRUE.,dummy)
      dummy => a_rrate
      CALL f2%NewField("rrate","Precipitation flux (positive downward)","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_evap
      CALL f2%NewField("evap","Net evap of rain-water","s^-1","tttt",.TRUE.,dummy)
      dummy => a_frc_prc
      CALL f2%NewField("frc_prc","Conditionally sampled rain fraction","-","ttmt",.TRUE.,dummy)
      dummy => a_prc_prc
      CALL f2%NewField("prc_prc","Conditionally sampled precipitation flux","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_frc_ran
      CALL f2%NewField("frc_ran","Rain water fraction","-","tttt",.TRUE.,dummy)
      dummy => a_hst_srf
      CALL f2%NewField("hst_srf","Histogram of surface rain rates","-","tttt",.TRUE.,dummy)
      dummy => a_sw_up
      CALL f2%NewField("sw_up","Upwelling shortwave radiation","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_sw_down
      CALL f2%NewField("sw_down","Downwelling shortwave radiation","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_lw_up
      CALL f2%NewField("lw_up","Upwelling longwave radiation","W/m^2","ttmt",.TRUE.,dummy)
      dummy => a_lw_down
      CALL f2%NewField("lw_down","Downwelling longwave radiation","W/m^2","ttmt",.TRUE.,dummy)

      IF (level >= 4) THEN

         IF (IsUsed(prtcl,'SO4')) lso4 = .TRUE.

         IF (IsUsed(prtcl,'OC')) loc = .TRUE.

         IF (IsUsed(prtcl,'BC')) lbc = .TRUE.

         IF (IsUsed(prtcl,'DU')) ldu = .TRUE.

         IF (IsUsed(prtcl,'SS')) lss = .TRUE.

         IF (IsUsed(prtcl,'NH')) lnh3 = .TRUE.

         IF (IsUsed(prtcl,'NO')) lno3 = .TRUE.
       ! ** BULK TEMPORAL STATISTICS FOR SALSA **
         ! Number concentrations
         dummy => a_Nc_ic
         CALL f1sbulk%NewField("Nc_ic","In-cloud CDNC","kg^-1","time",.TRUE.,dummy)
         dummy => a_Na_int
         CALL f1sbulk%NewField("Na_int","In-cloud interstitial aerosol number concentration","kg^-1","time",.TRUE.,dummy)
         dummy => a_Na_oc
         CALL f1sbulk%NewField("Na_oc","Aerosol number concentration outside clous","kg^-1","time",.TRUE.,dummy)

         dummy => a_SO4_ic
         CALL f1sbulk%NewField("SO4_ic","Cloud droplet SO4 mass mixing ratio","kg/kg","time",lso4,dummy)
         dummy => a_SO4_int
         CALL f1sbulk%NewField("SO4_int","SO4 mass mixing ratio in interstitial aerosols","kg/kg","time",lso4,dummy)
         dummy => a_SO4_oc
         CALL f1sbulk%NewField("SO4_oc","Aerosol SO4 mass mixing ratio outside clouds","kg/kg","time",lso4,dummy)

         dummy => a_OC_ic
         CALL f1sbulk%NewField("OC_ic","Cloud droplet OC mass mixing ratio","kg/kg","time",loc,dummy)
         dummy => a_OC_int
         CALL f1sbulk%NewField("OC_int","OC mass mixing ratio in interstitial aerosols","kg/kg","time",loc,dummy)
         dummy => a_OC_oc
         CALL f1sbulk%NewField("OC_oc","Aerosol OC mass mixing ratio outside clouds","kg/kg","time",loc,dummy)

         dummy => a_BC_ic
         CALL f1sbulk%NewField("BC_ic","Cloud droplet BC mass mixing ratio","kg/kg","time",lbc,dummy)
         dummy => a_BC_int
         CALL f1sbulk%NewField("BC_int","BC mass mixing ratio in interstitial aerosols","kg/kg","time",lbc,dummy)
         dummy => a_BC_oc
         CALL f1sbulk%NewField("BC_oc","Aerosol BC mass mixing ratio outside clouds","kg/kg","time",lbc,dummy)

         dummy => a_DU_ic
         CALL f1sbulk%NewField("DU_ic","Cloud droplet DU mass mixing ratio","kg/kg","time",ldu,dummy)
         dummy => a_DU_int
         CALL f1sbulk%NewField("DU_int","DU mass mixing ratio in interstitial aerosols","kg/kg","time",ldu,dummy)
         dummy => a_DU_oc
         CALL f1sbulk%NewField("DU_oc","Aerosol DU mass mixing ratio outside clouds","kg/kg","time",ldu,dummy)

         dummy => a_SS_ic
         CALL f1sbulk%NewField("SS_ic","Cloud droplet SS mass mixing ratio","kg/kg","time",lss,dummy)
         dummy => a_SS_int
         CALL f1sbulk%NewField("SS_int","SS mass mixing ratio in interstitial aerosols","kg/kg","time",lss,dummy)
         dummy => a_SS_oc
         CALL f1sbulk%NewField("SS_oc","Aerosol SS mass mixing ratio outside clouds","kg/kg","time",lss,dummy)

         dummy => a_NO_ic
         CALL f1sbulk%NewField("NO_ic","Cloud droplet NO3 mass mixing ratio","kg/kg","time",lno3,dummy)
         dummy => a_NO_int
         CALL f1sbulk%NewField("NO_int","NO3 mass mixing ratio in interstitial aerosols","kg/kg","time",lno3,dummy)
         dummy => a_NO_oc
         CALL f1sbulk%NewField("NO_oc","Aerosol NO3 mass mixing ratio outside clouds","kg/kg","time",lno3,dummy)

         dummy => a_NH_ic
         CALL f1sbulk%NewField("NH_ic","Cloud droplet NH3 mass mixing ratio","kg/kg","time",lnh3,dummy)
         dummy => a_NH_int
         CALL f1sbulk%NewField("NH_int","NH3 mass mixing ratio in interstitial aerosols","kg/kg","time",lnh3,dummy)
         dummy => a_NH_oc
         CALL f1sbulk%NewField("NH_oc","Aerosol NH3 mass mixing ratio outside clouds","kg/kg","time",lnh3,dummy)

         dummy => a_rmSO4dr
         CALL f1sbulk%NewField("rmSO4dr","Aerosol deposition of SO4","kg/m^2/s","time",lso4,dummy)
         dummy => a_rmSO4cl
         CALL f1sbulk%NewField("rmSO4cl","Cloud deposition of SO4","kg/^2/s","time",lso4,dummy)
         dummy => a_rmSO4pr
         CALL f1sbulk%NewField("rmSO4pr","Precipitation deposition of SO4","kg/m^2/s","time",lso4,dummy)

         dummy => a_rmOCdr
         CALL f1sbulk%NewField("rmOCdr","Aerosol deposition of OC","kg/m^2/s","time",loc,dummy)
         dummy => a_rmOCcl
         CALL f1sbulk%NewField("rmOCcl","Cloud deposition of OC","kg/^2/s","time",loc,dummy)
         dummy => a_rmOCpr
         CALL f1sbulk%NewField("rmOCpr","Precipitation deposition of OC","kg/m^2/s","time",loc,dummy)

         dummy => a_rmBCdr
         CALL f1sbulk%NewField("rmBCdr","Aerosol deposition of BC","kg/m^2/s","time",lbc,dummy)
         dummy => a_rmBCcl
         CALL f1sbulk%NewField("rmBCcl","Cloud deposition of BC","kg/^2/s","time",lbc,dummy)
         dummy => a_rmBCpr
         CALL f1sbulk%NewField("rmBCpr","Precipitation deposition of BC","kg/m^2/s","time",lbc,dummy)

         dummy => a_rmDUdr
         CALL f1sbulk%NewField("rmDUdr","Aerosol deposition of DU","kg/m^2/s","time",ldu,dummy)
         dummy => a_rmDUcl
         CALL f1sbulk%NewField("rmDUcl","Cloud deposition of DU","kg/^2/s","time",ldu,dummy)
         dummy => a_rmDUpr
         CALL f1sbulk%NewField("rmDUpr","Precipitation deposition of DU","kg/m^2/s","time",ldu,dummy)

         dummy => a_rmSSdr
         CALL f1sbulk%NewField("rmSSdr","Aerosol deposition of SS","kg/m^2/s","time",lss,dummy)
         dummy => a_rmSScl
         CALL f1sbulk%NewField("rmSScl","Cloud deposition of SS","kg/^2/s","time",lss,dummy)
         dummy => a_rmSSpr
         CALL f1sbulk%NewField("rmSSpr","Precipitation deposition of SS","kg/m^2/s","time",lss,dummy)

         dummy => a_rmNOdr
         CALL f1sbulk%NewField("rmNOdr","Aerosol deposition of NO3","kg/m^2/s","time",lno3,dummy)
         dummy => a_rmNOcl
         CALL f1sbulk%NewField("rmNOcl","Cloud deposition of NO3","kg/^2/s","time",lno3,dummy)
         dummy => a_rmNOpr
         CALL f1sbulk%NewField("rmNOpr","Precipitation deposition of NO3","kg/m^2/s","time",lno3,dummy)

         dummy => a_rmNHdr
         CALL f1sbulk%NewField("rmNHdr","Aerosol deposition of NH3","kg/m^2/s","time",lnh3,dummy)
         dummy => a_rmNHcl
         CALL f1sbulk%NewField("rmNHcl","Cloud deposition of NH3","kg/^2/s","time",lnh3,dummy)
         dummy => a_rmNHpr
         CALL f1sbulk%NewField("rmNHpr","Precipitation deposition of NH3","kg/m^2/s","time",lnh3,dummy)

            ! Water removal temporal statistics
         dummy => a_rmH2Oae
         CALL f1sbulk%NewField("rmH2Oae","Deposition of H2O with aerosols","kg/m^2/s","time",.TRUE.,dummy)
         dummy => a_rmH2Ocl
         CALL f1sbulk%NewField("rmH2Ocl","Deposition of H2O with cloud droplets","kg/m^2/s","time",.TRUE.,dummy)
         dummy => a_rmH2Opr
         CALL f1sbulk%NewField("rmH2Opr","Deposition of H2O with rain","kg/m^2/s","time",.TRUE.,dummy)

         ! ** BULK PROFILE OUTPUT FOR SALSA

         ! Bin dimensions, number concentrations and radius
         dummy => a_aea
         CALL f2sbulk%NewField("aea","Aerosol size bins, regime a","m","aea",.TRUE.,dummy)
         dummy => a_aeb
         CALL f2sbulk%NewField("aeb","Aerosol size bins, regime b","m","aeb",salsa_b_bins,dummy)
         dummy => a_cla
         CALL f2sbulk%NewField("cla","Cloud droplet size bins, regime a","m","cla",.TRUE.,dummy)
         dummy => a_clb
         CALL f2sbulk%NewField("clb","Cloud droplet size bins, regime b","m","clb",salsa_b_bins,dummy)
         dummy => a_prc
         CALL f2sbulk%NewField("prc","Precipitation size bins","m","prc",.TRUE.,dummy)
         dummy => a_Naa
         CALL f2sbulk%NewField("P_Naa","SALSA aerosol number concentration in regime A","kg^-1","tttt",.TRUE.,dummy)
         dummy => a_Nab
         CALL f2sbulk%NewField("P_Nab","SALSA aerosol number concentration in regime B","kg^-1","tttt",salsa_b_bins,dummy)
         dummy => a_Nca
         CALL f2sbulk%NewField("P_Nca","SALSA CDNC in regime A","kg^-1","tttt",.TRUE.,dummy)
         dummy => a_Ncb
         CALL f2sbulk%NewField("P_Ncb","SALSA CDNC in regime B","kg^-1","tttt",salsa_b_bins,dummy)
         dummy => a_Np
         CALL f2sbulk%NewField("P_Np","SALSA rdnc","kg^-1","tttt",.TRUE.,dummy)
         dummy => a_Rwaa
         CALL f2sbulk%NewField("P_Rwaa","SALSA mean aerosol wet radius, regime A","m","tttt",.TRUE.,dummy)
         dummy => a_Rwab
         CALL f2sbulk%NewField("P_Rwab","SALSA mean aerosol wet radius, regime B","m","tttt",salsa_b_bins,dummy)
         dummy => a_Rwca
         CALL f2sbulk%NewField("P_Rwca","SALSA mean cloud droplet radius, regime A","m","tttt",.TRUE.,dummy)
         dummy => a_Rwcb
         CALL f2sbulk%NewField("P_Rwcb","SALSA mean cloud droplet radius, regime B","m","tttt",salsa_b_bins,dummy)
         dummy => a_Rwp
         CALL f2sbulk%NewField("P_Rwp","SALSA mean drizzle drop radius","m","tttt",.TRUE.,dummy)

         dummy => a_cSO4a
         CALL f2sbulk%NewField("P_cSO4a","SALSA total mass mixing ratio of SO4 in aerosols","kg/kg","tttt",lso4,dummy)
         dummy => a_cSO4c
         CALL f2sbulk%NewField("P_cSO4c","SALSA total mass mixing ratio of SO4 in cloud droplets","kg/kg","tttt",lso4,dummy)
         dummy => a_cSO4p
         CALL f2sbulk%NewField("P_cSO4p","SALSA total mass mixing ratio of SO4 in drizzle drops","kg/kg","tttt",lso4,dummy)

         dummy => a_cOCa
         CALL f2sbulk%NewField("P_cOCa","SALSA total mass mixing ratio of OC in aerosols","kg/kg","tttt",loc,dummy)
         dummy => a_cOCc
         CALL f2sbulk%NewField("P_cOCc","SALSA total mass mixing ratio of OC in cloud droplets","kg/kg","tttt",loc,dummy)
         dummy => a_cOCp
         CALL f2sbulk%NewField("P_cOCp","SALSA total mass mixing ratio of OC in drizzle drops","kg/kg","tttt",loc,dummy)

         dummy => a_cBCa
         CALL f2sbulk%NewField("P_cBCa","SALSA total mass mixing ratio of BC in aerosols","kg/kg","tttt",lbc,dummy)
         dummy => a_cBCc
         CALL f2sbulk%NewField("P_cBCc","SALSA total mass mixing ratio of BC in cloud droplets","kg/kg","tttt",lbc,dummy)
         dummy => a_cBCp
         CALL f2sbulk%NewField("P_cBCp","SALSA total mass mixing ratio of BC in drizzle drops","kg/kg","tttt",lbc,dummy)

         dummy => a_cDUa
         CALL f2sbulk%NewField("P_cDUa","SALSA total mass mixing ratio of DU in aerosols","kg/kg","tttt",ldu,dummy)
         dummy => a_cDUc
         CALL f2sbulk%NewField("P_cDUc","SALSA total mass mixing ratio of DU in cloud droplets","kg/kg","tttt",ldu,dummy)
         dummy => a_cDUp
         CALL f2sbulk%NewField("P_cDUp","SALSA total mass mixing ratio of DU in drizzle drops","kg/kg","tttt",ldu,dummy)

         dummy => a_cSSa
         CALL f2sbulk%NewField("P_cSSa","SALSA total mass mixing ratio of SS in aerosols","kg/kg","tttt",lss,dummy)
         dummy => a_cSSc
         CALL f2sbulk%NewField("P_cSSc","SALSA total mass mixing ratio of SS in cloud droplets","kg/kg","tttt",lss,dummy)
         dummy => a_cSSp
         CALL f2sbulk%NewField("P_cSSp","SALSA total mass mixing ratio of SS in drizzle drops","kg/kg","tttt",lss,dummy)

         dummy => a_cNOa
         CALL f2sbulk%NewField("P_cNOa","SALSA total mass mixing ratio of NO3 in aerosols","kg/kg","tttt",lno3,dummy)
         dummy => a_cNOc
         CALL f2sbulk%NewField("P_cNOc","SALSA total mass mixing ratio of NO3 in cloud droplets","kg/kg","tttt",lno3,dummy)
         dummy => a_cNOp
         CALL f2sbulk%NewField("P_cNOp","SALSA total mass mixing ratio of NO3 in drizzle drops","kg/kg","tttt",lno3,dummy)

         dummy => a_cNHa
         CALL f2sbulk%NewField("P_cNHa","SALSA total mass mixing ratio of NH3 in aerosols","kg/kg","tttt",lnh3,dummy)
         dummy => a_cNHc
         CALL f2sbulk%NewField("P_cNHc","SALSA total mass mixing ratio of NH3 in cloud droplets","kg/kg","tttt",lnh3,dummy)
         dummy => a_cNHp
         CALL f2sbulk%NewField("P_cNHp","SALSA total mass mixing ratio of NH3 in drizzle drops","kg/kg","tttt",lnh3,dummy)

         dummy => a_cH2Oa
         CALL f2sbulk%NewField("P_cH2Oa","SALSA total mass mixing ratio of H2O in aerosols","kg/kg","tttt",.TRUE.,dummy)
         dummy => a_cH2Oc
         CALL f2sbulk%NewField("P_cH2Oc","SALSA total mass mixing ratio of H2O in cloud droplets","kg/kg","tttt",.TRUE.,dummy)
         dummy => a_cH2Op
         CALL f2sbulk%NewField("P_cH2Op","SALSA total mass mixing ratio of H2O in drizzle drops","kg/kg","tttt",.TRUE.,dummy)


         dummy => a_prl
         CALL f2sbulk%NewField("P_rl","Level 4 cloud water mixing ratio","kg/kg","tttt",.TRUE.,dummy)
         dummy => a_prr
         CALL f2sbulk%NewField("P_rr","Level 4 precipitation mixing ratio","kg/kg","tttt",.TRUE.,dummy)
         dummy => a_prv
         CALL f2sbulk%NewField("P_rv","Level 4 water vapor mixing ratio","kg/kg","tttt",.TRUE.,dummy)
         dummy => a_pRH
         CALL f2sbulk%NewField("P_RH","Level 4 relative humidity","%","tttt",.TRUE.,dummy)

         ! Stats for cloudy columns
         dummy => a_Na_c
         CALL f2sbulk%NewField("P_Na_c","Aerosol number concentration in cloudy columns","kg^-1","tttt",cloudy_col_stats,dummy)
         dummy => a_Nc_c
         CALL f2sbulk%NewField("P_Nc_c","Cloud droplet number concentration in cloudy columns","kg^-1","tttt",&
                               cloudy_col_stats,dummy)
         dummy => a_Np_c
         CALL f2sbulk%NewField("P_Np_c","Rain drop number concentration in cloudy columns","kg^-1","tttt",&
                               cloudy_col_stats,dummy)
         dummy => a_pcfrac
         CALL f2sbulk%NewField("P_cfrac","Fraction of cloudy columns","","tttt",cloudy_col_stats,dummy)
         dummy => a_clw_c
         CALL f2sbulk%NewField("P_clw_c","Cloud liquid water in cloudy columns","kg/kg","tttt",cloudy_col_stats,dummy)
         dummy => a_thl_c
         CALL f2sbulk%NewField("P_thl_c","Liquid water potential temperature in cloudy columns","K","tttt",&
                               cloudy_col_stats,dummy)


         ! ** BINNED PROFILE OUTPUT FOR SALSA **

         ! Aerosols
         ! Regime A
         dummy2 => a_Naba
         CALL f2aea%NewField("P_Naba","Aerosol number concentration in size bins A","kg^-1","ttztaea",lbinanl,dummy2)
         dummy2 => a_SO4aa
         CALL f2aea%NewField("P_SO4aa","Mass mixing ratio of SO4 in aerosol bins A","kg/kg","ttztaea",&
                              (lso4 .AND. lbinanl),dummy2)
         dummy2 => a_OCaa
         CALL f2aea%NewField("P_OCaa","Mass mixing ratio of OC in aerosol bins A","kg/kg","ttztaea",&
                              (lbinanl .AND. loc),dummy2)
         dummy2 => a_BCaa
         CALL f2aea%NewField("P_BCaa","Mass mixing ratio of BC in aerosol bins A","kg/kg","ttztaea",&
                              (lbinanl .AND. lbc),dummy2)
         dummy2 => a_DUaa
         CALL f2aea%NewField("P_DUaa","Mass mixing ratio of DU in aerosol bins A","kg/kg","ttztaea",&
                              (lbinanl .AND. ldu),dummy2)
         dummy2 => a_SSaa
         CALL f2aea%NewField("P_SSaa","Mass mixing ratio of SS in aerosol bins A","kg/kg","ttztaea",&
                              (lbinanl .AND. lss),dummy2)
         dummy2 => a_NHaa
         CALL f2aea%NewField("P_NHaa","Mass mixing ratio of NH3 in aerosol bins A","kg/kg","ttztaea",&
                              (lbinanl .AND. lnh3),dummy2)
         dummy2 => a_NOaa
         CALL f2aea%NewField("P_NOaa","Mass mixing ratio of NO3 in aerosol bins A","kg/kg","ttztaea",&
                              (lbinanl .AND. lno3),dummy2)

         ! Regime B
         dummy2 => a_Nabb
         CALL f2aeb%NewField("P_Nabb","Aerosol number concentration in size bins B","kg^-1","ttztaeb",&
                             (lbinanl .AND. salsa_b_bins),dummy2)
         dummy2 => a_SO4ab
         CALL f2aeb%NewField("P_SO4ab","Mass mixing ratio of SO4 in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. lso4),dummy2)
         dummy2 => a_OCab
         CALL f2aeb%NewField("P_OCab","Mass mixing ratio of OC in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. loc),dummy2)
         dummy2 => a_BCab
         CALL f2aeb%NewField("P_BCab","Mass mixing ratio of BC in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. lbc),dummy2)
         dummy2 => a_DUab
         CALL f2aeb%NewField("P_DUab","Mass mixing ratio of DU in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. ldu),dummy2)
         dummy2 => a_SSab
         CALL f2aeb%NewField("P_SSab","Mass mixing ratio of SS in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. lss),dummy2)
         dummy2 => a_NHab
         CALL f2aeb%NewField("P_NHab","Mass mixing ratio of NH3 in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. lnh3),dummy2)
         dummy2 => a_NOab
         CALL f2aeb%NewField("P_NOab","Mass mixing ratio of NO3 in aerosol bins B","kg/kg","ttztaeb",&
                              (lbinanl .AND. salsa_b_bins .AND. lno3),dummy2)

         ! Clouds
         ! Regime A
         dummy2 => a_Ncba
         CALL f2cla%NewField("P_Ncba","Cloud droplet number concentration in size bins A","kg^-1","ttztcla",lbinanl,dummy2)
         dummy2 => a_SO4ca
         CALL f2cla%NewField("P_SO4ca","Mass mixing ratio of SO4 in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. lso4),dummy2)
         dummy2 => a_OCca
         CALL f2cla%NewField("P_OCca","Mass mixing ratio of OC in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. loc),dummy2)
         dummy2 => a_BCca
         CALL f2cla%NewField("P_BCca","Mass mixing ratio of BC in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. lbc),dummy2)
         dummy2 => a_DUca
         CALL f2cla%NewField("P_DUca","Mass mixing ratio of DU in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. ldu),dummy2)
         dummy2 => a_SSca
         CALL f2cla%NewField("P_SSca","Mass mixing ratio of SS in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. lss),dummy2)
         dummy2 => a_NHca
         CALL f2cla%NewField("P_NHca","Mass mixing ratio of NH3 in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. lnh3),dummy2)
         dummy2 => a_NOca
         CALL f2cla%NewField("P_NOca","Mass mixing ratio of NO3 in cloud bins A","kg/kg","ttztcla",&
                              (lbinanl .AND. lno3),dummy2)

         ! Regime B
         dummy2 => a_Ncbb
         CALL f2clb%NewField("P_Ncbb","Cloud droplet number concentration in size bins B","kg^-1","ttztclb",&
                             (lbinanl .AND. salsa_b_bins),dummy2)
         dummy2 => a_SO4cb
         CALL f2clb%NewField("P_SO4cb","Mass mixing ratio of SO4 in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. lso4),dummy2)
         dummy2 => a_OCcb
         CALL f2clb%NewField("P_OCcb","Mass mixing ratio of OC in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. loc),dummy2)
         dummy2 =>  a_BCcb
         CALL f2clb%NewField("P_BCcb","Mass mixing ratio of BC in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. lbc),dummy2)
         dummy2 => a_DUcb
         CALL f2clb%NewField("P_DUcb","Mass mixing ratio of DU in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. ldu),dummy2)
         dummy2 => a_SScb
         CALL f2clb%NewField("P_SScb","Mass mixing ratio of SS in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. lss),dummy2)
         dummy2 => a_NHcb
         CALL f2clb%NewField("P_NHcb","Mass mixing ratio of NH3 in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. lnh3),dummy2)
         dummy2 => a_NOcb
         CALL f2clb%NewField("P_NOcb","Mass mixing ratio of NO3 in cloud bins B","kg/kg","ttztclb",&
                              (lbinanl .AND. salsa_b_bins .AND. lno3),dummy2)

         ! Precipitation
         dummy2 => a_Npb
         CALL f2prc%NewField("P_Npb","Number concentration of drizzle bins","kg^-1","ttztprc",lbinanl,dummy2)
         dummy2 => a_SO4pb
         CALL f2prc%NewField("P_SO4pb","Mass mixing ratio of SO4 in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. lso4),dummy2)
         dummy2 => a_OCpb
         CALL f2prc%NewField("P_OCpb","Mass mixing ratio of OC in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. loc),dummy2)
         dummy2 => a_BCpb
         CALL f2prc%NewField("P_BCpb","Mass mixing ratio of BC in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. lbc),dummy2)
         dummy2 => a_DUpb
         CALL f2prc%NewField("P_DUpb","Mass mixing ratio of DU in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. ldu),dummy2)
         dummy2 => a_SSpb
         CALL f2prc%NewField("P_SSpb","Mass mixing ratio of SS in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. lss),dummy2)
         dummy2 => a_NHpb
         CALL f2prc%NewField("P_NHpb","Mass mixing ratio of NH3 in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. lnh3),dummy2)
         dummy2 => a_NOpb
         CALL f2prc%NewField("P_NOpb","Mass mixing ratio of NO3 in drizzle bins","kg/kg","ttztprc",&
                              (lbinanl .AND. lno3),dummy2)

         IF (level >= 5) THEN

         dummy => a_Ni_ic
         CALL f1lvl5%NewField("Ni_ic","Ice number concentration in liquid clouds","kg^-1","time",.TRUE.,dummy)

         dummy => a_Ni_ii
         CALL f1lvl5%NewField("Ni_ii","Ice number concentration in icy regions","kg^-1","time",.TRUE.,dummy)

         dummy => a_Ni_is
         CALL f1lvl5%NewField("Ni_is","Ice number concentration in snowy regions","kg^-1","time",.TRUE.,dummy)

         dummy => a_Ns_ic
         CALL f1lvl5%NewField("Ns_ic","Snow number concentration in liquid clouds","kg^-1","time",.TRUE.,dummy)

         dummy => a_Ns_ii
         CALL f1lvl5%NewField("Ns_ii","Snow number concentration in icy regions","kg^-1","time",.TRUE.,dummy)

         dummy => a_Ns_is
         CALL f1lvl5%NewField("Ns_is","Snow number concentration in snowy regions","kg^-1","time",.TRUE.,dummy)

         dummy => a_Ri_ii
         CALL f1lvl5%NewField("Ri_ii","Mean ice radius in icy regions","m","time",.TRUE.,dummy)

         dummy => a_iwp_bar
         CALL f1lvl5%NewField("iwp_bar","Ice-water path","kg/m^2","time",.TRUE.,dummy)

         dummy => a_imax
         CALL f1lvl5%NewField("imax","Maximum ice water mixing ratio","kg/kg","time",.TRUE.,dummy)

         dummy => a_nicnt
         CALL f1lvl5%NewField("nicnt","Ice cell counts","#","time",.TRUE.,dummy)

         dummy => a_Rs_is
         CALL f1lvl5%NewField("Rs_is","Mean snow radius in snowy regions","m","time",.TRUE.,dummy)

         dummy => a_swp_bar
         CALL f1lvl5%NewField("swp_bar","Snow-water path","kg/m^2","time",.TRUE.,dummy)

         dummy => a_smax
         CALL f1lvl5%NewField("smax","Maximum snow water mixing ratio","kg/kg","time",.TRUE.,dummy)

         dummy => a_nscnt
         CALL f1lvl5%NewField("nscnt","Snow cell counts","#","time",.TRUE.,dummy)

         dummy => a_rmSO4ic
         CALL f1lvl5%NewField("rmSO4ic","Deposition of SO4 with ice","kg/m^2/s","time",lso4,dummy)

         dummy => a_rmSO4sn
         CALL f1lvl5%NewField("rmSO4sn","Deposition of SO4 with snow","kg/m^2/s","time",lso4,dummy)

         dummy => a_rmOCic
         CALL f1lvl5%NewField("rmOCic","Deposition of OC with ice","kg/m^2/s","time",loc,dummy)

         dummy => a_rmOCsn
         CALL f1lvl5%NewField("rmOCsn","Deposition of OC with snow","kg/m^2/s","time",loc,dummy)

         dummy => a_rmBCic
         CALL f1lvl5%NewField("rmBCic","Deposition of BC with ice","kg/m^2/s","time",lbc,dummy)

         dummy => a_rmBCsn
         CALL f1lvl5%NewField("rmBCsn","Deposition of BC with snow","kg/m^2/s","time",lbc,dummy)

         dummy => a_rmDUic
         CALL f1lvl5%NewField("rmDUic","Deposition of DU with ice","kg/m^2/s","time",ldu,dummy)

         dummy => a_rmDUsn
         CALL f1lvl5%NewField("rmDUsn","Deposition of DU with snow","kg/m^2/s","time",ldu,dummy)

         dummy => a_rmNOic
         CALL f1lvl5%NewField("rmNOic","Deposition of NO3 with ice","kg/m^2/s","time",lno3,dummy)

         dummy => a_rmNOsn
         CALL f1lvl5%NewField("rmNOsn","Deposition of NO3 with snow","kg/m^2/s","time",lno3,dummy)

         dummy => a_rmNHic
         CALL f1lvl5%NewField("rmNHic","Deposition of NH3 with ice","kg/m^2/s","time",lnh3,dummy)

         dummy => a_rmNHsn
         CALL f1lvl5%NewField("rmNHsn","Deposition of NH3 with snow","kg/m^2/s","time",lnh3,dummy)

         dummy => a_rmSSic
         CALL f1lvl5%NewField("rmSSic","Deposition of SS with ice","kg/m^2/s","time",lss,dummy)

         dummy => a_rmSSsn
         CALL f1lvl5%NewField("rmSSsn","Deposition of SS with snow","kg/m^2/s","time",lss,dummy)

         dummy => a_rmH2Oic
         CALL f1lvl5%NewField("rmH2Oic","Deposition of water with ice","kg/m^2/s","time",.TRUE.,dummy)

         dummy => a_rmH2Osn
         CALL f1lvl5%NewField("rmH2Osn","Deposition of water with snow","kg/m^2/s","time",.TRUE.,dummy)

         dummy => a_sfrac
         CALL f1lvl5%NewField("sfrac","Surface snow precipitation fraction","-","time",.TRUE.,dummy)

         dummy => a_sprcp
         CALL f1lvl5%NewField("sprcp","Surface snow precipitation rate","W/m^2","time",.TRUE.,dummy)

         dummy => a_ica
         CALL f2lvl5%NewField("ica","Ice size bins, regime a","m","ica",.TRUE.,dummy)

         dummy => a_icb
         CALL f2lvl5%NewField("icb","Ice size bins, regime b","m","icb",salsa_b_bins,dummy)

         dummy => a_snw
         CALL f2lvl5%NewField("snw","Snow size bins","m","snw",.TRUE.,dummy)

         dummy => a_Nia
         CALL f2lvl5%NewField("P_Nia","SALSA ice number concentration in regime A","kg^-1","tttt",.TRUE.,dummy)

         dummy => a_Nib
         CALL f2lvl5%NewField("P_Nib","SALSA ice number concentration in regime B","kg^-1","tttt",salsa_b_bins,dummy)

         dummy => a_Ns
         CALL f2lvl5%NewField("P_Ns","SALSA snow number concentration","m","tttt",.TRUE.,dummy)

         dummy => a_Rwia
         CALL f2lvl5%NewField("P_Rwia","SALSA mean ice radius, regime a","m","tttt",.TRUE.,dummy)

         dummy => a_Rwib
         CALL f2lvl5%NewField("P_Rwib","SALSA mean ice radius, regime b","m","tttt",salsa_b_bins,dummy)

         dummy => a_Rws
         CALL f2lvl5%NewField("P_Rws","SALSA mean snow radius","m","tttt",.TRUE.,dummy)

         dummy => a_cSO4i
         CALL f2lvl5%NewField("P_cSO4i","SALSA total mass mixing ratio of SO4 in ice","kg/kg","tttt",lso4,dummy)

         dummy => a_cSO4s
         CALL f2lvl5%NewField("P_cSO4s","SALSA total mass mixing ratio of SO4 in snow","kg/kg","tttt",lso4,dummy)

         dummy => a_cOCi
         CALL f2lvl5%NewField("P_cOCi","SALSA total mass mixing ratio of OC in ice","kg/kg","tttt",loc,dummy)

         dummy => a_cOCs
         CALL f2lvl5%NewField("P_cOCs","SALSA total mass mixing ratio of OC in snow","kg/kg","tttt",loc,dummy)

         dummy => a_cBCi
         CALL f2lvl5%NewField("P_cBCi","SALSA total mass mixing ratio of BC in ice","kg/kg","tttt",lbc,dummy)

         dummy => a_cBCs
         CALL f2lvl5%NewField("P_cBCs","SALSA total mass mixing ratio of BC in snow","kg/kg","tttt",lbc,dummy)

         dummy => a_cDUi
         CALL f2lvl5%NewField("P_cDUi","SALSA total mass mixing ratio of DU in ice","kg/kg","tttt",ldu,dummy)

         dummy => a_cDUs
         CALL f2lvl5%NewField("P_cDUs","SALSA total mass mixing ratio of DU in snow","kg/kg","tttt",ldu,dummy)

         dummy => a_cSSi
         CALL f2lvl5%NewField("P_cSSi","SALSA total mass mixing ratio of SS in ice","kg/kg","tttt",lss,dummy)

         dummy => a_cSSc
         CALL f2lvl5%NewField("P_cSSs","SALSA total mass mixing ratio of SS in snow","kg/kg","tttt",lss,dummy)

         dummy => a_cNOi
         CALL f2lvl5%NewField("P_cNOi","SALSA total mass mixing ratio of NO3 in ice","kg/kg","tttt",lno3,dummy)

         dummy => a_cNOs
         CALL f2lvl5%NewField("P_cNOs","SALSA total mass mixing ratio of NO3 in snow","kg/kg","tttt",lno3,dummy)

         dummy => a_cNHi
         CALL f2lvl5%NewField("P_cNHi","SALSA total mass mixing ratio of NH3 in ice","kg/kg","tttt",lnh3,dummy)

         dummy => a_cNHs
         CALL f2lvl5%NewField("P_cNHs","SALSA total mass mixing ratio of NH3 in snow","kg/kg","tttt",lnh3,dummy)

         dummy => a_cH2Oi
         CALL f2lvl5%NewField("P_cH2Oi","SALSA total mass mixing ratio of H2O in ice","kg/kg","tttt",.TRUE.,dummy)

         dummy => a_cH2Os
         CALL f2lvl5%NewField("P_cH2Os","SALSA total mass mixing ratio of H2O in snow","kg/kg","tttt",.TRUE.,dummy)

         dummy => a_pri
         CALL f2lvl5%NewField("P_ri","Ice water mixing ratio","kg/kg","tttt",.TRUE.,dummy)

         dummy => a_prs
         CALL f2lvl5%NewField("P_rs","Snow water mixing ratio","kg/kg","tttt",.TRUE.,dummy)

         dummy => a_prhi
         CALL f2lvl5%NewField("P_RHi","Relative humidity over ice","%","tttt",.TRUE.,dummy)

         dummy => a_srate
         CALL f2lvl5%NewField("srate","Snow deposition flux (positive downward)","W/m^2","ttmt",.TRUE.,dummy)


         dummy2 => a_Niba
         CALL f2ica%NewField("P_Niba","Ice number concentration in size bins A","kg^-1","ttztica",lbinanl,dummy2)

         dummy2 => a_SO4ia
         CALL f2ica%NewField("P_SO4ia","Mass mixing ratio of SO4 in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. lso4),dummy2)

         dummy2 => a_OCia
         CALL f2ica%NewField("P_OCia","Mass mixing ratio of OC in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. loc),dummy2)

         dummy2 => a_BCia
         CALL f2ica%NewField("P_BCia","Mass mixing ratio of BC in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. lbc),dummy2)

         dummy2 => a_DUia
         CALL f2ica%NewField("P_DUia","Mass mixing ratio of DU in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. ldu),dummy2)

         dummy2 => a_SSia
         CALL f2ica%NewField("P_SSia","Mass mixing ratio of SS in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. lss),dummy2)

         dummy2 => a_NHia
         CALL f2ica%NewField("P_NOia","Mass mixing ratio of NO3 in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. lno3),dummy2)

         dummy2 => a_NHia
         CALL f2ica%NewField("P_NHia","Mass mixing ratio of NH3 in ice bins A","kg/kg",&
         "ttztica",(lbinanl .AND. lnh3),dummy2)

         dummy2 => a_Nibb
         CALL f2icb%NewField("P_Nibb","Ice number concentration in size bins B","kg^-1",&
         "ttzticb",(lbinanl .AND. salsa_b_bins),dummy2)

         dummy2 => a_SO4ib
         CALL f2icb%NewField("P_SO4ib","Mass mixing ratio of SO4 in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. lso4 .AND. salsa_b_bins),dummy2)

         dummy2 => a_OCib
         CALL f2icb%NewField("P_OCib","Mass mixing ratio of OC in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. loc .AND. salsa_b_bins),dummy2)

         dummy2 => a_BCib
         CALL f2icb%NewField("P_BCib","Mass mixing ratio of BC in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. lbc .AND. salsa_b_bins),dummy2)

         dummy2 => a_DUib
         CALL f2icb%NewField("P_DUib","Mass mixing ratio of DU in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. ldu .AND. salsa_b_bins),dummy2)

         dummy2 => a_SSib
         CALL f2icb%NewField("P_SSib","Mass mixing ratio of SS in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. lss .AND. salsa_b_bins),dummy2)

         dummy2 => a_NOib
         CALL f2icb%NewField("P_NOib","Mass mixing ratio of NO3 in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. lno3 .AND. salsa_b_bins),dummy2)

         dummy2 => a_NHib
         CALL f2icb%NewField("P_NHib","Mass mixing ratio of NH3 in ice bins B","kg/kg",&
         "ttzticb",(lbinanl .AND. lnh3 .AND. salsa_b_bins),dummy2)

         dummy2 => a_Nsb
         CALL f2snw%NewField("P_Nsb","Number concentration of snow","kg^-1","ttztsnw",lbinanl,dummy2)

         dummy2 => a_SO4sb
         CALL f2snw%NewField("P_SO4sb","Mass mixing ratio of SO4 in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. lso4),dummy2)

         dummy2 => a_OCsb
         CALL f2snw%NewField("P_OCsb","Mass mixing ratio of OC in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. loc),dummy2)

         dummy2 => a_BCsb
         CALL f2snw%NewField("P_BCsb","Mass mixing ratio of BC in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. lbc),dummy2)

         dummy2 => a_DUsb
         CALL f2snw%NewField("P_DUsb","Mass mixing ratio of DU in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. ldu),dummy2)

         dummy2 => a_SSsb
         CALL f2snw%NewField("P_SSsb","Mass mixing ratio of SS in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. lss),dummy2)

         dummy2 => a_NOsb
         CALL f2snw%NewField("P_NOsb","Mass mixing ratio of NO3 in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. lno3),dummy2)

         dummy2 => a_NHsb
         CALL f2snw%NewField("P_NHsb","Mass mixing ratio of NH3 in snow bins","kg/kg",&
         "ttztsnw",(lbinanl .AND. lnh3),dummy2)


         END IF

         END IF

   END SUBROUTINE

   !
   ! -------------------------------------------------------------------------
   !
   INTEGER FUNCTION close_stat()

      USE netcdf

      close_stat = nf90_close(ncid1) + nf90_close(ncid2)
      IF (csflg) close_stat = close_stat + nf90_close(ncid3)

   END FUNCTION close_stat
   !
   ! -------------------------------------------------------------------------
   !
   REAL FUNCTION get_zi (n1, n2, n3, itype, sx, xx, z, threshold)

      INTEGER, INTENT (in) :: n1, n2, n3, itype
      REAL, INTENT (in)    :: xx(n1), z(n1), sx(n1,n2,n3), threshold

      INTEGER :: i, j, k, kk
      REAL    :: zibar, sval, dmy, scr(n2,n3)

      get_zi = -999.
      SELECT CASE(itype)
         CASE (1)
            !
            ! find level at which sx=threshold (xx is one over grid spacing)
            !
            zibar = 0.
            DO j = 3, n3-2
               DO i = 3, n2-2
                  k = 2
                  DO WHILE (k < n1-2 .AND. sx(k,i,j) > threshold)
                     k = k+1
                  END DO
                  IF (k == n1-2) zibar = -999.
                  IF (zibar /= -999.) zibar = zibar + z(k-1) +  &
                     (threshold - sx(k-1,i,j))/xx(k-1)     /  &
                     (sx(k,i,j) - sx(k-1,i,j) + epsilon(1.))
               END DO
            END DO
            IF (zibar /= -999.) get_zi = zibar/REAL((n3-4)*(n2-4))

         CASE(2)
            !
            ! find level of maximum gradient (xx is one over grid spacing)
            !
            scr = 0.
            DO j = 3, n3-2
               DO i = 3, n2-2
                  sval = 0.
                  DO k = 2, n1-5
                     dmy = (sx(k+1,i,j)-sx(k,i,j))*xx(k)
                     IF (dmy > sval) THEN
                        sval = dmy
                        scr(i,j) = z(k)
                     END IF
                  END DO
               END DO
            END DO
            get_zi = get_avg2dh(n2,n3,scr)

         CASE(3)
            !
            ! find level where xx is a maximum
            !
            sval = -huge(1.)
            kk = 1
            DO k = 2, n1
               IF (xx(k) > sval) THEN
                  kk = k
                  sval = xx(k)
               END IF
            END DO
            get_zi = z(kk)

         CASE(4)
            !
            ! find level where xx is a minimum
            !
            sval = huge(1.)
            kk = 1
            DO k = 2, n1-2
               IF (xx(k) < sval) THEN
                  kk = k
                  sval = xx(k)
               END IF
            END DO
            get_zi = z(kk)
      END SELECT

   END FUNCTION get_zi

END MODULE stat

