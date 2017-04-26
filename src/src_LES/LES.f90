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
PROGRAM ucla_les

  USE mpi_interface, ONLY : myid

  IMPLICIT NONE

  REAL :: t1, t2

  CALL cpu_time(t1)
  CALL driver
  CALL cpu_time(t2)

  IF (myid == 0) THEN
    PRINT "(/,' ',49('-')/,' ',A16,F10.1,' s')", '  Execution time: ', t2-t1
    STOP ' ..... Normal termination'
  END IF

CONTAINS

  !----------------------------------------------------------------------
  ! Subroutine Driver:  This is the main program driver.  It calls routines
  ! to read the model initialization file, and configure memory and pointes.
  ! It also calls the routines which initialize the model and timestep it.
  !
  SUBROUTINE driver

    USE grid, ONLY          : define_grid, define_vars, level, nxp, nyp, nzp, nxpart
    USE init, ONLY          : initialize
    USE step, ONLY          : stepper
    USE mpi_interface, ONLY : init_mpi, define_decomp,                    &
         init_alltoall_reorder, appl_finalize

    ! Added for SALSA
    USE mo_salsa_init, ONLY : define_salsa, salsa_initialize

    IMPLICIT NONE

    INTEGER ierror

    CALL init_mpi

    CALL define_parm

    IF (level >= 4) CALL define_salsa ! Read SALSA namelist etc.

    IF (level >= 4) CALL salsa_initialize ! All salsa variables are now initialized

    CALL define_decomp(nxp, nyp, nxpart)

    CALL define_grid

    CALL init_alltoall_reorder(nxp, nyp, nzp)

    CALL define_vars

    CALL initialize ! Added initialization of aerosol size distributions here + a single call
                    ! for SALSA to set up cloud microphysics
    CALL stepper

    CALL appl_finalize(ierror)

    RETURN
  END SUBROUTINE driver

  !
  ! ----------------------------------------------------------------------
  ! Subroutine Read_nl: Driver for reading model namelist
  !
  SUBROUTINE define_parm

    USE util, ONLY : fftinix,fftiniy
    USE sgsm, ONLY : csx, prndtl
    USE srfc, ONLY : isfctyp, zrough, ubmin, dthcon, drtcon
    USE step, ONLY : timmax, istpfl, corflg, outflg, frqanl, frqhis,          &
         strtim, radfrq, cntlat, nudge_time, nudge_zmin, nudge_zmax, &
         nudge_theta, tau_theta, nudge_rv, tau_rv, nudge_u, tau_u, &
         nudge_v, tau_v, nudge_ccn, tau_ccn
    USE grid, ONLY : deltaz, deltay, deltax, nzp, nyp, nxp, nxpart, &
         dtlong, dzrat,dzmax, th00, umean, vmean, isgstyp, naddsc, level,     &
         filprf, expnme, iradtyp, igrdtyp, nfpt, distim, runtype, CCN,        &
         Tspinup,sst, lbinanl
    USE init, ONLY : us, vs, ts, rts, ps, hs, ipsflg, itsflg,iseed, hfilin,   &
         zrand
    USE stat, ONLY : ssam_intvl, savg_intvl, mcflg, csflg
    USE forc, ONLY : radsounding, &        ! Juha: added for radiation background profile
                     div, case_name, &     ! Divergence, forcing case name
                     sfc_albedo, &         ! Surface albedo
                     useMcICA,RadConstPress,RadPrecipBins
    USE mcrp, ONLY : sed_aero, sed_cloud, sed_precp, sed_ice, sed_snow
    USE mpi_interface, ONLY : myid, appl_abort, ver, author

    IMPLICIT NONE

    NAMELIST /model/  &
         expnme    ,       & ! experiment name
         nxpart    ,       & ! whether partition in x direction?
         naddsc    ,       & ! Number of additional scalars
         savg_intvl,       & ! output statistics frequency
         ssam_intvl,       & ! integral accumulate/ts PRINT frequency
         mcflg,            & ! Mass conservation stats flag
         csflg,            & ! Column statistics flag
         corflg , cntlat , & ! coriolis flag
         nfpt   , distim , & ! rayleigh friction points, dissipation time
         level  , CCN    , & ! Microphysical model Number of CCN per kg of air
         iseed  , zrand  , & ! random seed
         nxp    , nyp    , nzp   ,  & ! number of x, y, z points
         deltax , deltay , deltaz , & ! delta x, y, z (meters)
         dzrat  , dzmax  , igrdtyp, & ! stretched grid parameters
         timmax , dtlong , istpfl , & ! timestep control
         runtype, hfilin , filprf , & ! type of run (INITIAL or HISTORY)
         frqhis , frqanl , outflg , & ! freq of history/anal writes, output flg
         iradtyp, radfrq , strtim , & ! radiation type flag
         isfctyp, ubmin  , zrough , & ! surface parameterization type
         sst    , dthcon , drtcon , & ! SSTs, surface flx parameters
         isgstyp, csx    , prndtl , & ! SGS model type, parameters
         ipsflg , itsflg ,          & ! sounding flags
         hs     , ps     , ts    ,  & ! sounding heights, pressure, temperature
         us     , vs     , rts   ,  & ! sounding E/W winds, water vapor
         umean  , vmean  , th00,    & ! gallilean E/W wind, basic state
         Tspinup, lbinanl,          & ! Length of spinup period in seconds
         nudge_time,                & ! Total nudging time (independent of spin-up)
         nudge_zmin, nudge_zmax, & ! Altitude (m) range for nudging
         nudge_theta, tau_theta,   & ! Temperature nudging
         nudge_rv, tau_rv,   & ! Water vapor mixing ratio nudging
         nudge_u, tau_u, nudge_v, tau_v,  & ! Horozontal wind nudging
         nudge_ccn, tau_ccn,   & ! Aerosol number concentration nudging
         radsounding, div, case_name, & ! Name of the radiation sounding file, divergence for LEVEL 4
         sfc_albedo,                  & ! Surface albedo
         useMcICA,           & ! use the Monte Carlo Independent Column Approximation method (T/F)
         RadConstPress,      & ! keep constant pressure levels (T/F),
         RadPrecipBins,      & ! add precipitation bins cloud water (0, 1, 2, 3,...)
         sed_aero, sed_cloud, sed_precp, sed_ice, sed_snow ! Sedimentation (T/F)

    NAMELIST /version/  &
         ver, author        ! Information about UCLALES-SALSA version and author

    ps       = 0.
    ts       = th00
    !
    ! these are for initializing the temp variables used in ffts in x and y
    ! directions.
    !
      fftinix=1
      fftiniy=1
    !
    ! read namelist from specified file
    !
    OPEN  (1,status='old',file='NAMELIST')
    READ  (1, nml=version)
    READ  (1, nml=model)
    CLOSE(1)

    !
    ! write file variable control to standard output
    !
    IF (myid == 0) THEN
       IF (runtype == 'HISTORY') THEN
          WRITE(*,601) expnme, hfilin, timmax
       ELSE
          WRITE(*,600) expnme, timmax
       END IF
       IF (outflg) WRITE(*,602) filprf, frqhis, frqanl, Tspinup
       !
       ! Do some cursory error checking in namelist variables
       !

       IF (min(nxp,nyp) < 5) THEN
          IF (myid == 0) PRINT *, '  ABORTING: min(nxp,nyp) must be > 4.'
          CALL appl_abort(0)
       END IF

       IF (nzp < 3 ) THEN
          IF (myid == 0) PRINT *, '  ABORTING: nzp must be > 2 '
          CALL appl_abort(0)
       END IF

       IF (cntlat < -90. .OR. cntlat > 90.) THEN
          IF (myid == 0) PRINT *, '  ABORTING: central latitude out of bounds.'
          CALL appl_abort(0)
       END IF
    END IF

600 FORMAT(//' ',49('-')/,' ',/,'  Initial Experiment: ',A50 &
         /,'  Final Time:         ',F8.1,' s'              )
601 FORMAT(//' ',49('-')/,' ',/,'  Restart Experiment: ',A50 &
         /,'  Restart File: ',A30,                           &
         /,'  Final Time: ',F10.1,' s'              )
602 FORMAT('  Output File Stem:   ',A50                      &
         /,'  History Frequency:  ',F7.1,                    &
         /,'  Analysis Frequency: ',F7.1,                    &
         /,'  Model spinup period: ',F7.1)

    RETURN
  END SUBROUTINE define_parm

END PROGRAM ucla_les

