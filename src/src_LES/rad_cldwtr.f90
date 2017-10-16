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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE cldwtr

  USE defs, ONLY : nv, mb
  IMPLICIT NONE
  INTEGER, SAVE :: nsizes
  LOGICAL, SAVE :: Initialized = .FALSE.
  LOGICAL, SAVE :: iceInitialized = .FALSE.
  LOGICAL, SAVE :: grpInitialized = .FALSE.
  INTEGER, SAVE :: mbs,mbir

  REAL, ALLOCATABLE :: re(:), fl(:), bz(:,:), wz(:,:), gz(:,:)
  REAL, ALLOCATABLE :: ap(:,:), bp(:,:), cps(:,:,:), dps(:,:), cpir(:,:)
  REAL, ALLOCATABLE :: bg(:), wgf(:), gg(:)
  REAL :: gwc

CONTAINS
  !
  !---------------------------------------------------------------------------
  ! Subroutine cloud_init initialize data arrays for the cloud model,
  ! checking for consistency between band structure of cloud model and CKD
  !
  SUBROUTINE init_cldwtr

    USE ckd, ONLY : band, center
    INTEGER, PARAMETER  :: nrec = 21600

    REAL, DIMENSION(mb) :: cntrs

    INTEGER             :: ib, i, nbands
    CHARACTER (len=12)  :: frmt

    OPEN( unit = 71, file = 'datafiles/cldwtr.dat', status = 'old', recl=nrec)
    READ(71,'(2I3)') nsizes, nbands
    IF (nbands /= mb .OR. nsizes*nbands*15 > nrec) &
         STOP 'TERMINATING: incompatible cldwtr.dat file'

    ALLOCATE (re(nsizes),fl(nsizes),bz(nsizes,mb),wz(nsizes,mb),gz(nsizes,mb))
    WRITE(frmt,'(A1,I2.2,A8)') '(',mb,'E15.7)    '
    READ(71,frmt) (cntrs(i), i=1,mb)
    DO i = 1, mb
       IF (spacing(1.) < abs(cntrs(i)- center(band(i))) ) &
            STOP 'TERMINATING: cloud properties not matched to band structure'
    END DO

    WRITE(frmt,'(A1,I2.2,A9)') '(',nsizes,'E15.7)   '
    READ(71,frmt) (re(i), i=1,nsizes)
    READ(71,frmt) (fl(i), i=1,nsizes)

    WRITE(frmt,'(A1,I4.4,A7)') '(',nsizes*mb,'E15.7) '
    READ(71,frmt) ((bz(i,ib), i=1,nsizes), ib=1,mb)
    READ(71,frmt) ((wz(i,ib), i=1,nsizes), ib=1,mb)
    READ(71,frmt) ((gz(i,ib), i=1,nsizes), ib=1,mb)
    CLOSE(71)

    IF (minval((/bz,wz,gz/)) < 0.) &
         STOP 'TERMINATING: cloud properties out of bounds'

    Initialized = .TRUE.

  END SUBROUTINE init_cldwtr

  !
  !---------------------------------------------------------------------------
  ! Subroutine cloud_init initialize data arrays for the cloud model,
  ! checking for consistency between band structure of cloud model and CKD
  !
  SUBROUTINE init_cldice

    USE ckd, ONLY : center
    INTEGER, PARAMETER :: nrec = 21600


    INTEGER            :: i, j

    mbs  = 6
    mbir = 12

    ALLOCATE (ap(3,mb),bp(4,mb),cps(4,4,mbs),dps(4,mbs),cpir(4,mbir))
    OPEN( unit = 71, file = 'datafiles/cldice.dat', status = 'old', recl=nrec)
    DO i = 1, mb
       READ(71,'(3E10.3)') ap(1,i), ap(2,i), ap(3,i)
    END DO
    DO i = 1, mb
       READ(71,'(4E12.5)') bp(1,i), bp(2,i), bp(3,i), bp(4,i)
    END DO
    DO i = 1, mbs
       DO j = 1, 4
          READ(71,'(4E12.5)') cps(1,j,i), cps(2,j,i), cps(3,j,i), cps(4,j,i)
       END DO
    END DO
    DO i = 1, mbs
       READ(71,'(4E12.5)') dps(1,i), dps(2,i), dps(3,i), dps(4,i)
    END DO
    DO i = 1, mbir
       READ(71,'(4E13.4)') cpir(1,i), cpir(2,i), cpir(3,i), cpir(4,i)
    END DO

    CLOSE(71)

    iceInitialized = .TRUE.

  END SUBROUTINE init_cldice

  !
  !---------------------------------------------------------------------------
  ! Subroutine cloud_init initialize data arrays for the cloud model,
  ! checking for consistency between band structure of cloud model and CKD
  !
  SUBROUTINE init_cldgrp

    INTEGER, PARAMETER :: nrec = 21600
    INTEGER            :: i
    CHARACTER (len=12) :: frmt

    gwc = 1.5e10

    ALLOCATE (bg(mb),wgf(mb),gg(mb))
    OPEN( unit = 71, file = 'datafiles/cldgrp.dat', status = 'old', recl=nrec)
    WRITE(frmt,'(A1,I2.2,A8)') '(',mb,'E10.3)    '
    READ(71,frmt) (bg(i), i=1,mb)
    WRITE(frmt,'(A1,I2.2,A8)') '(',mb,'F7.4)    '
    READ(71,frmt) (wgf(i), i=1,mb)
    READ(71,frmt) (gg(i), i=1,mb)

    CLOSE(71)

    grpInitialized = .TRUE.

  END SUBROUTINE init_cldgrp

  ! -----------------------------------------------------------------------
  ! Subroutine cloud_water:  calculates the optical depth (tw), single 
  ! scattering albedo (ww), and phase function (www(4)) given the cloud
  ! water [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes  
  !
  SUBROUTINE cloud_water ( ib, pre, pcw, dz, tw, ww, www )

    IMPLICIT NONE

    INTEGER, INTENT (in) :: ib
    REAL, DIMENSION (nv), INTENT (in) :: pre, pcw, dz
    REAL, INTENT (out) :: tw(nv), ww(nv), www(nv,4)

    INTEGER :: k, j, j0, j1
    REAL    :: gg, wght, cwmks

    IF (.NOT. Initialized) STOP 'TERMINATING: Cloud not Initialized'

    DO k = 1, nv
       cwmks = pcw(k)*1.e-3
       IF ( cwmks >= 1.e-8) THEN
          j = 0
          DO WHILE (j < nsizes .AND. pre(k) > re(MIN(j+1,nsizes))) ! Juha: purkkafix for (too) large pre
             j = j + 1
          END DO
          IF (j >= 1 .AND. j < nsizes) THEN
             j1 = j+1
             wght = (pre(k)-re(j))/(re(j1)-re(j))
             tw(k) = dz(k) * cwmks * ( bz(j,ib) / fl(j) +            &
                     ( bz(j1,ib) / fl(j1) - bz(j,ib) / fl(j) ) /     &
                     ( 1.0 / re(j1) - 1.0 / re(j) ) * ( 1.0 / pre(k) &
                     - 1.0 / re(j) ) )
             ww(k) = wz(j,ib) + (wz(j1,ib) - wz(j,ib) ) * wght
             gg    = gz(j,ib) + (gz(j1,ib) - gz(j,ib) ) * wght
          ELSE
             j0 = max(j,1)
             tw(k) = dz(k) * cwmks * (bz(j0,ib)/fl(j0))
             ww(k) = wz(j0,ib)
             gg    = gz(j0,ib)
          END IF
          www(k,1) = 3.0 * gg
          DO j = 2, 4
             wght = REAL(2*j+1)/REAL(2*j-1)
             www(k,j) = www(k,j-1) * gg * wght
          END DO
       ELSE
          www(k,:) = 0.0
          tw(k) = 0.0
          ww(k) = 0.0
          gg    = 0.
       END IF
       IF (ww(k) < 0.) PRINT*,'bad ww, ',ww(k),ib,k,cwmks
       IF (tw(k) < 0.) PRINT*,'bad tw, ',tw(k),ib,k,cwmks
    END DO

    RETURN
  END SUBROUTINE cloud_water

  ! -----------------------------------------------------------------------
  ! Subroutine cloud_ice:  calculates the optical depth (ti), single
  ! scattering albedo (wi), and phase function (wwi(4)) given the cloud
  ! ice [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes
  !
  SUBROUTINE cloud_ice ( ib, pde, pci, dz, ti, wi, wwi )

    IMPLICIT NONE

    INTEGER, INTENT (in) :: ib
    REAL, DIMENSION (nv), INTENT (in) :: pde, pci, dz
    REAL, INTENT (out)   :: ti(nv), wi(nv), wwi(nv,4)

    INTEGER :: ibr,k
    REAL    :: gg, cwmks
    REAL    :: fw1, fw2, fw3, wf1, wf2, wf3, wf4, x1, x2, x3, x4,  fd

    IF (.NOT. iceInitialized) STOP 'TERMINATING: Ice not Initialized'

    DO k = 1, nv
       cwmks = pci(k)!here we don't need the factor 1000
       IF ( (cwmks >= 1.e-5) .AND. (pde(k) > 0.) ) THEN
         fw1 = pde(k)
         fw2 = fw1 * pde(k)
         fw3 = fw2 * pde(k)
         ti(k) = dz(k) * cwmks * ( ap(1,ib) + &
                 ap(2,ib) / fw1 + ap(3,ib) / fw2 )
         wi(k) = 1.0 - ( bp(1,ib) + bp(2,ib) * fw1 + &
                 bp(3,ib) * fw2 + bp(4,ib) * fw3 )

         IF ( wi(k) < 0. ) PRINT*,'bad wi, ',wi(k),ib,k,bp(1,ib),bp(2,ib),bp(3,ib),bp(4,ib),fw1,fw2,fw3
         IF ( ti(k) < 0. ) PRINT*,'bad ti, ',ti(k),ib,k,cwmks,dz(k),ap(1,ib),ap(2,ib),ap(3,ib),fw1,fw2

         IF ( ib <= mbs ) THEN ! shortwave
           fd = dps(1,ib) + dps(2,ib) * fw1 + &
                dps(3,ib) * fw2 + dps(4,ib) * fw3
           wf1 = cps(1,1,ib) + cps(2,1,ib) * fw1 + &
                 cps(3,1,ib) * fw2 + cps(4,1,ib) * fw3
           wwi(k,1) = ( 1.0 - fd ) * wf1 + 3.0 * fd
           wf2 = cps(1,2,ib) + cps(2,2,ib) * fw1 + &
                 cps(3,2,ib) * fw2 + cps(4,2,ib) * fw3
           wwi(k,2) = ( 1.0 - fd ) * wf2 + 5.0 * fd
           wf3 = cps(1,3,ib) + cps(2,3,ib) * fw1 + &
                 cps(3,3,ib) * fw2 + cps(4,3,ib) * fw3
           wwi(k,3) = ( 1.0 - fd ) * wf3 + 7.0 * fd
           wf4 = cps(1,4,ib) + cps(2,4,ib) * fw1 + &
                 cps(3,4,ib) * fw2 + cps(4,4,ib) * fw3
           wwi(k,4) = ( 1.0 - fd ) * wf4 + 9.0 * fd
         ELSE ! longwave
           ibr = ib - mbs
           gg = cpir(1,ibr) + cpir(2,ibr) * fw1 + &
                cpir(3,ibr) * fw2 + cpir(4,ibr) * fw3
           x1 = gg
           x2 = x1 * gg
           x3 = x2 * gg
           x4 = x3 * gg
           wwi(k,1) = 3.0 * x1
           wwi(k,2) = 5.0 * x2
           wwi(k,3) = 7.0 * x3
           wwi(k,4) = 9.0 * x4
         END IF
       ELSE
          wwi(k,:) = 0.0
          ti(k) = 0.0
          wi(k) = 0.0
          gg    = 0.
       END IF
    END DO

    RETURN
  END SUBROUTINE cloud_ice


  ! -----------------------------------------------------------------------
  ! Subroutine cloud_grp:  calculates the optical depth (ti), single
  ! scattering albedo (wi), and phase function (wwi(4)) given the cloud
  ! ice [g/m^3] and effective radius [microns] by interpolating based on
  ! known optical properties at predefined sizes
  !
  ! tgr, wgr, and wwgr are the optical depth, single scattering albedo,
  ! and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
  ! due to the Mie scattering of graupel for a given layer.
  !                        Jan. 19, 1993
  ! -----------------------------------------------------------------------

  SUBROUTINE cloud_grp ( ib, pcg, dz, tgr, wgr, wwgr )

    IMPLICIT NONE

    INTEGER, INTENT (in) :: ib
    REAL, DIMENSION (nv), INTENT (in) :: pcg, dz
    REAL, INTENT (out)   :: tgr(nv), wgr(nv), wwgr(nv,4)

    INTEGER :: i
    REAL    :: cwmks
    REAL    :: y1, y2, y3, y4, x1, x2, x3, x4
    IF (.NOT. grpInitialized) STOP 'TERMINATING: Ice not Initialized'

    x1 = gg(ib)
    x2 = x1 * gg(ib)
    x3 = x2 * gg(ib)
    x4 = x3 * gg(ib)
    y1 = 3.0 * x1
    y2 = 5.0 * x2
    y3 = 7.0 * x3
    y4 = 9.0 * x4

    DO i = 1, nv
       cwmks = pcg(i)*1.e-3! convert to km
       IF ( cwmks < 1.0e-8 ) THEN
          tgr(i) = 0.0
          wgr(i) = 0.0
          wwgr(i,1) = 0.0
          wwgr(i,2) = 0.0
          wwgr(i,3) = 0.0
          wwgr(i,4) = 0.0
       ELSE
          tgr(i) = dz(i) * cwmks * bg(ib) / gwc
          wgr(i) = wgf(ib)
          wwgr(i,1) = y1
          wwgr(i,2) = y2
          wwgr(i,3) = y3
          wwgr(i,4) = y4
       END IF
    END DO

  END SUBROUTINE cloud_grp

  ! ---------------------------------------------------------------------------
  ! linear interpolation between two points, returns indicies of the
  ! interpolation points and weights
  !
  SUBROUTINE interpolate(x,ny,y,i1,i2,alpha)

    INTEGER, INTENT (in) :: ny
    REAL, INTENT (in)    :: x, y(ny)

    INTEGER, INTENT (out) :: i1, i2
    REAL, INTENT (out)    :: alpha

    IF (y(1) < y(2)) STOP 'TERMINATING: band centers increasing'

    i2 = 1
    DO WHILE (x < y(i2) .AND. i2 < ny)
       i2 = i2+1
    END DO
    i1 = max(1,i2-1)
    alpha = 1.

    IF (i2 /= i1) alpha = (x-y(i1))/(y(i2)-y(i1))
    IF (alpha < 0 .OR. alpha >1) PRINT 600, x, y(1), y(ny), alpha

    RETURN

600 FORMAT(/'CLOUD_INIT WARNING:  Extrapolating because data out of range', &
           /1x,'x = ',F8.1,', ymax = ',F7.0,', ymin =',F7.0,', alpha = ',F6.3)
  END SUBROUTINE interpolate

END MODULE cldwtr
