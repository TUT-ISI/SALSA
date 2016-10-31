!****************************************************************
!*                                                              *
!*   module MO_HAM_SALSA_UPDATE                                     *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************
!----------------------------------------------------------------------------
!!$   Copyright 2014 Atmospheric Research Centre of Eastern Finland,
!!$         Finnish Meteorological Institute, Kuopio, Finland
!!$
!!$   Licensed under the Apache License, Version 2.0 (the "License");
!!$   you may not use this file except in compliance with the License.
!!$   You may obtain a copy of the License at
!!$
!!$       http://www.apache.org/licenses/LICENSE-2.0
!!$
!!$   Unless required by applicable law or agreed to in writing, software
!!$   distributed under the License is distributed on an "AS IS" BASIS,
!!$   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!$   See the License for the specific language governing permissions and
!!$   limitations under the License.
!------------------------------------------------------------------------
MODULE mo_ham_salsa_update

CONTAINS

  SUBROUTINE distr_update(kproma, kbdim, klev, &
                          pnaero, pvols)

    USE mo_ham_submctl
    USE mo_kind, ONLY : dp

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels

    REAL(dp), INTENT(INOUT) ::      &
         pnaero(kbdim,klev,fn2b),   & ! number concentration of particles in each bin [#/m3]
         pvols(kbdim,klev,fn2b,5)     ! volume concentration of each chem. compound in regimes 1&2 [fxm]

    !-- Local variables ----------------------
    INTEGER :: ii, jj, kk, mm
    REAL(dp) :: zvpart, znfrac, zvfrac, zVrat, zVilo, zVihi, zVexc
    LOGICAL  :: within_bins

    DO jj = 1,klev
       DO ii = 1,kproma
          
          within_bins = .FALSE.

          !-- Check if the volume of the bin is within bin limits after update

          DO WHILE(.NOT.within_bins)
             within_bins = .TRUE.

          DO kk = fn2b-1,in1a,-1

             IF (pnaero(ii,jj,kk) > nlim) THEN
             IF (kk == fn2a) CYCLE

                zvpart = sum(pvols(ii,jj,kk,:))/pnaero(ii,jj,kk)

                IF (vlolim(kk) > zvpart .AND. kk == in1a) CYCLE

                IF(zvpart < vlolim(kk)) THEN
                   mm = kk - 1
                   IF(kk == in2b) mm = fn1a

                   pnaero(ii,jj,mm) = pnaero(ii,jj,mm) + pnaero(ii,jj,kk) 
                   pnaero(ii,jj,kk) = 0._dp
                   pvols(ii,jj,mm,:) = pvols(ii,jj,mm,:) + pvols(ii,jj,kk,:)
                   pvols(ii,jj,kk,:) = 0._dp
                   
                END IF

                !-- If size bin has not grown, cycle
                IF(zvpart <= pi6*dpmid(kk)**3) CYCLE

                !-- volume ratio of the size bin
                zVrat = vhilim(kk)/vlolim(kk)
                
                !-- particle volume at the low end of the bin
                zVilo = 2._dp*zvpart/(1._dp + zVrat)

                !-- particle volume at the high end of the bin
                zVihi = zVrat * zVilo

                !-- volume in the grown bin which exceeds 
                !   the bin upper limit
                zVexc = 1._dp/2._dp*(zVihi + vhilim(kk))

                !-- number fraction to be moved to the larger bin
                znfrac = min(1._dp,(zVihi-vhilim(kk)) / (zVihi - zVilo))
          
                !-- volume fraction to be moved to the larger bin
                zvfrac = znfrac * zVexc / (1._dp/2._dp*(zVihi+zVilo))

                !-- update bin
                mm = kk+1
                !-- volume
                pvols(ii,jj,mm,:) = pvols(ii,jj,mm,:) &
                     + znfrac * pnaero(ii,jj,kk) * zVexc * pvols(ii,jj,kk,:) / sum(pvols(ii,jj,kk,:))

                pvols(ii,jj,kk,:) = pvols(ii,jj,kk,:) &
                     - znfrac * pnaero(ii,jj,kk) * zVexc * pvols(ii,jj,kk,:) / sum(pvols(ii,jj,kk,:))

                !-- number
                pnaero(ii,jj,mm) = pnaero(ii,jj,mm) + znfrac * pnaero(ii,jj,kk)

                pnaero(ii,jj,kk) = pnaero(ii,jj,kk) * (1._dp - znfrac)


             END IF

             IF ( pnaero(ii,jj,kk) > nlim ) THEN
                zvpart = sum(pvols(ii,jj,kk,:))/pnaero(ii,jj,kk)

                IF(zvpart > vhilim(kk)) within_bins = .FALSE.

             END IF

          END DO ! - kk
          END DO ! - within_bins
       END DO    ! - ii
    END DO       ! - jj

665 FORMAT(999(E11.4,1X))

  END SUBROUTINE distr_update

END MODULE mo_ham_salsa_update
