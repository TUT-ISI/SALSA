!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsa_update
!!
!! \brief
!! Contains subroutines and functions that are used  
!! to calculate aerosol sizedistribution update
!!
!! \author Harri Kokkola (FZ Juelich)
!!
!! \responsible_coder
!! Harri Kokkola, Harri.Kokkola@fmi.fi
!!
!! \revision_history
!!   -# Martin G. Schultz (FZ Juelich) - original code (2009-10)
!!   -# Harri Kokkola (FMI) - Implementation of SALSA aerosol microphysics model (2014)
!!   
!! \limitations
!! None
!!
!! \details
!! None
!!
!! \bibliographic_references
!! None
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !****************************************************************
!*                                                              *
!*   module MO_HAM_SALSA_UPDATE                                     *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************

MODULE mo_ham_salsa_update

CONTAINS

  SUBROUTINE distr_update(kproma, kbdim, klev, &
                          pnaero, pvols, ptemp)

    USE mo_ham_salsactl
    USE mo_kind, ONLY : dp
    USE mo_ham, ONLY: subm_naerospec_nowat
    USE mo_math_constants, ONLY : pi, pi
    USE mo_physical_constants, ONLY : argas
    USE mo_species, ONLY :                      &
         speclist

    USE mo_ham_species, ONLY:                   &
         id_so4, id_wat


    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels

    REAL(dp), INTENT(IN) ::         &
         ptemp(kbdim,klev)

    REAL(dp), INTENT(INOUT) ::      &
         pnaero(kbdim,klev,fn2b),   & ! number concentration of particles in each bin [#/m3]
         pvols(kbdim,klev,fn2b,subm_naerospec_nowat)     ! volume concentration of each chem. compound in regimes 1&2 [fxm]

    !-- Local variables ----------------------
    INTEGER :: ii, jj, kk, mm
    REAL(dp) :: zvpart, znfrac, zvfrac, zVrat, zVilo, zVihi, zVexc, zaa, zbb, zepsv, zepspart
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
                   IF(zvpart <= pi/6.0*dpmid(kk)**3) CYCLE
   
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


    IF (lsol2b) THEN !eehol: included flag for moving of soluble fraction

       !-- Move soluble fraction of particles in subregime 2b to subregime 2a -
       DO kk = in2b,fn2b
          DO jj = 1,klev
             DO ii = 1,kproma
                
                IF (pnaero(ii,jj,kk) > nlim) THEN

                   !-- corresponding bin in regime 2a
                   
                   mm = kk - in2b + in2a  
                   
                   !-- soluble volume fraction of particles
                   
                   zepspart = (pvols(ii,jj,kk,1)+epsoc*pvols(ii,jj,kk,2))/sum(pvols(ii,jj,kk,:))
                   
                   zaa = 4.*speclist(id_wat)%moleweight/1000._dp*surfw0/(argas*speclist(id_wat)%density)   ! curvature (Kelvin) effect
                   zbb = 6.*speclist(id_wat)%moleweight/1000._dp*ions*speclist(id_so4)%density/ & 
                        (pi*speclist(id_wat)%density*speclist(id_so4)%moleweight/1000.) ! solute (Raoult) effect                
                   
                   zepsv = 4._dp*zaa**3/(27._dp*zbb*(log(slim))**2*(pi/6.0*dpmid(kk)**3))                
                   
                   IF (zepsv/ptemp(ii,jj)**3 < zepspart) THEN                   
                      
                      pnaero(ii,jj,mm) = pnaero(ii,jj,mm) + pnaero(ii,jj,kk)
                      pnaero(ii,jj,kk) = 0._dp
                      
                      pvols(ii,jj,mm,:) = pvols(ii,jj,mm,:) + pvols(ii,jj,kk,:)   
                      pvols(ii,jj,kk,:) = 0._dp                 
                      
                   END IF
                END IF
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE distr_update

END MODULE mo_ham_salsa_update
