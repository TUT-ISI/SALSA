!****************************************************************
!*	                                                        *
!*   module MO_HAM_SALSA_PROPERTIES                                 *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate particle properties during simulation         *
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
!----------------------------------------------------------------------------
MODULE mo_ham_salsa_properties


CONTAINS

  ! fxm: should sea salt form a solid particle when prh is very low
  !  (even though it could be mixed with e.g. sulphate)?
  ! fxm: crashes if no sulphate or sea salt
  ! fxm: do we really need to consider Kelvin effect for regime 2
  !********************************************************************
  !
  ! subroutine WETSIZE()
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates ambient sizes of particles by equilibrating
  !  soluble fraction of particles with water
  !
  !
  ! Method:
  ! ------- 
  ! Following chemical components are assumed water-soluble
  ! - (ammonium) sulphate (100%)
  ! - sea salt (100 %)
  ! - organic carbon (epsoc * 100%)
  !
  ! Exact thermodynamic considerations neglected
  ! - If particles contain no sea salt, calculation according to
  !  sulphate properties
  ! - If contain sea salt but no sulphate, calculation according to
  !  sea salt properties
  ! - If contain both sulphate and sea salt
  !  -> the molar fraction of these compounds determines
  !     which one of them is used as the basis of calculation
  !     
  ! If sulphate and sea salt coexist in a particle,
  !   it is assumed that the Cl is replaced by sulphate;
  !   thus only either sulphate + organics or sea salt + organics
  !   is included in the calculation of soluble fraction.
  !
  ! Molality parameterizations taken from table 1 of
  !  Tang: Mixed-salt aerosols of atmospheric importance,
  !   JGR, 102 (D2), 1883-1893 (1997)
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  ! 
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  ! Harri Kokkola (FMI) 2006
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  !
  !---------------------------------------------------------------------

  SUBROUTINE equilibration(kproma, kbdim, klev,                    &
                           pnaero, pvols, prh, ptemp, pcore, pdwet)

    USE mo_ham_submctl, ONLY : &
         pi6,          & ! pi/6
         in1a, fn1a,   &
         in2a, fn2a,   &
         in2b, fn2b,   &
         boltz,        & ! Boltzmann constant [J/K]
         nlim,         & ! lowest possible particle conc. in a bin [#/m3]
    
                         ! molar masses [kg/mol]
         mwa,          & ! water
                         ! molecular volumes [m3]
         mvsu,         & ! sulphate 
                         ! density [kg/m3]
         rhowa,        & ! water
         
         surfw0,       & ! surface tension of water [J/m2]
         dpmid,        & ! mid dry diameter of each bin [m]
         epsoc,        & ! fxm

         rhosu, msu,   & ! properties of compounds
         rhooc, moc,   &
         rhobc, mbc,   &
         rhoss, mss


    USE mo_kind, ONLY : dp

    IMPLICIT NONE

    !-- input variables -------------
    INTEGER, INTENT(in) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels 

    REAL(dp), INTENT(in) ::        &     
         pnaero(kbdim,klev,fn2b),  & ! particle concentration [#/m3]
         pvols(kbdim,klev,fn2b,5), & ! total volume concentrations of each
                                     ! chem. compound in a size bin [fxm]
         prh(kbdim,klev),          & ! relative humidity [0-1]
         ptemp(kbdim,klev)           ! temperature [K]


    !-- output variables -------------
    REAL(dp), INTENT(out) ::       &
         pcore(kbdim,klev,fn2b),   & ! particle dry volume [fxm]
         pdwet(kbdim,klev,fn2b)      ! particle ambient diameter [m]


    !-- local variables --------------
    INTEGER :: ii, jj, kk             ! loop indices

    REAL(dp) ::      &
         zbinmol(4), &   ! binary molality of individual components [mol/kg]
         zvpart(5),  &   ! volume of chem. compounds in one particle [fxm]
         zke,        &   ! Kelvin term
         zaw,        &   ! water activity [0-1]

         zlwc,       &   ! liquid water content [kg/m3-air]
         zdold,      &   !
         zrh

    !----------------------------------------------------------------------
    !-- 1) Regime 1: sulphate and partly water-soluble OC -----------------
    zke = 1.001_dp
    pdwet = 0._dp
    DO kk = in1a,fn1a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid

             !-- initialize
             zbinmol = 0._dp
             zdold = 1._dp

             IF ((pnaero(ii,jj,kk) > nlim)) THEN

                !-- volume of sulphate and OC in one particle [fxm]

                zvpart(1:2) = pvols(ii,jj,kk,1:2)/pnaero(ii,jj,kk)

                !-- total volume of one dry particle [fxm] 
                pcore(ii,jj,kk)   = sum(zvpart(1:2))

                ! Relative Humidity:
                zrh = prh(ii,jj)
                zrh = MAX(zrh , 0.05_dp)
                zrh = MIN(zrh , 0.95_dp)

                DO WHILE(abs(pdwet(ii,jj,kk)/zdold-1.) > 1.e-2_dp) 
                   zdold = max(pdwet(ii,jj,kk),1.e-20_dp)

                   zaw = zrh/zke

                   !-- binary molalities [mol/kg]
                   zbinmol(1) =                      & ! sulphate
                        + 1.1065495e+2_dp            & 
                        - 3.6759197e+2_dp * zaw      &  
                        + 5.0462934e+2_dp * zaw**2   &
                        - 3.1543839e+2_dp * zaw**3   &
                        + 6.770824e+1_dp  * zaw**4 

                   zbinmol(2) = 1./(zaw*mwa)-1./mwa ! organic carbon

                   !
                   ! Calculate the liquid water content (kg/m3-air) using ZSR
                   ! (see e.g. equation (9.98) in Seinfeld and Pandis (1998))
                   !
                   zlwc = (pvols(ii,jj,kk,1)*(rhosu/msu))/zbinmol(1) + &
                          epsoc * pvols(ii,jj,kk,2)*(rhooc/moc)/zbinmol(2)
                   
                   !-- particle wet radius [m] 
                   pdwet(ii,jj,kk) = (zlwc/pnaero(ii,jj,kk)/rhowa/pi6 + &
                                     pcore(ii,jj,kk)/pi6)**(1._dp/3._dp)

                   zke = exp(2._dp*surfw0*mvsu/(boltz*ptemp(ii,jj)*pdwet(ii,jj,kk)))
                   !-- Kelvin effect 

                END DO

             ELSE
                !-- 1.2) empty bins given bin average values ----------------- 
                pdwet(ii,jj,kk) = dpmid(kk)
                pcore(ii,jj,kk) = pi6*dpmid(kk)**3
             END IF
          END DO
       END DO
    END DO

    !-- 2) Regime 2a: sulphate, OC, BC and sea salt ----------------------------

    ! loops over:


    
    DO kk = in2a,fn2a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid

             zke = 1.02_dp
            
             !-- initialize
             zbinmol = 0._dp

             !-- 1) particle properties calculated for non-empty bins ---------
             IF ((pnaero(ii,jj,kk) > nlim)) THEN

                !-- volume in one particle [fxm]
                zvpart = 0._dp
                zvpart = pvols(ii,jj,kk,:)/pnaero(ii,jj,kk)

                !-- total volume of one dry particle [fxm] 
                pcore(ii,jj,kk) = sum(zvpart)

                ! Relative Humidity:
                zrh = prh(ii,jj)
                zrh = MAX(zrh , 0.37_dp)
                zrh = MIN(zrh , 0.95_dp)

                DO WHILE(abs(pdwet(ii,jj,kk)/zdold-1.) > 1.e-12_dp)
                   zdold = max(pdwet(ii,jj,kk),1.e-20_dp)

                zaw = zrh/zke

                !-- binary molalities [mol/kg]
                zbinmol(1) =                     &  ! sulphate
                     + 1.1065495e+2_dp           & 
                     - 3.6759197e+2_dp * zaw     &  
                     + 5.0462934e+2_dp * zaw**2  &
                     - 3.1543839e+2_dp * zaw**3  &
                     + 6.770824e+1_dp  * zaw**4 

                zbinmol(2) = 1._dp/(zaw*mwa)-1._dp/mwa ! organic carbon

                zbinmol(4) =                     &  ! sea salt (NaCl)
                     + 5.875248e+1_dp            &  ! 
                     - 1.8781997e+2_dp * zaw     &  
                     + 2.7211377e+2_dp * zaw**2  &
                     - 1.8458287e+2_dp * zaw**3  &
                     + 4.153689e+1_dp  * zaw**4 

                !-- calculate the liquid water content (kg/m3-air)
                zlwc = (pvols(ii,jj,kk,1)*(rhosu/msu))/zbinmol(1) +                 &
                       epsoc * (pvols(ii,jj,kk,2)*(rhooc/moc))/zbinmol(2) +         &
                       (pvols(ii,jj,kk,4)*(rhoss/mss))/zbinmol(4)

                !-- particle wet radius [m] 
                pdwet(ii,jj,kk) = (zlwc/pnaero(ii,jj,kk)/rhowa/pi6 +  &
                     pcore(ii,jj,kk)/pi6)**(1._dp/3._dp)

                !-- Kelvin effect 

                zke = exp(2._dp*surfw0*mvsu/(boltz*ptemp(ii,jj)*pdwet(ii,jj,kk)))

                END DO


             ELSE
                !-- 2.2) empty bins given bin average values ------------------------- 
                pdwet(ii,jj,kk) = dpmid(kk)
                pcore(ii,jj,kk) = pi6*dpmid(kk)**3
             END IF

             !-- calculate the wet and dry radius for subregime 2b
             IF (kk == in2a) THEN
                pdwet(ii,jj,in2b:fn2b) = dpmid(in2b:fn2b)
                pcore(ii,jj,in2b:fn2b) = pi6*dpmid(in2b:fn2b)**3
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE equilibration

END MODULE mo_ham_salsa_properties
