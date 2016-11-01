

!****************************************************************
!*                                                              *
!*   module MO_HAM_SALSA_DYNAMICS                               *
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
!----------------------------------------------------------------------------
MODULE mo_salsa_dynamics


CONTAINS

  ! this calculated for empty bins too!!!
  ! fxm: test well, esp. self-coagulation (but other bits too!)
  ! AL_note: Diagnostic variables of cond and nucl mass
  !********************************************************************
  !
  ! subroutine COAGULATION(kproma,kbdim,klev, &
  !       pnaero,pvols,pdwet, &
  !       pcore, ptstep)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates particle loss and change in size distribution
  !  due to (Brownian) coagulation
  !
  !
  ! Method:
  ! -------  
  ! Semi-implicit, non-iterative method:
  !  Volume concentrations of the smaller colliding particles
  !  added to the bin of the larger colliding particles.
  !  Start from first bin and use the updated number and volume
  !  for calculation of following bins. NB! Our bin numbering
  !  does not follow particle size in regime 2.
  !
  !Schematic for bin numbers in different regimes:
  !        	 1             			2             
  !    +-------------------------------------------+
  !  a | 1 | 2 | 3 || 4 | 5 | 6 | 7 |  8 |  9 | 10||
  !  b |           ||11 |12 |13 |14 | 15 | 16 | 17||
  !    +-------------------------------------------+
  !
  ! Exact coagulation coefficients for each pressure level
  !  are calculated in subroutine SET_COAGC (in mo_salsa_init) 
  !  which is called once at the beginning of the simulation 
  !  from model driver. In subroutine COAGULATION, these exact 
  !  coefficients are scaled according to current particle wet size
  !  (linear scaling).
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
  ! Tommi Bergman (FMI) 2012
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  !
  !---------------------------------------------------------------------


  SUBROUTINE coagulation(kproma, kbdim,  klev,        &
                         pnaero, pvols,  pdwet,       &
                         pcore,  ptstep, ptemp, ppres)

    USE mo_submctl, ONLY:        &
         in1a,                       & ! size bin indices
         in2a, fn2a,                 &
         in2b, fn2b,                 &
         dpmid ,&                      ! mid dry diameter of size bins [m]
         pi6
    !USE mo_aero_mem_salsa, ONLY :    &
    !     d_cond_so4,                 & ! diagnostic variable for condensated mass of so4
    !     d_nuc_so4                     !diagnostic variable for nucleated mass of so4
    USE mo_salsa_init, only: coagc
    USE mo_kind, ONLY : dp

    IMPLICIT NONE


    !-- Input and output variables -------------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical klev 

    REAL(dp), INTENT(INOUT) ::      &     
         pnaero(kbdim,klev,fn2b),   & ! particle concentration [#/m3]
         pvols(kbdim,klev,fn2b,5)     ! total volume concentrations of each
                                      !   chem. compound in a size bin [fxm]

    REAL(dp), INTENT(IN) ::         &
         pdwet(kbdim,klev,fn2b),    & ! particle wet radius [m]
         pcore(kbdim,klev,fn2b),    & ! particle dry volume [fxm]
         ptstep,                    & ! time step [s]
         ptemp(kbdim,klev),         &
         ppres(kbdim,klev)
    !-- Local variables ------------------------
    INTEGER ::                      &
         ii, jj, kk, ll, mm, nn,    & ! loop indices 
         index_2a, index_2b           ! corresponding bin in regime 2a/2b

    REAL(dp) ::                     &
         zntemp(kbdim,klev),        & ! variable for saving pnaero(fn2b) temporarily
         zcc(fn2b,fn2b), & ! updated coagulation coefficients [m3/s]
         zminusterm,                & ! coagulation loss in a bin [1/s] 
         zplusterm(5)                 ! coagulation gain in a bin [fxm/s]
                                      ! (for each chemical compound)
         !zcc(kbdim,klev,fn2b,fn2b), & ! updated coagulation coefficients [m3/s]

    REAL(dp) :: &
         zmpart(fn2b)   ! approximate mass of particles [kg]
    REAL(dp) :: &
         temppi,pressi,pdmm,pdnn

    !-----------------------------------------------------------------------------
    !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
    !      CoagSink ~ Dp in continuum regime, thus we calculate
    !      'effective' number concentration of coarse particles

    zntemp = pnaero(:,:,fn2b)

    !-- 2) Updating coagulation coefficients -------------------------------------

    !-- particle mass; density of 1500 kg/m3 assumed [kg] 
    zmpart = pi6*dpmid(1:fn2b)**3*1500.


     DO jj = 1,klev      ! vertical grid
        DO ii = 1,kproma ! horizontal kproma in the slab
           temppi=ptemp(ii,jj)
           pressi=ppres(ii,jj)

        
           DO mm = 1,fn2b         ! smaller colliding particle
              DO nn = mm,fn2b            ! larger colliding particle 
              pdmm=pdwet(ii,jj,mm)
              pdnn=pdwet(ii,jj,nn)
              zcc(mm,nn) = coagc(dpmid(mm),dpmid(nn),zmpart(mm),zmpart(nn),temppi,pressi)&
                     *dpmid(mm)*pdnn/(dpmid(nn)*pdmm)
              zcc(nn,mm)=zcc(mm,nn)
             END DO
          END DO
    
	! Original version:
    ! loops over
    !DO nn = 1,fn2b            ! larger colliding particle
    !   DO mm = 1,fn2b         ! smaller colliding particle
    !      DO jj = 1,klev      ! vertical grid
    !         DO ii = 1,kproma ! horizontal kproma in the slab

               !zcc(ii,jj,mm,nn) = coagtable(ii,jj,mm,nn)*dpmid(mm)*pdwet(ii,jj,nn)/(dpmid(nn)*pdwet(ii,jj,mm))
    !            zcc(ii,jj,mm,nn) = coagc(dpmid(mm),dpmid(nn),zmpart(mm),zmpart(nn),ptemp(ii,jj),ppres(ii,jj))*dpmid(mm)*pdwet(ii,jj,nn) &
    !                 /(dpmid(nn)*pdwet(ii,jj,mm))

    !         END DO
    !      END DO
    !   END DO
    !END DO

    !-- 3) New particle and volume concentrations after coagulation -------------

    ! loops over
    DO kk = in1a,fn2b      ! bins that we want to update
       !DO jj = 1,klev      ! vertical grid
          !DO ii = 1,kproma ! horizontal kproma in the slab

             !-- 3.1) Bin in regime 1 -----------------------------------

             IF (kk < in2a) THEN

                !-- Particles lost by coagulation with larger particles
                !zminusterm = sum(zcc(ii,jj,kk,kk+1:fn2b)*pnaero(ii,jj,kk+1:fn2b))
                zminusterm = sum(zcc(kk,kk+1:fn2b)*pnaero(ii,jj,kk+1:fn2b))

                !-- Particle volume gained by coagulation with smaller particles

                zplusterm = 0._dp

                IF(kk > in1a) THEN

                   DO ll = in1a, kk-1
                      
                      !zplusterm(1:2) = zplusterm(1:2) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:2)
                      zplusterm(1:2) = zplusterm(1:2) + zcc(ll,kk)*pvols(ii,jj,ll,1:2)                      
                   END DO

                END IF

                !-- Volume and number concentrations after coagulation update [fxm]
                !pvols(ii,jj,kk,1:2) = pvols(ii,jj,kk,1:2)/(1._dp + ptstep*(zminusterm - 1./pcore(ii,jj,kk)*zplusterm(1:2)))
                pvols(ii,jj,kk,1:2) = (pvols(ii,jj,kk,1:2)+ptstep*zplusterm(1:2)*pnaero(ii,jj,kk))/ &
                                    (1._dp + ptstep*zminusterm) 
                !pnaero(ii,jj,kk)    = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                !                      0.5_dp*ptstep*zcc(ii,jj,kk,kk)*pnaero(ii,jj,kk))
                pnaero(ii,jj,kk)    = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                                      0.5_dp*ptstep*zcc(kk,kk)*pnaero(ii,jj,kk))
             ELSE 

                !-- 3.2) Bin in regime 2a -----------------------------------

                IF (kk < in2b) THEN

                   zplusterm = 0._dp

                   !-- Find corresponding size bin in subregime 2b
                   index_2b = kk - in2a + in2b

                   !-- Particles lost to larger particles in regimes 2a, 2b and 3
                   zminusterm = 0._dp

                   DO ll = kk+1, fn2a

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll) ! 2a
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll) ! 2a
                   END DO

                   DO ll = index_2b+1, fn2b

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll) ! 2b 
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll) ! 2b + 3

                   END DO

                   !-- Particle volume gained from smaller particles in regimes 1, 2a and 2b

                   ! sulphate

                   DO ll = in1a, kk-1

                      !zplusterm(1:2) = zplusterm(1:2) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:2) ! 1 + 2a
					  zplusterm(1:2) = zplusterm(1:2) + zcc(ll,kk)*pvols(ii,jj,ll,1:2) ! 1 + 2a

                   END DO
                   
                   DO ll = in2a, kk-1

                      !zplusterm(3:5) = zplusterm(3:5) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,3:5)
                      zplusterm(3:5) = zplusterm(3:5) + zcc(ll,kk)*pvols(ii,jj,ll,3:5)
                   
                   END DO
                      
                   DO ll = in2b, index_2b

                      !zplusterm = zplusterm + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,:) ! 2b
                      zplusterm = zplusterm + zcc(ll,kk)*pvols(ii,jj,ll,:) ! 2b

                   END DO
                   
                   !-- 3.2) Bin in regime 2b -----------------------------------

                ELSE

                   !-- Find corresponding size bin in subregime 2a
                   index_2a = kk - in2b + in2a

                   !-- Particles lost to larger particles in regimes 2b and 3, 
                   !  as well as to 'same sized' and larger particles in regime 2a             

                   zminusterm = 0._dp
                   zplusterm  = 0._dp

                   DO ll = kk+1, fn2b

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll) ! 2b + 3
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll) ! 2b + 3

                   END DO

                   DO ll = index_2a, fn2a

                      !zminusterm = zminusterm + zcc(ii,jj,kk,ll)*pnaero(ii,jj,ll)
                      zminusterm = zminusterm + zcc(kk,ll)*pnaero(ii,jj,ll)

                   END DO

                   !-- Particle volume gained from smaller particles in regimes 1, 2a and 2b
                   !  (except not corresponding bin in regime 2a even if it is slightly smaller)

                   ! sulphate

                   DO ll = in1a, index_2a-1

                      !zplusterm(1:2) = zplusterm(1:2) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:2)
                      zplusterm(1:2) = zplusterm(1:2) + zcc(ll,kk)*pvols(ii,jj,ll,1:2)
                   END DO


                   DO ll = in2a, index_2a-1

                      !zplusterm(3:5) = zplusterm(3:5) + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,3:5)
                      zplusterm(3:5) = zplusterm(3:5) + zcc(ll,kk)*pvols(ii,jj,ll,3:5)

                   END DO

                   DO ll = in2b, kk-1

                      !zplusterm = zplusterm + zcc(ii,jj,ll,kk)*pvols(ii,jj,ll,1:5)
                      zplusterm = zplusterm + zcc(ll,kk)*pvols(ii,jj,ll,1:5)

                   END DO

                END IF

                
                !-- Volume and number concentrations after coagulation update [fxm]
                !pvols(ii,jj,kk,:) = pvols(ii,jj,kk,:)/ &
                !                    (1._dp + ptstep*(zminusterm - 1./pcore(ii,jj,kk)*zplusterm))
                pvols(ii,jj,kk,:) = (pvols(ii,jj,kk,:)+ptstep*zplusterm*pnaero(ii,jj,kk))/ &
                                    (1._dp + ptstep*zminusterm)

                !pnaero(ii,jj,kk) = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                !                   0.5_dp*ptstep*zcc(ii,jj,kk,kk)*pnaero(ii,jj,kk))
                pnaero(ii,jj,kk) = pnaero(ii,jj,kk)/(1. + ptstep*zminusterm  + &
                                   0.5_dp*ptstep*zcc(kk,kk)*pnaero(ii,jj,kk))

             END IF
          END DO
       END DO
    END DO

    ! fxm: here we approximate that the sea salt regime 2b particles have
    ! gained by coagulation can be treated as sulphate
    pvols(:,:,in2b:fn2b,1) = pvols(:,:,in2b:fn2b,1) + pvols(:,:,in2b:fn2b,4)
    pvols(:,:,in2b:fn2b,4) = 0._dp

    pvols = MAX(pvols, 0._dp)

  END SUBROUTINE coagulation


  ! fxm: calculated for empty bins too
  ! fxm: same diffusion coefficients and mean free paths used for sulphuric acid
  !      and organic vapours (average values? 'real' values for each?)
  !********************************************************************
  !
  ! subroutine CONDENSATION(kproma, kbdim, klev, &
  !                         pnaero, pvols, pdwet,     &
  !                         pcsa, pcocnv, pcocsv,     &
  !                         ptemp, ppres, ptstep)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates the increase in particle volume and 
  !  decrease in gas phase concentrations due to condensation 
  !  of sulphuric acid and two organic compounds (non-volatile
  !  and semivolatile)
  !
  !
  ! Method:
  ! -------
  ! Regime 3 particles only act as a sink for condensing vapours
  !  while their size and composition does not change.
  ! Exception: Soluble fraction of regime 3c particles can change
  !  and thus they can be moved to regime 3b 
  !
  ! New gas and aerosol phase concentrations calculated according
  !  to Jacobson (1997): Numerical techniques to solve 
  !  condensational and dissolutional growth equations 
  !  when growth is coupled to reversible reactions, 
  !  Aerosol Sci. Tech., 27, pp 491-498.
  !
  ! fxm: one should really couple with vapour production and loss terms as well
  !      should nucleation be coupled here as well????
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
  !
  !---------------------------------------------------------------
  !
  ! Following parameterization has been used:
  ! ------------------------------------------
  !
  ! Molecular diffusion coefficient of condensing vapour [m2/s]
  !  (Reid et al. (1987): Properties of gases and liquids,
  !   McGraw-Hill, New York.)
  !
  ! D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / &
  !  {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
  !   
  ! M_air = 28.965 : molar mass of air [g/mol]
  ! d_air = 19.70  : diffusion volume of air
  ! M_h2so4 = 98.08  : molar mass of h2so4 [g/mol]
  ! d_h2so4 = 51.96  : diffusion volume of h2so4
  !
  !---------------------------------------------------------------

  SUBROUTINE condensation(kproma, kbdim,  klev,   krow,              &
                          pnaero, pvols,  pdwet,                     &
                          pcsa,   pcocnv, pcocsv,                    &
                          prh,    ptemp,  ppres,  ptstep, ppbl)

    USE mo_salsa_nucleation

    USE mo_submctl,    ONLY :   &
         pi,                        & 
         pi6,                       & ! pi/6 
         in1a, in2a,                & ! size bin indices
         fn2b,                      & 
         nbin,                      & ! number of size bins in each regime
         
         boltz,                     & ! Boltzmann constant [J/K]
         rg,                        & ! molar gas constant [J/(mol K)]
         pstand,                    & ! standard pressure [Pa]
         msu,                       & ! molar mass of sulphate [kg/mol]
         mvsu, mvoc,                & ! molecular volumes of sulphate and OC [kg]
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         
         epsoc,                     & ! soluble fraction of organics (scaled to sulphate)
         massacc,                   & ! mass accomodation coefficients in each bin
         dpmid,                     & ! mid dry diameter of each bin [m]
         n3,                        & ! number of molecules in one 3 nm particle [1]
         nsnucl,                   & ! nucleation
         rhosu
    USE mo_kind, ONLY : dp
    !USE mo_aero_mem_salsa, only: d_nuc_so4,d_cond_so4
    USE mo_time_control,   ONLY: delta_time,time_step_len
    USE mo_constants,      ONLY: g
   IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev,                      & ! number of vertical klev 
         krow

    REAL(dp), INTENT(IN) ::        &  
         pdwet(kbdim,klev,fn2b),   & ! wet diameter of particles in each bin [m]
         prh(kbdim,klev),          & ! relative humidity [0-1]
         ptemp(kbdim,klev),        & ! ambient temperature [K]
         ppres(kbdim,klev),        & ! ambient pressure [Pa]
         ptstep                     ! timestep [s]


    INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level

    REAL(dp), INTENT(INOUT) ::     &
         pnaero(kbdim,klev,fn2b),  & ! number concentration of particles in each bin [#/m3]
         pvols(kbdim,klev,fn2b,5), & ! volume concentration of each chem. compound in regimes 1&2 [fxm]
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev)          ! semivolatile organic concentration [#/m3]


    !-- Local variables ----------------------
    INTEGER :: ii, jj                ! loop indices

    REAL(dp) ::                    &
         zvisc,                    & ! viscosity of air [kg/(m s)]
         zdfvap,                   & ! air diffusion coefficient [m2/s]
         zmfp,                     & ! mean free path of condensing vapour [m]
         zcs_tot,                  & ! total condensation sink [1/s]
         zcs_ocsv,                 & ! condensation sink for semivolatile organics [1/s]
         zcs_su,                   & ! condensation sink for sulfate [1/s]
         zcs_ocnv,                 & ! condensation sink for nonvolatile organics [1/s]
                                     ! vapour concentration after time step [#/m3]
         zcvap_new1,               & ! sulphuric acid
         zcvap_new2,               & ! nonvolatile organics
         zcvap_new3,               & ! semivolatile organics
                                     ! change in vapour concentration [#/m3]
         zdvap1,                   & ! sulphuric acid
         zdvap2,                   & ! nonvolatile organics
         zdvap3,                   & ! semivolatile organics
         
         zdfpart(in1a+1),          & ! particle diffusion coefficient
         zknud(fn2b),              & ! particle Knudsen number
         zbeta(fn2b),              & ! transitional correction factor
         zcolrate(fn2b),           & ! collision rate of molecules to particles [1/s]
         zcolrate_ocnv(fn2b),      & ! collision rate of organic molecules to particles [1/s]
         zdvolsa(fn2b),            & ! change of sulphate volume in each bin [fxm]
         zdvoloc(fn2b),            & !    - " - organics 
         zj3n3(kbdim,klev,2),      & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zn_vs_c,                  & ! ratio of nucleation of all mass transfer in the smallest bin
         zxsa(kbdim,klev),         & ! ratio of sulphuric acid and organic vapor in 3nm particles 
         zxocnv(kbdim,klev)         

    real(dp)::ztmst,zqtmst
    !-- Initializations
    ztmst = time_step_len
    zqtmst = 1.0_dp/ztmst

    zj3n3 = 0._dp
    
    !------------------------------------------------------------------------------
    
    IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,                     &
                                 pnaero, pdwet,                                    &
                                 ptemp,  prh,    ppres,                            &
                                 pcsa,   pcocnv, ptstep, zj3n3, zxsa, zxocnv, ppbl)
    zdvolsa=0._dp
    zn_vs_c=0._dp
    DO jj = 1,klev
       DO ii = 1,kproma

          zdvoloc = 0._dp

          !-- 1) Properties of air and condensing gases --------------------

          zvisc  = (7.44523e-3_dp*ptemp(ii,jj)**1.5_dp)/(5093._dp*(ptemp(ii,jj)+110.4_dp))! viscosity of air [kg/(m s)] 
          zdfvap = 5.1111e-10_dp*ptemp(ii,jj)**1.75_dp*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]
          zmfp   = 3._dp*zdfvap*sqrt(pi*msu/(8._dp*rg*ptemp(ii,jj)))                      ! mean free path [m]


          !-- 2) Transition regime correction factor for particles ---------
          !  
          !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
          !  Topics in current aerosol research, Pergamon.  
          !
          !  Size of condensing molecule considered only for 
          !  nucleation mode (3 - 20 nm) 
          !

          !-- particle Knudsen number
          zknud(in1a:in1a+1) = 2.*zmfp/(pdwet(ii,jj,in1a:in1a+1)+d_sa) 
          zknud(in1a+2:fn2b) = 2.*zmfp/pdwet(ii,jj,in1a+2:fn2b)

          !-- transitional correction factor
          zbeta = (zknud + 1.)/(0.377_dp*zknud+1._dp+4._dp/ &
                  (3._dp*massacc)*(zknud+zknud**2))  

          !-- 3) Collision rate of molecules to particles -------------------
          !
          !  Particle diffusion coefficient considered only for 
          !  nucleation mode (3 - 20 nm) 
          !

          !-- particle diffusion coefficient [m2/s]
          zdfpart = boltz*ptemp(ii,jj)*zbeta(in1a:in1a+1)/ &    
                    (3._dp*pi*zvisc*pdwet(ii,jj,in1a:in1a+1))  

          !-- collision rate [1/s]
          zcolrate(in1a:in1a+1) = 2._dp*pi*(pdwet(ii,jj,in1a:in1a+1)+d_sa)* & 
                                  (zdfvap+zdfpart)*zbeta(in1a:in1a+1)* &
                                  pnaero(ii,jj,in1a:in1a+1)
         
          zcolrate(in1a+2:fn2b) = 2._dp*pi*pdwet(ii,jj,in1a+2:fn2b)*zdfvap* &
                                  zbeta(in1a+2:fn2b)*pnaero(ii,jj,in1a+2:fn2b)

          !-- 4) Condensation sink [1/s] -------------------------------------

          zcs_tot = sum(zcolrate)             ! total sink    
          
          !-- 5) Changes in gas-phase concentrations and particle volume -----
          !
          !  Sulphuric acid and non-volatile organic compound
          !  condense onto all sized particles. Semivolatile
          !  organic compound condenses only onto particles
          !  in regimes 2 and 3.
          !
          !  Regime 3 particles act only as a sink for vapours
          !  and their size and composition are not updated.
          !  (except for soluble fraction of subregime 3c) 
          !
          
          
          !--- 5.1) Organic vapours ------------------------
          
          !---- 5.1.1) Non-volatile organic compound: condenses onto all bins 

          IF(pcocnv(ii,jj) > 1.e-10_dp .and. zcs_tot > 1.e-30_dp) THEN

             zn_vs_c = 0._dp

             IF(zj3n3(ii,jj,2) > 1._dp) zn_vs_c = (zj3n3(ii,jj,2))/(zj3n3(ii,jj,2) + &
                                                  pcocnv(ii,jj) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrate_ocnv = zcolrate
             zcolrate_ocnv(in1a) = zcolrate_ocnv(in1a) + zj3n3(ii,jj,2)/pcocnv(ii,jj)

             zcs_ocnv = zcs_tot + zj3n3(ii,jj,2)/pcocnv(ii,jj) ! total sink for organic vapor

             zcvap_new2 = pcocnv(ii,jj)/(1.+ptstep*zcs_ocnv)   ! new gas phase concentration [#/m3]
             zdvap2 = pcocnv(ii,jj) - zcvap_new2               ! change in gas concentration [#/m3]
             pcocnv(ii,jj) = zcvap_new2                        ! updating vapour concentration [#/m3]
             
             zdvoloc = zcolrate_ocnv(in1a:fn2b)/zcs_ocnv*mvoc*zdvap2 ! volume change of particles 
                                                                     !  [m3(OC)/m3(air)]

             pvols(ii,jj,in1a:fn2b,2) = pvols(ii,jj,in1a:fn2b,2) + & !-- change of volume
                                        zdvoloc                      !   due to condensation in 1a-2b

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
             ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
             IF (zxocnv(ii,jj) > 0._dp) THEN 
                pnaero(ii,jj,in1a) = pnaero(ii,jj,in1a) + &
                     zn_vs_c * zdvoloc(in1a)/mvoc/(n3*zxocnv(ii,jj))
             END IF
             
          END IF
          !---- 5.1.2) Semivolatile organic compound: regimes 1, 2 and 3
          
          IF(pcocsv(ii,jj) > 1.e-10_dp .and. sum(zcolrate(in2a:fn2b)) > 1.e-30_dp) THEN

             zcs_ocsv = sum(zcolrate(in2a:fn2b))               ! sink for semivolatile organics
             zcvap_new3 = pcocsv(ii,jj)/(1.+ptstep*zcs_ocsv)   ! new gas phase concentration [#/m3]
             zdvap3 = pcocsv(ii,jj) - zcvap_new3               ! change in gas concentration [#/m3]
             pcocsv(ii,jj) = zcvap_new3                        ! updating gas concentration [#/m3]
             
             zdvoloc(in2a:fn2b) = zdvoloc(in2a:fn2b) + &       ! volume change of particles 
                  zcolrate(in2a:fn2b)/zcs_ocsv*mvoc*zdvap3     !  [m3(OC)/m3(air)]

             pvols(ii,jj,in1a:fn2b,2) = &                      !-- change of volume due
                  pvols(ii,jj,in1a:fn2b,2) + zdvoloc           !   due to condensation in 1a-2b
             
          END IF
          ! nucleh2so4=nucleh2so4+pcsa(ii,jj)
          IF(pcsa(ii,jj) > 1.e-10_dp .and. zcs_tot > 1.e-30_dp) THEN

             !-- Ratio of mass transfer betwwn nucleation and condensation

             zn_vs_c = 0._dp

             IF(zj3n3(ii,jj,1) > 1._dp) zn_vs_c = (zj3n3(ii,jj,1)) / &
                                              (zj3n3(ii,jj,1) +  &
                                              pcsa(ii,jj) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrate(in1a) = zcolrate(in1a) + zj3n3(ii,jj,1) / pcsa(ii,jj)

             zcs_su = zcs_tot + zj3n3(ii,jj,1) / pcsa(ii,jj)      ! total sink for sulfate

             !--- 5.2) Sulphuric acid -------------------------
             !
             zcvap_new1 = pcsa(ii,jj) /(1.+ptstep*zcs_su)         ! new gas phase concentration [#/m3]
             zdvap1 = pcsa(ii,jj) - zcvap_new1                    ! change in gas concentration [#/m3]
             pcsa(ii,jj) = zcvap_new1                             ! updating vapour concentration [#/m3]
             
             zdvolsa = zcolrate(in1a:fn2b)/zcs_su*mvsu*zdvap1     ! volume change of particles
             ! [m3(SO4)/m3(air)] by condensation

             !-- Change of volume concentration of sulphate in aerosol [fxm]
             pvols(ii,jj,in1a:fn2b,1) = pvols(ii,jj,in1a:fn2b,1) + zdvolsa

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             IF (zxsa(ii,jj) > 0._dp) THEN
                pnaero(ii,jj,in1a) = pnaero(ii,jj,in1a) +          &
                     zn_vs_c * zdvolsa(in1a)/mvsu/(n3*zxsa(ii,jj))
             END IF

             ! Diagnostic for nucleation of so4

             !d_nuc_so4(ii,krow)=d_nuc_so4(ii,krow)+zn_vs_c*zdvolsa(in1a)*pdpg(ii,jj)*rhosu*zqtmst*delta_time/3.0_dp
             ! Diagnostic for condensation of so4
             !d_cond_so4(ii,krow)=d_cond_so4(ii,krow)+(sum(zdvolsa)-zn_vs_c*zdvolsa(in1a))*zqtmst*delta_time*rhosu*pdpg(ii,jj)/3.0_dp
             
          END IF

        END DO
     END DO


  END SUBROUTINE condensation


END MODULE mo_salsa_dynamics
