!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_salsa_properties.f90
!!
!! \brief        
!!   Contains subroutines and functions that are used           
!!   to calculate particle properties during simulation
!!
!! \author Harri Kokkola (FMI)
!!
!! \responsible_coder
!! Harri Kokkola, harri.kokkola@fmi.fi
!!
!! \revision_history
!!   -# H. Korhonen (FMI) - original code (2005)
!!   -# H. Kokkola (FMI) - original code (2006-2014)
!!   -# M. Niskanen(FMI) - change from 3 subregions to 2 subregions (2012)
!!   -# A. Laakso (FMI) - ECHAM6-HAMMOZ implementation (2013)
!!
!! \limitations
!! None
!!
!! \details
!! This module contain equilibration subroutine, which calculates ambient sizes of particles by equilibrating
!! soluble fraction of particles with water. Equilibration subroutine is called from mo_ham_salsa.
!! Parameter lists and flags are defined in mo_ham_salsactl.
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

  ! >> thk: VBS
  FUNCTION zbinmol_inv(id_vols, zaw) RESULT(res)
    ! Function to calculate inverse of the binary molalities of the 
    ! different substances the parameter epsoc has been included 
    ! here to make the calculation more portable.
    ! If VBS (nsoa == 2) is on the same function as for OC is used

    USE mo_kind, ONLY : &
         dp                      ! double precision

    USE mo_species,      ONLY: &
         speclist

    USE mo_ham, ONLY: &
         nsoa                    ! which soa algorithm is used

    USE mo_ham_species, ONLY: &
         id_so4, id_oc, id_ss, id_wat

    USE mo_ham_salsactl, ONLY : &
         epsoc           

#ifdef HAMMOZ
    USE mo_ham_vbsctl, ONLY: &
         vbs_ngroup,           & ! number of VBS bins
         vbs_set                 ! VBS
#endif
    
    IMPLICIT NONE

    INTEGER, INTENT(in) :: id_vols
    REAL(dp), INTENT(in):: zaw
    REAL(dp)            :: res !-- binary molalities [mol/kg]
 
    
    ! for easier code readability
    INTEGER, PARAMETER :: salsa_su = 1
    INTEGER, PARAMETER :: salsa_oc = 2
    INTEGER, PARAMETER :: salsa_bc = 3
    INTEGER, PARAMETER :: salsa_ss = 4
    INTEGER, PARAMETER :: salsa_du = 5

    ! for iteration
    INTEGER :: jg

    IF (id_vols == salsa_su) THEN ! sulphate
       res = 1._dp/(                    &
            + 1.1065495e+2_dp           & 
            - 3.6759197e+2_dp * zaw     &  
            + 5.0462934e+2_dp * zaw**2  &
            - 3.1543839e+2_dp * zaw**3  &
            + 6.770824e+1_dp  * zaw**4  &
       )

    ELSE IF (id_vols == salsa_oc) THEN ! organic carbon
       !epsoc: only a fraction of oc is soluble
       res = epsoc/(                                         & 
            +1._dp/(zaw*(speclist(id_wat)%moleweight*1e-3_dp)) &
            -1._dp/(speclist(id_wat)%moleweight*1e-3_dp)       &
       )

    ELSE IF (id_vols == salsa_ss) THEN ! sea salt (NaCl)
       res =  1._dp/(                   &
            + 5.875248e+1_dp            &
            - 1.8781997e+2_dp * zaw     &  
            + 2.7211377e+2_dp * zaw**2  &
            - 1.8458287e+2_dp * zaw**3  &
            + 4.153689e+1_dp  * zaw**4  &
       )
    ELSE
       res = 0.0
    END IF
    
#ifdef HAMMOZ
    !the VBS -- using same function as for OC:
    IF (nsoa == 2) THEN
       !if OC is part of the VBS, res will be set here again, but as
       !the same mathematical function is used everywhere, the effect is the same
       DO jg = 1,vbs_ngroup
          ! easier access
          IF (vbs_set(jg)%id_vols == id_vols) THEN
             !epsoc: only a fraction of oc is soluble
             res = epsoc/(                                         & 
                  +1._dp/(zaw*(speclist(id_wat)%moleweight*1e-3_dp)) &
                  -1._dp/(speclist(id_wat)%moleweight*1e-3_dp)       &
             )
          END IF
       END DO
    END IF
#endif
  END FUNCTION zbinmol_inv
  ! << thk
    

  SUBROUTINE equilibration(kproma, kbdim, klev,                    &
                           pnaero, pvols, prh, ptemp, pcore, pdwet)

    USE mo_species, ONLY :                      &
         speclist

    USE mo_ham_species, ONLY:                   &
         id_so4, id_oc, id_ss, id_wat

    USE mo_ham_salsactl, ONLY : &
         in1a, fn1a,   &
         in2a, fn2a,   &
         in2b, fn2b,   &
         nlim,         & ! lowest possible particle conc. in a bin [#/m3]
         surfw0,       & ! surface tension of water [J/m2]
         dpmid,        & ! mid dry diameter of each bin [m]
         lhydr2b,      & ! hydration calculation for insoluble region
         epsoc           ! thk: should this be made dynamic?

    USE mo_physical_constants, ONLY: ak, avo

    USE mo_math_constants, ONLY: pi

    USE mo_kind, ONLY : dp

    USE mo_ham, ONLY:&         
         sizeclass,                  & ! list of bin properties
         subm_naerospec_nowat,       & ! number of aerosol species (excluding water)
         subm_aerospec_nowat,       & ! map from salsa species to speclist
         nsoa                          ! switch to turn on VBS

#ifdef HAMMOZ
    USe mo_ham_vbsctl, ONLY:    &
         vbs_ngroup,            & ! number of VBS groups
         vbs_set,               & ! list of VBS grouops
         laqsoa,                & ! switch to turn on wet SOA
         aqsoa_ngroup,          & ! number of wet SOA groups 
         aqsoa_set                ! list of wet SOA grouops
    ! << thk
#endif

    IMPLICIT NONE

    !-- input variables -------------
    INTEGER, INTENT(in) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels 

    REAL(dp), INTENT(in) ::        &     
         pnaero(kbdim,klev,fn2b),  & ! particle concentration [#/m3]
         pvols(kbdim,klev,fn2b,subm_naerospec_nowat), & ! total volume concentrations of each
                                ! chem. compound in a size bin [fxm]
         prh(kbdim,klev),          & ! relative humidity [0-1]
         ptemp(kbdim,klev)           ! temperature [K]


    !-- output variables -------------
    REAL(dp), INTENT(out) ::       &
         pcore(kbdim,klev,fn2b),   & ! particle dry volume [fxm]
         pdwet(kbdim,klev,fn2b)      ! particle ambient diameter [m]


    ! >> thk: VBS
    !-- local variables --------------
    INTEGER :: ii, jj, kk, js        ! loop indices

    ! -- speclist indexing
    INTEGER :: id_spec
    ! << thk`
    
    REAL(dp) ::      &
         ! >> thk: VBS -- issue ???
         zbinmol(subm_naerospec_nowat), &   ! binary molality of individual components [mol/kg]
         ! << thk
         zvpart(subm_naerospec_nowat), &   ! volume of chem. compounds in one particle [fxm]         
         zke,        &   ! Kelvin term
         zaw,        &   ! water activity [0-1]         
         zlwc,       &   ! liquid water content [kg/m3-air]
         zdold,      &   !
         zrh,        &   ! 
         zmvsu           ! molar volume


    zmvsu = (speclist(id_so4)%moleweight/1000.)/avo/speclist(id_so4)%density

    !----------------------------------------------------------------------
    !-- 1) Regime 1: sulphate and partly water-soluble OC -----------------

    pdwet(1:kproma,:,:) = 0._dp
    
    ! >> thk: VBS -- issue ???
    DO kk = in1a,fn2b      ! size bin, eehol: expanded the loop from fn1a to fn2b
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid

             !-- initialize
             zke = 1.001_dp !<--eehol: Kelvin effect for all bins changed to reduce IF statements
             
             zbinmol = 0._dp
             zdold = 1._dp

             IF ((pnaero(ii,jj,kk) > nlim)) THEN

                !-- volume of sulphate and OC in one particle [fxm]

                zvpart(:) = pvols(ii,jj,kk,:)/pnaero(ii,jj,kk)

                !-- total volume of one dry particle [fxm] 
                ! ??? thk: can we just sum over all species???
                !pcore(ii,jj,kk)   = sum(zvpart(1:2))
                pcore(ii,jj,kk)   = sum(zvpart)

                ! Relative Humidity:
                zrh = prh(ii,jj)
                zrh = MAX(zrh , 0.05_dp) !eehol: RH changed to reduce IF statements
                zrh = MIN(zrh , 0.95_dp)
                
                !<--eehol: hydration calculation according to lsoluble flag
                IF (sizeclass(kk)%lsoluble) THEN
                   !DO WHILE(abs(pdwet(ii,jj,kk)/zdold-1.) > 1.e-2_dp)
                   DO WHILE(abs(pdwet(ii,jj,kk)-zdold) > 1.e-2_dp*abs(zdold)) ! avoiding one division, eehol: 1e-12 changed to 1e-2 as there is no need for such precision
                      zdold = max(pdwet(ii,jj,kk),1.e-20_dp)
                      
                      zaw = zrh/zke
                      
                      zlwc = 0.0_dp
                      DO js = 1,subm_naerospec_nowat ! looping over species
                         
                         id_spec = subm_aerospec_nowat(js)
                         ! Calculate the liquid water content (kg/m3-air) using ZSR
                         ! (see e.g. equation (9.98) in Seinfeld and Pandis (1998))
                         ! thk: for insoluble substances sol_frac = 0, so basically a +0 calculation
                         ! thk: the 1e3 converts the mole mass from g/mol to kg/mol
                         ! thk: the factor epsoc that accounts for OC being only partly soluble
                         !      is included in the function zbinmol_inv
                         zlwc = zlwc+pvols(ii,jj,kk,js)*speclist(id_spec)%density/&
                              speclist(id_spec)%moleweight*1e3_dp*zbinmol_inv(js,zaw)
                         ! 
                      END DO !js
                      
                      !-- particle wet radius [m] 
                      pdwet(ii,jj,kk) = (zlwc/pnaero(ii,jj,kk)/speclist(id_wat)%density/(pi/6.0) + &
                           pcore(ii,jj,kk)/(pi/6.0))**(1._dp/3._dp)
                      
                      !-- Kelvin effect 
                      zke = exp(2._dp*surfw0*zmvsu/(ak*ptemp(ii,jj)*pdwet(ii,jj,kk)))
                      
                   END DO
                ELSE !eehol: if not soluble then dwet = ddry
                   pdwet(ii,jj,kk) = (pcore(ii,jj,kk)/(pi/6.0))**(1._dp/3._dp)
                END IF
                !-->eehol
             ELSE
                !-- 1.2) empty bins given bin average values ----------------- 
                pdwet(ii,jj,kk) = dpmid(kk)
                pcore(ii,jj,kk) = (pi/6.0)*dpmid(kk)**3
             END IF
          END DO
       END DO
    END DO
    
  END SUBROUTINE equilibration

END MODULE mo_ham_salsa_properties
