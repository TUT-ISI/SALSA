!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_submctl.f90
!!
!! \brief
!! mo_ham_submctl contains parameters, switches and initialization routines for the m7 aerosol module.
!!
!! \author Elisabetta Vignatti (JRC/EI)
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# E. Vignati and J. Wilson (JRC/EI) - original code (2000)
!!   -# P. Stier (MPI-Met) (2001/2002)
!!   -# J. Kazil (MPI-Met) (2008)
!!   -# D. O'Donnell (MPI-Met) (2007-2007)
!!   -# M.G. Schultz (FZ Juelich) - new module struture (2009)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!!   -# A. Laakso (FMI) - salsa_initialize (2013-02-25)
!! 
!! \limitations
!! Currently, there are two index lists for aerosol species: aero_idx in mo_species
!! and subm_aerospec in this module. I hope these are identical for the current model set-up 
!! in preparation for CMIP5. Later, one may wish to distinguish between the two: aero_idx
!! could contain additional aerosol species (e.g. from MOZART or climatologies), and this could
!! mess up the M7 code. If this can be generalized: fine. if not we should keep the two 
!! lists separate. mo_ham_rad (for example) works on aero_idx to be independent of M7 specifics.
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

MODULE mo_ham_submctl

  USE mo_kind,             ONLY: dp

  IMPLICIT NONE

  PRIVATE

  !M7 and SALSA
  PUBLIC :: lscoag, lscond, nsnucl
  !SALSA:
  PUBLIC :: locgas, lsol2b, act_coeff,nj3

  PUBLIC :: in1a,in2a,in2b,fn1a,fn2a,fn2b,nbins
  PUBLIC :: nbin, nbin2, nbin3,reglim,nlim,nreg
  PUBLIC :: epsv,vhilim,vlolim,vratiohi,vratiolo,dpmid,sigma
  PUBLIC :: pi, pi6, rg, avog, boltz, cpa, mair
  PUBLIC :: rhosu,rhooc, rhobc,rhoss, rhodu, rhowa  
  PUBLIC :: msu,mdu,n3,massacc,d_sa,pstand,mss,mbc,moc,epsoc,mwa,slim,ions,mvsu,mvoc,mvss,surfw0
  PUBLIC :: recalc

  !--- 1) Define and pre-set switches for the processes of M7: -----------------------

  !--- Physical:
  !Switches for both M7 and SALSA aerosol microphysical processes
  LOGICAL :: lscoag     = .TRUE.    ! Coagulation
  LOGICAL :: lscond     = .TRUE.    ! Condensation of H2SO4
  
  
  ! 1) Switches for M7 aerosol microphysical processes ------------------------
  INTEGER :: nwater     = 1         ! Aerosol water uptake scheme:
                                    !
                                    ! nwater = 0 Jacobson et al., JGR 1996
                                    !        = 1 Kappa-Koehler theory based approach (Petters and Kreidenweis, ACP 2007)

  INTEGER :: nsnucl     = 2         ! Choice of the H2SO4/H2O nucleation scheme:
                                    ! M7:
                                    ! nsnucl = 0 off
                                    !        = 1 Vehkamaeki et al., JGR 2002
                                    !        = 2 Kazil and Lovejoy, ACP 2007
                                    ! SALSA:
                                    ! 0 = off   
                                    ! 1 = binary nucleation
                                    ! 2 = activation type nucleation
                                    ! 3 = kinetic nucleation
                                    ! 4 = ternary nucleation
                                    ! 5 = nucleation with ORGANICs
                                    ! 6 = activation type of nucleation with H2SO4+ORG
                                    ! 7 = heteromolecular nucleation with H2SO4*ORG
                                    ! 8 = homomolecular nucleation of  H2SO4 + 
                                    !           heteromolecular nucleation with H2SO4*ORG
                                    ! 9 = homomolecular nucleation of  H2SO4 and ORG + 
                                    !           heteromolecular nucleation with H2SO4*ORG

  INTEGER :: nonucl     = 1         ! Choice of the organic nucleation scheme:
                                    ! 
                                    ! nonucl = 0 off
                                    !        = 1 Activation nucleation, Kulmala et al., ACP 2006
                                    !        = 2 Activation nucleation, Laakso et al., ACP 2004
  
  LOGICAL :: lgcr       = .TRUE.    ! Calculate ionization due to galactic cosmic rays
  
  REAL(dp):: nsolact    = -99.99_dp ! Solar activity parameter [-1,1]; if outside of
                                    ! this range (as per default), then the model will
                                    ! determine the solar activity based on the model
                                    ! calendar date; otherwise, it will use the user
                                    ! set solar activity parameter throughout the run.
                                    ! -1 is solar minimum, 1 solar maximum.
  
  ! 1) Switches for SALSA aerosol microphysical processes ------------------------ 

  LOGICAL :: locgas = .FALSE.,&   ! emission of organic carbon in gas phase
             lsol2b = .FALSE.     ! repartitioning of insoluble material in 
                                  ! case of increase in solubility 

  LOGICAL :: recalc   = .FALSE.   ! recalculation of wet diameter between
                                  ! calculation of microphysical processes

  INTEGER ::                    & ! J3 parametrization
             nj3 = 1              ! 1 = condensational sink (Kerminen&Kulmala, 2002)
                                  ! 2 = coagulational sink (Lehtinen et al. 2007)
                                  ! 3 = coagS+self-coagulation (Anttila et al. 2010)
  REAL(dp) :: act_coeff=1.e-7_dp  ! activation coefficient

  !--Indices corresponding to number concentration, radius
  !    and chemical components in each subregime

  INTEGER, PARAMETER ::            &
   nreg = 2                        ! number of main size regimes

  REAL(dp), PARAMETER ::                       &
   reglim(nreg+2) =                            & ! low/high diameter limits
    (/ 3.e-9_dp, 5.e-8_dp, 7.e-7_dp, 1.e-5_dp /) ! of main size regimes [m]

   INTEGER, PARAMETER :: &
   nbin(nreg)   = (/ 3, 7 /)   ! number of bins in each main regime

  INTEGER, PARAMETER ::      &
   nbin2 = 4,                & ! number of bins in former 2-region
   nbin3 = nbin(2) - nbin2     ! number of bins in former 3-region

  INTEGER, PARAMETER ::      & ! number/radius: start index
   in1a = 1,                 & ! regime 1a
   in2a = in1a + nbin(1),    & ! regime 2a
   in2b = in2a + nbin(2),    & ! regime 2b

                               ! number/radius: last index
   fn1a = in2a - 1,          & ! regime 1a
   fn2a = fn1a + nbin(2),    & ! regime 2a
   fn2b = fn2a + nbin(2),    & ! regime 2b

   nbins = fn2b                ! total number of size bins


  REAL(dp) ::               &
       epsv(nbins),         &
       vhilim(fn2b),        &
       vlolim(fn2b),        &
       vratiohi(fn2b),      &
       vratiolo(fn2b),      &
       dpmid(fn2b),         &
       sigma(fn2b),         &
       csr_strat_wat(fn2b), &
       csr_strat_mix(fn2b), &
       csr_strat_ice(fn2b), &
       csr_conv(fn2b),      &
       zbcr(fn2b)

  REAL(dp), PARAMETER ::     &
   avog   = 6.0221e+23_dp,   & ! Avogadro number (#/mol)
   boltz  = 1.3807e-23_dp,   & ! Boltzmann constant (J/K)
   grav   = 9.81_dp,         & ! gravitational acceleration (m/s^2)
   pstand = 1.01325e+5_dp,   & ! standard pressure (Pa)
   rg     = 8.314_dp,        & ! molar gas constant (J/(mol K)) 
   pi     = 3.1415927_dp,    & ! self explanatory
   pi6    = 0.5235988_dp,    & ! pi/6
   cpa    = 1010._dp,        & ! specific heat of dry air, constant P (J/kg/K)
   mair   = 28.97e-3_dp,     & ! molar mass of air (mol/kg)
   deltav = 1.096e-7_dp,     & ! vapor jump length (m)
   deltaT = 2.16e-7_dp,      & ! thermal jump length (m)
   alphaT = 0.96_dp,         & ! thermal accomodation coefficient
   alphac = 1.0_dp             ! condensation coefficient

  REAL(dp), PARAMETER ::     & ! molar mass [kg/mol]
   msu = 98.08e-3_dp,        & ! sulphate   
   moc = 150.e-3_dp,         & ! organic carbon
   mbc = 12.e-3_dp,          & ! black carbon
   mss = 58.44e-3_dp,        & ! sea salt (NaCl)
   mdu = 100.e-3_dp,         & ! mineral dust
   mwa = 18.016e-3_dp,       & ! water
                               !
                               ! densities [kg/m3]
   rhosu = 1830._dp,         & ! sulphate
   rhooc = 2000._dp,         & ! organic carbon
   rhobc = 2000._dp,         & ! black carbon
   rhoss = 2165._dp,         & ! sea salt (NaCl)
   rhodu = 2650._dp,         & ! mineral dust
   rhowa = 1000._dp,         & ! water 
                               !
                               ! volume of molecule [kg/#]
   mvsu = msu/avog/rhosu,    & ! sulphate
   mvoc = moc/avog/rhooc,    & ! organic carbon
   mvss = mss/avog/rhoss,    & ! sea salt
                               !
   volratio =                & ! ratio of molecular volumes for
    (msu*rhoss)/(rhosu*mss), & ! sulphate and sea salt
                               !
   n3 = 158.79_dp               ! number of H2SO4 molecules in 3 nm cluster 
                               !  assuming d_sa = 5.54 ???     
  !-- 4.3) Properties of condensing vapours

  REAL(dp), PARAMETER ::                               & ! diameter of condensing molecule [m]
      d_sa   = 5.539376964394570e-10_dp,               &

      d_oc   = 6.195906936656752e-10_dp,               &

      d_h2o  = 3.851565216195334e-10_dp

  REAL(dp), PARAMETER :: &
       slim = 1.005_dp,  & ! water saturation used as limit
       ions = 3.0_dp,    & ! van't Hoff factor (ions produced upon dissociation)
       surfw0 = 0.073_dp,& ! surface tension of pure water @ ~ 293 K [J/m2]
       epsoc = 0.15_dp     ! water uptake of organic material

  !-- 7) Parameters for cloud activation

  REAL(dp), PARAMETER :: crcut=0.035*1E-6_dp ! Assumed lower cut-off of the
                                             ! aerosol size distribution [m]

  !--- Ulrike: included for activation in convective clouds
  REAL(dp), PARAMETER :: crcut_cv=0.025*1E-6_dp ! Assumed lower cut-off of the
  
  
  REAL(dp) :: cfracn(fn2b)
  
  REAL(dp) :: zfracn(fn2b),   &
              zfracn_cv(fn2b)
  REAL(dp), PARAMETER :: &
   massacc(nbins) = 1._dp

  REAL(dp), PARAMETER :: &
   nlim = 1._dp,         & ! number conc. limit below which bin empty  [#/m3] 
   m3_2_um3 = 1.e+18_dp    ! conversion factor for volume from m3 to um3

  INTEGER, PUBLIC :: iso4b(fn2b), iocb(fn2b), ibcb(fn2b), idub(fn2b), issb(fn2b)
  !INTEGER(fn2b), PUBLIC :: iso4b, iocb
  !INTEGER(fn2b-fn1a), PUBLIC :: ibcb, idub
  !INTEGER(fn2a-fn1a), PUBLIC :: issb
 
  !--- 12) Service routines for initialization and auxiliary computations ----------

END MODULE mo_ham_submctl
