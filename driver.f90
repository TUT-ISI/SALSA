PROGRAM driver
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
  USE driver_input
  USE mo_ham_salsa
  USE mo_ham_submctl
  USE mo_ham_salsa_init
  USE mo_ham_salsa_sizedist
  USE mo_kind
  USE mo_ham_salsa_cloud
  USE mo_time_control

  IMPLICIT NONE

  INTEGER :: ii, jj, kk

  ! variables that can be deleted upon model coupling:
  REAL(dp) :: &
       core(fn2b)

  !**********************************************************
  !*                                                        *
  !* I) Tracers/subroutines to be coupled with host model   *
  !*                                                        *
  !*  NB: When coupling, check units carefully!!!!          *
  !*                                                        *
  !**********************************************************

  !-----------------------------------------------------------------
  !-- Tracers (provided by/to host model) --------------------------
  !-----------------------------------------------------------------

  !-- aerosol tracers ----------------- 
  INTEGER, PARAMETER :: nmod = 7
  REAL(dp) :: &
       
       zaernl(kbdim,klev,fn2b), & ! number concentration of aerosol particles
                                  ! for each grid point (kbdim,klev) and each
                                  ! size bin (fn2b) [#/m3]
       v_aero(kbdim,klev,45),& ! for output layout...
       zaerml(kbdim,klev,2*fn1a+5*(fn2a-in2a+1)+4*(fn2b-in2b+1)), &
       vols(kbdim,klev,fn2b,5),&  ! volume concentration of chem. compounds
       zm6rp(kbdim,klev,fn2b), &  ! dumb variables included for ECHAM
       zm6dry(kbdim,klev,fn2b),&  !
       zrhop(kbdim,klev,fn2b), &  ! 
       zww(kbdim,klev,fn2b),   &  !
       zw(kbdim,klev),          & ! updraft velocity (m/s)
       zcd(kbdim,klev),         & ! number of activated droplets [#/m3]
       n(nmod),sigmag(nmod),dpg(nmod) 


  !-- gas compound tracers ------------
  REAL(dp) :: &
       
       zgso4(kbdim,klev), & ! sulphuric acid concentration in gas phase
                                ! for each grid point (kbdim,klev) [#/m3]
       
       zgocnv(kbdim,klev),  & ! non-volatile organic vapour concentration
                                ! for each grid point (kbdim,klev) [#/m3]
       
       zgocsv(kbdim,klev)     ! semivolatile organic vapour concentration
                                ! for each grid point (kbdim,klev) [#/m3]


  !-- atmospheric conditions --------------
  REAL(dp) :: &
       
       pap(kbdim,klev),  & ! atmospheric pressure at time t+1 
                                ! for each grid point (kbdim,klev) [Pa]
       
       zrh(kbdim,klev), & ! atmospheric relative humidity at time t+1 
                                ! for each grid point (kbdim,klev) [0-1]
       
       pt(kbdim,klev),  &       ! atmospheric temperature at time t+1 
                                ! for each grid point (kbdim,klev) [K]
       zpbl(kbdim)              ! boundary layer height

  !-- local latitude index
  INTEGER :: krow

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! FOR TESTING PURPOSE ONLY!!!
  !
  REAL(dp) :: timein, timeout, ztmst
  INTEGER :: i

  OPEN(14,FILE='output.dat',STATUS='unknown')

  ! Set values to 'ambient tracers'. 
  ! NB! Later these variables provided by host model

  ztmst=time_step_len
  zpbl=1
  ! Ambient pressure (Pa)
  pap = 101325.
  ! Ambient temperature (K)
  pt = 298._dp
  do i = 1,kbdim
     ! Saturation ratio of gas phase water
     zrh(i,1:klev) = 0.3_dp!(/0.1, 0.2, 0.3, 0.4, 0.5/)!0.6! [fxm]
     ! Updraft velocity (m/s)
     zw(i,1:klev) = 1._dp
  end do
  
  zgso4 = 0._dp!5.E14_dp
  zgocnv = 0._dp!5.E14_dp
  zgocsv = 0._dp!1.E14_dp
  vols = 0._dp

  !--------------------------------------------------------------------------------
  !-- Routines that must be called before aerosol calculations started ------------
  !--------------------------------------------------------------------------------
  !
  !  Set size bin spacing for particles
  !  The output (i.e. low and high boundary limits of bins, 
  !  bin mid-diameter etc.) is stored in global variables in mod_fxm
  ! 
  !  NB: This must be done somewhere in the host model -
  !  but only for one time i.e. before any aerosol calculations are started
  !
  
  call set_sizebins()

  call actcurve()

  !--------------------------------------------------------------------------------
  !
  !  Calculate coagulation coefficients for particles:
  !  different set for each vertical (pressure) level
  !
  !  The values are calculated for bin mid-diameters and scaled each
  !  time step according to actual particle wet size 
  !
  !  NB: This must be done somewhere in the host model -
  !  but only for one time i.e. before any aerosol calculations are started
  !

  call set_coagc(klev,pap,pt)

 !***********************************************************
  !
  !  FOR TESTING PURPOSE ONLY
  !
  !***********************************************************
  !              
  ! II) Initial particle size distribution, seven log-normal modes
  !              

  ! NB! This only for easier calculation of initial volumes

  core(in1a:fn2b) = pi6 * dpmid(in1a:fn2b)**3


  ! Stdev of modes
  sigmag = (/ 1.59, 1.59, 1.59, 1.59, 2.0, 2.0, 2.0 /)
  ! Mean diameter of the modes (m)
  dpg = (/0.03, 0.2, 0.2, 0.2, 6.0, 6.0, 6.0/)*1.e-6_dp
  ! Number concentration of modes (#/cm3)
  n = (/1000.0, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0/)*1.e6_dp

  CALL size_distribution(kproma, kbdim,  klev,   &
       n, dpg, sigmag, zaernl)

  ! Volume concentrations (m3/m3)

  DO jj = 1, klev
     DO ii = 1, kproma


        vols(ii,jj,in1a:fn1a,1) = 0.2_dp*zaernl(ii,jj,in1a:fn1a)*core(in1a:fn1a)
        vols(ii,jj,in1a:fn1a,2) = 0.8_dp*zaernl(ii,jj,in1a:fn1a)*core(in1a:fn1a)

        vols(ii,jj,in2a:fn2a,1) = 0.5_dp*zaernl(ii,jj,in2a:fn2a)*core(in2a:fn2a)
        vols(ii,jj,in2a:fn2a,2) = 0.05_dp*zaernl(ii,jj,in2a:fn2a)*core(in2a:fn2a)
        vols(ii,jj,in2a:fn2a,3) = 0.05_dp*zaernl(ii,jj,in2a:fn2a)*core(in2a:fn2a)
        vols(ii,jj,in2a:fn2a,4) = 0.1_dp*zaernl(ii,jj,in2a:fn2a)*core(in2a:fn2a)
        vols(ii,jj,in2a:fn2a,5) = 0.3_dp*zaernl(ii,jj,in2a:fn2a)*core(in2a:fn2a)

        vols(ii,jj,in2b:fn2b,1) = 0.01*zaernl(ii,jj,in2b:fn2b)*core(in2b:fn2b)
        vols(ii,jj,in2b:fn2b,2) = 0.01_dp*zaernl(ii,jj,in2b:fn2b)*core(in2b:fn2b)
        vols(ii,jj,in2b:fn2b,3) = 0.01_dp*zaernl(ii,jj,in2b:fn2b)*core(in2b:fn2b)
        vols(ii,jj,in2b:fn2b,5) = 0.97_dp*zaernl(ii,jj,in2b:fn2b)*core(in2b:fn2b)

     END DO
  END DO

  ! Set up mass concentrations (for compatibility with ECHAM6)

  jj = 0

    DO ii = in1a, fn1a

       jj = jj + 1

       iso4b(ii) = jj
       !--- Sulfate mass
       zaerml(1:kproma,:,iso4b(ii)) = vols(1:kproma,:,ii,1) /(1.e6_dp *mvsu)

    END DO

    DO ii = in1a, fn1a

       jj = jj + 1

       iocb(ii) = jj
       !--- Organic carbon
       zaerml(1:kproma,:,iocb(ii))  = vols(1:kproma,:,ii,2) * rhooc / 1.e-9_dp

    END DO

    DO ii = in2a, fn2a

       jj = jj + 1

       iso4b(ii) = jj
       !--- Sulfate volume
       zaerml(1:kproma,:,iso4b(ii)) = vols(1:kproma,:,ii,1) / (1.e6_dp *mvsu)

    END DO

    DO ii = in2a, fn2a

       jj = jj + 1

       iocb(ii) = jj
       !--- Organic carbon
       zaerml(1:kproma,:,iocb(ii)) =  vols(1:kproma,:,ii,2) * rhooc / 1.e-9_dp

    END DO

    DO ii = in2a, fn2a

       jj = jj + 1

       ibcb(ii) = jj
       !--- Black carbon
       zaerml(1:kproma,:,ibcb(ii)) =  vols(1:kproma,:,ii,3) * rhobc / 1.e-9_dp

    END DO

    DO ii = in2a, fn2a

       jj = jj + 1

       issb(ii) = jj
       !--- Sea salt
       zaerml(1:kproma,:,issb(ii)) =  vols(1:kproma,:,ii,4) * rhoss / 1.e-9_dp

    END DO

    DO ii = in2a, fn2a

       jj = jj + 1

       idub(ii) = jj
       !--- Mineral dust
       zaerml(1:kproma,:,idub(ii)) =  vols(1:kproma,:,ii,5) * rhodu / 1.e-9_dp

    END DO

    DO ii = in2b, fn2b

       jj = jj + 1

       iso4b(ii) = jj

       !--- Sulfate volume
       zaerml(1:kproma,:,iso4b(ii)) = vols(1:kproma,:,ii,1) / (1.e6_dp *mvsu)
       
    END DO

    DO ii = in2b, fn2b

       jj = jj + 1

       iocb(ii) = jj

       !--- Organic carbon
       zaerml(1:kproma,:,iocb(ii))  = vols(1:kproma,:,ii,2) * rhooc / 1.e-9_dp
    END DO

    DO ii = in2b, fn2b

       jj = jj + 1

       ibcb(ii) = jj

       !--- Black carbon
       zaerml(1:kproma,:,ibcb(ii))  = vols(1:kproma,:,ii,3) * rhobc / 1.e-9_dp

    END DO

    DO ii = in2b, fn2b

       jj = jj + 1

       idub(ii) = jj

       !--- Mineral dust
       zaerml(1:kproma,:,idub(ii))  = vols(1:kproma,:,ii,5) * rhodu / 1.e-9_dp

    END DO

  !*************************************************
  !*                                               *
  !*  III) REAL STUFF BEGINS HERE !!!!             *
  !*                                               *
  !*************************************************

  CALL CPU_TIME(timein)
  
  DO ii = 1, 5000
 
     IF(ii > 200 .AND. ii < 1000) THEN
        zgso4 = 5.E14_dp
     END IF

     zgocnv=0._dp
     zgocsv=0._dp

     CALL salsa(kproma,  kbdim,   klev,    krow, &  ! ECHAM indices
                pap,     zrh,     pt,      ztmst,&  ! Pressure, RH, temperature, time step length
                zgso4,   zgocnv,  zgocsv,        &  ! [H2SO4(g)], [OCNV(g)], [OCSV(g)]
                zaerml,  zaernl,                 &  ! Aerosol volume and number
                zm6rp,   zm6dry,  zrhop,   zww,  &  ! Aerosol properties
                zpbl                              ) ! Planetary boundary layer top level

!!$     CALL salsa(kproma, kbdim,  klev, krow,            &
!!$          papp1, prelhum, ptp1, time_step_len,         &
!!$          c_h2so4, c_ocnv,  c_ocsv,                    &
!!$          vols, zaernl, zm6rp,  zm6dry, zrhop,  zww, &
!!$          zw, zcd, ppbl)

    DO jj = in1a, fn1a
       !--- Sulfate volume
       vols(1:kproma,:,jj,1) = zaerml(1:kproma,:,iso4b(jj)) * 1.e6_dp *mvsu
       !--- Organic carbon
       vols(1:kproma,:,jj,2) = zaerml(1:kproma,:,iocb(jj))/rhooc * 1.e-9_dp
    END DO

    DO jj = in2a, fn2a
       !--- Sulfate volume
       vols(1:kproma,:,jj,1) = zaerml(1:kproma,:,iso4b(jj)) * 1.e6_dp *mvsu
       !--- Organic carbon
       vols(1:kproma,:,jj,2) = zaerml(1:kproma,:,iocb(jj))/rhooc * 1.e-9_dp
       !--- Black carbon
       vols(1:kproma,:,jj,3) = zaerml(1:kproma,:,ibcb(jj))/rhobc * 1.e-9_dp
       !--- Sea salt
       vols(1:kproma,:,jj,4) = zaerml(1:kproma,:,issb(jj))/rhoss * 1.e-9_dp
       !--- Mineral dust
       vols(1:kproma,:,jj,5) = zaerml(1:kproma,:,idub(jj))/rhodu * 1.e-9_dp
    END DO

    DO jj = in2b, fn2b
       !--- Sulfate volume
       vols(1:kproma,:,jj,1) = zaerml(1:kproma,:,iso4b(jj)) * 1.e6_dp *mvsu
       !--- Organic carbon
       vols(1:kproma,:,jj,2) = zaerml(1:kproma,:,iocb(jj))/rhooc * 1.e-9_dp
       !--- Black carbon
       vols(1:kproma,:,jj,3) = zaerml(1:kproma,:,ibcb(jj))/rhobc * 1.e-9_dp
       !--- Mineral dust
       vols(1:kproma,:,jj,5) = zaerml(1:kproma,:,idub(jj))/rhodu * 1.e-9_dp
    END DO

     CALL cloud_activation(kproma, kbdim, klev, &
          zaernl, vols, pt, pap, zcd, zw)
! Changing to the old v_aero layout:

     v_aero = 0._dp
     write(14,665) zm6rp(1,1,in1a:fn2a), &
          zaernl(1,1,in1a:fn2a)
  END DO
  
665 FORMAT(99(E11.4,1X))

  CALL CPU_TIME(timeout)

END PROGRAM driver

