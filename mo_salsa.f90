MODULE mo_salsa
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
  
  PRIVATE

  ! -- subroutines
  PUBLIC :: salsa

CONTAINS

  SUBROUTINE salsa(kproma,   kbdim,   klev,    krow,                &
       ppres,    prh,     ptemp,   ptstep,              &
       pc_h2so4, pc_ocnv, pc_ocsv, paerml,  pnaero,    &
       pm6rp,    pm6dry,  prhop,   pww,    ppbl)

    USE mo_salsa_properties
    USE mo_salsa_dynamics
    USE mo_salsa_nucleation
    USE mo_salsa_update
    USE mo_salsa_init

    USE mo_submctl, ONLY :      &
         pi6,                       & ! pi/6
         avog,                      &
         in1a,  in2a,               & ! size section and composition indices
         in2b,  fn1a,               &
         fn2a,  fn2b,               &
         lscond,                    & ! switches for aerosol processes
         lscoag,                    &
         dpmid,                     &
         rhosu, msu, mvsu,          & ! properties of compounds
         rhooc, moc,                &
         rhobc, mbc,                &
         rhoss, mss,                &
         rhowa,                     &
         rhodu, mdu,                &
         nlim,                      &
         recalc                       ! logical switch for recalculation of zdwet

    USE mo_submctl,only:&
         iso4b,&
         iocb,&
         ibcb,&
         idub,&
         issb

    USE mo_salsa_trac, ONLY:    &
         idt_ms4,                   &
         idt_moc,                   &
         idt_mbc,                   &
         idt_mss,                   &
         idt_mdu

    USE mo_kind, ONLY : dp


    IMPLICIT NONE

    INTEGER, PARAMETER  :: naerocomp = 10
    !-- Input parameters and variables --------------------------------------------------
    INTEGER, INTENT(in) ::          &
         kproma,                    & ! number of horiz. grid points in a slab (='kproma')
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev,                      & ! number of vertical levels (='klev')
         krow                         ! local latitude index


    REAL(dp), INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         prh  (kbdim,klev),            & ! relative humidity at each grid point [0-1]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         ptstep                          ! time step [s]


    !-- Input variables that are changed within --------------------------------------
    REAL(dp), INTENT(inout) ::      & ! gas phase concentrations at each grid point [#/m3]
         paerml(kbdim,klev,naerocomp), & ! aerosol mass for individual compounds [molec. cm-3 for sulfate and ug m-3 for bc, oc, ss, and dust] 
         pc_h2so4(kbdim,klev),      & ! sulphuric acid
         pc_ocnv (kbdim,klev),      & ! nonvolatile organic compound
         pc_ocsv (kbdim,klev),      & ! semivolatile organic compound
         
                                ! aerosol distribution properties at each grid point
         pnaero(kbdim,klev,fn2b)   ,& ! number concentration in each size bin [#/m3]
         pm6rp (kbdim,klev,fn2b)   ,& ! mean mode actual radius (wet radius for soluble modes
                                ! and dry radius for insoluble modes) [cm]
         pm6dry(kbdim,klev,fn2a)   ,& ! dry radius [cm]
         pww   (kbdim,klev,fn2b)   ,& ! aerosol water content for each mode [kg(water) m-3(air)]
         prhop (kbdim,klev,fn2b)   ,& ! mean mode particle density [g cm-3]
         ppbl(kbdim)
    

    INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level

    !-- Output variables -----------------------------------------------------------------

    !   in each bin

    !-- Local variables ------------------------------------------------------------------

    LOGICAL::moving_center
    INTEGER :: ii, jj, kk, nn

    REAL(dp) ::                       &
         zvols (kbdim,klev,fn2b,5) ,  & ! volume concentrations
         zcore   (kbdim,klev,fn2b),   & ! volume of the core in one size bin
         zdwet   (kbdim,klev,fn2b),   &
         zvq     (kbdim,klev,fn2b),   &
         zddry

    zpbl(:) = int(ppbl(:))

    moving_center=.false.

    zcore = 0._dp
    zddry = 0._dp
    zvols(1:kproma,:,:,:)=0.0_dp

    !>> convert mass concentrations to volume concentrations for SALSA

    DO ii = in1a, fn1a
       !--- Sulfate volume
       zvols(1:kproma,:,ii,1) = paerml(1:kproma,:,iso4b(ii)) * 1.e6_dp *mvsu
       !--- Organic carbon
       zvols(1:kproma,:,ii,2) = paerml(1:kproma,:,iocb(ii))/rhooc * 1.e-9_dp
    END DO

    DO ii = in2a, fn2a
       !--- Sulfate volume
       zvols(1:kproma,:,ii,1) = paerml(1:kproma,:,iso4b(ii)) * 1.e6_dp *mvsu
       !--- Organic carbon
       zvols(1:kproma,:,ii,2) = paerml(1:kproma,:,iocb(ii))/rhooc * 1.e-9_dp
       !--- Black carbon
       zvols(1:kproma,:,ii,3) = paerml(1:kproma,:,ibcb(ii))/rhobc * 1.e-9_dp
       !--- Sea salt
       zvols(1:kproma,:,ii,4) = paerml(1:kproma,:,issb(ii))/rhoss * 1.e-9_dp
       !--- Mineral dust
       zvols(1:kproma,:,ii,5) = paerml(1:kproma,:,idub(ii))/rhodu * 1.e-9_dp
    END DO

    DO ii = in2b, fn2b
       !--- Sulfate volume
       zvols(1:kproma,:,ii,1) = paerml(1:kproma,:,iso4b(ii)) * 1.e6_dp *mvsu
       !--- Organic carbon
       zvols(1:kproma,:,ii,2) = paerml(1:kproma,:,iocb(ii))/rhooc * 1.e-9_dp
       !--- Black carbon
       zvols(1:kproma,:,ii,3) = paerml(1:kproma,:,ibcb(ii))/rhobc * 1.e-9_dp
       !--- Mineral dust
       zvols(1:kproma,:,ii,5) = paerml(1:kproma,:,idub(ii))/rhodu * 1.e-9_dp
    END DO

    DO kk = in1a,fn1a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(pnaero(ii,jj,kk) < nlim) CYCLE
             zddry = (sum(zvols(ii,jj,kk,1:2))/pnaero(ii,jj,kk)/pi6)**(1._dp/3._dp)
             IF(zddry < 1.e-10_dp) THEN
                pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + zvols(ii,jj,kk,1) * rhosu / msu * avog
                pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + zvols(ii,jj,kk,2) * rhooc / moc * avog
                pnaero(ii,jj,kk)  = 0.0_dp
                zvols(ii,jj,kk,:) = 0.0_dp
             END IF
          END DO
       END DO
    END DO

    DO kk = in2a,fn2a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(pnaero(ii,jj,kk) < nlim) CYCLE
             zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi6)**(1._dp/3._dp)
             IF(zddry < 1.e-10_dp) THEN
                pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + zvols(ii,jj,kk,1) * rhosu / msu * avog
                pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + zvols(ii,jj,kk,2) * rhooc / moc * avog
                pnaero(ii,jj,kk)  = 0.0_dp
                zvols(ii,jj,kk,:) = 0.0_dp
             END IF
          END DO
       END DO
    END DO

    DO kk = in2b,fn2b      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(pnaero(ii,jj,kk) < nlim) CYCLE
             zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi6)**(1._dp/3._dp)
             IF(zddry < 1.e-10_dp) THEN
                pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + zvols(ii,jj,kk,1) * rhosu / msu * avog
                pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + zvols(ii,jj,kk,2) * rhooc / moc * avog
                pnaero(ii,jj,kk)  = 0.0_dp
                zvols(ii,jj,kk,:) = 0.0_dp
             END IF
          END DO
       END DO
    END DO

    CALL equilibration(kproma, kbdim, klev,        &
         pnaero, zvols, prh, ptemp,  &
         zcore,  zdwet) 

    IF (lscoag) THEN

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
       CALL coagulation(kproma, kbdim,  klev,        &
            pnaero, zvols,  zdwet,       &
            zcore,  ptstep, ptemp, ppres)

    END IF


    !- For more accurate condensation sink, wet diameter must be recalculated
    !  NOTE: this is much more time consuming
    IF (recalc) CALL equilibration(kproma, kbdim, klev,       &
         pnaero, zvols, prh, ptemp, &
         zcore,  zdwet) 

    IF (lscond) CALL condensation(kproma,   kbdim,   klev,    krow,               &
         pnaero,   zvols,   zdwet,                       &
         pc_h2so4, pc_ocnv, pc_ocsv,                     &
         prh,      ptemp,   ppres,   ptstep, zpbl)

    DO nn = 1, 2

       CALL distr_update(kproma, kbdim, klev, &
                         pnaero, zvols)

       DO kk = in1a,fn1a      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid
                IF(pnaero(ii,jj,kk) < nlim) CYCLE
                zddry = (sum(zvols(ii,jj,kk,1:2))/pnaero(ii,jj,kk)/pi6)**(1._dp/3._dp)
                IF(zddry < 1.e-10_dp) THEN
                   pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + zvols(ii,jj,kk,1) * rhosu / msu * avog
                   pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + zvols(ii,jj,kk,2) * rhooc / moc * avog
                   pnaero(ii,jj,kk)  = 0.0_dp
                   zvols(ii,jj,kk,:) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       DO kk = in2a,fn2a      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid
                IF(pnaero(ii,jj,kk) < nlim) CYCLE
                zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi6)**(1._dp/3._dp)
                IF(zddry < 1.e-10_dp) THEN
                   pc_h2so4(ii,jj)     = pc_h2so4(ii,jj) + zvols(ii,jj,kk,1) * rhosu / msu * avog
                   pc_ocsv(ii,jj)      = pc_ocsv(ii,jj) + zvols(ii,jj,kk,2) * rhooc / moc * avog
                   pnaero(ii,jj,kk)    = 0.0_dp
                   zvols(ii,jj,kk,1:5) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       DO kk = in2b,fn2b      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid
                IF(pnaero(ii,jj,kk) < nlim) CYCLE
                zddry = (sum(zvols(ii,jj,kk,:))/pnaero(ii,jj,kk)/pi6)**(1._dp/3._dp)
                IF(zddry < 1.e-10_dp) THEN
                   pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + zvols(ii,jj,kk,1) * rhosu / msu * avog
                   pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + zvols(ii,jj,kk,2) * rhooc / moc * avog
                   pnaero(ii,jj,kk)  = 0.0_dp
                   zvols(ii,jj,kk,:) = 0.0_dp
                END IF
             END DO
          END DO
       END DO

       CALL equilibration(kproma, kbdim, klev,        &
            pnaero, zvols, prh, ptemp,  &
            zcore,  zdwet) 

       pm6dry = (zcore(:,:,1:fn2a)/pi6)**(1._dp/3._dp)/2._dp * 100._dp
       pm6rp  = zdwet/2._dp * 100._dp
       pww    = (pi6*zdwet**3 - zcore)*rhowa*pnaero 

    END DO

    !---------------------------------------------------------------


    DO jj = 1,klev      ! vertical grid
       DO ii = 1,kproma ! horizontal kproma in the slab
          DO kk = in1a, fn1a
             zvq(ii,jj,kk) = (zvols(ii,jj,kk,1) + zvols(ii,jj,kk,2) + &
                  pww(ii,jj,kk)/ rhowa)
             
             ! check if enough aerosol mass present
             
             IF(zvq(ii,jj,kk) > 1.e-30_dp) THEN
                prhop(ii,jj,kk) = (zvols(ii,jj,kk,1)*rhosu            &
                     + zvols(ii,jj,kk,2)*rhooc + pww(ii,jj,kk)) /     &
                     zvq(ii,jj,kk)                                    &
                     ! conversion from kg/m3 to g/cm3
                     /1000._dp
             ELSE
                
                ! if not enough mass, assume density of water in g/cm3 to avoid NaN
                prhop(ii,jj,kk) = rhowa/1000._dp

             END IF

          END DO

          DO kk = in2a, fn2a
             zvq(ii,jj,kk) = (zvols(ii,jj,kk,1) +                     &
                              zvols(ii,jj,kk,2) +                     &
                              zvols(ii,jj,kk,3) +                     &
                              zvols(ii,jj,kk,4) +                     &
                              zvols(ii,jj,kk,5) +                     &
                              pww(ii,jj,kk)/rhowa)

             IF(zvq(ii,jj,kk) > 1.e-30_dp) THEN 

                prhop(ii,jj,kk) = (zvols(ii,jj,kk,1)*rhosu +          &
                                   zvols(ii,jj,kk,2)*rhooc +          &
                                   zvols(ii,jj,kk,3)*rhobc +          &
                                   zvols(ii,jj,kk,4)*rhoss +          &
                                   zvols(ii,jj,kk,5)*rhodu +          &
                                   pww(ii,jj,kk))/                    &
                                   zvq(ii,jj,kk) / 1000._dp

             ELSE

                prhop(ii,jj,kk) = rhowa/1000._dp

             END IF

          END DO

          DO kk = in2b, fn2b 

             zvq(ii,jj,kk) = (zvols(ii,jj,kk,1) +                      &
                              zvols(ii,jj,kk,2) +                      &
                              zvols(ii,jj,kk,3) +                      &
                              zvols(ii,jj,kk,5) +                      &
                              pww(ii,jj,kk) / rhowa)

             IF(zvq(ii,jj,kk) > 1.e-30_dp) THEN
                prhop(ii,jj,kk) = (zvols(ii,jj,kk,1)*rhosu +           &
                                   zvols(ii,jj,kk,2)*rhooc +           &
                                   zvols(ii,jj,kk,3)*rhobc +           &
                                   zvols(ii,jj,kk,5)*rhodu +           &
                                   pww(ii,jj,kk))/                     &
                                   zvq(ii,jj,kk)/1000._dp
             ELSE
                prhop(ii,jj,kk) = rhobc/1000._dp
             END IF

          END DO
       END DO
    END DO

    !>> convert volume concentrations to mass concentrations for HAM

    DO ii = in1a, fn1a
       !--- Sulfate volume
       paerml(1:kproma,:,iso4b(ii)) = zvols(1:kproma,:,ii,1) / (1.e6_dp * mvsu)
       !--- Organic carbon
       paerml(1:kproma,:,iocb(ii)) = zvols(1:kproma,:,ii,2) * rhooc * 1.e9_dp
    END DO

    DO ii = in2a, fn2a
       !--- Sulfate volume
       paerml(1:kproma,:,iso4b(ii)) = zvols(1:kproma,:,ii,1)  / (1.e6_dp * mvsu)
       !--- Organic carbon
       paerml(1:kproma,:,iocb(ii)) = zvols(1:kproma,:,ii,2) * rhooc * 1.e9_dp
       !--- Black carbon
       paerml(1:kproma,:,ibcb(ii)) = zvols(1:kproma,:,ii,3) * rhobc * 1.e9_dp
       !--- Sea salt
       paerml(1:kproma,:,issb(ii)) = zvols(1:kproma,:,ii,4) * rhoss * 1.e9_dp
       !--- Mineral dust
       paerml(1:kproma,:,idub(ii)) = zvols(1:kproma,:,ii,5) * rhodu * 1.e9_dp
    END DO

    DO ii = in2b, fn2b
       !--- Sulfate volume
       paerml(1:kproma,:,iso4b(ii)) = zvols(1:kproma,:,ii,1)  / (1.e6_dp * mvsu)
       !--- Organic carbon
       paerml(1:kproma,:,iocb(ii)) = zvols(1:kproma,:,ii,2) * rhooc * 1.e9_dp
       !--- Black carbon
       paerml(1:kproma,:,ibcb(ii)) = zvols(1:kproma,:,ii,3) * rhobc * 1.e9_dp
       !--- Mineral dust
       paerml(1:kproma,:,idub(ii)) = zvols(1:kproma,:,ii,5) * rhodu * 1.e9_dp
    END DO

  END SUBROUTINE salsa

END MODULE mo_salsa
