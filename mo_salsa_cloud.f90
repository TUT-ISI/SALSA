MODULE mo_salsa_cloud

  !*********************************************************
  !  MOD_AERO CLOUD
  !*********************************************************
  ! 
  ! Purpose: Calculates the number of activated cloud 
  ! droplets according to parameterizations by:
  !
  ! Abdul-Razzak et al: "A parameterization of aerosol activation - 
  !                      3. Sectional representation"
  !                      J. Geophys. Res. 107, 10.1029/2001JD000483, 2002. 
  !                      [Part 3]
  !
  ! Abdul Razzak et al: "A parameterization of aerosol activation - 
  !                      1. Single aerosol type"
  !                      J. Geophys. Res. 103, 6123-6130, 1998. 
  !                      [Part 1]
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  ! 
  !
  ! Coded by:
  ! ---------
  ! T. Anttila (FMI)     2007
  ! H. Kokkola (FMI)     2007
  ! A.-I. Partanen (FMI) 2007
  !

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

  !*********************************************************

CONTAINS 

  SUBROUTINE cloud_activation(kproma, kbdim, klev, &
       n_aero, vols, temp, pres, cd, w)

    USE mo_kind,        ONLY : dp

    USE mo_submctl, ONLY :                      &
         grav,                                      &
         rg,                                        & ! molar gas constant 
                                                      ! [J/(mol K)]
         slim,                                      &
         surfw0,                                    & ! surface tension 
                                                      ! of water [J/m2]
         nbin,                                      & ! number of size bins 
                                                      ! in subranges
         nlim,                                      & ! lowest possible particle conc. in a bin [#/m3]
         rhosu, msu,                                & ! properties of compounds
         rhooc, moc,                                &
         rhoss, mss,                                &
         rhowa, mwa,                                &
         pi,                                        &
         pi6,                                       &
         cpa,                                       &
         mair,                                      &
         in1a,in2a,in2b,fn1a,fn2a, fn2b,            & ! size regime bin indices
         vlolim, vhilim,                            & ! low and high volume
                                                      ! limits for bins
         vratiohi, vratiolo,                        & ! ratio of volume in the bin (high and low) boundary and 
                                                      ! the middle of the bin

         dpmid                                        ! diameter in the
                                                      !  middle of the bin

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::              &
             kproma,                    & ! number of horiz. grid points 
             kbdim,                     & ! dimension for arrays 
             klev                       ! number of vertical levels 

    REAL(dp), INTENT(INOUT)::             &
             n_aero(kbdim,klev,fn2b), & ! number concentration  [#/m3]
             pres(kbdim,klev),        & ! atmospheric pressure at grid point [Pa]
             temp(kbdim,klev),        & ! temperature [K]
             vols(kbdim,klev,fn2b,5), & ! total volume concentrations of each
                                          ! chem. compound in a size bin
             w(kbdim,klev)


    REAL(dp), INTENT(OUT) ::            &
             cd(kbdim,klev)             ! number of cloud droplets


    !-- local variables --------------
    INTEGER :: ii, jj, kk, mm             ! loop indices

    REAL(dp) ::                         &
             sil,                       & ! critical supersaturation 
                                          !     at the upper bound of the bin 
             siu,                       & !  "  at the lower bound of the bin
             scrit(fn2b),               & !  "  at the center of the bin
             aa,                        & ! curvature (Kelvin) effect [m]
             bb,                        & ! solute (Raoult) effect [m3]
             ns(fn2b),                  & ! number of moles of solute
             nshi,                      & !             " at the upper bound of the bin
             nslo,                      & !             " at the lower bound of the bin
             ratio,                     & ! volume ratio
             vmiddle(kbdim,klev,fn2b),  & ! volume in the middle of the bin [m3]
             s_max,                     & ! maximum supersaturation
             s_eff,                     & ! effective supersaturation
             x, x1, x2, x3, a1, sum1,   & ! technical variables
             ka1,                       & ! thermal conductivity
             dv1,                       & ! diffusion coefficient
             Gc,                        & ! growth coefficient
             alpha,                     & ! see Abdul-Razzak and Ghan, part 3
             gamma,                     & ! see Abdul-Razzak and Ghan, part 3
             L,                         & ! latent heat of evaporation
             ps,                        & ! saturation vapor pressure of water [Pa]
             khi,                       & ! see Abdul-Razzak and Ghan, part 3
             theta,                     & ! see Abdul-Razzak and Ghan, part 3
             frac(fn2b),                & ! fraction of activated droplets in a bin
             ntot,                      & ! total number conc of particles [#/m3]
             dinsol(fn2b),              & ! diameter of the insoluble fraction [m]
             dinsolhi,                  & !    "   at the upper bound of a bin [m]
             dinsollo,                  & !    "   at the lower bound of a bin [m]
             dcrit(kbdim,klev,fn2b),    & ! critical diameter [m]
             V,                         & ! updraft velocity [m/s]
             rref,                      & ! reference radius [m]
             A, B, D_p0, dmx,           & !
             vlo, k

    bb = 6._dp*mwa/(pi*rhowa)             ! Raoult effect [m3]
                                          ! NOTE!
                                          ! bb must be multiplied
                                          ! by the number of moles of
                                          ! solute

    cd = 0._dp

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          vmiddle(ii,jj,in1a:fn2b) = pi6*dpmid(in1a:fn2b)**3

       END DO
    END DO

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid
          V  = w(ii,jj)

          aa = 4._dp*mwa*surfw0/(rg*rhowa*temp(ii,jj)) ! Kelvin effect [m]
          A  = aa * 1.e6_dp                            !     "         [um] 
          x  = 4._dp*aa**3/(27._dp*bb)

          ntot = 0._dp
          sum1 = 0._dp

          !-- subrange 1a

          !-- calculation of critical superaturation in the middle of the bin

          !-- volume in the middle of the bin

          DO kk = in1a, fn1a

             IF (n_aero(ii,jj,kk) > nlim) THEN

                !-- number of moles of solute in one particle [mol]
                ns(kk) = (3._dp*vols(ii,jj,kk,1)*rhosu/msu   +                &
                     vols(ii,jj,kk,2)*rhooc/moc)/                             &
                     n_aero(ii,jj,kk)

                !-- critical supersaturation, Kohler equation
                scrit(kk) = exp(sqrt(x/ns(kk))) - 1._dp

                !-- sums in equation (8), part 3
                ntot = ntot + n_aero(ii,jj,kk)
                sum1 = sum1 + n_aero(ii,jj,kk)/scrit(kk)**(2._dp/3._dp)

             END IF

          END DO

          !-- subrange 2a

          DO kk = in2a, fn2a

             IF (n_aero(ii,jj,kk) > nlim .OR. (vols(ii,jj,kk,1) + vols(ii,jj,kk,2) + &
                                               vols(ii,jj,kk,4))/ (vols(ii,jj,kk,3) + &
                                               vols(ii,jj,kk,5)) > 1.e-20_dp) THEN

                ns(kk) = (3._dp*vols(ii,jj,kk,1)*rhosu/msu  +                 &
                        vols(ii,jj,kk,2) * rhooc/moc +                        &
                        vols(ii,jj,kk,4) * rhoss/mss)/                        &
                        n_aero(ii,jj,kk)  

                !-- critical supersaturation, Kohler equation
                scrit(kk) = exp(sqrt(x/ns(kk))) - 1._dp

                !-- sums in equation (8), part 3
                ntot = ntot + n_aero(ii,jj,kk)
                sum1 = sum1 + n_aero(ii,jj,kk)/scrit(kk)**(2._dp/3._dp)

             END IF

          END DO

          !-- subrange 2b

          DO kk = in2b, fn2b

             IF (n_aero(ii,jj,kk) > nlim .OR. (vols(ii,jj,kk,1) + vols(ii,jj,kk,2)  + &
                                               vols(ii,jj,kk,4))/ (vols(ii,jj,kk,3) + &
                                               vols(ii,jj,kk,5)) > 1.e-20_dp) THEN

                ns(kk) = (3._dp*vols(ii,jj,kk,1)*rhosu/msu  +                 &
                        vols(ii,jj,kk,2) * rhooc/moc +                        &
                        vols(ii,jj,kk,4) * rhoss/mss)/                        &
                        n_aero(ii,jj,kk)  

                !-- critical supersaturation, Kohler equation
                scrit(kk) = exp(sqrt(x/ns(kk))) - 1._dp

                !-- sums in equation (8), part 3
                ntot = ntot + n_aero(ii,jj,kk)
                sum1 = sum1 + n_aero(ii,jj,kk)/scrit(kk)**(2._dp/3._dp)

             END IF

          END DO

          IF(ntot < nlim) CYCLE

          !-- latent heat of evaporation [J/kg]
          L     = 2.501e6_dp-2370._dp*(temp(ii,jj)-273.15_dp)

          !-- saturation vapor pressure of water [Pa]
          a1    = 1._dp-(373.15_dp/temp(ii,jj))
          ps    = 101325._dp*                                                 &
               exp(13.3185_dp*a1-1.976_dp*a1**2-0.6445_dp*a1**3-0.1299_dp*a1**4)     

          !-- part 1, eq (11)
          alpha = grav*mwa*L/(cpa*rg*temp(ii,jj)**2)-                            &
               grav*mair/(rg*temp(ii,jj))

          !-- part 1, eq (12)
          gamma = rg*temp(ii,jj)/(ps*mwa) &
               + mwa*L**2/(cpa*pres(ii,jj)*mair*temp(ii,jj))

          !-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
          x1 = pres(ii,jj) / 101325._dp
          dv1= 1.e-4_dp * (0.211_dp/x1) * ((temp(ii,jj)/273._dp)**1.94_dp)

          rref = 10.e-9_dp
          !-- corrected diffusivity, part 1, eq (17)
          ! dv = dv1 / (rref/(rref + deltaV) + (dv1/(rref * alphac)) *        &
          !     SQRT(2.*pi*mwa/(rg*temp(ii,jj))))

          !-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
          ka1= 1.e-3_dp * (4.39_dp + 0.071_dp * temp(ii,jj))

          !-- growth coefficient, part 1, eq (16)
          !-- (note: here uncorrected diffusivities and conductivities are used
          !    based on personal communication with H. Abdul-Razzak, 2007)
          Gc = 1._dp/(rhowa*rg*temp(ii,jj)/(ps*dv1*mwa) +                      &
               L*rhowa/(ka1*temp(ii,jj)) * (L*mwa/(temp(ii,jj)*rg)-1._dp))

          !-- effective critical supersaturation: part 3, eq (8)
          s_eff = (ntot/sum1)**(3._dp/2._dp)

          !-- part 3, equation (5)
          theta = ((alpha*V/Gc)**(3._dp/2._dp))/(2._dp*pi*rhowa*gamma*ntot)

          !-- part 3, equation (6)
          khi = (2._dp/3._dp)*aa*SQRT(alpha*V/Gc)

          !-- maximum supersaturation of the air parcel: part 3, equation (9)
          s_max = s_eff / SQRT(0.5_dp*(khi/theta)**(3._dp/2._dp)              &
               + ((s_eff**2)/(theta+3._dp*khi))**(3._dp/4._dp))

          frac = 0._dp

          DO kk = in1a, fn2b

             frac(kk) = 0._dp

             IF (n_aero(ii,jj,kk) < nlim) CYCLE

             !-- moles of solute in particle at the upper bound of the bin
             nshi = ns(kk)*vratiohi(kk)

             !-- critical supersaturation
             sil = exp(sqrt(x/nshi)) - 1._dp

             IF(s_max < sil) CYCLE

             !-- moles of solute at the lower bound of the bin:
             nslo = ns(kk)*vratiolo(kk)

             !-- critical supersaturation
             siu = exp(sqrt(x/nslo)) - 1._dp

             !-- fraction of activated in a bin, eq (13), part 3
             frac(kk) = min(1._dp,log(s_max/sil)/log(siu/sil))

          END DO

          !-- number of cloud droplets in a grid cell
          
          cd(ii,jj) = sum(frac(in1a:fn2b)*n_aero(ii,jj,in1a:fn2b))

       END DO

    END DO

  END SUBROUTINE cloud_activation

END MODULE mo_salsa_cloud
