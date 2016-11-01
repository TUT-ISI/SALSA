!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_salsa_trac.f90
!!
!! \brief
!! mo_salsa_trac contains routines to requests tracers for ECHAM/HAM and 
!! prescribes their physical and chemical properties.
!! It controls the aerosol physics by providing the necessary switches.
!!
!! \author Philip Stier (MPI-Met)
!!
!! \responsible_coder
!! Martin G. Schultz, m.schultz@fz-juelich.de
!!
!! \revision_history
!!   -# P. Stier (MPI-Met) - original code (2001)
!!   -# D. O'Donnell (MPI-Met) - code generalization and changes for soa (2009-02-xx)
!!   -# K. Zhang (MPI-Met) - adaption for new species list and tracer defination (2009-08-11) 
!!   -# M.G. Schultz (FZ Juelich) - cleanup and adaptation to new structure (2009-11-20)
!!   -# T. Bergman (FMI) - nmod->nclass to facilitate new aerosol models (2013-02-05)
!! 	 -# A. Laakso (FMI) - Tracers for SALSA
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
!!### Questionable whether this module is needed with new emissions scheme.
!! idt_ seem to be used primarily for identiyfing tracers for emissions.


MODULE mo_salsa_trac

  ! Parameters:
  ! -----------
  ! User defined flags: density   density                    [kg m-3]
  !                     osm       osmotic coefficient        [???]
  !                     nion      number of ions the tracer 
  !                               dissolves into             [1]

  USE mo_kind,          ONLY: dp

  USE mo_submctl,   ONLY: in1a, fn1a,           &
                              in2a, fn2a,           &
                              in2b, fn2b, nbins
			      
                        
  
  IMPLICIT NONE

  !--- Public entities:

  PUBLIC :: idt_dms,   idt_so2,   idt_so4,   idt_ocnv,     &
            idt_ms4,   idt_moc,   idt_mbc,   idt_mss,      & ! SALSA indices
            idt_mdu,   idt_mws,   idt_n, idt_mwa
 
  !--- Module variables:
  !
  !    Tracer indices:
  !
  !    Legend: iABBCD
  !
  !            A:  m  = particle mass mixing ratio, n number mixing ratio
  !            BB: s4 = sulfate, bc/oc = black/organic carbon, du = dust, ss = seasalt
  !            C:  n  = nucleation , k = Aitken, a = accumulation, c = coarse mode
  !            D:  i  = insoluble,  s = soluble

  INTEGER :: idt_dms    ! mass mixing ratio dms
  INTEGER :: idt_so2    ! mass mixing ratio so2
  INTEGER :: idt_so4    ! mass mixing ratio so4
  INTEGER :: idt_ocnv   ! mass mixing ratio nonvolatile organic

  INTEGER :: idt_cdnc_ham   ! cloud droplet number concentration
  INTEGER :: idt_icnc_ham   ! ice   cristal number concentration

  INTEGER :: idt_ms4(fn2b)    ! mass mixing ratio sulfate
  INTEGER :: idt_moc(fn2b)    ! mass mixing ratio organic carbon
  INTEGER :: idt_mbc(fn2b)    ! mass mixing ratio black carbon
  INTEGER :: idt_mss(fn2a)    ! mass mixing ratio seasalt
  INTEGER :: idt_mdu(fn2b)    ! mass mixing ratio mineral dust
  INTEGER :: idt_mws(fn2b)    ! mass mixing ratio water soluble
  INTEGER :: idt_n(fn2b)      ! number mixing ratio
  INTEGER :: idt_mwa(nbins)   ! mass mixing ratio aerosol water

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> 
!! Define HAM tracers 
!! 
!! @author see module info 
!!
!! $Id: 1423$
!!
!! @par Revision History
!! see module info 
!!
!! @par This subroutine is called by
!! to_be_filled
!!
!! @par Externals:
!! <ol>
!! <li>none
!! </ol>
!!
!! @par Notes
!! 
!! @par Responsible coder
!! m.schultz@fz-juelich.de
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE ham_salsa_set_idt

  USE mo_submctl,     ONLY:  iso4b, ibcb,iocb,issb, idub,	&
											in1a, in2a, in2b,  &
    		                               fn1a, fn2a, fn2b

  INTEGER :: 			i, kk
  
  kk=0
  DO i = in1a,fn2b
     kk = kk + 1
     idt_ms4(i) = kk 
     idt_n(i) = i
  END DO
  
  DO i = in1a,fn2b
     kk = kk + 1
     idt_moc(i)= kk
  END DO
  
  DO i = in2a,fn2b
     kk = kk + 1
     idt_mbc(i) = kk
  END DO
     
  DO i = in2a,fn2b
     kk = kk + 1
     idt_mdu(i) = kk
  END DO
     
  DO i = in2a,fn2a
     kk = kk + 1
     idt_mss(i) = kk
  END DO

  DO i = in1a,fn2a
     kk = kk + 1
     idt_mwa(i) = kk
  END DO
  
END SUBROUTINE ham_salsa_set_idt 

END MODULE mo_salsa_trac
