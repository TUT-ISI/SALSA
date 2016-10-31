MODULE mo_time_control

  ! ----------------------------------------------------------------------------
  !+
  !
  ! ECHAM Date and Time Control Interface Module
  ! --------------------------------------------
  !
  ! * interface between the time control modules and ECHAM
  !
  ! * definitions of global date/time constants
  !
  ! * the date/time functions are controlled by three structures
  !
  !  TIME_EVENT ..... time dependend handling of operations
  !
  !  TIME_MANAGER ... time axis informations of the model
  !                   like the start date, present position (= the time step)
  !                   time interval between two steps
  !
  !  TIME_DAYS ...... specific dates of the model history
  !                   like start date, stop date, ...
  !
  ! Authors:
  !
  ! I. Kirchner, MPI, April 2000
  ! I. Kirchner, MPI, October/December 2000
  ! I. Kirchner, MPI, March 2001, revision
  ! S. Legutke,  MPI, Aug   2001, separate events for model coupling read/write
  ! I. Kirchner, MPI, Sep/Oct 2001, revision
  ! I. Kirchner, MPI, Aug 2002, add change_present_date
  ! I. Kirchner, FUB, February 2003, revision/code review
  ! L. Kornblueh MPI, April 2003, additional features for debugging and testing
  ! U. Schulzweida, MPI, March 2007, added weights for daily interpolation
  ! S. Rast, MPI, April 2010, added weights for interpolation with respect to 
  !                           radiation time step
  ! D. Klocke, MPI, Nov 2010, changed time step lenght for first time step in
  !                           NWP mode and reset switch after initial time step
  !
  ! external modules
  !
  ! ECHAM specific modules
  !   mo_kind
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
!-
  USE mo_kind,            ONLY: dp

  IMPLICIT NONE

  PUBLIC

  ! ***** definition of DATE/TIME variables/constants


  REAL(dp)           :: delta_time    = 1.0_dp ! distance of adjacent times
  REAL(dp)           :: time_step_len = 1.0_dp ! forecast time step,
                                               ! at beginning equal delta_time

CONTAINS

END MODULE mo_time_control
!
! ------------------------------------------------------------------------------
