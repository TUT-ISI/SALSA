MODULE driver_input
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
  !-------------------------------------------------------
  !-------------------------------------------------------
  !
  ! These are defined in ECHAM5 - here given for
  !	testing purposes	
  !
  !  NB: Check units carefully when coupling with ECHAM5
  !
  !-------------------------------------------------------

  USE mo_kind

  INTEGER, PARAMETER :: &
       
       kproma = 1, &  ! number of horiz. grid points in the slab
       kbdim = 1,  &  ! defines one dimension of matrices
                                ! (for regular grid kbdim = kproma ????)
       klev = 1       ! number of vertical levels

END MODULE driver_input
