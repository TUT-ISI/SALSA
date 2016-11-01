MODULE mo_time_control

  USE mo_kind,            ONLY: dp

  IMPLICIT NONE

  PUBLIC

  REAL(dp)           :: delta_time    = 1.0_dp ! distance of adjacent times
  REAL(dp)           :: time_step_len = 1.0_dp ! forecast time step,
                                               ! at beginning equal delta_time

CONTAINS

END MODULE mo_time_control
!
! ------------------------------------------------------------------------------
