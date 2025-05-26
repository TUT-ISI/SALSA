Module mo_time_control
  USE mo_kind, only : dp
  IMPLICIT NONE
  ! <-- thk: bug fix
  !REAL(dp), PUBLIC   :: time_step_len    = 720.0_dp
  REAL(dp), PUBLIC :: time_step_len
  !REAL(dp), PUBLIC :: delta_time !eehol: not needed
  ! this does not work:
  !ASSOCIATE(time_step_len=>YRPHY2%TSPHY) !TeMi

  CONTAINS
  ! --> thk
END MODULE mo_time_control

