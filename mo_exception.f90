MODULE mo_exception

  USE mo_doctor, ONLY: nerr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn

  INTEGER, PARAMETER :: em_none = 0 
  INTEGER, PARAMETER :: em_info = 1
  INTEGER, PARAMETER :: em_warn = 2

  CHARACTER(512) :: message_text = ''

CONTAINS

  SUBROUTINE finish (name, text, exit_no)

    CHARACTER(*) :: name
    CHARACTER(*), OPTIONAL :: text
    INTEGER, OPTIONAL :: exit_no
    INTEGER           :: iexit

    IF (PRESENT(text)) THEN
      WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      WRITE (nerr,'(1x,a,a)') TRIM(name), ': '
    ENDIF

    WRITE (nerr,'(/,80("-"),/)')

    WRITE (nerr,'(/,80("*"),/)')

    STOP

  END SUBROUTINE finish

  SUBROUTINE message (name, text, kout, klevel)

    CHARACTER (*) :: name, text
    INTEGER, INTENT(in), OPTIONAL :: kout
    INTEGER, INTENT(in), OPTIONAL :: klevel

    INTEGER :: iout
    INTEGER :: ilevel

    IF (PRESENT(kout)) THEN
       iout = kout
    ELSE
       iout = nerr
    END IF
    
    IF (PRESENT(klevel)) THEN
       ilevel = klevel
    ELSE
       ilevel = em_none
    END IF
    
    SELECT CASE (ilevel)
    CASE (em_none)
    CASE (em_info)
    CASE (em_warn)
    END SELECT
    
    IF (name == '') THEN
       WRITE(iout,'(1x,a)') TRIM(text)
    ELSE
       WRITE(iout,'(1x,a,": ",a)') TRIM(name), TRIM(text)
    END IF
    
  END SUBROUTINE message
  
END MODULE mo_exception
