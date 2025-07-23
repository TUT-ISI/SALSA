MODULE TM5M7_OPTICS_DATA

    ! hhalonen

    USE parkind1, ONLY: JPIM, JPRB

    IMPLICIT NONE

    ! wavelength type (to be used by methods using the optics)
    TYPE, PUBLIC :: WAVELENDEP 
        REAL  :: wl             ! user requested wavelength    unit = um (e.g. 0.550)
        REAL, DIMENSION(7) :: n = 0.0_JPRB ! SO4, BC, OC, SOA, SS, DU, WATER
        REAL, DIMENSION(7) :: k = 0.0_JPRB ! SO4, BC, OC, SOA, SS, DU, WATER
        LOGICAL :: split = .false.
        LOGICAL :: insitu = .false.
    END TYPE WAVELENDEP

    INTEGER(KIND=JPIM), PARAMETER :: NASWBAND = 14
    INTEGER(KIND=JPIM), PARAMETER  :: NALWBAND = 16
    REAL(KIND=JPRB),DIMENSION(NALWBAND)   ::  ALWWN1
    REAL(KIND=JPRB),DIMENSION(NALWBAND)   ::  ALWWN2

    TYPE(WAVELENDEP), PUBLIC :: ASWBAND(NASWBAND)

END MODULE TM5M7_OPTICS_DATA
