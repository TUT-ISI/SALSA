SUBROUTINE TM5M7_INIT
    ! Init routine for IFS tm5m7 aerosol
    ! TM5M7_INIT is called from mo_ham_rad.
    ! Author: Hermanni Halonen
    USE TM5M7_OPTICS_DATA, ONLY : ASWBAND, ALWWN1, ALWWN2
    USE parkind1,          ONLY : JPRB

    ! Wavelength values
    ASWBAND(14)%wl = 5.254_JPRB
    ASWBAND(13)%wl = 0.257_JPRB
    ASWBAND(12)%wl = 0.313_JPRB
    ASWBAND(11)%wl = 0.398_JPRB
    ASWBAND(10)%wl = 0.530_JPRB
    ASWBAND( 9)%wl = 0.697_JPRB
    ASWBAND( 8)%wl = 0.973_JPRB
    ASWBAND( 7)%wl = 1.269_JPRB
    ASWBAND( 6)%wl = 1.447_JPRB
    ASWBAND( 5)%wl = 1.767_JPRB
    ASWBAND( 4)%wl = 2.040_JPRB
    ASWBAND( 3)%wl = 2.308_JPRB
    ASWBAND( 2)%wl = 2.752_JPRB
    ASWBAND( 1)%wl = 3.407_JPRB

    !LW wavenumbers for ham optics
    ALWWN1 = (/ & !< Spectral band lower boundary in wavenumbers
    &   10._JPRB, 350._JPRB, 500._JPRB, 630._JPRB, 700._JPRB, 820._JPRB, &
    &  980._JPRB,1080._JPRB,1180._JPRB,1390._JPRB,1480._JPRB,1800._JPRB, &
    & 2080._JPRB,2250._JPRB,2380._JPRB,2600._JPRB/)
    ALWWN2 = (/ & !< Spectral band upper boundary in wavenumbers
    &  350._JPRB, 500._JPRB, 630._JPRB, 700._JPRB, 820._JPRB, 980._JPRB, &
    & 1080._JPRB,1180._JPRB,1390._JPRB,1480._JPRB,1800._JPRB,2080._JPRB, &
    & 2250._JPRB,2380._JPRB,2600._JPRB,3250._JPRB/)

END SUBROUTINE TM5M7_INIT