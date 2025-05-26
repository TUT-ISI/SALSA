MODULE driver_tracer
  USE mo_kind, ONLY:&
       dp
  USE mo_tracdef, ONLY:&
       ntrac, GAS, AEROSOL, GAS_OR_AEROSOL

  IMPLICIT NONE
  
  ! global variables
  REAL(dp) :: ztmst
  INTEGER, PARAMETER    :: nmaxtrac = 200   ! maximum amount of definable tracers
  REAL(dp), ALLOCATABLE :: pxtm1(:,:,:)     ! field to hold tracer values of last time step, units: [kg/kg] or [#/kg]
  REAL(dp), ALLOCATABLE :: pxtte(:,:,:)     ! field to hold tracer tendencies, units: [kg/kg] or [#/kg]
  
  ! SALSA-specific variables
  INTEGER               :: ngastrac         ! number of gas tracers in the model  -- replace with subm_ngastrac?
  INTEGER, ALLOCATABLE  :: gastrac(:)       ! field to map gas tracer ids



CONTAINS

  FUNCTION new_tracer()
    ! Simple function to count the amout of tracers in the code, which can be used
    ! for mapping tracers of different types
    INTEGER :: new_tracer

    IF (ntrac < nmaxtrac) THEN
       ntrac = ntrac + 1
       new_tracer = ntrac
    ELSE
       WRITE (*,*) 'Trying to define too many tracers -- if you need more; increase "nmaxtrac"'
       STOP
    END IF
  END FUNCTION new_tracer


  SUBROUTINE allocate_tracers(kbdim, kproma, klev, ntrac, pxtm1, pxtte)
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::&
         kbdim, & ! horizontal model dimension
         kproma,& ! horizontal dimension maximum iterator
         klev,  & ! vertical model dimension
         ntrac    ! number of model tracers
    REAL(dp), INTENT(inout), ALLOCATABLE :: pxtm1(:,:,:)     ! field to hold tracer variables
    REAL(dp), INTENT(inout), ALLOCATABLE :: pxtte(:,:,:)     ! field to hold tracer variables
    

    ! Subroutine to allocate memory for pxtm1 and pxtte. 
    ! IMPORTANT: only call init_pxtm1 AFTER all
    ! tracers have been registered with new_tracer
    IF (ALLOCATED(pxtm1)) THEN
       DEALLOCATE(pxtm1)
    END IF

    IF (ALLOCATED(pxtte)) THEN
       DEALLOCATE(pxtte)
    END IF

    ALLOCATE(pxtm1(kbdim, klev, ntrac))
    ALLOCATE(pxtte(kbdim, klev, ntrac))

    pxtm1(1:kproma,:,:) = 0.0_dp
    pxtte(1:kproma,:,:) = 0.0_dp

  END SUBROUTINE allocate_tracers



  SUBROUTINE time_step(kbdim, kproma, klev, ntrac, pxtm1, pxtte)
    ! advance the model by one time step
    ! this means adding the tracer tendency to the current tracer values
    ! and then setting the tracer tendencies to zero

    IMPLICIT NONE
    INTEGER, INTENT(IN) ::&
         kbdim, & ! horizontal model dimension
         kproma,& ! horizontal model dimension upper bound
         klev,  & ! vertical model dimension
         ntrac    ! number of model tracers
    REAL(dp), INTENT(INOUT) :: &
         pxtm1(kbdim,klev,ntrac), & !
         pxtte(kbdim,klev,ntrac)    !
    
    pxtm1(1:kproma,:,:) = pxtm1(1:kproma,:,:) + pxtte(1:kproma,:,:)*ztmst
    pxtte(1:kproma,:,:) = 0.0_dp

  END SUBROUTINE time_step












  SUBROUTINE tracer_to_salsa(kbdim, kproma, klev, pxtm1, pxtte, pcgas, pvols, pnaero, prhoa)

    USE mo_physical_constants, ONLY:&
         avo

    USE mo_species, ONLY:&
         speclist, nspec

    USE mo_ham, ONLY:&
         aerocomp, &
         sizeclass, nclass, &
         subm_naerospec_nowat, &
         nsoa

    USE mo_ham_salsactl, ONLY:&
         in1a, fn1a, in2a, fn2a, in2b, fn2b, &
         iso4b, iocb, ibcb, issb, idub
#ifdef HAMMOZ
    USE mo_ham_vbsctl, ONLY:&
         laqsoa, &
         t_vbs_group, vbs_set,   vbs_ngroup, &
         t_aq_soa,    aqsoa_set, aqsoa_ngroup
#endif

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: &
         kbdim,     & ! model horizontal dimension
         kproma,    & ! horizontal looping index
         klev!,      & ! model vertical dimension
    REAL(dp), INTENT(IN) :: &
         pxtm1(kbdim,klev,ntrac), & ! tracer of previous time step
         pxtte(kbdim,klev,ntrac), & ! tracer tendency
         prhoa(kbdim,klev)          ! air density
    REAL(dp), INTENT(OUT) :: &
         pcgas(kbdim,klev,ngastrac), & !
         pvols(kbdim,klev,nclass,subm_naerospec_nowat), & !
         pnaero(kbdim,klev,nclass)   !

    ! Local variables
#ifdef HAMMOZ
    TYPE(t_vbs_group), POINTER :: zgroup
    TYPE(t_aq_soa), POINTER :: zgroupaq
#endif
    INTEGER :: aero_idx
    INTEGER :: ii, jn, jt, jv, jspec


    ! initializing all concentrations to zero
    pcgas(:,:,:)   = 0.0_dp
    pvols(:,:,:,:) = 0.0_dp
    pnaero(:,:,:)  = 0.0_dp


    ! ----------------------------------------------------------------
    ! Converting gas phase mass mixing ratios to gas concentrations
    ! ----------------------------------------------------------------
    ! jn is the index in pcgas
    ! jt is the index in pxtm1 and pxtte
    ! jspec is the index in speclist
    DO jn = 1,ngastrac
       jt = gastrac(jn)
       ! have to loop through speclist -- introduce extra mapping?
       DO jspec = 1,nspec
          IF (speclist(jspec)%idt == jt) THEN       
             pcgas(1:kproma,:,jn) = avo/speclist(jspec)%moleweight*prhoa(1:kproma,:)*(&
                  pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst&
             )
             EXIT
          END IF
       END DO
    END DO


    ! ----------------------------------------------------------------
    ! Converting particle mass mixing ratios to volume concentrations
    ! ----------------------------------------------------------------
    ! 1a
    DO ii = in1a, fn1a
       !--- Sulfate volume
       jn = iso4b(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,1) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Organic carbon
       jn = iocb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,2) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

    END DO

    ! 2a
    DO ii = in2a, fn2a
       !--- Sulfate volume
       jn = iso4b(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,1) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Organic carbon
       jn = iocb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,2) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)


       !--- Black carbon
       jn = ibcb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,3) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Sea salt
       jn = issb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,4) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Mineral dust
       jn = idub(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,5) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       
    END DO

    ! 2b
    DO ii = in2b, fn2b
       !--- Sulfate volume
       jn = iso4b(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,1) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Organic carbon
       jn = iocb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,2) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Black carbon
       jn = ibcb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,3) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)

       !--- Mineral dust
       jn = idub(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       pvols(1:kproma,:,ii,5) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
      
    END DO

    ! VBS
#ifdef HAMMOZ
    IF (nsoa == 2) THEN
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, vbs_ngroup
                zgroup => vbs_set(jv)               ! pointer to VBS group description
                aero_idx = zgroup%idx(ii)           ! index for (redundant) zaerml
                jt = aerocomp(aero_idx)%idt         ! index for tracer fields
                pvols(1:kproma, :, ii, zgroup%id_vols) = &
                     prhoa(1:kproma,:) / speclist(zgroup%spid)%density*&
                     (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
             END DO
          END IF
       END DO
       IF (laqsoa) THEN    
          DO ii = 1,nclass
             IF ( sizeclass(ii)%lsoainclass ) THEN
                DO jv = 1, aqsoa_ngroup
                   zgroupaq => aqsoa_set(jv)           ! pointer to aqsoa group description
                   aero_idx = zgroupaq%idx(ii)         ! index for (redundant) zaerml
                   jt = aerocomp(aero_idx)%idt         ! index for tracer fields
                   pvols(1:kproma, :, ii, zgroupaq%id_aqsoa) = &
                        prhoa(1:kproma,:) / speclist(zgroupaq%spid)%density*&
                        (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
                END DO
             END IF
          END DO
       END IF
    END IF
#endif
    ! ----------------------------------------------------------------
    ! Converting particle number mixing ratios to number concentrations
    ! ----------------------------------------------------------------
    DO jn=1, nclass
       jt = sizeclass(jn)%idt_no     
       pnaero(1:kproma,:,jn) = prhoa(1:kproma,:)*(pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
    END DO

  END SUBROUTINE tracer_to_salsa





  SUBROUTINE salsa_to_tracer(kbdim, kproma, klev, pxtm1, pxtte, pcgas, pvols, pnaero, prhoa)
    USE mo_physical_constants, ONLY:&
         avo

    USE mo_species, ONLY:&
         speclist, nspec

    USE mo_ham, ONLY:&
         aerocomp, &
         sizeclass, nclass, &
         subm_naerospec_nowat, &
         nsoa

    USE mo_ham_salsactl, ONLY:&
         in1a, fn1a, in2a, fn2a, in2b, fn2b, &
         iso4b, iocb, ibcb, issb, idub
#ifdef HAMMOZ
    USE mo_ham_vbsctl, ONLY:&
         laqsoa, &
         t_vbs_group, vbs_set,   vbs_ngroup, &
         t_aq_soa,    aqsoa_set, aqsoa_ngroup
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: &
         kbdim,     & ! model horizontal dimension
         kproma,    & ! horizontal looping index
         klev!,      & ! model vertical dimension
    REAL(dp), INTENT(OUT) :: &
         pxtte(kbdim,klev,ntrac)    !
    REAL(dp), INTENT(IN) :: &
         pxtm1(kbdim,klev,ntrac), & !
         pcgas(kbdim,klev,ngastrac), & !
         pvols(kbdim,klev,nclass,subm_naerospec_nowat), & !
         pnaero(kbdim,klev,nclass), &   !
         prhoa(kbdim,klev)          ! air density

    ! Local variables
#ifdef HAMMOZ
    TYPE(t_vbs_group), POINTER :: zgroup
    TYPE(t_aq_soa), POINTER :: zgroupaq
#endif
    INTEGER :: aero_idx
    INTEGER :: ii, jn, jt, jv, jspec




    ! ----------------------------------------------------------------
    ! Converting molecular gas concentrations to gas phase mass mixing ratios 
    ! ----------------------------------------------------------------
    ! jn is the index in pcgas
    ! jt is the index in pxtm1 and pxtte
    ! jspec is the index in speclist
    DO jn = 1,ngastrac
       jt = gastrac(jn)
       ! have to loop through speclist -- introduce extra mapping?
       DO jspec = 1,nspec
          IF (speclist(jspec)%idt == jt) THEN
             pxtte(1:kproma,:,jt) = (&
                  speclist(jspec)%moleweight/(avo*prhoa(1:kproma,:))*pcgas(1:kproma,:,jn) - pxtm1(1:kproma,:,jt)&
             )/ztmst
             EXIT
          END IF
       END DO
    END DO


    ! ----------------------------------------------------------------
    ! Converting particle mass mixing ratios to volume concentrations
    ! ----------------------------------------------------------------
    ! 1a
    DO ii = in1a, fn1a
       !--- Sulfate volume
       jn = iso4b(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,1) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,1)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Organic carbon
       jn = iocb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,2) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,2)-pxtm1(1:kproma,:,jt)&
       )/ztmst

    END DO

    ! 2a
    DO ii = in2a, fn2a
       !--- Sulfate volume
       jn = iso4b(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,1) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,1)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Organic carbon
       jn = iocb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,2) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,2)-pxtm1(1:kproma,:,jt)&
       )/ztmst


       !--- Black carbon
       jn = ibcb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,3) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,3)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Sea salt
       jn = issb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,4) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,4)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Mineral dust
       jn = idub(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,5) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,5)-pxtm1(1:kproma,:,jt)&
       )/ztmst
       
    END DO

    ! 2b
    DO ii = in2b, fn2b
       !--- Sulfate volume
       jn = iso4b(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,1) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,1)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Organic carbon
       jn = iocb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,2) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,2)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Black carbon
       jn = ibcb(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,3) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,3)-pxtm1(1:kproma,:,jt)&
       )/ztmst

       !--- Mineral dust
       jn = idub(ii)
       jt = aerocomp(jn)%idt
       jspec = aerocomp(jn)%spid
       !pvols(1:kproma,:,ii,5) = prhoa(1:kproma,:) / speclist(jspec)%density * (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (&
            speclist(jspec)%density/prhoa(1:kproma,:)*pvols(1:kproma,:,ii,5)-pxtm1(1:kproma,:,jt)&
       )/ztmst
      
    END DO
#ifdef HAMMOZ
    ! VBS
    IF (nsoa == 2) THEN
       DO ii = 1,nclass
          IF ( sizeclass(ii)%lsoainclass ) THEN
             DO jv = 1, vbs_ngroup
                zgroup => vbs_set(jv)               ! pointer to VBS group description
                aero_idx = zgroup%idx(ii)           ! index for (redundant) zaerml
                jt = aerocomp(aero_idx)%idt         ! index for tracer fields
                !pvols(1:kproma, :, ii, zgroup%id_vols) = &
                !     prhoa(1:kproma,:) / speclist(zgroup%spid)%density*&
                !     (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
                pxtte(1:kproma,:,jt) = (&
                     speclist(zgroup%spid)%density/prhoa(1:kproma,:)*&
                     pvols(1:kproma,:,ii,zgroup%id_vols)-pxtm1(1:kproma,:,jt)&
                )/ztmst
             END DO
          END IF
       END DO
       IF (laqsoa) THEN    
          DO ii = 1,nclass
             IF ( sizeclass(ii)%lsoainclass ) THEN
                DO jv = 1, aqsoa_ngroup
                   zgroupaq => aqsoa_set(jv)           ! pointer to aqsoa group description
                   aero_idx = zgroupaq%idx(ii)         ! index for (redundant) zaerml
                   jt = aerocomp(aero_idx)%idt         ! index for tracer fields
                   !pvols(1:kproma, :, ii, zgroupaq%id_aqsoa) = &
                   !     prhoa(1:kproma,:) / speclist(zgroup%spid)%density*&
                   !     (pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
                   pxtte(1:kproma,:,jt) = (&
                        speclist(zgroupaq%spid)%density/prhoa(1:kproma,:)*&
                        pvols(1:kproma,:,ii,zgroupaq%id_aqsoa)-pxtm1(1:kproma,:,jt)&
                   )/ztmst
                END DO
             END IF
          END DO
       END IF
    END IF
#endif
    ! ----------------------------------------------------------------
    ! Converting number concentrations to particle number mixing ratios
    ! ----------------------------------------------------------------
    DO jn=1, nclass
       jt = sizeclass(jn)%idt_no     
       !pnaero(1:kproma,:,jn) = prhoa(1:kproma,:)*(pxtm1(1:kproma,:,jt)+pxtte(1:kproma,:,jt)*ztmst)
       pxtte(1:kproma,:,jt) = (pnaero(1:kproma,:,jn)/prhoa(1:kproma,:)-pxtm1(1:kproma,:,jt))/ztmst
    END DO



  END SUBROUTINE salsa_to_tracer


END MODULE driver_tracer
