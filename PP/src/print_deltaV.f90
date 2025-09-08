!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM print_deltaV
  !-----------------------------------------------------------------------
  !
  ! Print the matrix elements of e^iq.r on the KS basis
  ! Requires a preceeding bands calculation to get the KS states at k and k+q
  !
  USE parameters,       ONLY : npk
  USE constants,        ONLY : rytoev
  USE kinds,            ONLY : DP
  USE klist,            ONLY : ngk, igk_k, nks, xk
  USE io_files,         ONLY : prefix, tmp_dir
  USE io_files,         ONLY : iunwfc, nwordwfc
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE environment,      ONLY : environment_start, environment_end
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : mp_startup
  USE mp_images,        ONLY : intra_image_comm
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : noncolin, npol
  USE fft_base,         ONLY : dffts
  USE buffers,          ONLY : get_buffer, open_buffer
  USE control_flags,    ONLY : io_level
  USE cell_base,        ONLY : at
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: outdir
  INTEGER :: ios
  LOGICAL :: needwf = .TRUE.
  INTEGER, PARAMETER :: N_MAX_BOXES = 999
  !
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:)
  ! e^iqr|psi_k>
  COMPLEX (DP), ALLOCATABLE :: deltaV(:,:)
  ! \Delta^q V(k)_ij = <psi_k+q,i|e^iqr|psi_k,j> 
  COMPLEX(DP) , ALLOCATABLE ::  aux (:)
  ! work space
  COMPLEX(DP), ALLOCATABLE, TARGET :: evk (:,:)
  ! wfc at k 
  COMPLEX(DP), POINTER :: evq (:,:)
  ! wfc at k+q
  !
  INTEGER :: ik, npw, npwq, ikk, ikq
  !
  INTEGER :: ibnd, jbnd, ig, eff_ik
  ! counter on bands
  !
  REAL(DP) :: xq(3), xk_(3)
  !
  LOGICAL :: exst, exst_mem, lgamma
  CHARACTER (len=10):: file_ik
  CHARACTER (len=50):: file_q1
  CHARACTER (len=50):: file_q2
  CHARACTER (len=50):: file_q3
  CHARACTER (len=50):: file_name
  !
  NAMELIST / pp / outdir, prefix
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  !
  CALL environment_start ( 'PRINT_DV' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, pp, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, intra_image_comm )
  IF (ios /= 0) CALL errore ('print_deltaV', 'reading namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir,   ionode_id, intra_image_comm )
  CALL mp_bcast( prefix,    ionode_id, intra_image_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  !
  ! BODY of the Program
  !
  CALL open_buffer(iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst)
  !
  ALLOCATE (dvpsi ( npwx*npol , nbnd))
  ALLOCATE (deltaV (nbnd, nbnd))
  ALLOCATE( evk(npwx*npol,nbnd) )
  !
  xq(:)=xk(:,2)-xk(:,1)
  call cryst_to_cart (1, xq, at, - 1)
  WRITE(stdout,'(/,5X, "INFO xq    =", 3F8.4, " [Cart]"  )') xk(:,2)-xk(:,1)
  WRITE(stdout,'(  5X, "INFO xq    =", 3F8.4, " [Crys]",/)') xq(:)
  WRITE (file_q1, '(3F8.4)') xq(1)
  WRITE (file_q2, '(3F8.4)') xq(2)
  WRITE (file_q3, '(3F8.4)') xq(3)
  !
  lgamma = .FALSE.
  IF ( ALL( ABS(xq(:)) < 1.d-5 ) ) lgamma = .TRUE.
  !
  IF (lgamma) then
     WRITE(stdout,'(/,5X, "INFO evq is a pointer to evk"  )') 
     !  q=0  : evq is a pointer to evk
     evq  => evk
  ELSE
     !  q!=0 : evq is allocated and calculated at point k+q
     ALLOCATE (evq ( npwx*npol , nbnd))
  endif
  !
  DO ik = 1, nks, 2
    ! effective k point index not counting the k+q points
    eff_ik = ik/2+1
    !
    DeltaV = CMPLX(0,0,kind=DP)
    !
    WRITE(stdout,'(/, 5X, "INFO ik  =", I5)') ik
    WRITE(stdout,'(   5X, "INFO ikq =", I5)') ik+1
    !
    ALLOCATE (aux(dffts%nnr))
    !
    !  The unperturbed wavefunctions must be multiplied by e^{iqr}.
    !  This means that we have to order the coefficients with the mesh
    !  of k+q. We do this by distributing the vectors on the FFT mesh
    !  with the indices of k+G (igk) and then saving them with the mesh
    !  of k+q+G (igkq)
    !
    dvpsi(:,:) = (0.D0, 0.D0)
    !
    ikk = ik
    ikq = ik+1
    npw = ngk(ikk)
    npwq= ngk(ikq)
    !
    WRITE(stdout,'(5X, "INFO xk    =", 3F8.4)') xk(:,ik)
    WRITE(stdout,'(5X, "INFO xkq   =", 3F8.4)') xk(:,ikq)
    !
    WRITE(stdout, '(5X, "INFO npw  =", I5)') npw
    WRITE(stdout, '(5X, "INFO npwq =", I5)') npwq
    !
    WRITE(stdout, '(5X, "INFO iunwfc  =", I5)') iunwfc
    !
    CALL get_buffer(evk, nwordwfc, iunwfc, ikk)
    CALL get_buffer(evq, nwordwfc, iunwfc, ikq)
    !
    DO ibnd = 1, nbnd
       aux(:) = (0.D0, 0.D0)
       DO ig = 1, npw
          aux (dffts%nl (igk_k (ig,ikk) ) ) = evk (ig, ibnd)
       ENDDO
       !
       DO ig = 1, npwq
          dvpsi (ig, ibnd) = aux (dffts%nl (igk_k (ig,ikq) ) )
       ENDDO
       IF (noncolin) THEN
          aux(:) = (0.d0, 0.d0)
          DO ig = 1, npw
             aux (dffts%nl (igk_k (ig,ikk) ) ) = evk (ig+npwx, ibnd)
          ENDDO
          !
          DO ig = 1, npwq
             dvpsi (ig+npwx, ibnd) = aux (dffts%nl (igk_k (ig,ikq) ) )
          ENDDO
       END IF
       !
    ENDDO
    !
    CALL ZGEMM( 'C', 'N', nbnd, nbnd, npwx*npol, (1.d0,0.d0), &
                      evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), deltaV, nbnd )
    !          
    CALL mp_sum( deltaV, intra_bgrp_comm )
    !
    !WRITE(*,*) evk(1:3,1)
    !WRITE(*,*) evq(1:3,1)
    !WRITE(*,*) evq(1:3,1)/evk(1:3,1)
    !WRITE(*,*) SUM( CONJG(evq(:,1))*evk(:,1))
    !WRITE(*,*) ABS(SUM( CONJG(evq(:,1))*evk(:,1)))
    !
    ! Writing on file
    WRITE (file_ik,'(i0)') eff_ik
    file_name = 'DeltaV_q' // trim(adjustl(file_q1)) // '_'   &  
                           // trim(adjustl(file_q2)) // '_'   & 
                           // trim(adjustl(file_q3)) // '_ik' & 
                           // trim(adjustl(file_ik)) // '.dat'
    OPEN (file=trim(file_name), unit=101)
    write(*,*) file_name
    WRITE(101, '("# Wave vector of the perturbation xq = ", 3F12.8, " [Cart]")') xk(:,2)-xk(:,1)
    WRITE(101, '("# Wave vector of the perturbation xq = ", 3F12.8, " [Crys]")') xq
    WRITE(101, '("# ik     = ", I5)') eff_ik
    xk_(:) = xk(:,ik)
    call cryst_to_cart (1, xk_, at, - 1)
    WRITE(101, '("# xk(ik) = ", 3F12.8, " [Cart]")') xk(:,ik)
    WRITE(101, '("# xk(ik) = ", 3F12.8, " [Crys]")') xk_(:)
    WRITE(101, '("# ibnd jbnd DV(ik)_ij = <u_k+q,i |e^iq.r | u_k,j> ", 3F12.8)') 
    DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
      WRITE(101, '(2I5, 3X, 2E24.16)') ibnd, jbnd, deltaV(ibnd,jbnd)
     ENDDO
    ENDDO
    CLOSE(101)
    !
    DEALLOCATE (aux)
    !
  ENDDO
  !
  DEALLOCATE (dvpsi )
  DEALLOCATE (deltaV)
  DEALLOCATE (evk)
  IF (lgamma) THEN
     IF(associated(evq)) NULLIFY(evq)
  ELSE
     IF(associated(evq)) DEALLOCATE(evq)
  ENDIF
  !
  CALL environment_end ( 'PRINT_DV' )
  !
  CALL stop_pp
  !
END PROGRAM print_deltaV

