_chargedens_fromfileMODULE module_chargedens_fromfile
  !
  PRIVATE
  !
  ! Subroutines
  !
  PUBLIC :: print_symmetries,symmetrize_rhor
  !
CONTAINS
  !--------------------------------------------------------------------
  SUBROUTINE print_dmat(iun,fileout,dmat,nsym,l)
    !--------------------------------------------------------------------
    USE kinds,                ONLY : DP
    !
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: fileout
    INTEGER,      INTENT(IN) :: iun,nsym,l 
    REAL(DP),     INTENT(IN) :: dmat(:,:,:)
    INTEGER           :: m0,ldim,isym,m00,ios
    !
    ldim = 2*l+1
    open(unit=iun,file=trim(fileout),action='write',status='replace',iostat = ios)
    write(iun,"('#  D',i1,' FILE')") l
    write(iun,"('#  l nsym')") 
    write(iun,"('# ',i5,x,i5)") l, nsym
    do isym = 1,nsym
      write(iun,"('#  isym ',i5)") isym
      do m0 = 1,ldim
        write(iun,"('#  first-left (row) index',i5)") m0
        do m00 = 1,ldim
          write(iun,"(E23.16,x,E23.16)") dmat(m0,m00,isym)
        enddo
      enddo
    enddo
    write(iun,"('#  END D',i1,' FILE')") l
    close(iun)
    !
  END SUBROUTINE print_dmat 
  !
  !--------------------------------------------------------------------
  SUBROUTINE print_irt(iun,fileout,irt,nsym,nat)
    !--------------------------------------------------------------------
    USE kinds,                ONLY : DP
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN) :: fileout
    INTEGER,      INTENT(IN) :: iun,nsym,nat
    INTEGER,      INTENT(IN) :: irt(:,:)
    INTEGER           :: na,isym,ios
    !
    open(unit=iun,file=trim(fileout),action='write',status='replace',iostat = ios)
    write(iun,"('#  IRT FILE')") 
    do isym = 1,nsym
      write(iun,"('#  isym ',i5)") isym
      do na = 1,nat
          write(iun,"(i5)",advance='no') irt(isym,na) 
      enddo
      write(iun,"()")
    enddo
    write(iun,"('#  END IRT FILE')") 
    close(iun)
    !
  END SUBROUTINE print_irt 
  !
  !--------------------------------------------------------------------
  SUBROUTINE print_symmetries(fileout)
    !--------------------------------------------------------------------
    USE kinds,                ONLY : DP
    USE symm_base,            ONLY : nsym, irt, t_rev, d1, d2, d3
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp 
    USE io_files,             ONLY : prefix, tmp_dir
    USE io_global,            ONLY : stdout, ionode, ionode_id
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN) :: fileout
    INTEGER, EXTERNAL :: find_free_unit
    INTEGER           :: iun,ios
    INTEGER           :: n
    !
    CALL d_matrix (d1, d2, d3)
    IF ( ionode )  THEN
      iun = find_free_unit() 
      CALL print_dmat(iun,'d1_'//trim(fileout),d1,nsym,1) 
      CALL print_dmat(iun,'d2_'//trim(fileout),d2,nsym,2) 
      CALL print_dmat(iun,'d3_'//trim(fileout),d3,nsym,3) 
      CALL print_irt(iun,'irt_'//trim(fileout),irt,nsym,nat)
    ENDIF
    !
  END SUBROUTINE print_symmetries 
  !
  !--------------------------------------------------------------------
  SUBROUTINE symmetrize_rhor(rhor)
    !--------------------------------------------------------------------
    USE kinds,                ONLY : DP
    USE symme,                ONLY : sym_rho, sym_rho_init, sym_rho_deallocate
    USE control_flags,        ONLY : gamma_only
    USE fft_rho,              ONLY : rho_r2g,rho_g2r 
    USE gvect,                ONLY : ngm  
    USE lsda_mod,             ONLY : nspin
    USE fft_base,             ONLY : dfftp
    USE noncollin_module,     ONLY : nspin_mag
    !
    IMPLICIT NONE
    !
    REAL(DP),    INTENT(INOUT)   :: rhor(:,:)
    COMPLEX(DP), ALLOCATABLE     :: rhog(:,:)
    !
    ! ... bring rho(r) to G-space (use psic as work array)
    !
    allocate(rhog(ngm,nspin))
    rhog(:,:) = cmplx(0.d0,0.d0,kind=DP)
    call rho_r2g(dfftp, rhor, rhog)
    !
    ! symmetrize rhog
    !
    CALL sym_rho_init (gamma_only )
    CALL sym_rho ( nspin_mag, rhog )
    CALL sym_rho_deallocate()
    !
    ! retur to real space
    !
    CALL rho_g2r (dfftp, rhog, rhor)
    deallocate(rhog)
    !
  END SUBROUTINE symmetrize_rhor 
  !
END MODULE module_chargedens_fromfile
!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM do_chargedens_fromfile
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions,
  ! calculates Lowdin charges, spilling parameter, projected DOS
  ! or computes the LDOS in a volume given in input as function of energy
  !
  ! See files INPUT_PROJWFC.* in Doc/ directory for usage
  ! IMPORTANT: since v.5 namelist name is &projwfc and no longer &inputpp
  !
  USE io_files,             ONLY : prefix, tmp_dir
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE environment,          ONLY : environment_start, environment_end
  USE mp_global,            ONLY : mp_startup
  USE kinds,                ONLY : DP
  USE mp_images,            ONLY : intra_image_comm
  USE cell_base,            ONLY : omega
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE ener,                 ONLY : ef
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_rho,              ONLY : rho_r2g,rho_g2r 
  USE fft_interfaces,       ONLY : fwfft, invfft, fft_interpolate
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, wk, xk, &
                                   nkstot, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : sym_rho, sym_rho_init, sym_rho_deallocate
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions,        ONLY : evc, psic, psic_nc
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin, npol, lspinorb, nspin_mag
  USE upf_spinorb,          ONLY : fcoef
  USE io_files,             ONLY : restart_dir,iunmix
  USE pw_restart_new,       ONLY : read_collected_wfc
  USE becmod,               ONLY : calbec
  USE uspp_init,            ONLY : init_us_2
  USE mp_pools,             ONLY : me_pool, nproc_pool, my_pool_id, npool, &
                                   inter_pool_comm, intra_pool_comm, root_pool
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE io_rho_xml,           ONLY : write_scf,read_scf
  USE matrix_inversion
  USE symm_base,  ONLY : nsym, irt, t_rev, d1, d2, d3
  USE module_chargedens_fromfile, ONLY : print_symmetries,symmetrize_rhor 
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filename_rho,outdir,outdir_old,outdir_svd,filename_dmat
  INTEGER :: ios,ik,ig,npw,ibnd,ibnd_prime,ipol,is,iun
  LOGICAL :: needwf = .TRUE.
  LOGICAL :: mix_rho,change_rho 
  REAL(DP):: mixing_beta_chgdens,dens_real, dens_im
  REAL(DP), ALLOCATABLE :: rbecp(:,:), rhor_aux(:,:),lambda(:)
  COMPLEX(DP), ALLOCATABLE :: becp(:,:), becp_nc(:,:,:), be1(:,:), be2(:,:)
  COMPLEX(DP), ALLOCATABLE :: density_mat(:,:,:,:) !ib,ib',is,ik 
  COMPLEX(DP), ALLOCATABLE :: psic_nc_nbnd(:,:,:),psicpw_nbnd(:,:),evc_new(:,:),evc_dot_evc_new(:,:),evc_dot_evc_new_inv(:,:)
  COMPLEX(DP), ALLOCATABLE :: density_mat_ik(:,:),density_mat_ik_new(:,:) 
  LOGICAL :: exst
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  NAMELIST / chargedens_fromfile / outdir, prefix, change_rho, filename_rho, mixing_beta_chgdens, outdir_old, &
                                   mix_rho, filename_dmat
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  !
  CALL environment_start ( 'CHARGEDENS_FROMFILE' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filename_rho= ' '
  mix_rho = .false.
  mixing_beta_chgdens = 0.7d0
  outdir_old = ''
  change_rho = .false.
  filename_dmat = ''
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, chargedens_fromfile, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, intra_image_comm )
  IF (ios /= 0) CALL errore ('do_chargedens_fromfile', 'reading chargedens_fromfile namelist', abs (ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir,              ionode_id, intra_image_comm )
  CALL mp_bcast( prefix,               ionode_id, intra_image_comm )
  CALL mp_bcast( filename_rho,         ionode_id, intra_image_comm )
  CALL mp_bcast( mixing_beta_chgdens,  ionode_id, intra_image_comm )
  CALL mp_bcast( outdir_old,           ionode_id, intra_image_comm )
  CALL mp_bcast( mix_rho,              ionode_id, intra_image_comm )
  CALL mp_bcast( change_rho,           ionode_id, intra_image_comm )
  CALL mp_bcast( filename_dmat,        ionode_id, intra_image_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file_new ( needwf )
  !
  ! if needed print symmetry files
  if (len(trim(filename_dmat)) > 0) then
    call print_symmetries(filename_dmat) 
  endif
  !
  ! read density_mat from file
  !
  iun = find_free_unit() 
  ALLOCATE (rhor_aux(dfftp%nnr,nspin))
  rhor_aux(:,:) = 0.d0  
  if (change_rho) then
    ALLOCATE(density_mat(nbnd,nbnd,nks,nspin))
    ALLOCATE (psic_nc_nbnd(dfftp%nnr,nspin,nbnd))
    becsum(:,:,:) = 0.d0
    open(unit=iun,file=trim(filename_rho),action="read",status="old",iostat = ios)
    IF (ios /= 0) CALL errore ('do_chargedens_fromfile', 'reading densityMat.data', abs (ios) )
    do ik = 1,nks
      read(iun,*)
      do ibnd = 1, nbnd
        read(iun,*)
        do ibnd_prime = 1,nbnd
          read(iun,"(E23.16,x,E23.16)") dens_real, dens_im 
          density_mat(ibnd,ibnd_prime,ik,1) = cmplx(dens_real,dens_im,kind=DP)
        enddo
      enddo
    enddo
    density_mat(:,:,:,nspin) = density_mat(:,:,:,1)
    close(iun)
    ! here an up down reading is necessary if non collinear
    !
    if (nspin /= 1) CALL errore ('do_chargedens_fromfile', 'not implemented for 2 spins', abs (nspin) ) 
    IF (gamma_only) THEN
       ALLOCATE (rbecp(nkb,nbnd))
    ELSE
      IF (noncolin) THEN
         ALLOCATE (becp_nc(nkb,npol,nbnd))
         IF ( ANY(upf(1:ntyp)%has_so) ) THEN
           ALLOCATE(be1(nhm,2))
           ALLOCATE(be2(nhm,2))
         ENDIF
      ELSE
         ALLOCATE (becp(nkb,nbnd))
      ENDIF
    ENDIF
    current_spin = 1
    DO ik = 1, nks
      IF (lsda) current_spin = isk (ik)
      CALL read_collected_wfc ( restart_dir(), ik, evc )
      npw = ngk(ik)
      CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
      IF (gamma_only) THEN
         CALL calbec ( npw, vkb, evc, rbecp )
      ELSEIF (noncolin) THEN
         CALL calbec ( npw, vkb, evc, becp_nc )
      ELSE
         CALL calbec ( npw, vkb, evc, becp )
      ENDIF
      !
      DO ibnd = 1, nbnd
        IF (noncolin) THEN
          psic_nc = (0.d0,0.d0)
          DO ig = 1, npw
             psic_nc(dffts%nl(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
             psic_nc(dffts%nl(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
          ENDDO
          DO ipol=1,npol
             CALL invfft ('Wave', psic_nc(:,ipol), dffts)
             psic_nc_nbnd(:,ipol,ibnd) = psic_nc(:,ipol)
          ENDDO
        ELSE
          psic(1:dffts%nnr) = (0.d0,0.d0)
          DO ig = 1, npw
             psic (dffts%nl (igk_k(ig,ik) ) ) = evc (ig, ibnd)
          ENDDO
          IF (gamma_only) THEN
             DO ig = 1, npw
                psic (dffts%nlm(igk_k (ig,ik) ) ) = conjg(evc (ig, ibnd))
             ENDDO
          ENDIF
          CALL invfft ('Wave', psic, dffts)
          psic_nc_nbnd(:,1,ibnd) = psic
          psic_nc_nbnd(:,nspin,ibnd) = psic
        ENDIF
        !
        ! case of density_mat = identity
        !
        IF (.not.change_rho) THEN
          IF (noncolin) THEN
            DO ipol=1,npol
              rhor_aux(:,current_spin) = rhor_aux(:,current_spin) + & 
                                          dble(psic_nc(:,ipol) * wg (ibnd, ik) * conjg(psic_nc(:,ipol))) / omega
            ENDDO
          ELSE
            rhor_aux(:,current_spin) = rhor_aux(:,current_spin) + & 
                                          dble(psic * wg (ibnd, ik) * conjg(psic)) / omega 
          ENDIF
        ENDIF
      ENDDO
      ! change rho
      DO ibnd_prime = 1, nbnd
        DO ibnd = 1, nbnd
            DO is=1,nspin
              rhor_aux(:,current_spin) = rhor_aux(:,current_spin) + & 
                dble(psic_nc_nbnd(:,is,ibnd) * wk(ik) * density_mat(ibnd,ibnd_prime,ik,current_spin)  & 
                                                      * conjg(psic_nc_nbnd(:,is,ibnd_prime))) / omega 
               ! probably wrong treatment of spin
            ENDDO
        ENDDO
      ENDDO
    ENDDO
    IF (gamma_only) THEN
       DEALLOCATE(rbecp)
    ELSE
       IF (noncolin) THEN
          IF ( ANY(upf(1:ntyp)%has_so) ) THEN
             DEALLOCATE(be1)
             DEALLOCATE(be2)
          ENDIF
          DEALLOCATE(becp_nc)
       ELSE
          DEALLOCATE(becp)
       ENDIF
    ENDIF
    !
    deallocate(density_mat)
    !
    ! symmetrize rho
    !
    call symmetrize_rhor(rhor_aux)
    !
  endif
  ! mix if needed
  if (mix_rho) then
    if (.not. change_rho) CALL rho_g2r (dfftp, rho%of_g, rhor_aux)
    ! change outdir manually
    ! outdir_svd = tmp_dir
    ! tmp_dir = trimcheck(outdir_old)
    ! CALL read_scf ( rho, nspin, gamma_only )
    ! CALL rho_g2r (dfftp, rho%of_g, rho%of_r)
    ! rhor_aux(:,:) = (1.0d0 - mixing_beta_chgdens) * rho%of_r(:,:) + mixing_beta_chgdens * rhor_aux(:,:)
    ! tmp_dir = outdir_svd
    ! CALL read_scf ( rho, nspin, gamma_only )
    CALL open_mix_file( iunmix, 'mix_chargedens_fromfile', exst )
    IF ( my_pool_id == root_pool ) CALL mix_rho( rho, rhoin, &
            mixing_beta, dr2, tr2_min, iter, nmix, iunmix, conv_elec )
    !
  endif
  !
  ! ... bring rho(r) to G-space (use psic as work array)
  !
  rho%of_r(:,:) = rhor_aux(:,:)
  rho%of_g(:,:) = cmplx(0.d0,0.d0,kind=DP)
  call rho_r2g(dfftp, rho%of_r, rho%of_g )
  !
  if (change_rho .or. mix_rho) then
    call write_scf ( rho, nspin)
  endif
  deallocate(rhor_aux)
  !
  ! becsum term not computed (see local_dos.f90 line 255)
  !
  CALL environment_end ( 'CHARGEDENS_FROMFILE' )
  !
  CALL stop_pp
  !
  ! ... bring the unsymmetrized rho(r) to G-space (use psic as work array)
  !
  ! DO is = 1, nspin
  !    psic(1:dffts%nnr) = rho%of_r(1:dffts%nnr,is)
  !    psic(dffts%nnr+1:) = 0.0_dp
  !    CALL fwfft ('Rho', psic, dffts)
  !    rho%of_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
  !    rho%of_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
  ! END DO
  !
END PROGRAM do_chargedens_fromfile
