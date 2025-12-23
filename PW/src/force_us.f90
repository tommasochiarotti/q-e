!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE force_us( forcenl )
  !----------------------------------------------------------------------------
  !! The nonlocal potential contribution to forces.
  !
  USE io_global,            ONLY : ionode, ionode_id
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only, offload_type
  USE cell_base,            ONLY : tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk, igk_k, wk
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, vkb, qq_at, deeq, qq_so, deeq_nc, ofsbeta
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin, lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : calbec, becp, bec_type, &
                                   allocate_bec_type, deallocate_bec_type, &
                                   allocate_bec_type_acc, deallocate_bec_type_acc
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
  USE mp_world,         ONLY : world_comm
  USE mp,                   ONLY : mp_sum, mp_bcast, mp_barrier
  USE uspp_init,            ONLY : init_us_2
  USE pw_restart_new,       ONLY : read_collected_wfc
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: forcenl(3,nat)
  !! the nonlocal contribution
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)   ! contains g*|beta>
  TYPE(bec_type) :: becd                  ! contains <dbeta|psi>
  COMPLEX(DP) :: deff_nc
  REAL(DP) :: deff, fnl
  INTEGER :: npw, ik, ipol, ig, na, na_s, na_e, mykey
  INTEGER :: nt, ibnd, ibnd1, ibnd2, ibnd3, nhnt, ih, jh, ijkb0, ikb, jkb, is, js, ijs
  INTEGER :: ios, iun
  COMPLEX(DP), ALLOCATABLE :: density_mat(:,:,:)
  REAL(DP) :: dens_real, dens_im
  INTEGER, EXTERNAL :: find_free_unit
  COMPLEX(DP), ALLOCATABLE :: evc_old(:,:), overlap(:,:)


  !
  forcenl(:,:) = 0.D0

  ! READ OCCMAT
  ALLOCATE(density_mat(nbnd,nbnd,nks))
  density_mat = (0.0_dp, 0.0_dp)
  iun = find_free_unit()
  !if (ionode) then
  open(unit=iun,file=trim("densityMat.dat"),action="read",status="old",iostat = ios)
  do ik = 1,nks
    if (ios == 0) then
      read(iun,*)
      do ibnd = 1, nbnd
        read(iun,*)
        do ibnd1 = 1,nbnd
          read(iun,"(E23.16,x,E23.16)") dens_real, dens_im
          density_mat(ibnd,ibnd1,ik) = cmplx(dens_real,dens_im,kind=DP)
          ! density_mat(ibnd1,ibnd,ik) = conjg(cmplx(dens_real,dens_im,kind=DP))
          ! if (ibnd == ibnd1) then
          !   density_mat(ibnd, ibnd,ik) = cmplx(dens_real, 0)
          ! endif
          !density_mat(ibnd1,ibnd,ik) = cmplx(dens_real,-dens_im,kind=DP)
          !density_mat(ibnd, ibnd1, ik) = (1._dp, 0._dp)
        enddo
      enddo
    else
      do ibnd = 1, nbnd
        density_mat(ibnd, ibnd, ik) = wg(ibnd,ik) / wk(ik)
      enddo
    endif
  enddo
  close(iun)
  !endif
  !CALL mp_bcast( density_mat, ionode_id, intra_bgrp_comm )
  !CALL mp_bcast( density_mat, ionode_id, inter_pool_comm )
  !
  CALL allocate_bec_type_acc( nkb, nbnd, becp, intra_bgrp_comm )
  CALL allocate_bec_type_acc( nkb, nbnd, becd, intra_bgrp_comm )
  ALLOCATE( vkb1(npwx,nkb) )
  !$acc data create(vkb1)
  !
  ! ... the forces are summed over K-points
  !
  DO ik = 1, nks
     npw = ngk(ik)
     ALLOCATE(overlap(nbnd, nbnd))
     ALLOCATE(evc_old(npwx, nbnd))
     overlap = (0.0_dp, 0.0_dp)
     evc_old = (0.0_dp, 0.0_dp)
     CALL read_collected_wfc ( "./results_old/SrVO3.save/", ik, evc_old )
     !
     IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
     !$acc data present (evc, vkb, becp)
     CALL calbec( offload_type, npw, vkb, evc_old, becp )
     !$acc end data
     !
    ! density_mat(:,:,ik) = (density_mat(:,:,ik) + TRANSPOSE(CONJG(density_mat(:,:,ik))))/2.
     DO ipol = 1, 3
        DO jkb = 1, nkb
           DO ig = 1, npw
              ! vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk_k(ig,ik))
              vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * (g(ipol,igk_k(ig,ik)) + xk(ipol,ik))
           ENDDO
        ENDDO
        CALL calbec( offload_type, npw, vkb1, evc_old, becd )
        ! becp = <beta|psi>, becd = <dbeta/dG_ipol|psi>
        DO na = 1, nat
        fnl = 0.0_dp
        nt = ityp(na)
        nhnt = nh(nt)
        ijkb0 = ofsbeta(na)
        DO ibnd = 1, nbnd
          DO ih = 1, nhnt
            DO jh = 1, nhnt
              ikb = ijkb0 + ih
              jkb = ijkb0 + jh
              deff = deeq(ih,jh,na,current_spin) !- et(ibnd,ik) * qq_at(ih,jh,na)
              ! fnl = fnl + wg(ibnd,ik) * deff *  &
              !      DBLE(CONJG(becp%k(ikb,ibnd)) * becd%k(jkb,ibnd))
              DO ibnd1 = 1, nbnd
                fnl = fnl + wk(ik) * deff *  &
                DBLE(CONJG(becp%k(ikb,ibnd)) * &
                becd%k(jkb,ibnd1) * density_mat(ibnd1, ibnd, ik) + &
                CONJG(becd%k(ikb,ibnd)) * &
                becp%k(jkb,ibnd1) * density_mat(ibnd1, ibnd, ik))/2.
              END DO
            END DO
          END DO
        END DO
        ! factor 2 from Ry a.u. (e^2=2)? tpiba from k+G, minus sign
        forcenl(ipol,na) = forcenl(ipol,na) - 2.0_dp * tpiba* fnl
        END DO
     ENDDO
     DEALLOCATE(evc_old)
     DEALLOCATE(overlap)
  ENDDO
  !
  !$acc end data
  DEALLOCATE( vkb1 )
  CALL deallocate_bec_type_acc( becd )
  CALL deallocate_bec_type_acc( becp )
  !
  RETURN
  !
END SUBROUTINE force_us
