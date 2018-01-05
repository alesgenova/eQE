!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_cc (forcecc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE ions_base,            ONLY : atm
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE large_cell_base,      ONLY : alatl => alat, omegal => omega!, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_base,             ONLY : dfftl
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl
  USE gvecl,                ONLY : ngml => ngm, nll => nl, g_large => g!, gg, ngl, gl, igtongl
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  use fde,                  ONLY : linterlock, fde_dotsonlarge, rho_fde_large, &
                                    native_cell, reduced_cell, &
                                    rho_gauss_fde_large, rhog_gauss_fde_large, &
                                    do_fde, rho_gauss, rhog_gauss, &
                                    rho_gauss_fde, rhog_gauss_fde, &
                                    rho_core_fde, rhog_core_fde, &
                                    rho_core_fde_large, rhog_core_fde_large, &
                                    fde_nspin, fde_frag_nspin, fde_fake_nspin, rho_fde, &
                                    use_gaussians, fde_kin_is_nl, &
                                    nonlocalkernel, nonlocalkernel_fde
  use fde_routines
!
  !
  implicit none
  !
  !   first the dummy variable
  !
  real(DP) :: forcecc (3, nat)
  ! output: the local forces on atoms

  integer :: ipol, ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms
  integer :: is


  real(DP), allocatable :: vxc (:,:), rhocg (:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  real(DP)  ::  arg, fact

  real(DP), allocatable :: vts (:,:), vts_frag (:,:), aux_fde (:,:), auxl_fde(:,:), auxl(:,:)
  real(dp), allocatable :: gauss(:) !, erfvec
  logical :: conv_elec = .false.
  real(DP)  ::  ekin, ekin0, trash, alpha, sigma
  !
  integer, external :: atomic_number
  integer :: z_core
  logical :: do_cc

  !
  forcecc(:,:) = 0.d0
  if ( ANY ( upf(1:ntyp)%nlcc ) ) go to 15
  if ( do_fde .and. use_gaussians ) go to 15
  return
  !
15 continue
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc(dfftp%nnr,nspin) )
  !
    !
  ! If FDE, apply NLCC also to Ts
  if (do_fde ) then
    allocate ( vts(dfftp%nnr,nspin) )
    allocate ( vts_frag(dfftp%nnr,nspin) )
    vts = 0.d0
    vts_frag = 0.d0
    vxc = 0.d0
    !
    if (fde_frag_nspin /= fde_nspin) then
      allocate( aux_fde(dfftp%nnr, fde_nspin) )
      aux_fde(:,:) = 0.d0
      call fde_fake_nspin(.true.)
      if ( conv_elec .or. &
           fde_dotsonlarge .or. &
           fde_kin_is_nl) then
        allocate( auxl_fde(dfftl%nnr,fde_nspin) )
        !
        ! Ts 
        !
        call fde_kin(rho_fde_large, rho_gauss_fde_large, &
            rhog_gauss_fde_large, ekin, auxl_fde, native_cell, nonlocalkernel_fde) !dfftl, ngml, g_large, nll, omegal, .true.)
        do is = 1, fde_nspin
          call copy_pot_l2f(auxl_fde(:,is), aux_fde(:,is))
        enddo
        call v_join_spin (vts, aux_fde)
        !
        ! XC
        !
        aux_fde = 0.d0
        auxl_fde = 0.d0
        call v_xc_large(rho_fde_large, rho_core_fde_large, &
                      rhog_core_fde_large, etxc, trash, auxl_fde)
        do is = 1, fde_nspin
            call copy_pot_l2f(auxl_fde(:,is), aux_fde(:,is))
        enddo
        call v_join_spin (vxc, aux_fde)

        deallocate(auxl_fde)
      else
        !
        ! Ts
        !
        call fde_kin(rho_fde, rho_gauss_fde, rhog_gauss_fde, ekin, &
                        aux_fde, reduced_cell, nonlocalkernel_fde) !dfftp, ngm, g, nl, omega, .false.)
        call v_join_spin (vts, aux_fde)
        !
        ! XC
        !
        aux_fde = 0.d0
        call v_xc(rho_fde, rho_core_fde, rhog_core_fde, etxc, trash, aux_fde)
        call v_join_spin (vxc, aux_fde)
      endif
      call fde_fake_nspin(.false.)
      deallocate(aux_fde)
    else
      if ( conv_elec .or. & 
           fde_dotsonlarge .or. &
           fde_kin_is_nl) then
         allocate( auxl(dfftl%nnr,nspin) )
         !
         ! Ts
         !
         call fde_kin(rho_fde_large, rho_gauss_fde_large, &
                 rhog_gauss_fde_large, ekin, auxl, native_cell, nonlocalkernel_fde)!dfftl, ngml, g_large, nll, omegal, .true.)
         do is = 1, fde_nspin
           call copy_pot_l2f( auxl(:,is), vts(:,is) )
         enddo
         !
         ! XC
         !
         auxl = 0.d0
         call v_xc_large(rho_fde_large, rho_core_fde_large, &
                            rhog_core_fde_large, etxc, trash, auxl)
         do is = 1, fde_nspin
           call copy_pot_l2f( auxl(:,is), vxc(:,is) )
         enddo
         deallocate(auxl)
      else
         ! 
         ! Ts
         !
         call fde_kin(rho_fde, rho_gauss_fde, rhog_gauss_fde, ekin, &
                              vts, reduced_cell, nonlocalkernel_fde)!dfftp, ngm, g, nl, omega, .false.)
         !
         ! XC
         !
         call v_xc(rho_fde, rho_core_fde, rhog_core_fde, etxc, trash, vxc)
      endif
    endif
  
    ! V_Ts[\rho_I]
    vts_frag = 0.d0
    call fde_kin(rho, rho_gauss, rhog_gauss, ekin0, vts_frag, reduced_cell, nonlocalkernel)!&
                          !dfftp, ngm, g, nl, omega, .false.)
    !
    vxc = vxc + vts - vts_frag
    !
    deallocate( vts_frag, vts )
  
    else
      call v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
    endif
  
  !
  psic=(0.0_DP,0.0_DP)
  if (nspin == 1 .or. nspin == 4) then
     do ir = 1, dfftp%nnr
        psic (ir) = vxc (ir, 1)
     enddo
  else
     do ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     enddo
  endif
  deallocate (vxc)
  CALL fwfft ('Rho', psic, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  allocate ( rhocg(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then

        call drhoc (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, rgrid(nt)%r,&
             rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        do na = 1, nat
           if (nt.eq.ityp (na) ) then
              do ig = gstart, ngm
                 arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                      + g (3, ig) * tau (3, na) ) * tpi
                 do ipol = 1, 3
                    forcecc (ipol, na) = forcecc (ipol, na) + tpiba * omega * &
                         rhocg (igtongl (ig) ) * CONJG(psic (dfftp%nl (ig) ) ) * &
                         CMPLX( sin (arg), cos (arg) ,kind=DP) * g (ipol, ig) * fact
                 enddo
              enddo
           endif
        enddo
     endif
  enddo
  !
  call mp_sum(  forcecc, intra_bgrp_comm )
  !
  deallocate (rhocg)
  !
  return
end subroutine force_cc
