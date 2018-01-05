!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_rhoc
  !-----------------------------------------------------------------------
  !
  !    This routine computes the core charge on the real space 3D mesh
  !
  !
  USE io_global, ONLY : stdout
  USE io_global, ONLY : ionode
  USE kinds,     ONLY : DP
  USE atom,      ONLY : msh, rgrid
  USE uspp_param,ONLY : upf
  USE ions_base, ONLY : ntyp => nsp
  USE ions_base, ONLY : atm
  USE cell_base, ONLY : omega, tpiba2
  USE ener,      ONLY : etxcc
  USE fft_base,  ONLY : dfftp
  USE fft_base,  ONLY : dfftl, grid_gather, grid_scatter, grid_scatter_large
  USE fft_interfaces,ONLY : invfft
  USE fft_interfaces,ONLY : fwfft
  USE gvect,     ONLY : ngm, ngl, gl, igtongl
  USE scf,       ONLY : rho_core, rhog_core, scf_type
  USE lsda_mod,  ONLY : nspin
  USE vlocal,    ONLY : strf
  USE control_flags, ONLY : gamma_only
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : inter_fragment_comm
  USE constants, ONLY : tpi
  USE fde
  use fde_routines

  !
  implicit none
  !
  real(DP), parameter :: eps = 1.d-10

  complex(DP) , allocatable :: aux (:)
  ! used for the fft of the core charge

  real(DP) , allocatable ::  rhocg(:)
  ! the radial fourier transform
  real(DP) ::  rhoima, rhoneg, rhorea
  real(DP) ::  corecharge
  ! used to check the core charge
  real(DP) ::  vtxcc
  ! dummy xc energy term
  type(scf_type) :: dum

  real(dp), allocatable :: raux(:), gauss(:), rauxl(:)
  complex(dp) , allocatable :: gaux(:), gauxl(:)
  ! FDE auxiliary gathered grid array

  real(dp) :: sigma, alpha, rmatch, sigma2, arg
  integer :: z_core
  integer, external :: atomic_number
  real(DP), external :: qe_erf

  integer :: ir, nt, ng
  ! counter on mesh points
  ! counter on atomic types
  ! counter on g vectors

  etxcc = 0.0_DP
  if ( ANY( upf(1:ntyp)%nlcc ) ) goto 10
  
  rhog_core(:) = 0.0_DP
  rho_core(:)  = 0.0_DP

  if (do_fde) then
     rhog_gauss(:) = 0._DP
     rho_gauss(:)  = 0._DP
  endif


  return

10 continue
  allocate (aux( dfftp%nnr))    
  allocate (rhocg( ngl))    
  aux (:) = (0.0_DP, 0.0_DP)
  !
  !    the sum is on atom types
  !    1) the sum is on atom types
  !    2) strf ind g-space are on SMALL CELL
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then
        !
        !     drhoc compute the radial fourier transform for each shell of g vec
        !
        call drhoc (ngl, gl, omega, tpiba2, msh (nt), rgrid(nt)%r, &
             rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        !
        !     multiply by the structure factor and sum
        !
        do ng = 1, ngm
           aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + strf(ng,nt) * rhocg(igtongl(ng))
        enddo
     endif
  enddo
  if (gamma_only) then
     do ng = 1, ngm
        aux(dfftp%nlm(ng)) = CONJG(aux(dfftp%nl (ng)))
     end do
  end if
  !
  rhog_core(:) = aux(dfftp%nl(:))
  !
  !   the core charge in real space
  !
  CALL invfft ('Rho', aux, dfftp)
  !
  !    test on the charge and computation of the core energy
  !
  rhoneg = 0.d0
  rhoima = 0.d0
  do ir = 1, dfftp%nnr
     rhoneg = rhoneg + min (0.d0,  DBLE (aux (ir) ) )
     rhoima = rhoima + abs (AIMAG (aux (ir) ) )
     rho_core(ir) =  DBLE (aux(ir))
     !
     ! NOTE: Core charge is computed in reciprocal space and brought to real
     ! space by FFT. For non smooth core charges (or insufficient cut-off)
     ! this may result in negative values in some grid points.
     ! Up to October 1999 the core charge was forced to be positive definite.
     ! This induces an error in the force, and probably stress, calculation if
     ! the number of grid points where the core charge would be otherwise neg
     ! is large. The error disappears for sufficiently high cut-off, but may be
     ! rather large and it is better to leave the core charge as it is.
     ! If you insist to have it positive definite (with the possible problems
     ! mentioned above) uncomment the following lines.  SdG, Oct 15 1999
     !
     !         rhorea = max ( DBLE (aux (ir) ), eps)
     !         rho_core(ir) = rhorea
     !
  enddo
  rhoneg = rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  rhoima = rhoima / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  !
  call mp_sum(  rhoneg, intra_bgrp_comm )
  call mp_sum(  rhoima, intra_bgrp_comm )
  !
  IF (rhoneg < -1.0d-6 .OR. rhoima > 1.0d-6) &
       WRITE( stdout, '(/5x,"Check: negative/imaginary core charge=",2f12.6)')&
       rhoneg, rhoima
  !
  ! calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
  ! The term was present in previous versions of the code but it shouldn't
  !
  !   call create_scf_type(dum)
  !   dum%of_r(:,:) = 0.0_DP
  !   dum%of_g(:,:) = (0.0_DP, 0.0_DP)
  !   
  !   call v_xc( dum, rho_core, rhog_core, etxcc, vtxcc, aux )
  ! 
  !   call destroy_scf_type(dum)
  !   WRITE( stdout, 9000) etxcc
  !   WRITE( stdout,  * ) 'BEWARE it will be subtracted from total energy !'
  !

  if (do_fde .and. use_gaussians) then
     ! now do the same but for the compensating gaussian charges
     aux (:) = (0._DP, 0._DP)
     do nt = 1, ntyp

        allocate (gauss(1:msh(nt)))
        z_core = (atomic_number(atm(nt)) - upf(nt)%zp)
        sigma = 0.3d0
        alpha = z_core/(sigma*tpi)**(3.d0/2.d0)

        if (atomic_number(atm(nt)) == 1) alpha = 0.0      

        gauss(1:msh(nt)) = alpha*exp(-rgrid(nt)%r(1:msh(nt))**2.d0/(2.d0*sigma)) 

        ! see above for what this routine does
        call drhoc (ngl, gl, omega, tpiba2, msh (nt), rgrid(nt)%r, &
             rgrid(nt)%rab, gauss, rhocg)
        
        !if (nt == 1) then
        !do ir = 1 , msh(nt) - 1
        !   write(666,*) rgrid(nt)%r(ir), gauss(ir), rgrid(nt)%r(ir+1) - rgrid(nt)%r(ir)
        !enddo
        !endif

        deallocate (gauss)
        
        !     multiply by the structure factor and sum
        do ng = 1, ngm
           aux(nl(ng)) = aux(nl(ng)) + strf(ng,nt) * rhocg(igtongl(ng))
        enddo
     enddo

     if (gamma_only) then
        do ng = 1, ngm
           aux(nlm(ng)) = CONJG(aux(nl (ng)))
        end do
     end if

     rhog_gauss(:) = aux(nl(:))

     !   the core charge in real space
     CALL invfft ('Dense', aux, dfftp)


  !    test on the charge and computation of the core energy
  !
  rhoneg = 0.d0
  rhoima = 0.d0
  corecharge = 0.d0
  rho_gauss = 0.d0
  do ir = 1, dfftp%nnr
     rhoneg = rhoneg + min (0.d0,  DBLE (aux (ir) ) )
     rhoima = rhoima + abs (AIMAG (aux (ir) ) )
     rho_gauss(ir) =  DBLE (aux(ir))
     corecharge = corecharge + rho_gauss(ir)
     !
  enddo
  !
  call mp_sum(  rhoneg, intra_bgrp_comm )
  call mp_sum(  rhoima, intra_bgrp_comm )
  call mp_sum(  corecharge, intra_bgrp_comm )
  !
  rhoneg = rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  rhoima = rhoima / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  corecharge = omega*corecharge / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  !
  IF (rhoneg < -1.0d-6 .OR. rhoima > 1.0d-6 .or. corecharge > 0.d0) &
       WRITE(stdout,*) 'WARNING:'
       WRITE( stdout, '(/5x,"Check: negative/imaginary GAUSS core charge=",4f12.6 i8 )')&
       corecharge, rhoneg, rhoima, omega, (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     !do ir = 1, dfftp%nnr
     !   rho_gauss(ir) =  DBLE (aux(ir))
     !enddo
  endif ! fde

  deallocate (rhocg)
  deallocate (aux)

  if (do_fde ) then
     if(.not.use_gaussians) then
        ! copy core to gaussians otherwise the 
        ! nonadditive kin goes without core density
        ! if use_gaussians=.false.
        call copy_rho_core_to_rho_gauss
     endif
  endif

  ! generate rho_core_fde
  if (do_fde) then
     !
     rhog_core_fde(:) = rhog_core(:)
     rho_core_fde(:) = rho_core(:)
     !
     !
     !
     ! Allocate gathered arrays
     if (ionode) then
       allocate(raux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
       if (linterlock) allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x) )
     endif
     !
     !
     ! Allocate scattered array (distributed matrices)
     allocate( gaux(dfftp%nnr) )
     if (linterlock) allocate( gauxl(dfftl%nnr) )
     !
     !
     !
     ! gather rho_core_fde in the small grid in raux
     call grid_gather(rho_core_fde, raux)
     !
     !
     !
     !
     ! take gathered raux array on small grid 
     !   and copy it on large grid in rauxl
     if (linterlock) then
       ! gather on large grid
       if (ionode) then
         rauxl = 0.d0
         ! copy gathered array on large grid
         rauxl(f2l(:)) = raux(:)
         ! sum across subsystems to get total core density
         ! nota bene: small cells overlap - but each ionode
         !            contain only the core charge of its 
         !            associated subsystem.
         call mp_sum(rauxl, inter_fragment_comm)
         !
         ! all subsystems have rauxl, and can be 
         ! copied back to raux.
         !
         ! now raux has TOTAL core charge carved in the
         ! small cell
         raux = 0.d0
         raux(:) = rauxl(f2l(:))
       endif
       !
       !
       ! scatter rauxl on large grid, thus
       ! rho_core_fde_large is on large grid (scattered)
       call grid_scatter_large(rauxl, rho_core_fde_large)
       !
       !
       ! generate FFT of rho_core_fde_large
       ! scattered and on large grid g-space
       gauxl(:) = cmplx(rho_core_fde_large(:), 0.d0, kind=dp)
       call fwfft ('Custom', gauxl, dfftl)
       rhog_core_fde_large(1:ngml) = gauxl(nll(1:ngml))
       !
       ! scatter raux on small grid
       call grid_scatter(raux, rho_core_fde)
       !
       !
       ! generate FFT of rho_core_fde
       ! scattered and on small grid g-space
       gaux(:) = cmplx(rho_core_fde(:), 0.d0, kind=dp)
       call fwfft ('Dense', gaux, dfftp)
       rhog_core_fde(1:ngm) = gaux(nl(1:ngm))
       !
     else
       ! no linterlock
       if (ionode) call mp_sum(raux, inter_fragment_comm)
       call grid_scatter(raux, rho_core_fde)
       call c_grid_gather_sum_scatter(rhog_core_fde)
     endif
     !
     ! Same as above but with rho_gauss. 
     !
     rhog_gauss_fde(:) = rhog_gauss(:)
     rho_gauss_fde(:) = rho_gauss(:)
     !
     call grid_gather(rho_gauss_fde, raux)

     if (linterlock) then
       if (ionode) then
         rauxl = 0.d0
         rauxl(f2l(:)) = raux(:)
         call mp_sum(rauxl, inter_fragment_comm)
         raux = 0.d0
         raux(:) = rauxl(f2l(:))
       endif
       !
       call grid_scatter_large(rauxl, rho_gauss_fde_large)
       gauxl(:) = cmplx(rho_gauss_fde_large(:), 0.d0, kind=dp)
       call fwfft ('Custom', gauxl, dfftl)
       rhog_gauss_fde_large(1:ngml) = gauxl(nll(1:ngml))
       !
       call grid_scatter(raux, rho_gauss_fde)
       gaux(:) = cmplx(rho_gauss_fde(:), 0.d0, kind=dp)
       call fwfft ('Dense', gaux, dfftp)
       rhog_gauss_fde(1:ngm) = gaux(nl(1:ngm))
       !
     else
       if (ionode) call mp_sum(raux, inter_fragment_comm)
       call grid_scatter(raux, rho_gauss_fde)
       call c_grid_gather_sum_scatter(rhog_gauss_fde)
     endif
     !
     if (ionode) then
       deallocate(raux)
       if (linterlock) deallocate(rauxl)
     endif
     !
     deallocate( gaux )
     if (linterlock) deallocate(gauxl)
     !
  endif

  !
  return

  ! 9000 format (5x,'core-only xc energy         = ',f15.8,' Ry')

end subroutine set_rhoc

