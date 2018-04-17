! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This module contains the routines implementing the Frozen Density Embedding
! method [1,2]. The FDE method is implemented thanks to image parallelization,
! and a serial version compilation is not possible. The SCF is performed on
! each fragment in parallel, and the non-additive potentials and energies
! are calculated at every SCF step.
!
! Coded by D. Ceresoli and A. Genova, during the April 2013 PRG hackathon
! in Rutgers-Newark (http://www.michelepavanello.com/pdf/hackathon.pdf)
!
! [1] Genova, Ceresoli & Pavanello, J. Chem. Phys. (2014)
! [2] see: $QEDIR/FDE-Notes/fde-equations.tex

module fde_routines
implicit none

contains
!--------------------------------------------------------------------------
SUBROUTINE allocate_fde
  !--------------------------------------------------------------------------
  !
  ! allocate memory for FDE calculations
  !
  use ions_base,                only : nat, ityp, ntyp => nsp, if_pos, fixatom
  use lsda_mod,                 only : nspin
  use gvect,                    only : ngm
  use gvecl,                    only : ngml => ngm
  use mp,                       only : mp_sum, mp_barrier, mp_bcast
  use mp_images,                only : nimage, inter_fragment_comm, intra_image_comm
  use klist,                    only : nelec
  use scf,                      only : create_scf_type
  use scf_large,                only : create_scf_type_large => create_scf_type
  use fft_base,                 only : dfftp, dfftl
  use io_global,                only : ionode, ionode_id, stdout
  use fde
  implicit none
  integer, allocatable :: nat_(:)
  integer :: i
  integer :: ierr

  ! structure factor and density
  allocate( strf_fde(ngm,ntyp) , stat=ierr)
  call fde_fake_nspin(.true.)
  call create_scf_type(rho_fde)
  rho_fde%is_fde=.true.
  call create_scf_type_large(rho_fde_large)
  rho_fde_large%is_fde=.true.
  call fde_fake_nspin(.false.)

  ! the core charge (for NLCC)
  allocate( rho_core_fde(dfftp%nnr), rhog_core_fde(ngm) , stat=ierr)
  allocate( rho_gauss(dfftp%nnr), rhog_gauss(ngm) , stat=ierr)
  allocate( rho_gauss_fde(dfftp%nnr), rhog_gauss_fde(ngm) , stat=ierr)

  ! arrays mapping the space from the large grid to the fragment grid and viceversa
  if (linterlock) then
    allocate( f2l(dfftp%nr1*dfftp%nr2*dfftp%nr3) )! , l2f(dfftp%nr1*dfftp%nr2*dfftp%nr3) )
    allocate( strf_fde_large(ngml,ntyp) )
    allocate( strf_large(ngml,ntyp) )
    allocate( rho_core_fde_large(dfftl%nnr), rhog_core_fde_large(ngml) , stat=ierr)
    allocate( rho_gauss_fde_large(dfftl%nnr), rhog_gauss_fde_large(ngml) , stat=ierr)
    !allocate( tau_large(3,nat) )
    !allocate( pot_wall(dfftp%nnr) )
    !call potential_wall(dfftp, pot_wall)
  endif

  ! get total number of atoms and number of fixed atoms in md/relax
  nat_fde = nat
  if ( ionode ) call mp_sum(nat_fde, inter_fragment_comm)
  call mp_bcast( nat_fde, ionode_id, intra_image_comm )
  allocate(tau_fde(3,nat_fde), ityp_fde(nat_fde), nat_(nfragments), atm_offset(nfragments), stat=ierr)
  allocate(force_fde(3,nat_fde), if_pos_fde(3,nat_fde), stat=ierr)
  !
  fixatom_fde = fixatom
  if ( ionode ) call mp_sum(fixatom_fde, inter_fragment_comm)
  call mp_bcast( fixatom_fde, ionode_id, intra_image_comm )

  ! get the atom offset for each fragment
  nat_(:) = 0
  nat_(currfrag) = nat
  if (ionode) call mp_sum(nat_, inter_fragment_comm)
  call mp_bcast( nat_, ionode_id, intra_image_comm )

  atm_offset(1) = 1
  do i = 2, nfragments
     atm_offset(i) = atm_offset(i-1) + nat_(i-1)
  enddo

  ! scatter ityp
  i = atm_offset(currfrag)
  ityp_fde(:) = 0
  ityp_fde(i:i+nat-1) = ityp(1:nat)
  if (ionode) call mp_sum(ityp_fde, inter_fragment_comm)
  call mp_bcast( ityp_fde, ionode_id, intra_image_comm )

  ! scatter if_pos
  i = atm_offset(currfrag)
  if_pos_fde(1:3,:) = 0
  if_pos_fde(1:3,i:i+nat-1) = if_pos(1:3,1:nat)
  if (ionode)  call mp_sum(if_pos_fde, inter_fragment_comm)
  call mp_bcast( if_pos_fde, ionode_id, intra_image_comm )

  ! get total number of electrons
  nelec_fde = nelec
  if (ionode) call mp_sum(nelec_fde, inter_fragment_comm)
  call mp_bcast( nelec_fde, ionode_id, intra_image_comm )

END SUBROUTINE allocate_fde

!
!--------------------------------------------------------------------------
SUBROUTINE gen_simulation_cells()
  !--------------------------------------------------------------------------
  use fde_types, only : new => new_cell
  use fde, only : tau_large, native_cell, reduced_cell, &
                  strf_large, strf_fde, strf_fde_large
  use ions_base, only : tau, nat, nsp, ityp
  use cell_base, only : at, &
                        bg, &
                        ibrav, &
                        alat, &
                        omega, &
                        celldm, &
                        tpiba, &
                        tpiba2
  use large_cell_base, only : l_at => at, &
                              l_bg =>  bg, &
                              l_ibrav => ibrav, &
                              l_alat =>  alat, &
                              l_omega =>  omega, &
                              l_celldm =>  celldm, &
                              l_tpiba => tpiba, &
                              l_tpiba2 => tpiba2
  use fft_base, only : dfftp, dfftl
  use vlocal, only : strf
  use gvect, only : ngm, &
                    ngm_g, &
                    ngl, &
                    ngmx, &
                    ecutrho, &
                    gcutm, &
                    gstart, &
                    gg, &
                    gl, &
                    g, &
                    igtongl, &
                    mill, &
                    ig_l2g, &
                    !sortedig_l2g, &
                    mill_g, &
                    eigts1, &
                    eigts2, &
                    eigts3
  use gvecl, only : l_ngm => ngm, &
                    l_ngm_g => ngm_g, &
                    l_ngl => ngl, &
                    l_ngmx => ngmx, &
                    l_ecutrho => ecutrho, &
                    l_gcutm => gcutm, &
                    l_gstart => gstart, &
                    l_gg => gg, &
                    l_gl => gl, &
                    l_g => g, &
                    l_igtongl => igtongl, &
                    l_mill => mill, &
                    l_ig_l2g => ig_l2g, &
                    !l_sortedig_l2g => sortedig_l2g, &
                    l_mill_g => mill_g, &
                    l_eigts1 => eigts1, &
                    l_eigts2 => eigts2, &
                    l_eigts3 => eigts3

! make the subsystem-based cell. use strf=strf, strf_fde=strf_fde
  call new(reduced_cell, &
           at=at, bg=bg, omega=omega, alat=alat, celldm=celldm, &
           tpiba=tpiba, tpiba2=tpiba2, &
           nat=nat, tau=tau, nsp=nsp, ityp=ityp, &
           ngm=ngm, ngm_g=ngm_g, ngl=ngl, ngmx=ngmx, &
           ecutrho=ecutrho, gcutm=gcutm, &
           strf=strf, strf_fde=strf_fde, &
           dfftp=dfftp, nl=dfftp%nl, nlm=dfftp%nlm, gstart=gstart, &
           gg=gg, gl=gl, g=g, igtongl=igtongl, mill=mill, &
           ig_l2g=ig_l2g, mill_g=mill_g, &
           eigts1=eigts1, eigts2=eigts2, eigts3=eigts3, &
           is_native_cell=.false.)

! make the subsystem-based cell. use strf=strf_large, strf_fde=strf_fde_large
  call new(native_cell, &
           at=l_at, bg=l_bg, omega=l_omega, alat=l_alat, celldm=l_celldm, &
           tpiba=l_tpiba, tpiba2=l_tpiba2, &
           nat=nat, tau=tau_large, nsp=nsp, ityp=ityp, &
           ngm=l_ngm, ngm_g=l_ngm_g, ngl=l_ngl, ngmx=l_ngmx, &
           ecutrho=l_ecutrho, gcutm=l_gcutm, &
           strf=strf_large, strf_fde=strf_fde_large, &
           dfftp=dfftl, nl=dfftl%nl, nlm=dfftl%nlm, gstart=l_gstart, &
           gg=l_gg, gl=l_gl, g=l_g, igtongl=l_igtongl, mill=l_mill, &
           ig_l2g=l_ig_l2g, mill_g=l_mill_g, &
           eigts1=l_eigts1, eigts2=l_eigts2, eigts3=l_eigts3, &
           is_native_cell=.true.)

END SUBROUTINE gen_simulation_cells

!--------------------------------------------------------------------------
SUBROUTINE gather_coordinates_nobcast(vec_frag, vec)
  !--------------------------------------------------------------------------
  !
  ! gather coordinates/force from fragments
  !
  use ions_base,                only : nat
  use mp,                       only : mp_sum, mp_bcast
  use mp_images,                only : nimage, inter_fragment_comm, intra_image_comm
  use io_global,                only : ionode, ionode_id
  use fde
  implicit none
  real(dp), intent(in) :: vec_frag(3,nat)
  real(dp), intent(out) :: vec(3,nat_fde)
  integer :: i

  i = atm_offset(currfrag)
  vec(:,:) = 0.d0
  vec(1:3,i:i+nat-1) = vec_frag(1:3,1:nat)
  if (ionode)  call mp_sum(vec, inter_fragment_comm)
!  call mp_bcast( vec, ionode_id, intra_image_comm )

END SUBROUTINE gather_coordinates_nobcast


!--------------------------------------------------------------------------
SUBROUTINE gather_coordinates(vec_frag, vec)
  !--------------------------------------------------------------------------
  !
  ! gather coordinates/force from fragments
  !
  use ions_base,                only : nat
  use mp,                       only : mp_sum, mp_bcast
  use mp_images,                only : nimage, inter_fragment_comm, intra_image_comm
  use io_global,                only : ionode, ionode_id
  use fde
  implicit none
  real(dp), intent(in) :: vec_frag(3,nat)
  real(dp), intent(out) :: vec(3,nat_fde)
  integer :: i

  i = atm_offset(currfrag)
  vec(:,:) = 0.d0
  vec(1:3,i:i+nat-1) = vec_frag(1:3,1:nat)
  if (ionode)  call mp_sum(vec, inter_fragment_comm)
  call mp_bcast( vec, ionode_id, intra_image_comm )

END SUBROUTINE gather_coordinates


!--------------------------------------------------------------------------
SUBROUTINE scatter_coordinates(vec, vec_frag)
  !--------------------------------------------------------------------------
  !
  ! scatter coordinates/forces to fragments
  !
  use ions_base,                only : nat
  use fde
  implicit none
  real(dp), intent(in) :: vec(3,nat_fde)
  real(dp), intent(out) :: vec_frag(3,nat)
  integer :: i

  i = atm_offset(currfrag)
  vec_frag(1:3,1:nat) = vec(1:3,i:i+nat-1)

END SUBROUTINE scatter_coordinates


!----------------------------------------------------------------------------
SUBROUTINE fde_summary
  !--------------------------------------------------------------------------
  !
  ! ... Print a short summary of FDE fragments
  !
  use ions_base,                only : tau, atm, nat
  use large_cell_base,          only : atl => at, bgl => bg
  use io_global,                only : stdout
  use fde
  use input_parameters,         only : fde_kin_funct
  implicit none
  integer :: ifrag, i, na, ipol, na_

  call gather_coordinates(tau_large, tau_fde)
!  if (linterlock) call reg_tau(tau, nat, atl, bgl)
  write(stdout,*)
  write(stdout,'(5X,''===== FROZEN DENSITY EMBEDDING METHOD ======================'')')
  write(stdout,'(5X,''Number of fragments  :'',I4)') nfragments
  write(stdout,'(5X,''Total number of atoms:'',I4)') nat_fde
  write(stdout,'(5X,''Current fragment     :'',I4)') currfrag
  if ( fde_fat .and. .false. ) then
       write(stdout,'(5X,"FDE SCF Type         : Freeze and Thaw")')
  elseif ( .false. ) then
       write(stdout,'(5X,"FDE SCF Type         : Regular Scf")')
  endif
  if ( fde_init_rho ) then
       write(stdout,'(5X,"Initial FDE Density  : Isolated KS")')
  else
       write(stdout,'(5X,"Initial FDE Density  : Atomic")')
  endif
  !write(stdout,'(5X,''Use gaussians        : '',L1)') use_gaussians
  write(stdout,'(5X,''Kinetic energy funct.: '',A10)') fde_kin_funct
  write(stdout,'(5X,''Fragment Charge      :'',F6.2)') fde_frag_charge
  write(stdout,'(5X,''Total    Charge      :'',F6.2)') fde_total_charge
  if ( fde_frag_nspin == 1 ) then
     write(stdout,'(5X,''Fragment Spin        : Closed Shell'')')
  else
     write(stdout,'(5X,''Fragment Spin        : Open Shell'')')
  endif
  if ( fde_nspin == 1 ) then
     write(stdout,'(5X,''Total    Spin        : Closed Shell'')')
  else
     write(stdout,'(5X,''Total    Spin        : Open Shell'')')
  endif
  do ifrag = 1, nfragments
     write(stdout,'(5X,''Fragment #'',I4,'':'')') ifrag
     i = atm_offset(ifrag)
     if (ifrag < nfragments) then
        na_ = atm_offset(ifrag+1)-atm_offset(ifrag)
     else
        na_ = nat_fde - atm_offset(ifrag) + 1
     endif
     write(stdout,'(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
             (na, atm(ityp_fde(na+i-1)), na+i-1, (tau_fde(ipol,na+i-1), ipol=1,3), na=1,na_)
  enddo
  write(stdout,'(5X,''===== FROZEN DENSITY EMBEDDING METHOD ======================'')')
  write(stdout,*)

END SUBROUTINE fde_summary


!----------------------------------------------------------------------------
SUBROUTINE fde_nonadditive(rho, rho_fde, rho_core, rhog_core, rho_core_fde, &
                           rhog_core_fde, rho_gauss, rhog_gauss, rho_gauss_fde, &
                           rhog_gauss_fde, ekin, ekin0, &
                           etxc, etxc0, v)
  !----------------------------------------------------------------------------
  !
  ! ... Add non-additive potential, kinetic and XC, and caculate energies
  !
  use kinds,                    only : dp
  use constants,                only : pi
  use fft_base,                 only : dfftp
  use fft_base,                 only : dfftl
  USE fft_interfaces,           ONLY : fwfft, invfft
  use gvect,                    only : gcutm, gg, g
  use gvecl,                    only : gcutml => gcutm, ggl =>gg, gl => g
  !use gvect,                    only : ngm, gcutm, gg, g, nl
  !use gvecl,                    only : ngml => ngm, gcutml => gcutm, ggl =>gg, gl => g, nll => nl
  use input_parameters,         only : fde_xc_funct,fde_kin_funct
  use scf,                      only : scf_type, create_scf_type, &
                                       destroy_scf_type, scf_type_copy
  use mp_global,                only : intra_pool_comm
  use mp_bands,                 only : intra_bgrp_comm
  use mp_large,                 only : intra_lgrp_comm
  use mp,                       only : mp_sum
  use cell_base,                only : omega
  use large_cell_base,          only : omegal => omega
  use fde,                      only : fde_nspin, fde_frag_nspin, fde_fake_nspin, &
                                       fde_regrho, currfrag, &
                                       rho_fde_large, linterlock, &
                                       rho_gauss_fde_large, rhog_gauss_fde_large, &
                                       rho_core_fde_large, rhog_core_fde_large, &
                                       fde_dotsonlarge, calc_ener_in_v_of_rho, &
                                       native_cell, reduced_cell, &
                                       NonlocalKernel_fde, NonlocalKernel, &
                                       fde_kin_is_nl, fde_overlap,     &
                                       pot_wall, saop_add, saop_nadd
  USE ions_base,                ONLY : tau, nat
  use io_global,                only : stdout, ionode
  use lsda_mod,                 only : nspin
  use spin_orb,                 only : domag
  use funct ,                   only : set_iexch, set_icorr, set_igcx, &
                                       set_igcc, set_inlc, set_dft_from_name, &
                                       get_iexch, get_icorr, get_igcx, &
                                       get_igcc, get_inlc
  USE control_flags,            ONLY : conv_elec
#if defined(__SAOP)
  !USE saop, only: v_rho_saop
#endif


  implicit none
  !-- parameters --------------------------------------------------------------
  type(scf_type), intent(inout) :: rho
  type(scf_type), intent(inout) :: rho_fde
  real(dp), intent(in)       :: rho_core(dfftp%nnr)
  complex(dp), intent(in)    :: rhog_core(:)
  real(dp), intent(in)       :: rho_core_fde(dfftp%nnr)
  complex(dp), intent(in)    :: rhog_core_fde(:)
  real(dp), intent(in)       :: rho_gauss(dfftp%nnr)
  complex(dp), intent(in)    :: rhog_gauss(:)
  real(dp), intent(in)       :: rho_gauss_fde(dfftp%nnr)
  complex(dp), intent(in)    :: rhog_gauss_fde(:)
  real(dp), intent(out)      :: ekin, etxc, ekin0, etxc0
  real(dp), intent(inout)    :: v(dfftp%nnr,nspin)

  !-- local variables ---------------------------------------------------------

  type(scf_type)           :: rhot, rhot_fde ! regularized density

  real(dp),    allocatable :: aux(:,:), aux_fde(:,:)
  real(dp),    allocatable :: auxl(:,:), auxl_fde(:,:), auxl2(:,:), auxl3(:,:)
  real(dp),    allocatable :: aux_frag(:,:)
  real(dp)                 :: trash, trash2
  integer                  :: save_iexch, save_icorr, save_igcx, save_igcc, save_inlc
  integer                  :: ipnt, iu, is, ir
  real(dp),    allocatable :: RegularDensity(:)
  complex(dp), allocatable :: RegularDensityG(:)
  integer :: innercount=0

#if defined(__SAOP)
  logical :: fde_saop
#endif

!----------------------------------------------------------------------------
!   -- BEGIN REG RHO -- to be moved to a wrapper subroutine
!            (possibly make a  module for regrho)

!  call RegularizeRho(rho)

 if (fde_regrho) then
  call errore('fde_regrho', 'keyword not yet released', 1)
     iu = 100*currfrag

    call create_scf_type (rhot,.true.)
    call scf_type_COPY (rho,rhot)

    call create_scf_type (rhot_fde,.true.)
    call scf_type_COPY (rho_fde,rhot_fde)

    allocate( RegularDensity(dfftp%nnr), RegularDensityG(dfftp%nnr) )
    RegularDensityG(1:dfftp%nnr) = (0.0d0,0.0d0)
    RegularDensity(1:dfftp%nnr) = 0.0d0
    !call PartitionWeight(RegularDensity,RegularDensityG,.false.)
    !call GlueRegRhoToRho(rho,RegularDensity,RegularDensityG)

! verbose
!  if (ionode) then
!   do ipnt=1,dfftp%nr1
!     write(iu+1,'(I3,1X,2f16.12)') ipnt, RegularDensity(ipnt), rho%of_r(ipnt,1)
!   enddo
!   do ipnt=3*ngm/4,ngm
!     write(iu+2,'(I3,1X,2f16.12)') ipnt, RegularDensityG(ipnt), real(rho%of_g(ipnt,1))
!   enddo
!  endif
! verbose end

    deallocate( RegularDensityG)
    deallocate( RegularDensity)

  call update_rho_fde (rho, .true.)
  ! Now rho_fde is the regularized one, and rhot_fde is the original rho_fde
  !   word choice a little unfortunate ...
  ! Same with the subsystem rho and rhot
 endif

!   -- END REG RHO --
!----------------------------------------------------------------------------


  ! it is assumed that the total density has been already calculated
  allocate( aux(dfftp%nnr,nspin) )
  allocate( aux_frag(dfftp%nnr,nspin) )
  if( linterlock .and. (conv_elec .or. fde_dotsonlarge .or. &
                                    calc_ener_in_v_of_rho) ) then
    allocate( auxl(dfftl%nnr, nspin) )
    if (fde_frag_nspin /= fde_nspin) allocate( auxl_fde(dfftl%nnr,fde_nspin) )
  endif
  if (fde_frag_nspin /= fde_nspin) allocate( aux_fde(dfftp%nnr,fde_nspin) )

  !----------------------------------------------------------------------
  ! non-additive kinetic based on Thomas-Fermi
  !----------------------------------------------------------------------
  ! the full system
  ekin = 0.d0
  aux(:,:) = 0.d0
  if (fde_frag_nspin /= fde_nspin) then
     aux_fde(:,:) = 0.d0
     call fde_fake_nspin(.true.)
     if ( linterlock .and. ( conv_elec .or. fde_dotsonlarge .or. &
                                    calc_ener_in_v_of_rho .or. fde_kin_is_nl )) then
        call fde_kin(rho_fde_large, rho_gauss_fde_large, &
           rhog_gauss_fde_large, ekin, auxl_fde, native_cell, NonlocalKernel_fde) ! dfftl, ngml, gl, nll, omegal, .true.)
        do is = 1, fde_nspin
          call copy_pot_l2f(auxl_fde(:,is), aux_fde(:,is))
        enddo
     else
     call fde_kin(rho_fde, rho_gauss_fde, rhog_gauss_fde, ekin, &
                        aux_fde, reduced_cell, NonlocalKernel_fde) !dfftp, ngm, g, nl, omega, .false.)
     endif
     call v_join_spin (aux, aux_fde)
     call fde_fake_nspin(.false.)
  else
    if ( linterlock .and. ( conv_elec .or. fde_dotsonlarge .or. &
                                    calc_ener_in_v_of_rho .or. fde_kin_is_nl)) then
       call fde_kin(rho_fde_large, rho_gauss_fde_large, &
               rhog_gauss_fde_large, ekin, auxl, native_cell, NonlocalKernel_fde) !dfftl, ngml, gl, nll, omegal, .true.)
       do is = 1, fde_nspin
         call copy_pot_l2f( auxl(:,is), aux(:,is) )
       enddo
    else
       call fde_kin(rho_fde, rho_gauss_fde, rhog_gauss_fde, ekin, &
                            aux, reduced_cell, NonlocalKernel_fde) ! dfftp, ngm, g, nl, omega, .false.)
    endif
  endif
  write(stdout,*) "EKIN ", ekin
  v(:,:) = v(:,:) + aux(:,:)

  ! the fragment
  ekin0 = 0.d0
  aux_frag(:,:) = 0.d0
  call fde_kin(rho, rho_gauss, rhog_gauss, ekin0, aux_frag, &
                         reduced_cell, NonlocalKernel)
!                         dfftp, ngm, g, nl, omega, .false.)
write(stdout,*) "EKIN_FRAG ", ekin0
  v(:,:) = v(:,:) - aux_frag(:,:)
!  aux = aux-aux_frag
! verbose
!   do ipnt=1,dfftp%nr1
!     write(iu+3,'(I3,1X,2f16.12)') ipnt, aux(ipnt,1)
!   enddo
! verbose end

  !----------------------------------------------------------------------
  ! Add ad-hoc repulsive potential given by the core/total densities 
  ! of the environment
  !----------------------------------------------------------------------

  if (fde_overlap) then
      aux(:,:) = 0.d0
      call fde_compute_overlap(rho, rho_fde, rho_core, rhog_core, rho_core_fde, &
                           rhog_core_fde, ekin0, aux, reduced_cell)
      v = v + aux
  endif

  !----------------------------------------------------------------------
  ! non-additive XC
  !----------------------------------------------------------------------

  ! Generally the XC functional for the nonadditive part could be different from the one used within the fragments. fde_xc_funct needs to be a (semi)local functional

  if ( trim(fde_xc_funct) /= 'SAME' ) then
    save_iexch = get_iexch()
    save_icorr = get_icorr()
    save_igcx  = get_igcx()
    save_igcc  = get_igcc()
    save_inlc  = get_inlc()
    ! Unset the functional
    call set_iexch( -1 )
    call set_icorr( -1 )
    call set_igcx( -1 )
    call set_igcc( -1 )
    call set_inlc( -1 )
    call set_dft_from_name( fde_xc_funct)
  endif

  ! the full system
  etxc = 0.d0
  aux(:,:) = 0.d0
  trash = 0.d0
  if (fde_frag_nspin /= fde_nspin) then
    aux_fde(:,:) = 0.d0
    call fde_fake_nspin(.true.)
    if ( linterlock .and. ( conv_elec .or. fde_dotsonlarge .or. &
                                  calc_ener_in_v_of_rho)) then
      if (saop_nadd) then
        write(stdout,*) "SAOP NADD RHO_FDE_LARGE FAKESPIN"
#if defined(__SAOP)
        !call v_rho_saop(rho_fde_large, rho_core_fde_large, rhog_core_fde_large, etxc, trash, auxl_fde, native_cell, fde_nspin)
#else
      call errore('fde_nonadditive', "For SAOP, must compiled eQE with -D__SAOP flag", 1)
#endif
      else
        call v_xc_large(rho_fde_large, rho_core_fde_large, &
                      rhog_core_fde_large, etxc, trash, auxl_fde)
      endif
      do is = 1, fde_nspin
        call copy_pot_l2f(auxl_fde(:,is), aux_fde(:,is))
      enddo
    else
      if (saop_nadd) then
        write(stdout,*) "SAOP NADD RHO_FDE FAKESPIN"
#if defined(__SAOP)
        !call v_rho_saop(rho_fde, rho_core_fde, rhog_core_fde, etxc, trash, aux_fde, reduced_cell, fde_nspin)
#else
      call errore('fde_nonadditive', "For SAOP, must compiled eQE with -D__SAOP flag", 1)
#endif
      else
        call v_xc(rho_fde, rho_core_fde, rhog_core_fde, etxc, trash, aux_fde)
      endif
    endif
    call v_join_spin (aux, aux_fde)
    call fde_fake_nspin(.false.)
  else
    if ( linterlock .and. ( conv_elec .or. fde_dotsonlarge .or. &
                                  calc_ener_in_v_of_rho)) then
      if (saop_nadd) then
        write(stdout,*) "SAOP NADD RHO_FDE_LARGE"
#if defined(__SAOP)
        !plot the potential, for debugging
        if (conv_elec) then
          allocate( auxl2(dfftl%nnr, nspin) )
          allocate( auxl3(dfftl%nnr, nspin) )
          !call v_rho_saop(rho_fde_large, rho_core_fde_large, rhog_core_fde_large, &
                          etxc, trash, auxl, native_cell, nspin, auxl2, auxl3)
          call plot_large(auxl(:,1),'vsaop_conv.pp       ')
          call plot_large(auxl2(:,1),'vpbe_conv.pp       ')
          call plot_large(auxl3(:,1),'vlb_conv.pp       ')
          deallocate(auxl2)
          deallocate(auxl3)
        else
          !call v_rho_saop(rho_fde_large, rho_core_fde_large, rhog_core_fde_large, etxc, trash, auxl, native_cell, nspin)
        endif
#else
      call errore('fde_nonadditive', "For SAOP, must compiled eQE with -D__SAOP flag", 1)
#endif
      else
        call v_xc_large(rho_fde_large, rho_core_fde_large, &
                          rhog_core_fde_large, etxc, trash, auxl)
      endif
      do is = 1, fde_nspin
        call copy_pot_l2f( auxl(:,is), aux(:,is) )
      enddo
    else
      if (saop_nadd) then
        write(stdout,*) "SAOP NADD RHO_FDE"
#if defined(__SAOP)
        !call v_rho_saop(rho_fde, rho_core_fde, rhog_core_fde, etxc, trash, aux, reduced_cell, nspin)
#else
      call errore('fde_nonadditive', "For SAOP, must compiled eQE with -D__SAOP flag", 1)
#endif
      else
        call v_xc(rho_fde, rho_core_fde, rhog_core_fde, etxc, trash, aux)
      endif
    endif
  endif

  !write(stdout, *) "etxc = ", etxc
  !write(stdout, *) "vtxc = ", trash
  !write(stdout, *) "vsum = ", sum(aux(:,:))
  write(stdout,*) "EXC ", etxc
  v(:,:) = v(:,:) + aux(:,:)

  ! the fragment
  etxc0 = 0.d0

#if defined(__SAOP)

  if (.false.) then
    write(stdout,*) "SAOP DEBUG"
    !call v_rho_saop(rho, rho_core, rhog_core, trash, trash2, aux, reduced_cell, fde_frag_nspin)
  endif

  !if ( innercount > 0 ) then
  if ( .false. ) then
    ! only calculate saop after the first iteration, since the evc array is not initialized before
    if (fde_saop) then
      !call v_rho_saop(rho, rho_core, rhog_core, trash, trash2, aux, reduced_cell, fde_frag_nspin)
      v(:,:) = v(:,:) + aux(:,:)
    endif

    ! subtract the GGA potential (which was added earlier), so we are effectively using saop.
    aux(:,:) = 0.d0
    call v_xc(rho, rho_core, rhog_core, etxc0, trash, aux)
    !if (fde_saop) v(:,:) = v(:,:) - aux(:,:)
    write(*,*) "E_XC: ", etxc0


    open(unit=668, file='density.dat', status="REPLACE")
    write(668,"(E15.7)") (rho%of_r(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    open(unit=668, file='vxc.dat', status="REPLACE")
    write(668,"(E15.7)") (aux(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

  endif
#endif

  if ( trim(fde_xc_funct) /= 'SAME' ) then
    aux_frag(:,:) = 0.d0
    if (saop_nadd) then
    write(stdout,*) "SAOP NADD RHO"
#if defined(__SAOP)
      !call v_rho_saop(rho, rho_core, rhog_core, etxc0, trash, aux_frag, reduced_cell, nspin)
#else
      call errore('fde_nonadditive', "For SAOP, must compiled eQE with -D__SAOP flag", 1)
#endif
    else
      call v_xc(rho, rho_core, rhog_core, etxc0, trash, aux_frag)
    endif
    
    open(unit=668, file='vxc.dat', status="REPLACE")
    write(668,"(E15.7)") (aux_frag(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)
    v(:,:) = v(:,:) - aux_frag(:,:)

    ! Revert back to the intrafragment XC functional
    call set_iexch( save_iexch )
    call set_icorr( save_icorr )
    call set_igcx( save_igcx )
    call set_igcc( save_igcc )
    call set_inlc( save_inlc )
  endif

!    aux = aux-aux_frag
! verbose
!   do ipnt=1,dfftp%nr1
!     write(iu+4,'(I3,1X,2f16.12)') ipnt, aux(ipnt,1)
!   enddo
! verbose end

  ! Add potential wall on the edges of the small cell
  !do is = 1, fde_frag_nspin
  !   v(:,is) = v(:,is) + pot_wall(:)
  !enddo

  deallocate (aux_frag, aux)
  if( linterlock .and. (conv_elec .or. fde_dotsonlarge .or. &
                                    calc_ener_in_v_of_rho) ) then
    deallocate( auxl )
    if (fde_frag_nspin /= fde_nspin) deallocate( auxl_fde )
  endif

  if (fde_frag_nspin /= fde_nspin) deallocate( aux_fde )

 if (fde_regrho) then
  call scf_type_COPY (rhot_fde,rho_fde)
  call scf_type_COPY (rhot,rho)
  call destroy_scf_type (rhot_fde)
  call destroy_scf_type (rhot)
 endif

 innercount = innercount + 1

  return

END SUBROUTINE fde_nonadditive

!----------------------------------------------------------------------------
SUBROUTINE fde_kin(rho, rho_core, rhog_core, ene, v, cell, NonlocalKernel)!dfftp, ngm, g, nl, omega, llarge)
  !----------------------------------------------------------------------------
  !
  ! Depending on the KEDF used, call either the GGA routine or the non-local subroutine
  !
  use kinds,                    only : dp
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use input_parameters,         only : fde_kin_funct
  use fde_types,                only : simulation_cell
  use kernel_module,            only : KernelType
  implicit none
   !-- parameters --------------------------------------------------------------
  type(scf_type), intent(inout)     :: rho
  type(simulation_cell), intent(in) :: cell
  type(KernelType), intent(inout)   :: NonlocalKernel
  real(dp), intent(in)              :: rho_core(cell%dfftp%nnr)
  complex(dp), intent(in)           :: rhog_core(cell%ngm)
  real(dp), intent(out)             :: ene 
  real(dp), intent(out)             :: v(cell%dfftp%nnr,nspin)


  select case (trim(fde_kin_funct))
    case('TF','VW','DGE2','LLP','PW86','LC94','APBEK','revAPBEK','TFP')
       call fde_kin_gga(rho, rho_core, rhog_core, ene, v, trim(fde_kin_funct) &
                                          , cell )
!                                          , dfftp, ngm, g, nl, omega, llarge)
    case('LDA')
       call fde_kin_lda(rho, rho_core, ene, v, cell%dfftp, cell%omega)
    case('GPRHO0')
       call fde_kin_gp_rho0(rho, rho_core, rhog_core, ene, v, cell, NonlocalKernel)
    case('GPLDA')
       call fde_kin_gp_lda(rho, rho_core, rhog_core, ene, v)
    case('NONE')
       return
    case default
       call errore('fde_kin_funct', 'unknown functional'//trim(fde_kin_funct), 1)
  end select

END SUBROUTINE fde_kin


!----------------------------------------------------------------------------
SUBROUTINE fde_kin_lda(rho, rho_core, ene, v, dfftp, omega)
  !----------------------------------------------------------------------------
  !
  ! ... generalized Thomas-Fermi
  !
  use kinds,                    only : dp
  use constants,                only : pi
!  use fft_base,                 only : dfftp
!  use cell_base,                only : omega
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use spin_orb,                 only : domag
  use mp_global,                only : intra_bgrp_comm
  use mp,                       only : mp_sum
  USE fft_types,  ONLY : fft_type_descriptor
  implicit none
  !-- parameters --------------------------------------------------------------
  type(scf_type), intent(in) :: rho
  type(fft_type_descriptor) , intent(in) :: dfftp
  real(dp), intent(in)       :: rho_core(dfftp%nnr)
  real(dp), intent(out)      :: ene
  real(dp), intent(out)      :: v(dfftp%nnr,nspin)
  real(dp), intent(in)       :: omega


  !-- local variables ---------------------------------------------------------
  real(dp), allocatable :: rhoout(:)
  real(dp), parameter :: hart_to_ryd = 2.d0
  real(dp), parameter :: Cf = (3.d0/10.d0) * (3.d0*pi*pi) ** (2.d0/3.d0)
  real(dp), parameter :: vanishing_rho = 1.d-10
  real(dp) :: e0, domega, fac
  integer :: i, ispin, nspin0

  ! TODO: non-collinear
  nspin0 = nspin
  if (nspin == 4) nspin0 = 1
  if (nspin==4 .and. domag) nspin0 = 2
  fac = 1.d0 / dble(nspin0)


  allocate(rhoout(dfftp%nnr))

  domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  ene = 0.d0
  v(:,:) = 0.d0
  do ispin = 1, nspin0

     rhoout(:)  = fac * rho_core(:)  + rho%of_r(:,ispin)

     do i = 1, dfftp%nnr
        if ( rhoout(i) > vanishing_rho ) then
          v(i,ispin) = (5.d0/3.d0) * Cf * rhoout(i)**(2.d0/3.d0) * hart_to_ryd
          e0 = Cf * rhoout(i)**(5.0/3.0) * hart_to_ryd
          ene = ene + e0
        endif
     enddo
  enddo
  ene = ene * domega
!  vrho = vrho * domega

  ! TODO: check if those are needed when using more than 1 CPU/fragment
!  call mp_sum(ene, intra_bgrp_comm)
  call mp_sum(ene, dfftp%comm)
  !!call mp_sum(vrho, intra_bgrp_comm)
  deallocate(rhoout)
  return

END SUBROUTINE fde_kin_lda






!----------------------------------------------------------------------------
SUBROUTINE fde_compute_overlap(rho, rho_fde, rho_core, rhog_core, rho_core_fde, &
                       rhog_core_fde, ene, v, cell)
  !----------------------------------------------------------------------------
  !
  ! Subroutine to compute a repulsive potential E_rep \propto (rho | rho_fde - rho )
  !
  use kinds,                    only : dp
  use constants,                only : pi
  use io_global,                only : stdout
  USE fft_types,                ONLY : fft_type_descriptor
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum
  use mp_global,                only : intra_bgrp_comm
  use cell_base,                only : alat
  use spin_orb,                 only : domag
  use fde,                      only : fde_regrho, fde_overlap_c
  use fde_types,                only : simulation_cell
  implicit none
  !-- input ----------------------------------------------------

  type(simulation_cell), intent(in) :: cell
  type(scf_type),        intent(in) :: rho
  type(scf_type),        intent(in) :: rho_fde
  real(dp),              intent(in) :: rho_core    (cell%dfftp%nnr)
  real(dp),              intent(in) :: rho_core_fde(cell%dfftp%nnr)
  complex(dp),           intent(in) :: rhog_core    (cell%ngm)
  complex(dp),           intent(in) :: rhog_core_fde(cell%ngm)
  !-- outputs ------------------------------------------------------
  real(dp),    intent(out)    :: ene
  real(dp),    intent(out)    :: v(cell%dfftp%nnr,nspin)
  !-- local parameters ---------------------------------------------------------
  real(dp), parameter      :: hart_to_ryd = 2.d0
  real(dp), parameter      :: epsr = 1d-8, epsg = 1d-10
  !-- local variables ---------------------------------------------------------
  real(dp),    allocatable :: rho_frag(:), rho_total(:)
  integer      :: ngm, ipnt
  integer      :: nspin0, ispin
  real(dp)     :: domega, omega, fac, coefficient
  logical      :: llarge
  !-- pointers to cell%  ---------------------------------------------------------
  type(fft_type_descriptor) , pointer :: dfftp
  real(dp),                   pointer :: g(:,:)
  integer,                    pointer :: nl(:)

  ! TODO: non-collinear
  nspin0 = nspin
  if (nspin == 4)             nspin0 = 1
  if (nspin == 4 .and. domag) nspin0 = 2

  fac = 1.d0 / dble(nspin0)

 ! associate the cell% params
  dfftp => cell%dfftp
  g => cell%g
  nl => cell%nl
  ngm = cell%ngm
  omega = cell%omega
  llarge = cell%is_native_cell

  ! allocate memory
  allocate( rho_frag(dfftp%nnr), rho_total(dfftp%nnr) )

 ! compute volume elements
  domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

  ene    = 0.d0
  v(:,:) = 0.d0

 !
 !
  ispin_: do ispin = 1, nspin0
        if (fde_overlap_c .le. 10000.0) then
          rho_frag(:)  = fac * rho_core(:)
          rho_total(:) = fac * rho_core_fde(:)
          coefficient  = fde_overlap_c*hart_to_ryd
          on_grid1: do ipnt = 1, dfftp%nnr
            v(ipnt,ispin) = coefficient*max( (rho_total(ipnt) - rho_frag(ipnt)), 0.0d0  )
            ene           = ene + max( (rho_total(ipnt) - rho_frag(ipnt)), 0.0d0  ) * abs(rho%of_r(ipnt,ispin))
          enddo on_grid1
        elseif (fde_overlap_c .gt. 10000.0 .and. fde_overlap_c .le. 20000.0 ) then
          rho_frag(:)  = fac * rho_core(:)     + rho%of_r(:,ispin)
          rho_total(:) = fac * rho_core_fde(:) + rho_fde%of_r(:,ispin)
          coefficient  = ( fde_overlap_c - 10000.0d0 )*hart_to_ryd
          on_grid2: do ipnt = 1, dfftp%nnr
            v(ipnt,ispin) = coefficient*max( (rho_total(ipnt) - rho_frag(ipnt)), 0.0d0  )
            ene           = ene + max( (rho_total(ipnt) - rho_frag(ipnt)), 0.0d0  ) * abs(rho%of_r(ipnt,ispin))
          enddo on_grid2
        elseif (fde_overlap_c .gt. 20000.0 .and. fde_overlap_c .le. 30000.0 ) then
          rho_total(:) = fac * rho_core_fde(:)
          rho_frag(:)  = fac * rho_core(:)
          coefficient  = ( fde_overlap_c - 20000.0d0 )*hart_to_ryd
          on_grid3: do ipnt = 1, dfftp%nnr
            v(ipnt,ispin) = -coefficient*max( rho_frag(ipnt), 0.0d0  )
            ene           = ene + max( rho_total(ipnt) - rho_frag(ipnt), 0.0d0  ) * abs(rho%of_r(ipnt,ispin))
          enddo on_grid3
        endif
  enddo ispin_

  call mp_sum(ene, intra_bgrp_comm)
  write(stdout,'(a,f8.3)') " ==> Overlap: ", ene*domega
  deallocate(rho_frag,rho_total)

END SUBROUTINE fde_compute_overlap



!----------------------------------------------------------------------------
SUBROUTINE fde_kin_gga(rho, rho_core, rhog_core, ene, v, funct, cell)!dfftp, ngm, g, nl, omega, llarge)
  !----------------------------------------------------------------------------
  !
  ! Subroutine to compute kinetic energy and potential using several GGA functionals
  !
  use kinds,                    only : dp
  use constants,                only : pi
  use io_global,                only : stdout
  USE fft_types,                ONLY : fft_type_descriptor
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
  use cell_base,                only : alat
  use spin_orb,                 only : domag
  use input_parameters,         only : fde_kin_funct
  use fde,                      only : fde_regrho
  use fde_types,                only : simulation_cell
  use fde_functionals
  implicit none
  !-- parameters --------------------------------------------------------------
  type(scf_type), intent(in) :: rho
  type(simulation_cell), intent(in) :: cell
  real(dp), intent(in)       :: rho_core(cell%dfftp%nnr)
  complex(dp), intent(in)    :: rhog_core(cell%ngm)
  real(dp), intent(out)      :: ene!, vrho
  real(dp), intent(out)      :: v(cell%dfftp%nnr,nspin)
  character(len=*)              :: funct
  !-- local variables ---------------------------------------------------------
  real(dp), parameter :: hart_to_ryd = 2.d0
  real(dp), parameter :: Cf = (3.d0/10.d0) * (3.d0*pi*pi) ** (2.d0/3.d0)
  real(dp), parameter :: Cs = 1.d0 / (2.d0 * (3.d0*pi*pi) ** (1.d0/3.d0))
  real(dp), parameter :: epsr = 1d-8, epsg = 1d-10
  !----------------------------------------------------------------------------
  !type(fft_type_descriptor) , pointer :: dfftp
  !real(dp), pointer :: g(:,:)
  !integer, pointer :: nl(:)
  !integer :: ngm
  !real(dp) :: omega
  !logical :: llarge
  !
  real(dp), allocatable :: s(:), Fs(:), dF_dS(:), grho(:,:), rhoout(:)
  real(dp), allocatable :: prodotto(:,:), dprodotto(:), mod_grho(:)
  complex(dp), allocatable :: rhogsum(:)
  real(dp) :: domega, fac
  integer :: i, ispin, nspin0


  associate( &
    dfftp => cell%dfftp, &
    ngm => cell%ngm, &
    nl => cell%nl, &
    g => cell%g, &
    omega => cell%omega, &
    tpiba => cell%tpiba, &
    llarge => cell%is_native_cell &
  )

  ! TODO: non-collinear
  nspin0 = nspin
  if (nspin == 4) nspin0 = 1
  if (nspin==4 .and. domag) nspin0 = 2
  fac = 1.d0 / dble(nspin0)

  ! allocate memory
  allocate( s(dfftp%nnr), Fs(dfftp%nnr), dF_ds(dfftp%nnr) )
  allocate( grho(3,dfftp%nnr), rhoout(dfftp%nnr), rhogsum(ngm) )
  allocate( mod_grho(dfftp%nnr), prodotto(3,dfftp%nnr), dprodotto(dfftp%nnr) )

  domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  ene = 0.d0
  v(:,:) = 0.d0

  ! loop over spins
  do ispin = 1, nspin0

     ! add core charge
     rhoout(:)  = fac * rho_core(:)  + rho%of_r(:,ispin)
     rhogsum(:) = fac * rhog_core(:) + rho%of_g(:,ispin)

     ! fix for the spin-polarized case (ref: DFT of Atoms and
     !                                  molecules, Parr & Yang, pp. 174)

     ! compute \nabla\rho

     !if (llarge) then
     !  call gradrho_large(dfftp%nnr, rhogsum, ngm, g, nl, grho)
     !else
     !  call gradrho(dfftp%nnr, rhogsum, ngm, g, nl, grho)
     !endif
     call fft_gradient_g2r(dfftp, rhogsum, g, grho)


     !if ( fde_regrho ) call GradRegularRho(rhoout,grho,rho%is_fde)

     mod_grho(:) = sqrt(grho(1,:)**2 + grho(2,:)**2 + grho(3,:)**2)

     ! compute s
     s = 0.d0
!$omp parallel do
     do i = 1, dfftp%nnr
        if (abs(rhoout(i)) > epsr .and. mod_grho(i) > epsg) &
            s(i) = fac**(1.d0/3.d0) * Cs * mod_grho(i) / abs(rhoout(i))**(4.d0/3.d0)
     enddo
!$omp end parallel do

     ! compute F(s) and dF(s)/ds
     select case (trim(funct))
        case('TF')
           call fde_kin_tf(dfftp%nnr, s, Fs, dF_ds)
        case('VW')
           call fde_kin_vw(dfftp%nnr, s, Fs, dF_ds)
        case('DVW')
           call fde_kin_dvw(dfftp%nnr, s, Fs, dF_ds)
        case('DGE2')
           call fde_kin_dge2(dfftp%nnr, s, Fs, dF_ds)
        case('LLP')
           call fde_kin_llp(dfftp%nnr, s, Fs, dF_ds)
        case('PW86')
           call fde_kin_pw86(dfftp%nnr, s, Fs, dF_ds)
        case('LC94')
           call fde_kin_lc94(dfftp%nnr, s, Fs, dF_ds)
        case('APBEK')
           call fde_kin_apbek(dfftp%nnr, s, Fs, dF_ds)
        case('revAPBEK')
           call fde_kin_rapbek(dfftp%nnr, s, Fs, dF_ds)
        case('TFP')
           call fde_kin_tfp(dfftp%nnr, s, Fs, dF_ds)
        case default
           call errore('fde_kin_funct', 'unknown functional'//trim(funct), 1)
     end select

     ! compute (dG/d\nabla\rho) * (\nabla\rho/|\nabla\rho|)
     prodotto = 0.d0
!$omp parallel do
     do i = 1, dfftp%nnr
        if (abs(rhoout(i)) > epsr .and. mod_grho(i) > epsg) &
        prodotto(1:3,i) = fac**(-4.d0/3.d0) * Cs * abs(rhoout(i))**(1.d0/3.d0) * dF_ds(i) * &
                          grho(1:3,i) / mod_grho(i)
     enddo
!$omp end parallel do

     !if (llarge) then
     !  call grad_dot_large(dfftp%nnr, prodotto, ngm, g, nl, alat, dprodotto)
     !else
     !  call grad_dot(dfftp%nnr, prodotto, ngm, g, nl, alat, dprodotto)
     !endif

     call fft_graddot(dfftp, prodotto, g, dprodotto)

     v(:,ispin) = fac**(-5.d0/3.d0) * (5.d0/3.d0) * Cf * abs(rhoout)**(2.d0/3.d0) * Fs - &
                  fac**(-5.d0/3.d0) * (4.d0/3.d0) * Cf * abs(rhoout)**(2.d0/3.d0) * s * dF_ds - &
                  Cf * dprodotto
     v(:,ispin) = fac * v(:,ispin) ! CHECK!!!

     ! subtract rho_core
     ! Why removing rho core from the evaluation of the energy?

     !   rhoout(:) = rhoout(:) - fac*rho_core(:)

     !vrho = vrho + sum(rhoout(:)*v(:,ispin)) ! Check if 0.5 factor is needed
     ene = ene + fac * Cf * sum(abs(fac**(-1.d0) * rhoout(:))**(5.0/3.0) * Fs(:))
  enddo  ! ispin

  ! convert to Hartree
  ene = ene * domega*hart_to_ryd
  !vrho = vrho * domega*hart_to_ryd
  v = v*hart_to_ryd

!  call mp_sum(ene, intra_bgrp_comm)
  call mp_sum(ene, dfftp%comm)
  !call mp_sum(vrho, intra_bgrp_comm)

  end associate

  return

END SUBROUTINE fde_kin_gga


!----------------------------------------------------------------------------
SUBROUTINE fde_kin_gp_rho0(rho, rho_core, rhog_core, ene, v, cell, Kernel)
  !----------------------------------------------------------------------------
  !
  ! Subroutine to compute the Genova-Pavanello nonlocal kinetic potential
  !
  use kinds,                    only : dp
  use constants,                only : pi
  !use fft_base,                 only : dfftp
  use fft_types,                only : fft_type_descriptor
  USE fft_interfaces,           ONLY : fwfft, invfft
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
  !use cell_base,                only : at, alat, tpiba, tpiba2, omega
  !use gvect,                    only : ngm, nl, g, gg, gstart
  use spin_orb,                 only : domag
  use input_parameters,         only : fde_kin_funct
  use fde,                      only : fde_r0, fde_gp_rhot, fde_gp_alpha
  use fde_types,                only : simulation_cell
  use io_global,                only : stdout
  use klist,                    only : nelec
  use io_global,                only : ionode
  use kernel_module
  implicit none
  !-- parameters --------------------------------------------------------------
  type(scf_type), intent(inout)   :: rho
  type(simulation_cell), intent(in) :: cell
  type(KernelType), intent(inout) :: Kernel
  real(dp), intent(in)       :: rho_core(cell%dfftp%nnr)
  complex(dp), intent(in)    :: rhog_core(cell%ngm)
  real(dp), intent(out)      :: ene
  real(dp), intent(out)      :: v(cell%dfftp%nnr,nspin)
  !-- local variables ---------------------------------------------------------
  real(dp), parameter :: hart_to_ryd = 2.d0
  real(dp), parameter :: Cf = (3.d0/10.d0) * (3.d0*pi*pi) ** (2.d0/3.d0)
  real(dp), parameter :: Cs = 1.d0 / (2.d0 * (3.d0*pi*pi) ** (1.d0/3.d0))
  real(dp), parameter :: epsr = 1d-8, epsg = 1d-10
  !----------------------------------------------------------------------------
  real(dp)       :: domega, fac, ene0, trash
  !real(dp)       :: rho_star, kf_star, charge, kf_tmp
  real(dp)       :: deltat, t_tmp
  real(dp)       :: aux_lob(cell%dfftp%nnr,nspin), eta(cell%ngm)

  complex(dp), allocatable :: aux_cplx(:)
  !
  ! pointers
  type(fft_type_descriptor), pointer :: dfftp
  integer, pointer :: ngm, gstart, nl(:)
  real(dp), pointer :: at(:,:), alat, omega, tpiba, tpiba2
  real(dp), pointer :: g(:,:), gg(:)
  !
  integer :: i, ispin, nspin0, index0
  integer :: ipnt, MaxPoints
  !
  dfftp => cell%dfftp
  ngm => cell%ngm
  gstart => cell%gstart
  nl => cell%nl
  at => cell%at
  alat => cell%alat
  tpiba => cell%tpiba
  tpiba2 => cell%tpiba2
  omega => cell%omega
  g => cell%g
  gg => cell%gg
  !
  ! TODO: non-collinear
  nspin0 = nspin
  if (nspin == 4) nspin0 = 1
  if (nspin==4 .and. domag) nspin0 = 2

  fac = 1.d0 / dble(nspin0)

  ! allocate memory
  allocate( aux_cplx(dfftp%nnr))

  domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

  ! get rho_T(G) - the "kinetic lump"
  aux_cplx(1:ngm) = rho%of_g(1:ngm,1)*Kernel%KernelFunction(1:ngm,1)

  ! Calculate the Potential
  v(:,:) = 0.d0
  ene0 = 0.d0
  if ( cell%is_native_cell) then
    call v_h_large( aux_cplx , ene0, trash, v )
  else
    call v_h( aux_cplx , ene0, trash, v )
  endif

! Calculate the contribution to v and E of revAPBEK
  ene0 = 0.d0

  ene0 = 0.5d0*sum( v(:,1) * rho%of_r(:,1) )*domega
  call mp_sum(ene0, intra_bgrp_comm)

  aux_lob(:,:) = 0.d0
  call fde_kin_gga(rho, rho_core, rhog_core, ene0, aux_lob, 'TFP', cell)
  ene = ene0

! ****** DEBUG PRINT *********
!  do i=1,dfftp%nr1
!   write(101,*) i*alat/dfftp%nr1, v(i,1)
!  enddo
! **** END DEBUG PRINT *******

  v = v + aux_lob

! ****** DEBUG PRINT *********
!  do i=1,dfftp%nr1
!   write(102,*) i*alat/dfftp%nr1, v(i,1)
!  enddo
! **** END DEBUG PRINT *******

  deallocate( aux_cplx)

END SUBROUTINE fde_kin_gp_rho0
!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
SUBROUTINE fde_kin_gp_lda(rho, rho_core, rhog_core, ene, v)
  !----------------------------------------------------------------------------
  !
  ! Subroutine to compute LDA vaersion of
  !  the Genova-Pavanello nonlocal kinetic potential
  !
  use kinds,                    only : dp
  use constants,                only : pi
  use fft_base,                 only : dfftp
  USE fft_interfaces,           ONLY : fwfft, invfft
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum, mp_max
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
  use cell_base,                only : at, alat, tpiba, tpiba2, omega
  use gvect,                    only : ngm, g, gg, gstart
  use spin_orb,                 only : domag
  use input_parameters,         only : fde_kin_funct
  use fde,                      only : fde_r0, fde_gp_rhot, fde_gp_alpha, &
                                       reduced_cell
  use io_global,                only : stdout
  use klist,                    only : nelec
  use io_global,                only : ionode, stdout
  use splinelib , only : spline , splint
  implicit none
  !-- parameters --------------------------------------------------------------
  type(scf_type), intent(in) :: rho
  real(dp), intent(in)       :: rho_core(dfftp%nnr)
  complex(dp), intent(in)    :: rhog_core(ngm)
  real(dp), intent(out)      :: ene
  real(dp), intent(out)      :: v(dfftp%nnr,nspin)
  !-- local variables ---------------------------------------------------------
  real(dp), parameter :: hart_to_ryd = 2.d0
  real(dp), parameter :: Cf = (3.d0/10.d0) * (3.d0*pi*pi) ** (2.d0/3.d0)
  real(dp), parameter :: Cs = 1.d0 / (2.d0 * (3.d0*pi*pi) ** (1.d0/3.d0))
  real(dp), parameter :: epsr = 1d-8, epsg = 1d-10
  !----------------------------------------------------------------------------
  real(dp)       :: domega, fac, ene0
  real(dp)       :: rho_star, kf_star, charge, kf_tmp, rho_max, kf_max
  real(dp)       :: deltat, t_tmp
  real(dp)       :: trash
  real(dp)       :: aux_lob(dfftp%nnr,nspin), eta(ngm)
!  complex(dp), external :: g_of_eta

  complex(dp), allocatable :: aux_cplx(:), w(:)
  real(dp) , allocatable :: v_of_kf(:,:,:), kfs(:), ydata(:), d2y(:)

  integer :: i, ispin, nspin0, index0, ikf, nkf, ir, ir_ofmax
  integer :: ipnt, MaxPoints

  ! TODO: non-collinear
  nspin0 = nspin
  if (nspin == 4) nspin0 = 1
  if (nspin==4 .and. domag) nspin0 = 2

  fac = 1.d0 / dble(nspin0)


  ! allocate memory
  allocate( aux_cplx(dfftp%nnr), w(ngm))
  nkf = 10! number of points in the kf grid
  allocate( kfs(nkf) , v_of_kf(dfftp%nnr, nspin0, nkf), ydata(nkf), d2y(nkf) )

  domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

! ***** DISABLED ********
  ! loop over spins
!  rho_star = 0.d0
!  do ispin = 1, nspin0
!     rho_star = rho_star + sum(rho%of_r(:,ispin)**2)
!  end do
! ***********************


  ! calculate charge of fragment
   charge = 0.d0
   if ( gstart == 2 ) then
      charge = omega*real( rho%of_g(1,1) )
      if ( nspin0 == 2 ) charge = charge + omega*real( rho%of_g(1,2) )
   end if

  call mp_sum(charge, intra_bgrp_comm)
  call mp_sum(rho_star, intra_bgrp_comm)

  if ( nkf == 1 ) then
    !rho_star = rho_star*domega/charge
    rho_max = fde_r0 ! the one parameter here
  else
  !
  ! find the maximum value of the density
  rho_max = 0.d0
    do ir = 1 , dfftp%nnr
      !rho_max = max( rho_max , rho%of_r(ir,1) )
      if (rho%of_r(ir,1) > rho_max ) then
         rho_max = rho%of_r(ir,1)
         ir_ofmax = ir
      endif
      if ( nspin0 == 2 ) rho_max = max( rho_max , rho%of_r(ir,2) )
    enddo
    call mp_max(rho_max, intra_bgrp_comm )
  endif

  kf_max  = (3.d0*pi*pi*rho_max)**(1.d0/3.d0)
  ! to get rid of the dependency on the rho_star parameter, and in order to obtain a position dependent kf(r), we calculate the potential several times, according to a discre number of values from kf = 0 to kf = (3.d0*pi*pi*rho_max)**(1.d0/3.d0)

  kf_loop : do ikf = 1 , nkf  ! will probably have to include ikf = 0 later on (large q limit?)

    kf_star = kf_max * dble(ikf) / dble(nkf)
    kfs(ikf) = kf_star
  ! Calculate the Kernel in reciprocal space v
  w(:) = 0.0d0
  MaxPoints = 100

  deltat = 1.d0/dble(MaxPoints)
  t_tmp  = 1.d0/dble(MaxPoints)

  hypercorrelation_: do ipnt=1,MaxPoints
      kf_tmp = kf_star*(t_tmp**(1.d0/3.d0))
      eta(:) = sqrt(gg(:)*tpiba2)/(2.d0*kf_tmp)
       do i=1,ngm
         w(i) = w(i) + g_of_eta(eta(i),gg(i)*tpiba2) * (pi*pi/kf_tmp)
       enddo
       t_tmp = t_tmp + deltat
  enddo hypercorrelation_

! regular GP here
  w(1:ngm) = w(1:ngm)/4.d0/pi/dble(MaxPoints)

! CZQL correction
  w(1:ngm) = w(1:ngm) + (fde_gp_rhot/charge)*exp(-fde_gp_alpha*gg(1:ngm)*tpiba2)/4.d0/pi

! ****** DEBUG PRINT *********
!  if (ionode) then
!  open (unit=99, file='gp_kernel_czql.dat')
!  ene0 = -1.d0
!  do i=1,ngm
!   if (gg(i) == ene0 ) cycle
!   write(99,*) sqrt(gg(i)*tpiba2), w(i)
!   ene0 = gg(i)
!  enddo
!  close(99)
!  endif
! **** END DEBUG PRINT *******


  ! get rho_T(G) - the "kinetic lump"
  aux_cplx(1:ngm) = rho%of_g(1:ngm,1)*w(1:ngm)


  ! Calculate the Potential
  v_of_kf(:,:,ikf) = 0.d0
  aux_lob(:,:) = 0.d0
  ene0 = 0.d0
  call v_h( aux_cplx , ene0, trash, aux_lob )
  v_of_kf(:,:,ikf) = aux_lob(:,:)

  end do kf_loop

  ! for each point calculate the value of kf and interpolate the non local contribution of the potential accordingly from the discrete copies calculated before
  do ir = 1 , dfftp%nnr
    do ispin = 1, nspin0
      kf_tmp = (3.d0*pi*pi*abs(rho%of_r(ir,ispin)))**(1.d0/3.d0)
      if ( nkf > 1 .and. kf_tmp < kfs(1) ) then ! we don't want to extrapolate...
        v(ir,ispin) = v_of_kf(ir,ispin,1) !0.d0
      else
        ydata(:) = v_of_kf(ir,ispin,1:nkf)
        d2y = 0.d0
        !
        CALL spline( kfs, ydata(:), 0.0_DP, 0.0_DP, d2y  )
        !
        v(ir,ispin) = splint( kfs, ydata(:), d2y, kf_tmp )
        ! DEBUG BLOCK
        !if (ir == ir_ofmax ) then
        !  write(stdout,*) ydata , v(ir,1)
        !  flush(stdout)
        !endif
        ! !!!!
      endif
    enddo
  enddo


! ***** DISABLED ********
!   Calculate E of TF to be used as
!       the value for the NAKE
! with no contribution to the potential
!
!  ene0 = 0.d0
!  call fde_kin_lda(rho, rho_core, ene0, aux_lob)
!  ene =  ene0
! ***********************

! Calculate the contribution to E of MyGGA
  ene0 = 0.d0
  aux_lob(:,:) = 0.d0
  call fde_kin_gga(rho, rho_core, rhog_core, ene0, aux_lob, 'TFP', reduced_cell)!&
                                   !dfftp, ngm, g, nl, omega, .false.)
  ene = ene0

! ****** DEBUG PRINT *********
!  do i=1,dfftp%nr1
!   write(101,*) i*alat/dfftp%nr1, v(i,1)
!  enddo
! **** END DEBUG PRINT *******

  v = v + aux_lob

! ****** DEBUG PRINT *********
!  do i=1,dfftp%nr1
!   write(102,*) i*alat/dfftp%nr1, v(i,1)
!  enddo
! **** END DEBUG PRINT *******

  deallocate( aux_cplx, w)
  deallocate( kfs, v_of_kf, ydata )

END SUBROUTINE fde_kin_gp_lda


!--------------------------------
FUNCTION g_of_eta(eta,gg)
!--------------------------------
!   The -1/Chi_Lindhard routine
!   w/out the pi^2/k_f prefactor (included later)
!--------------------------------
  use kinds,                    only : dp
  implicit none

  complex(dp)       :: g_of_eta

  real(dp)       :: eta
  real(dp)       :: gg
  real(dp)       :: g_of_eta_real

! from Y. A. Wang and E. A. Carter, "Orbital-Free Kinetic Energy Density
! Functional Theory," in "Theoretical Methods in Condensed Phase Chemistry," S.
! D. Schwartz, Ed., within the series "Progress in Theoretical Chemistry and
! Physics," Kluwer, 117-84 (2000).
  real(dp) :: pa1=1.d0, pa2=8.d0,   pa3=104.d0, pa4=1048.d0,    pa5=24536.d0
  real(dp) :: qa1=3.d0, qa2=45.d0,  qa3=945.d0, qa4=14175.d0,   qa5=467775.d0
  real(dp) :: pb1=1.d0, pb2=8.d0,   pb3=8.d0,   pb4=12728.d0,   pb5=551416.d0
  real(dp) :: qb1=5.d0, qb2=175.d0, qb3=375.d0, qb4=1010625.d0, qb5=65690625.d0
  real(dp) :: pb6=41587384.d0, qb6=6897515625.d0


     if (abs(eta-1.0d0) <= 1.0d-3 ) then
      g_of_eta_real =  (2.0d0 - 3.d0 + 3.d0/5.d0)            * gg
     elseif (abs(eta) < (0.8) ) then
      g_of_eta_real =  (1.0d0 + 3.d0/5.d0                  + &
                             (pa1/qa1 - 3.d0) * eta**2.0d0 + &
                             (pa2/qa2)        * eta**4.0d0 + &
                             (pa3/qa3)        * eta**6.0d0 + &
                             (pa4/qa4)        * eta**8.0d0 + &
                             (pa5/qa5)        * eta**10.0d0) * gg
     elseif (abs(eta) > 2.0d0 ) then
      g_of_eta_real = - 3.0d0 * ( &
                             (pb2/qb2) * eta**(-2.0d0) + &
                             (pb3/qb3) * eta**(-4.0d0) + &
                             (pb4/qb4) * eta**(-6.0d0) + &
                             (pb5/qb5) * eta**(-8.0d0) + &
                             (pb6/qb6) * eta**(-10.0d0) )    * gg
     else
  g_of_eta_real = (1.0d0 / (0.5d0 + ((1.d0-eta*eta)/(4.d0*eta)) * &
                   log(abs((1.d0+eta)/(1.d0-eta)))) - 3.d0 * eta * eta + 3.d0/5.d0 ) &
                   * gg
     endif

 g_of_eta = cmplx(g_of_eta_real,0.d0, kind=dp)
  return

END FUNCTION g_of_eta


!----------------------------------------------------------------------------
SUBROUTINE v_h_fde( rhog, rhor, rhog_fde, ehart, charge, v )
  !----------------------------------------------------------------------------
  !
  ! ... Hartree potential V_H[rho_fde](r) and Hartree energy as
  ! ... int d^3 rho(r) V_H[rho_fde](r)
  !
  use constants,                only : fpi, e2
  use kinds,                    only : dp
  use fft_base,                 only : dfftp, dffts, dfftl
  use fft_interfaces,           only : fwfft, invfft
  use gvect,                    only : ngm, gg, gstart
  use gvecl,                    only : ngml => ngm
  use lsda_mod,                 only : nspin
  use cell_base,                only : omega, tpiba2
  use control_flags,            only : gamma_only, iverbosity
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
  use mp,                       only : mp_sum, mp_barrier
  use mp_images,                only : intra_image_comm, inter_fragment_comm, my_image_id
  use mp_world, only : world_comm
  use fde,                      only : fde_nspin, fde_frag_nspin, &
                                       rho_fde, fde_si, fde_si_all2all, &
                                       fde_si_alpha, fde_fake_nspin, &
                                       rho_fde_large, linterlock, f2l, &
                                       fde_si_vec, fde_true_ehart
  use mp_images,                only : my_image_id
  use wavefunctions_module,     only : evc
  use klist,                    only : nelec
  use wvfct,                    only : npw, npwx, wg
  use io_global,                only : stdout, ionode, ionode_id
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  !-- parameters --------------------------------------------------------------
  complex(dp), intent(in) :: rhog(ngm,nspin)
  complex(dp), intent(in) :: rhog_fde(ngm,fde_nspin)
  real(dp), intent(in) ::    rhor(dfftp%nnr,nspin)
  real(dp), intent(inout) :: v(dfftp%nnr,nspin)
  real(dp), intent(out)   :: ehart, charge  ! charge of single fragment
  !-- local variables ---------------------------------------------------------
  real(dp), allocatable    :: vh(:,:)
  real(dp), allocatable    :: aux_fde(:,:), auxl_fde(:,:)
  real(dp), allocatable    :: v_si(:), raux1(:), rauxl(:)
  real(dp), allocatable    :: auxl(:,:)
  complex(dp), allocatable :: gaux(:), tmp_gaux(:), gauxl(:), tmp_gauxl(:)
  complex(dp), allocatable :: mgaux(:)
  real(dp)                 :: tot_charge, domega
  integer                  :: is,myband
  real(dp)                 :: esi
  !character(len=16)            :: image_label, filename
  !character(len=6), external :: int_to_char


  ! calculate total Hartree potential and add it to XC potential
  allocate( vh(dfftp%nnr,nspin) )
  if (linterlock) allocate(auxl(dfftl%nnr,fde_nspin))
  vh(:,:) = v(:,:)
  if (fde_frag_nspin /= fde_nspin) then
     allocate( aux_fde(dfftp%nnr,fde_nspin) )
     aux_fde = 0.d0
     call fde_fake_nspin(.true.)
     if ( linterlock ) then
        allocate( auxl_fde(dfftl%nnr, fde_nspin) )
        auxl_fde = 0.d0
        call v_h_large( rho_fde_large%of_g, ehart, tot_charge, auxl_fde )
        fde_true_ehart = ehart
!        if (iverbosity>10) write(stdout,*) 'TRUE ehart = ' , ehart
        do is = 1 , fde_nspin
           call copy_pot_l2f( auxl_fde(:,is), aux_fde(:,is))
        enddo
        deallocate( auxl_fde )
     else
        call v_h( rhog_fde, ehart, tot_charge, aux_fde )
     endif
     call v_join_spin( v, aux_fde )
     call fde_fake_nspin(.false.)
     deallocate(aux_fde)
     v = v + vh
  else
     if (linterlock) then
       auxl = 0.d0
       call v_h_large( rho_fde_large%of_g, ehart, tot_charge, auxl )
       fde_true_ehart = ehart
!        if (iverbosity>10) write(stdout,*) 'TRUE ehart = ' , ehart
       do is = 1 , fde_nspin
         call copy_pot_l2f( auxl(:,is), v(:,is))
       enddo
       if (my_image_id == 0 ) then
!       call plot_large(auxl(:,1),'vh0_large.pp    ')
!       call plot_dense(v(:,1),'vh0_dense.pp    ')
       else
!       call plot_large(auxl(:,1),'vh1_large.pp    ')
!       call plot_dense(v(:,1),'vh1_dense.pp    ')
       endif
!       call mp_barrier(world_comm)
!       stop
       v = v + vh
     else
       call v_h( rhog_fde, ehart, tot_charge, v )
     endif
  endif

! Include vandevondele--sprik SIE correction
!

  if ( sum(fde_si_vec(:)) > 0 .or. fde_si_all2all ) then
     allocate(v_si(dfftp%nnr))
     if (linterlock)  then
        allocate( auxl_fde(dfftl%nnr, 1) )
        auxl_fde = 0.d0
     endif

     v_si = 0.d0
     if ( sum(fde_si_vec(:)) > 0 ) then
        allocate(gaux(dfftp%nnr), tmp_gaux(ngm))

        if ( linterlock ) then

          !if ( fde_frag_nspin == 2 ) then
             if (iverbosity>-1) write(stdout,*) 'SI correction on the spin density'
             allocate( gauxl(dfftl%nnr), tmp_gauxl(ngml) )
             gauxl(:) = cmplx( rho_fde_large%of_r(:,1) - &
                               rho_fde_large%of_r(:,2), &
                               0.d0, kind=dp )
             call fwfft ('Rho', gauxl, dfftl)
             tmp_gauxl(1:ngml) = gauxl(dfftl%nl(1:ngml))
             call fde_fake_nspin( .true. , 1 )
             call v_h_large( tmp_gauxl, esi, charge, auxl_fde )
             call fde_fake_nspin( .false. )
             write(stdout,*) 'FDE_SI ehart = ' , esi
             !do is = 1 , fde_nspin
             call copy_pot_l2f( auxl_fde(:,1), v_si)
             !enddo
             deallocate( gauxl, tmp_gauxl )
          !endif

        elseif (fde_si) then

        if (fde_frag_nspin == 2 ) then
           if (iverbosity>-1) write(stdout,*) 'SI correction on the spin density'
           gaux(:) = CMPLX(abs(rhor(:,1)-rhor(:,2)),0.d0,kind=dp)
           call fwfft ('Rho', gaux, dfftp)
        elseif (fde_frag_nspin == 1) then
           if (iverbosity>-1) write(stdout,*) 'SI correction on the HOMO'
           if (iverbosity>-1) write(stdout,*) 'WARNING: SI correction only available for Norm Conserving PS'
           if (iverbosity>-1) write(stdout,*) 'WARNING: SI correction only available single k-point calculations at the moment'
           allocate(mgaux(dffts%nnr))
           myband = NINT(nelec/2.0d0)
           mgaux(dffts%nl(1:npw)) = evc(1:npw,myband)
           call invfft ( 'Wave', mgaux(:), dffts )
           mgaux(1:dffts%nnr) = CONJG(mgaux(1:dffts%nnr)) * mgaux(1:dffts%nnr)
           mgaux(:) = mgaux(:)*wg(myband,1)/omega/2.0d0
           domega = omega / dble(dffts%nr1 * dffts%nr2 * dffts%nr3)
           charge = sum(real(mgaux(1:dffts%nnr)))*domega
           call mp_sum(charge, intra_bgrp_comm)
           if (iverbosity>-1) write(stdout,*) 'Charge in SI correction = ',charge
           call fft_interpolate_complex(dffts, mgaux, dfftp, gaux,1)
           deallocate(mgaux)
           call fwfft('Rho', gaux, dfftp)
        endif

        tmp_gaux(1:ngm) = gaux(dfftp%nl(1:ngm))

        call fde_fake_nspin( .true. , 1 )
        call v_h(tmp_gaux, ehart, charge, v_si)
        call fde_fake_nspin( .false. )

        deallocate(gaux, tmp_gaux)
        endif

     endif

     if (fde_si_all2all .and. .not. linterlock) then
        if (iverbosity >-1) write(stdout,*) 'SI correction broadcasted across all the fragments'
        if (ionode) allocate(raux1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
        if (linterlock .and. .false.) then
           if (ionode) allocate(rauxl(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
           call gather_grid( dfftl, auxl_fde(:,1), rauxl )
           if (ionode) call mp_sum(rauxl, inter_fragment_comm)
           raux1(:) = rauxl(f2l(:))
           if (ionode) deallocate(rauxl)
        else
           call gather_grid( dfftp,v_si, raux1)
           if (ionode) call mp_sum(raux1, inter_fragment_comm)
        endif
        call scatter_grid( dfftp,raux1, v_si)
        if (ionode) deallocate(raux1)
     endif
     !   if (linterlock .and. ) then
     !      call copy_pot_l2f(auxl_fde(:,1), v_si(:))
     !   endif
     !endif

     if ( linterlock .and. &
          .not. fde_si_all2all .and. &
          .not. fde_si ) v_si = 0.d0

     !image_label = '_' // int_to_char(my_image_id)
     !filename = 'v_si' // trim(image_label) // '.pp       '
     !call plot_dense( v_si , filename )
     !call plot_large( rho_fde_large%of_r(:,1) - &
     !                 rho_fde_large%of_r(:,2) , 'rho_spin.pp     ' )

     if ( fde_frag_nspin == 1 ) then
        v(:,1) = v(:,1) - fde_si_alpha * v_si(:)
     elseif ( fde_frag_nspin == 2 ) then
        !if ( intsirho > 0.d0 ) then
           v(:,1) = v(:,1) - fde_si_alpha * v_si(:)
        !else
           v(:,2) = v(:,2) - fde_si_alpha * v_si(:)
        !endif
     endif

     deallocate(v_si)

     !stop

  endif

  vh(:,:) = v(:,:) - vh(:,:)

  ! calculate charge of fragment
  charge = 0.d0
  if ( gstart == 2 ) then
     charge = omega*real( rhog(1,1) )
     if ( nspin == 2 ) charge = charge + omega*real( rhog(1,2) )
  end if
  call mp_sum(charge, intra_bgrp_comm)

  ! calculate Hartree energy of fragment
  ! TODO: non-collinear
  domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  ehart = sum(rhor(:,1) * vh(:,1))
  if (nspin == 2)  ehart = ehart + sum(rhor(:,2) * vh(:,2))
  ehart = 0.5d0 * ehart * domega

  CALL mp_sum(  ehart , intra_bgrp_comm )

  deallocate (vh)
  return

END SUBROUTINE v_h_fde

!----------------------------------------------------------------------------
SUBROUTINE update_rho_fde (density, total)
  !----------------------------------------------------------------------------
  !
  !  Subroutine to update rho_fde. If total is .true. rho_fde is the sum of the densities
  !  of the fragments. If total is .false. rho_fde of each image is only updated with the
  !  contribution of the new rhoin of the single fragment
  !
  use mp,                       only : mp_sum, mp_bcast
  use mp_images,                only : inter_fragment_comm, intra_image_comm, my_image_id
  use fde
  use scf,                      only : scf_type, scf_type_copy
  use io_global,                only : stdout, ionode, ionode_id
  use fft_interfaces,           only : fwfft, invfft
  use fft_base,                 only : dfftp, dfftl
  use gvect ,                   only : ngm
  use gvecl ,                   only : ngml => ngm
  use lsda_mod,                 only : nspin
  use control_flags,            only : iverbosity
  use command_line_options,    only : fancy_parallel_
  use scatter_mod,          only : gather_grid, scatter_grid

  implicit none
  type(scf_type), intent(in) :: density ! scf_type, either rho or rhoin
  logical, intent(in) :: total
  integer             :: is
  real(dp), allocatable :: raux(:)
  real(dp), allocatable :: rauxl(:), rho_old_large(:,:), density_large(:,:)
  complex(dp), allocatable :: gaux(:)
  complex(dp), allocatable :: gauxl(:)

   call start_clock('update_rho')
     if (total) then
        if (iverbosity>10) write(stdout,*) 'TOTAL FDE DENSITY UPDATE'
        if (fde_frag_nspin == fde_nspin) then
           call scf_type_copy(density, rho_fde)
        else
           rho_fde%of_r(:,1) = 0.5_dp * density%of_r(:,1)
           rho_fde%of_r(:,2) = rho_fde%of_r(:,1)
           rho_fde%of_g(:,1) = 0.5_dp * density%of_g(:,1)
           rho_fde%of_g(:,2) = rho_fde%of_g(:,1)
        endif
        if (ionode) allocate(raux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
        if (ionode .and. linterlock) allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
        if (linterlock) allocate( gaux(dfftp%nnr), gauxl(dfftl%nnr) )
        !
        do is=1,fde_nspin
           !
           call gather_grid( dfftp,rho_fde%of_r(:,is), raux)
            if (ionode) then
              rauxl = 0.d0
              rauxl(f2l(:)) = raux(:)
              call mp_sum(rauxl, inter_fragment_comm)
              if (.not. fde_dotsonlarge .or. fde_overlap) then
                raux = 0.d0
                raux(:) = rauxl(f2l(:))
              endif
            endif
            !
            call scatter_grid( dfftl,rauxl, rho_fde_large%of_r(:,is))

            gauxl(:) = cmplx(rho_fde_large%of_r(:,is), 0.d0, kind=dp)
            call fwfft ('Rho', gauxl, dfftl)
            rho_fde_large%of_g(1:ngml,is) = gauxl(dfftl%nl(1:ngml))
            !
            if (.not. fde_dotsonlarge .or. fde_overlap) then
              call scatter_grid( dfftp,raux, rho_fde%of_r(:,is))
              gaux(:) = cmplx(rho_fde%of_r(:,is), 0.d0, kind=dp)
              call fwfft ('Rho', gaux, dfftp)
              rho_fde%of_g(1:ngm,is) = gaux(dfftp%nl(1:ngm))
            endif
            !
           !
        end do
        !
        if (ionode) deallocate(raux)
        !
        deallocate(gaux, gauxl)
        if (ionode) deallocate(rauxl)
     else
        if (iverbosity>10) write(stdout,*) 'PARTIAL FDE DENSITY UPDATE'
        if (fancy_parallel_) call errore('update_rho_fde','partial density update can only be combined with -nfp',1)
        allocate(rho_old_large(dfftl%nnr,fde_frag_nspin)) ! to hold rho_old on large grid
        allocate(density_large(dfftl%nnr,fde_frag_nspin)) ! to hold density on large grid
        allocate(gauxl(dfftl%nnr))
        rho_old_large(:,:)=0.d0
        density_large(:,:)=0.d0
        do is=1,nspin
           call copy_pot_f2l(rho_old%of_r(:,is),rho_old_large(:,is))
           call copy_pot_f2l(density%of_r(:,is),density_large(:,is))
        enddo
        if (fde_frag_nspin == fde_nspin) then
           rho_fde_large%of_r = rho_fde_large%of_r - rho_old_large + density_large
        else
           rho_fde_large%of_r(:,1) = rho_fde_large%of_r(:,1) - 0.5_dp * (rho_old_large(:,1) - density_large(:,1))
           rho_fde_large%of_r(:,2) = rho_fde_large%of_r(:,2) - 0.5_dp * (rho_old_large(:,1) - density_large(:,1))
        endif
        do is=1,fde_nspin
           gauxl(:)=cmplx(rho_fde_large%of_r(:,is),0.d0,kind=dp)
           call fwfft('Rho',gauxl,dfftl)
           rho_fde_large%of_g(1:ngml,is)=gauxl(dfftl%nl(1:ngml))
        enddo
        deallocate(rho_old_large,density_large,gauxl)
     end if

   call stop_clock('update_rho')
   return

END SUBROUTINE update_rho_fde

!----------------------------------------------------------------------------
SUBROUTINE fde_energies ( conv )
  !----------------------------------------------------------------------------
  !
  !  Subroutine to calculate the FDE total and non-additive energies
  !
  use ener,                     only : etot, ewld, ehart, vtxc, etxc, eband
  use mp,                       only : mp_sum, mp_barrier, mp_bcast
  use mp_images,                only : inter_fragment_comm, intra_image_comm
  use io_global,                only : ionode, ionode_id, stdout
  use mp_world,                 only : world_comm
  use fde
  use input_parameters,         only : fde_xc_funct

  implicit none
  logical , intent(in) :: conv
  etot_sum = etot
  call mp_barrier(world_comm) ! Some extra barriers here and there, just in case

  if ( conv .and. trim(fde_xc_funct) == 'SAME' ) then
    etxc0_nadd = etxc
  endif

  if (ionode) call mp_sum(etot_sum, inter_fragment_comm)
  call mp_bcast(etot_sum, ionode_id, intra_image_comm)
  if (ionode) call mp_sum(ekin0_nadd, inter_fragment_comm)
  call mp_bcast(ekin0_nadd, ionode_id, intra_image_comm)
  if (ionode) call mp_sum(etxc0_nadd, inter_fragment_comm)
  call mp_bcast(etxc0_nadd, ionode_id, intra_image_comm)

  call mp_barrier(world_comm)
!  if ( fde_fat ) etot_fde_old = etot_fde
  ekin_nadd = ekin_nadd - ekin0_nadd
  etxc_nadd = etxc_nadd - etxc0_nadd
  edc_fde = -dble(nfragments-1)*(ewld)
!  etot_fde = etot_sum + edc_fde + ekin_nadd + etxc_nadd
  etot_fde = etot_sum + edc_fde + ekin_nadd + etxc_nadd

END SUBROUTINE fde_energies


SUBROUTINE fat_accuracy (rho1, rho2 )
  !----------------------------------------------------------------------------
  !
  !  Subroutine to calculate the scf accuracy between two consecutive freeze and thaw cycles
  !
  use scf,             only : scf_type, mix_type, create_mix_type, destroy_mix_type, &
                              assign_scf_to_mix_type, rho_ddot, create_scf_type, destroy_scf_type
  use gvecs,           only : ngms
  use fde,             only : dr2_fat, fde_fake_nspin

  implicit none

  type(scf_type), intent(in) :: rho1, rho2
  type(mix_type)            :: drho_m!, rho2_m
  type(scf_type)            :: drho_s

  call fde_fake_nspin(.true.)

  call create_mix_type(drho_m)
  call create_scf_type(drho_s)
!  call create_mix_type(rho2_m)

  drho_s%of_r = rho1%of_r - rho2%of_r
  drho_s%of_g = rho1%of_g - rho2%of_g

  call assign_scf_to_mix_type(drho_s, drho_m)
!  call assign_scf_to_mix_type(rho2, rho2_m)

  dr2_fat = rho_ddot(drho_m, drho_m, ngms)

  call destroy_scf_type(drho_s)
  call destroy_mix_type(drho_m)

  call fde_fake_nspin(.false.)

  return

END SUBROUTINE fat_accuracy


SUBROUTINE v_split_spin ( v, aux )
  !----------------------------------------------------------------------------
  !
  !  Subroutine to split a non-spinpolarized potential into alpha/beta channels
  !  (if needed)
  !  v ----> aux
  !
  use kinds,                    only : dp
  use fde,                      only : fde_nspin, fde_frag_nspin
  use fft_base,                 only : dfftp

  real(dp), intent(inout)   :: v(dfftp%nnr, fde_frag_nspin), aux(dfftp%nnr, fde_nspin)

  if ( fde_nspin == fde_frag_nspin ) then
     aux = v
  else
!     aux(:,1) = 0.5_dp * v(:,1)
     aux(:,1) = v(:,1)
     aux(:,2) = aux(:,1)
  endif

  return

END SUBROUTINE v_split_spin

SUBROUTINE v_join_spin ( v, aux )
  !----------------------------------------------------------------------------
  !
  !  Subroutine to combine a spinpolarized potential
  !  (if needed)
  !  aux ----> v
  !
  use kinds,                    only : dp
  use fde,                      only : fde_nspin, fde_frag_nspin
  use fft_base,                 only : dfftp

  real(dp), intent(inout)   :: v(dfftp%nnr, fde_frag_nspin), aux(dfftp%nnr, fde_nspin)

  if ( fde_nspin == fde_frag_nspin ) then
     v = aux
  else
     v(:,1) = 0.5_dp * ( aux(:,1) + aux(:,2) )
  endif

  return

END SUBROUTINE v_join_spin


!----------------------------------------------------------------------------
SUBROUTINE fde_plot_electrostatics
  !----------------------------------------------------------------------------
  !
  !  Plot v_bare + v_h on the large cell, useful to get the vacuum level for work function evaluation   !
  USE fft_base,                 ONLY : dfftp, dfftl
  USE cell_base,                ONLY : celldm, at, ibrav
  USE large_cell_base,          ONLY : celldml => celldm, atl => at
  USE ions_base,                ONLY : ntyp => nsp, atm, zv
  USE gvect,                    ONLY : gcutm
  use gvecl,                    only : gcutml => gcutm
  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode, stdout
  USE io_files,                 ONLY : prefix
  USE fde
  use scf_large,          only : vltot_large => vltot
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  integer, external :: find_free_unit
  integer :: unit, is
  logical :: exst
  real(dp), allocatable :: rauxl1(:), rauxl2(:), auxl(:,:)
  real(dp) :: trash0, trash1
  character(20) :: charfrag

  !if ( currfrag /= 1) return

  allocate(auxl(dfftl%nnr, fde_nspin))
  if (ionode.and.&
      currfrag == 1 ) allocate( rauxl1(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x), &
                                rauxl2(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x) )

  ! V_bare
  call gather_grid( dfftl, vltot_large, rauxl1 )

  ! V_H
  auxl(:,:) = 0.d0
  call v_h_large( rho_fde_large%of_g, trash0, trash1, auxl )

  call gather_grid( dfftl, auxl(:,1), rauxl2 )

  ! sum and plot
  if (ionode.and.&
      currfrag == 1 ) then

      rauxl1 = rauxl1 + rauxl2
      !rauxl1 = rauxl2

      call plot_io (trim(prefix)//'_electro.pp', 'VH and Vbare', dfftl%nr1x, dfftl%nr2x, dfftl%nr3x, &
             dfftl%nr1,  dfftl%nr2,  dfftl%nr3, nat_fde, ntyp, ibrav, celldml, atl, &
             gcutml, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_large, rauxl1, +1)

      deallocate( rauxl1, rauxl2 )

  endif




END SUBROUTINE fde_plot_electrostatics


!----------------------------------------------------------------------------
SUBROUTINE fde_plot_density
  !----------------------------------------------------------------------------
  !
  !  Plot FDE total density and potential
  !
  USE fft_base,                 ONLY : dfftp, dfftl
  USE cell_base,                ONLY : celldm, at, ibrav
  USE large_cell_base,          ONLY : celldml => celldm, atl => at
  USE ions_base,                ONLY : ntyp => nsp, atm, zv, tau, ityp, nat
  USE gvect,                    ONLY : gcutm
  use gvecl,                    only : gcutml => gcutm
  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode, stdout
  USE io_files,                 ONLY : prefix
  USE fde
  use scf, only : rho
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  integer, external :: find_free_unit
  integer :: unit, is
  logical :: exst
  real(dp), allocatable :: raux(:), raux1(:), rauxl(:), auxl(:)
  character(20) :: charfrag

  !if ( currfrag /= 1) return

  allocate(raux(dfftp%nnr))
  if (linterlock) allocate(auxl(dfftl%nnr))
  if (ionode) allocate(raux1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

  ! first the charge density
!  raux = rho_fde%of_r(:,1)
!  if (fde_nspin == 2) raux = raux + rho_fde%of_r(:,2)

  if (sum(fde_printdensity_vec(:)) > 0) then
  flush(stdout)
  if (linterlock) then
    CALL dcopy (dfftl%nnr, rho_fde_large%of_r (1, 1), 1, auxl, 1)
    DO is = 2, fde_nspin
      CALL daxpy (dfftl%nnr, 1.d0, rho_fde_large%of_r (1, is), 1, auxl, 1)
    ENDDO

      if ( ionode .and. currfrag==1) allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
      call gather_grid( dfftl,auxl, rauxl)
      if ( ionode .and. currfrag==1) &
      call plot_io (trim(prefix)//'_fde_rho.pp', 'FDE rho', dfftl%nr1x, dfftl%nr2x, dfftl%nr3x, &
             dfftl%nr1,  dfftl%nr2,  dfftl%nr3, nat_fde, ntyp, ibrav, celldml, atl, &
             gcutml, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, rauxl, +1)

      flush(stdout)
      if (fde_nspin == 2) then
         ! the the magnetization
!         raux = rho_fde%of_r(:,1) - rho_fde%of_r(:,2)
         CALL dcopy (dfftl%nnr, rho_fde_large%of_r (1, 1), 1, auxl, 1)
         CALL daxpy (dfftl%nnr, - 1.d0, rho_fde_large%of_r (1, 2), 1, auxl, 1)
         call gather_grid( dfftl,auxl, rauxl)
         if ( ionode .and. currfrag==1 ) &
         call plot_io (trim(prefix)//'_fde_magn.pp', 'FDE magn', dfftl%nr1x, dfftl%nr2x, dfftl%nr3x, &
                dfftl%nr1,  dfftl%nr2,  dfftl%nr3, nat_fde, ntyp, ibrav, celldml, atl, &
                gcutml, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, rauxl, +1)
      endif
      if ( ionode .and. currfrag==1) deallocate(rauxl)
  else
    CALL dcopy (dfftp%nnr, rho_fde%of_r (1, 1), 1, raux, 1)
    DO is = 2, fde_nspin
      CALL daxpy (dfftp%nnr, 1.d0, rho_fde%of_r (1, is), 1, raux, 1)
    ENDDO

      call gather_grid( dfftp,raux, raux1)
      if ( ionode ) &
      call plot_io (trim(prefix)//'_fde_rho.pp', 'FDE rho', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
             dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav, celldm, at, &
             gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)

      flush(stdout)
      if (fde_nspin == 2) then
         ! the the magnetization
!         raux = rho_fde%of_r(:,1) - rho_fde%of_r(:,2)
         CALL dcopy (dfftp%nnr, rho_fde%of_r (1, 1), 1, raux, 1)
         CALL daxpy (dfftp%nnr, - 1.d0, rho_fde%of_r (1, 2), 1, raux, 1)
         call gather_grid( dfftp,raux, raux1)
         if ( ionode) &
         call plot_io (trim(prefix)//'_fde_magn.pp', 'FDE magn', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav, celldm, at, &
                gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
      endif
 endif ! linterlock
endif
      flush(stdout)

  if (fde_print_density_frag) then
      CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
      DO is = 2, fde_frag_nspin
        CALL daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
      ENDDO

      call gather_grid( dfftp,raux, raux1)
      if (currfrag<10) write(charfrag,'(i1)')currfrag
      if (currfrag>9) write(charfrag,'(i2)')currfrag
      if ( ionode ) &
      call plot_io (trim(prefix)//'_fde_rho_'//trim(charfrag)//'.pp', 'FDE frag', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
             dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
             gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux1, +1)
      if (fde_frag_nspin == 2) then
         ! the the magnetization
!         raux = rho_fde%of_r(:,1) - rho_fde%of_r(:,2)
         CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
         CALL daxpy (dfftp%nnr, - 1.d0, rho%of_r (1, 2), 1, raux, 1)
         call gather_grid( dfftp,raux, raux1)
         if ( ionode) &
         call plot_io (trim(prefix)//'_fde_magn_'//trim(charfrag)//'.pp', 'FDE magn', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
                gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux1, +1)
      endif
   endif

  if (fde_print_density_frag_large.and.linterlock) then
! saving fragment density on the large grid, needed for frozen calculation in TDDFT
      CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
      DO is = 2, fde_frag_nspin
        CALL daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, raux, 1)
      ENDDO

      call gather_grid( dfftp,raux, raux1)

      if (currfrag<10) write(charfrag,'(i1)')currfrag
      if (currfrag>9) write(charfrag,'(i2)')currfrag
      if ( ionode ) then
         allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
         rauxl=0.d0
         rauxl(f2l(:))=raux1(:)

         call plot_io (trim(prefix)//'_fde_rho_large_'//trim(charfrag)//'.pp', 'FDE frag', dfftl%nr1x, dfftl%nr2x, dfftl%nr3x, &
             dfftl%nr1,  dfftl%nr2,  dfftl%nr3, nat, ntyp, ibrav, celldml, atl, &
             gcutml, dual, ecutwfc, 0, atm, ityp, zv, tau_large, rauxl, +1)
      endif
      if (fde_frag_nspin == 2) then
         ! the the magnetization
!         raux = rho_fde%of_r(:,1) - rho_fde%of_r(:,2)
         CALL dcopy (dfftp%nnr, rho%of_r (1, 1), 1, raux, 1)
         CALL daxpy (dfftp%nnr, - 1.d0, rho%of_r (1, 2), 1, raux, 1)
         call gather_grid( dfftp,raux, raux1)
         if ( ionode) then
            rauxl=0.d0
            rauxl(f2l(:))=raux1(:)
            call plot_io (trim(prefix)//'_fde_magn_large_'//trim(charfrag)//'.pp', 'FDE magn', dfftl%nr1x, dfftl%nr2x, dfftl%nr3x, &
                dfftl%nr1,  dfftl%nr2,  dfftl%nr3, nat, ntyp, ibrav, celldml, atl, &
                gcutml, dual, ecutwfc, 0, atm, ityp, zv, tau_large, rauxl, +1)
         endif
      endif
      if (ionode) deallocate(rauxl)
   endif

  deallocate(raux)
  if (ionode) deallocate(raux1)
  if (linterlock) deallocate(auxl)

  return
END SUBROUTINE fde_plot_density


!----------------------------------------------------------------------------
SUBROUTINE fde_plot_embedpot
  !----------------------------------------------------------------------------
  !
  !  Plot FDE total density and potential
  !
  USE fft_base,                 ONLY : dfftp, dfftl
  USE fft_interfaces,           ONLY : fwfft
  USE cell_base,                ONLY : celldm, at, ibrav
  USE large_cell_base,          ONLY : atl => at, bgl => bg
  USE ions_base,                ONLY : ntyp => nsp, atm, zv, ityp, tau, nat
  USE gvect,                    ONLY : gcutm, ngm
  use gvecl,                    only : ngml => ngm, gcutml => gcutm, ggl =>gg, gl => g, &
                                       eigts1l => eigts1, &
                                       eigts2l => eigts2, &
                                       eigts3l => eigts3

  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode, ionode_id
  USE io_files,                 ONLY : prefix
  USE mp,                       ONLY : mp_bcast
  USE mp_images,                ONLY : my_image_id, intra_image_comm, inter_fragment_comm
  USE scf,                      ONLY : rho, rho_core, rhog_core
  use input_parameters,         only : fde_xc_funct
  use lsda_mod,                 only : nspin
  USE fde
  use vlocal ,                  only : strf
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  integer, external :: find_free_unit
  integer :: unit, is
  logical :: exst
  real(dp), allocatable :: raux(:), rauxl(:)
  real(dp), allocatable :: vemb(:,:), aux(:,:)
  real(dp), allocatable :: vembl(:,:), auxl(:,:)
  integer :: nat_emb
  integer, allocatable :: ityp_emb(:)
  real(dp), allocatable :: tau_frag_emb(:,:)
  complex(dp), allocatable :: environ_rhog(:,:)!, environ_rhog_core(:,:)
  complex(dp), allocatable :: strf_frag_large(:,:)
  complex(dp), allocatable :: eigts1l_emb(:,:) , eigts2l_emb(:,:) ,eigts3l_emb(:,:)
  real(dp), allocatable :: environ_rhor(:,:)
  complex(dp), allocatable :: gauxl(:)
  real(dp) :: trash0, trash1, trash2, trash3
  character(len=7)  :: channel, image_label
  character(len=6), external :: int_to_char
  integer :: frag_to_plot, ifrag



  ! because of fancy parallelization, only the embedding potential of a single fragment can ble plotted. Make sure all the processors are on the same page about it ...
  frag_to_plot = -1

  do ifrag = 1 , nfragments
    if (fde_plotemb_vec(ifrag) == 1) then
       frag_to_plot = ifrag
       exit
    endif
  enddo

  if ( frag_to_plot < 0 ) return

  allocate( vemb(dfftp%nnr,fde_frag_nspin), &
            environ_rhor(dfftl%nnr,fde_frag_nspin), &
            environ_rhog(ngml,fde_frag_nspin), &
            vembl(dfftl%nnr, fde_frag_nspin), &
            aux(dfftp%nnr,fde_frag_nspin) , &
            auxl(dfftl%nnr,fde_frag_nspin) )
  if (linterlock) allocate( gauxl(dfftl%nnr) )
  if (ionode) allocate(raux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  vemb(:,:) = 0.d0


  call fde_nonadditive(rho, rho_fde, rho_core, rhog_core, &
                          rho_core_fde, rhog_core_fde, &
                          rho_gauss, rhog_gauss, &
                          rho_gauss_fde, rhog_gauss_fde, &
                          trash0, &
                          trash1, &
                          trash2, &
                          trash3, vemb)


  if ( trim(fde_xc_funct) == 'SAME' ) then
    aux(:,:) = 0.d0
    call v_xc(rho, rho_core, rhog_core, trash0, trash1, aux)
    vemb(:,:) = vemb(:,:) - aux(:,:)
  endif

  if (fde_frag_nspin /= fde_nspin) then
    !TODO: what when only some fragments are openshell?
  else
    ! Calculate the density of the environment. Code is hard to read because of fancy parallelization

    if (linterlock) then
       do is = 1, nspin
          call gather_grid( dfftp,rho%of_r(:,is), raux)
          if (ionode) then
             allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
             rauxl = 0.d0
             rauxl(f2l(:)) = raux(:)
             call mp_bcast(rauxl, frag_to_plot -1, inter_fragment_comm)
          endif
          call scatter_grid( dfftl,rauxl, auxl(:,is)) ! Change it
          if (ionode) deallocate( rauxl )
          environ_rhor(:,is) = rho_fde_large%of_r(:,is) - auxl(:,is)
          gauxl(:) = cmplx(environ_rhor(:,is), 0.d0, kind=dp)
          call fwfft ('Rho', gauxl, dfftl)
          environ_rhog(1:ngml,is) = gauxl(dfftl%nl(1:ngml))
       enddo

       auxl(:,:) = 0.d0
       aux(:,:) = 0.d0

       call v_h_large(environ_rhog, trash0, trash1, auxl )
       do is = 1 , nspin
          call copy_pot_l2f( auxl(:,is), aux(:,is))
       enddo
       vemb(:,:) = vemb(:,:) + aux(:,:)

    endif ! linterlock


  endif


  ! calculate the strf of the fragment on the large grid and the local potential
  nat_emb = nat
  if (ionode) call mp_bcast( nat_emb, frag_to_plot -1, inter_fragment_comm )
  call mp_bcast( nat_emb, ionode_id, intra_image_comm )
  allocate( tau_frag_emb(3,nat_emb), ityp_emb(nat_emb))
  allocate( strf_frag_large(ngml, ntyp) )
  allocate( eigts1l_emb( -dfftl%nr1:dfftl%nr1, nat_emb ), &
            eigts2l_emb( -dfftl%nr2:dfftl%nr2, nat_emb ), &
            eigts3l_emb( -dfftl%nr3:dfftl%nr3, nat_emb ) )
  if (currfrag == frag_to_plot) then
     call scatter_coordinates(tau_fde, tau_frag_emb)
     ityp_emb = ityp
  endif

  if (ionode) call mp_bcast( ityp_emb, frag_to_plot -1, inter_fragment_comm )
  call mp_bcast( ityp_emb, ionode_id, intra_image_comm )

  if (ionode) call mp_bcast( tau_frag_emb, frag_to_plot -1, inter_fragment_comm )
  call mp_bcast( tau_frag_emb, ionode_id, intra_image_comm )
  ! nat, tau_frag, nsp, ityp
  CALL struc_fact( nat_emb, tau_frag_emb, ntyp, ityp_emb, ngml, gl, bgl, &
                  dfftl%nr1, dfftl%nr2, dfftl%nr3, strf_frag_large, eigts1l_emb, eigts2l_emb, eigts3l_emb )
  !
  call setlocal_fde_large(auxl, strf_fde_large - strf_frag_large)
  !
  call copy_pot_l2f( auxl(:,1), aux(:,1))
  !do is = 1 , nspin
     !call copy_pot_l2f( auxl(:,is), aux(:,is))
  !enddo



  !call setlocal_fde(aux(:,1), strf_fde(:,:) - strf(:,:))

  vemb(:,1) = vemb(:,1) + aux(:,1)
  if (fde_frag_nspin == 2) vemb(:,2) = vemb(:,2) + aux(:,1)

!  allocate(raux(dfftp%nnr))
  if ( currfrag == frag_to_plot ) then

  image_label = '_' // int_to_char(my_image_id)
  channel = '_alpha'

  do is = 1, fde_frag_nspin
    if (is == 2) channel = '_beta'
    call gather_grid( dfftp,vemb(:,is), raux)
    if ( ionode ) &
      call plot_io (trim(prefix)//'_embedpot'//trim(image_label)//trim(channel)//'.pp', &
              'Embedding Potential', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
              gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux, +1)
  enddo

!  deallocate(raux)
  endif

  deallocate(vemb, vembl, environ_rhor, environ_rhog)
  deallocate(aux, auxl)
  deallocate( strf_frag_large, eigts1l_emb, eigts2l_emb, eigts3l_emb )
  if (ionode) deallocate(raux)
  if (linterlock) deallocate(gauxl)


  return
END SUBROUTINE fde_plot_embedpot


!----------------------------------------------------------------------------
SUBROUTINE fde_plot_gradrho
  !----------------------------------------------------------------------------
  !
  !  Plot the gradient of the electron density (both frag and total) so that we can calculate the embedding potentials
  !  with other programs. Only closed shell case implemented for now.
  !
  USE fft_base,                 ONLY : dfftp, dfftl
  USE fft_interfaces,           ONLY : fwfft
  USE cell_base,                ONLY : celldm, at, ibrav
  USE large_cell_base,          ONLY : atl => at, bgl => bg
  USE ions_base,                ONLY : ntyp => nsp, atm, zv, ityp, tau, nat
  USE gvect,                    ONLY : gcutm, ngm, g
  USE gvecl,                    ONLY : ngm_l => ngm, g_l => g
  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode, ionode_id
  USE io_files,                 ONLY : prefix
  USE mp,                       ONLY : mp_bcast
  USE mp_images,                ONLY : my_image_id, intra_image_comm, inter_fragment_comm
  USE scf,                      ONLY : rho, rho_core, rhog_core
  use input_parameters,         only : fde_xc_funct
  use lsda_mod,                 only : nspin
  USE fde
  use vlocal ,                  only : strf
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  integer, external :: find_free_unit
  integer :: unit, is, ipol
  logical :: exst
  real(dp), allocatable :: raux(:), rauxl(:)
  real(dp), allocatable :: grho(:,:), aux(:,:)
  real(dp), allocatable :: grho_large(:,:), auxl(:,:)
  real(dp) :: trash0, trash1, trash2, trash3
  character(len=7)  :: channel, image_label, idir
  character(len=6), external :: int_to_char
  integer :: frag_to_plot, ifrag



  ! because of fancy parallelization, only the embedding potential of a single fragment can ble plotted. Make sure all the processors are on the same page about it ...
  frag_to_plot = -1

  do ifrag = 1 , nfragments
    if (fde_plotgrad_vec(ifrag) == 1) then
       frag_to_plot = ifrag
       exit
    endif
  enddo

  if ( frag_to_plot < 0 ) return

  if ( currfrag == frag_to_plot ) then
     allocate( grho(3,dfftp%nnr) )
     if (ionode) allocate(raux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  endif

  if (linterlock) allocate( grho_large(3,dfftl%nnr) )

  image_label = '_' // int_to_char(my_image_id)
  channel = '_alpha'

  ! Grad of rho_frag
  if ( currfrag == frag_to_plot ) then

      grho = 0.d0

      !call gradrho( dfftp%nnr, rho%of_g(:,1), ngm, g, dfftp%nl, grho )
      call fft_gradient_g2r(dfftp, rho%of_g(:,1), g, grho)
      

      do ipol = 1, 3
        if ( ipol == 1 ) then
           idir = 'x'
        elseif ( ipol == 2 ) then
           idir = 'y'
        elseif ( ipol == 3 ) then
           idir = 'z'
        endif

        call gather_grid( dfftp, grho(ipol,:), raux )


        if (ionode) then
           call plot_io (trim(prefix)//'_grad'//trim(idir)//trim(image_label)//trim(channel)//'.pp', &
                  'Gradrho', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                  dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
                  gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux, +1)
        endif

      enddo
  endif


 ! Grad rho_fde

  if (linterlock) then
     allocate( grho_large(3,dfftl%nnr) )
     grho_large = 0.d0
     !call gradrho_large( dfftl%nnr, rho_fde_large%of_g(:,1), ngm_l, g_l, dfftl%nl, grho_large )
     call fft_gradient_g2r(dfftl, rho_fde_large%of_g(:,1), g_l, grho_large)
     if (ionode) allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
     do ipol = 1, 3
        if ( ipol == 1 ) then
           idir = 'x'
        elseif ( ipol == 2 ) then
           idir = 'y'
        elseif ( ipol == 3 ) then
           idir = 'z'
        endif
        call gather_grid( dfftl, grho_large(ipol,:), rauxl )
        if (ionode) call mp_bcast(rauxl, frag_to_plot -1, inter_fragment_comm)
        if ( currfrag == frag_to_plot ) then
           if (ionode) then
              raux(:) = rauxl(f2l(:))
              call plot_io (trim(prefix)//'_grad'//trim(idir)//'_fde'//trim(channel)//'.pp', &
                 'Gradrho', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                 dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
                 gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux, +1)
           endif
        endif

     enddo
     if (ionode) deallocate(rauxl)

  else
      if ( currfrag == frag_to_plot ) then
          grho = 0.d0

          !call gradrho( dfftp%nnr, rho_fde%of_g(:,1), ngm, g, dfftp%nl, grho )
          call fft_gradient_g2r(dfftp, rho_fde%of_g(:,1), g, grho)

          do ipol = 1, 3
            if ( ipol == 1 ) then
               idir = 'x'
            elseif ( ipol == 2 ) then
               idir = 'y'
            elseif ( ipol == 3 ) then
               idir = 'z'
            endif

            call gather_grid( dfftp, grho(ipol,:), raux )


            if (ionode) then
               call plot_io (trim(prefix)//'_grad'//trim(idir)//'_fde'//trim(channel)//'.pp', &
                      'Gradrho', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
                      dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
                      gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux, +1)
            endif

          enddo
      endif


  endif

  if ( currfrag == frag_to_plot ) then
     deallocate( grho )
     if (ionode) deallocate(raux)
  endif

  if (linterlock) deallocate( grho_large )


  return

END SUBROUTINE fde_plot_gradrho


!----------------------------------------------------------------------------
SUBROUTINE fde_plot_allpot
  !----------------------------------------------------------------------------
  !
  !  Plot FDE total density and potential
  !
  USE fft_base,                 ONLY : dfftp
  USE cell_base,                ONLY : celldm, at, ibrav
  USE ions_base,                ONLY : ntyp => nsp, atm, zv
  USE gvect,                    ONLY : gcutm, ngm, g
  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode, stdout
  USE io_files,                 ONLY : prefix
  USE mp_images,                ONLY : my_image_id
  USE scf,                      ONLY : rho, rho_core, rhog_core
  use input_parameters,         only : fde_xc_funct
  USE fde
  use vlocal,                   only : strf
  use lsda_mod,                 only : nspin
  use cell_base,                only : omega
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  integer, external :: find_free_unit
  integer :: unit, is, i
  logical :: exst
  real(dp), allocatable :: raux(:), raux1(:)
  real(dp), allocatable :: vemb(:,:), aux(:,:)
!  complex(dp), allocatable :: environ_rhog(:,:)!, environ_rhog_core(:,:)
  real(dp) :: trash0, trash1, trash2, trash3
  character(len=7)  :: channel, image_label
  character(len=6), external :: int_to_char
  character(len=5) :: potlabels(4), spinlabels(2)

  potlabels(1) = 'vloc_'
  potlabels(2) = 'vxc_'
  potlabels(3) = 'vh_'
  potlabels(4) = 'vt_'

  spinlabels(1) = '_alpha'
  spinlabels(2) = '_beta'

  write(stdout,*) potlabels

!  allocate( vemb(dfftp%nnr,fde_frag_nspin), environ_rhog(ngm,fde_frag_nspin) )
  if (ionode) allocate(raux1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  allocate( aux(dfftp%nnr,fde_frag_nspin) )
  image_label = 'rho' // int_to_char(my_image_id)

  do i = 1 , 4

    select case (i)
      case(1)
        aux(:,:) = 0.d0
        call setlocal_fde(aux(:,1), strf)
        call gather_grid( dfftp,aux(:,1),raux1)
        if ( ionode ) &
        call plot_io(trim(prefix)//'_'//trim(potlabels(i))//trim(image_label)//'.pp', &
              trim(potlabels(i))//trim(image_label), dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at, &
              gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)

        if (currfrag==1) then
          aux(:,:) = 0.d0
          call setlocal_fde(aux(:,1), strf_fde)
          call gather_grid( dfftp,aux(:,1),raux1)
          if ( ionode ) &
          call plot_io (trim(prefix)//'_'//trim(potlabels(i))//'rhotot'//'.pp',&
              trim(potlabels(i))//'rhotot', dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at,&
              gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
        endif


      case(2)

        aux(:,:) = 0.d0
        call v_xc(rho, rho_core, rhog_core, trash0, trash1, aux)
        call gather_grid( dfftp,aux(:,1),raux1)
        if ( ionode ) then
           do is = 1 , nspin
              call plot_io(trim(prefix)//'_'//trim(potlabels(i))//trim(image_label)//trim(spinlabels(is))//'.pp', &
                 trim(potlabels(i))//trim(image_label), dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
                 dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at, &
                 gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
           enddo
        endif

        if (currfrag==1) then
          aux(:,:) = 0.d0
          call v_xc(rho_fde, rho_core_fde, rhog_core_fde, trash0, trash1, aux)
          call gather_grid( dfftp,aux(:,1),raux1)
          if ( ionode ) then
             do is = 1 , fde_nspin
                call plot_io (trim(prefix)//'_'//trim(potlabels(i))//'rhotot'//trim(spinlabels(is))//'.pp',&
                    trim(potlabels(i))//'rhotot', dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
                    dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at,&
                    gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
             enddo
          endif
        endif

      case(3)

        aux(:,:) = 0.d0
        call v_h( rho%of_g, trash0, trash1, aux )
        call gather_grid( dfftp,aux(:,1),raux1)
        if ( ionode ) then
           do is = 1 , nspin
              call plot_io(trim(prefix)//'_'//trim(potlabels(i))//trim(image_label)//trim(spinlabels(is))//'.pp', &
                 trim(potlabels(i))//trim(image_label), dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
                 dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at, &
                 gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
           enddo
        endif

        if (currfrag==1) then
          aux(:,:) = 0.d0
          call v_h( rho_fde%of_g, trash0, trash1, aux )
          call gather_grid( dfftp,aux(:,1),raux1)
          if ( ionode ) then
             do is = 1 , fde_nspin
                call plot_io (trim(prefix)//'_'//trim(potlabels(i))//'rhotot'//trim(spinlabels(is))//'.pp',&
                    trim(potlabels(i))//'rhotot', dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
                    dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at,&
                    gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
             enddo
          endif
        endif

      case(4)

        aux(:,:) = 0.d0
        call fde_kin( rho, rho_gauss, rhog_gauss, trash0, aux, reduced_cell, NonlocalKernel)!dfftp, ngm, g, nl, omega, .false. )
        call gather_grid( dfftp,aux(:,1),raux1)
        if ( ionode ) then
             do is = 1 , nspin
              call plot_io(trim(prefix)//'_'//trim(potlabels(i))//trim(image_label)//trim(spinlabels(is))//'.pp', &
                    trim(potlabels(i))//trim(image_label), dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
                    dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at, &
                    gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
             enddo
        endif

        if (currfrag==1) then
          aux(:,:) = 0.d0
          call fde_kin( rho_fde, rho_gauss_fde, rhog_gauss_fde, trash0, aux, reduced_cell, NonlocalKernel)!dfftp, ngm, g,nl, omega, .false. )
          call gather_grid( dfftp,aux(:,1),raux1)
          if ( ionode ) then
             do is = 1 , fde_nspin
                call plot_io (trim(prefix)//'_'//trim(potlabels(i))//'rhotot'//trim(spinlabels(is))//'.pp',&
                    trim(potlabels(i))//'rhotot', dfftp%nr1x, dfftp%nr2x,dfftp%nr3x, &
                    dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav,celldm,at,&
                    gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)
             enddo
          endif
        endif

    end select

  enddo



  deallocate(aux)
  if (ionode) deallocate(raux1)
!  deallocate(vemb,environ_rhog)

  return
END SUBROUTINE fde_plot_allpot


!----------------------------------------------------------------------------
SUBROUTINE plot_dense( invec, label )
  !----------------------------------------------------------------------------
  !
  !  Plot FDE total density and potential
  !
  use kinds, only : dp
  USE fft_base,                 ONLY : dfftp
  USE cell_base,                ONLY : celldm, at, ibrav
  USE ions_base,                ONLY : ntyp => nsp, atm, zv, tau, ityp, nat
  USE gvect,                    ONLY : gcutm, ngm
  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode
  USE io_files,                 ONLY : prefix
  USE mp_global,                ONLY : my_image_id
  use fde,                      only : ityp_fde, tau_fde, nat_fde, currfrag
  use scatter_mod,          only : gather_grid, scatter_grid

  implicit none
  real(dp) , intent(in) :: invec(dfftp%nnr)
  character(len=16) , intent(in) :: label
  integer, external :: find_free_unit
  character(20) :: charfrag
  integer :: unit, is
  logical :: exst
  real(dp), allocatable :: raux(:), raux1(:)

  write(charfrag,'(i1)')currfrag

!  allocate(raux(dfftp%nnr))
  if (ionode) allocate(raux1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

  call gather_grid( dfftp,invec, raux1)
  if ( ionode ) &
      call plot_io (trim(label), &
              'Dense Grid', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat, ntyp, ibrav, celldm, at, &
              gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, raux1, +1)

!  deallocate(raux)
  if (ionode) deallocate(raux1)

  return
END SUBROUTINE plot_dense

!----------------------------------------------------------------------------
SUBROUTINE plot_large( invec, label )
  !----------------------------------------------------------------------------
  !
  !  Plot FDE total density and potential
  !
  use kinds, only : dp
  USE fft_base,                 ONLY : dfftp => dfftl
  USE large_cell_base,                ONLY : celldm, at, ibrav
  USE ions_base,                ONLY : ntyp => nsp, atm, zv
  USE gvecl,                    ONLY : gcutm, ngm
  USE gvecs,                    ONLY : dual
  USE gvecw,                    ONLY : ecutwfc
  USE io_global,                ONLY : ionode
  USE io_files,                 ONLY : prefix
  USE mp_global,                ONLY : my_image_id
  use fde,                      only : ityp_fde, tau_fde, nat_fde , currfrag
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none
  real(dp) , intent(in) :: invec(dfftp%nnr)
  character(len=16) , intent(in) :: label
  character(len=1)  :: charfrag
  integer, external :: find_free_unit
  integer :: unit, is
  logical :: exst
  real(dp), allocatable :: raux(:), raux1(:)

  write(charfrag,'(i1)')currfrag

!  allocate(raux(dfftp%nnr))
  if (ionode) allocate(raux1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))

  call gather_grid( dfftp,invec, raux1)
  if ( ionode ) &
      call plot_io (trim(label), &
              'Large Grid', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, &
              dfftp%nr1,  dfftp%nr2,  dfftp%nr3, nat_fde, ntyp, ibrav, celldm, at, &
              gcutm, dual, ecutwfc, 0, atm, ityp_fde, zv, tau_fde, raux1, +1)

!  deallocate(raux)
  if (ionode) deallocate(raux1)

  return
END SUBROUTINE plot_large

!----------------------------------------------------------------------------
SUBROUTINE c_grid_gather_sum_scatter( vec )
  !----------------------------------------------------------------------------
  !
  ! ... gathers nproc distributed data on the first processor of every pool,
  ! ....... sum the gathered quantity across all the fragments,
  ! ....... scatter the gathered and summed data to the processors of every pool
  !
  ! ... COMPLEX*16  vec  = distributed variable (ngm)
  !
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : ngm, ngm_g, ig_l2g
  USE fft_base,  ONLY : dfftp
  USE mp,        ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_images, ONLY : inter_fragment_comm
  USE io_global, only : ionode, ionode_id
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), intent(inout) :: vec(ngm)
  COMPLEX(DP), allocatable :: gathered_vec( : )
  !
#if defined (__MPI)
  !
  INTEGER :: i, idx
  !
  ! Gathering

  allocate(gathered_vec(ngm_g))

  gathered_vec(:) = cmplx(0.d0, 0.d0, kind=dp)

  do i = 1, ngm
    idx = ig_l2g(i)
    gathered_vec(idx) = vec(i)
  enddo

  call mp_sum(gathered_vec, dfftp%comm)

  ! Summing across fragments

  if (ionode) call mp_sum(gathered_vec, inter_fragment_comm)

  ! Scattering
  call mp_bcast(gathered_vec, ionode_id, dfftp%comm)
  do i = 1, ngm
    idx = ig_l2g(i)
    vec(i) = gathered_vec(idx)
  end do

  deallocate(gathered_vec)

#else
  CALL errore('c_grid_gather_sum_scatter', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE c_grid_gather_sum_scatter

!----------------------------------------------------------------------
SUBROUTINE setlocal_fde(pot, strf)
  !----------------------------------------------------------------------
  !
  !    This routine computes the local potential in real space vltot(ir)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  USE ions_base, ONLY : zv, ntyp => nsp
  USE cell_base, ONLY : omega
  USE extfield,  ONLY : tefield, dipfield, etotefield
  USE gvect,     ONLY : igtongl, gg
!  USE scf,       ONLY : rho, v_of_0, vltot
  USE vlocal,    ONLY : vloc !, strf
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : invfft
  USE gvect,     ONLY : ngm
  USE control_flags, ONLY : gamma_only
!  USE mp_bands,  ONLY : intra_bgrp_comm
!  USE mp_pools,  ONLY : intra_pool_comm
!  USE mp,        ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_loc, do_comp_mt
  USE esm,       ONLY : esm_local, esm_bc, do_comp_esm
!  USE qmmm,      ONLY : qmmm_add_mm_field
!  USE fde,       ONLY : strf_fde

  !
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT)  :: pot(dfftp%nnr)
  complex(dp) , intent(in) :: strf(ngm,ntyp)
  COMPLEX(DP), ALLOCATABLE :: aux (:), v_corr(:)!, strf_embed(:,:)
  ! auxiliary variable
  INTEGER :: nt, ng
  ! counter on atom types
  ! counter on g vectors
  !
  ALLOCATE (aux( dfftp%nnr))
  !ALLOCATE ( strf_embed(ngm, ntyp) )
  aux(:)=(0.d0,0.d0)
  !strf_embed(:,:) = strf_fde(:,:)
  !strf_embed(:,:) = strf_fde(:,:) - strf(:,:)
  !
  IF (do_comp_mt) THEN
      ALLOCATE(v_corr(ngm))
      CALL wg_corr_loc(omega,ntyp,ngm,zv,strf,v_corr)
      aux(dfftp%nl(:)) = v_corr(:)
      DEALLOCATE(v_corr)
  END IF
  !
  do nt = 1, ntyp
     do ng = 1, ngm
        aux (dfftp%nl(ng))=aux(dfftp%nl(ng)) + vloc (igtongl (ng), nt) * strf (ng, nt)
     enddo
  enddo
  IF (gamma_only) THEN
      DO ng = 1, ngm
          aux (dfftp%nlm(ng)) = CONJG(aux (dfftp%nl(ng)))
      END DO
  END IF
  !
  IF ( do_comp_esm .AND. ( esm_bc .NE. 'pbc' ) ) THEN
     !
     ! ... Perform ESM correction to local potential
     !
      CALL esm_local ( aux )
  ENDIF
  !
  ! ... v_of_0 is (Vloc)(G=0)
  !
!  v_of_0=0.0_DP
!  IF (gg(1) < eps8) v_of_0 = DBLE ( aux (nl(1)) )
  !
!  CALL mp_sum( v_of_0, intra_bgrp_comm )
  !
  ! ... aux = potential in G-space . FFT to real space
  !
  CALL invfft ('Rho', aux, dfftp)
  !
  pot (:) =  DBLE (aux (:) )
  !
  ! ... If required add an electric field to the local potential
  !
!  IF ( tefield .AND. ( .NOT. dipfield ) )  &
!      CALL add_efield(vltot,etotefield,rho%of_r,.TRUE.)
  !
  !  ... Add the electrostatic field generated by MM atoms
  !  in a QM/MM calculation to the local potential
  !
!  CALL qmmm_add_mm_field()
  !
  DEALLOCATE(aux)
  !DEALLOCATE(strf_embed)
  !
  RETURN
END SUBROUTINE setlocal_fde

SUBROUTINE copy_rho_core_to_rho_gauss
USE fde
USE scf, only: rho_core, rhog_core
USE io_global, only: stdout
IMPLICIT NONE
  rho_gauss = rho_core
  rhog_gauss = rhog_core
  rho_gauss_fde = rho_core_fde
  rhog_gauss_fde= rhog_core_fde
END SUBROUTINE copy_rho_core_to_rho_gauss

SUBROUTINE copy_rho_gauss_to_rho_core
USE fde
USE scf, only: rho_core, rhog_core
USE io_global, only: stdout
IMPLICIT NONE
  rho_core = rho_gauss
  rhog_core = rhog_gauss
  rho_core_fde = rho_gauss_fde
  rhog_core_fde = rhog_gauss_fde
END SUBROUTINE copy_rho_gauss_to_rho_core

!----------------------------------------------------------------------
SUBROUTINE setlocal_fde_large(pot, strf)
  !----------------------------------------------------------------------
  !
  !    This routine computes the local potential in real space vltot(ir)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  USE ions_base, ONLY : zv, ntyp => nsp
  USE large_cell_base, ONLY : omega
  USE extfield,  ONLY : tefield, dipfield, etotefield
  USE gvecl,     ONLY : igtongl, gg
!  USE scf,       ONLY : rho, v_of_0, vltot
  USE vlocal_large,    ONLY : vloc => vloc_large !, strf
  USE fft_base,  ONLY : dfftl, dfftp => dfftl
  USE fft_interfaces,ONLY : invfft
  USE gvecl,     ONLY : ngm
  USE control_flags, ONLY : gamma_only
!  USE mp_bands,  ONLY : intra_bgrp_comm
!  USE mp_pools,  ONLY : intra_pool_comm
!  USE mp,        ONLY : mp_sum
  USE martyna_tuckerman, ONLY : wg_corr_loc, do_comp_mt
!  USE esm,       ONLY : esm_local, esm_bc, do_comp_esm
!  USE qmmm,      ONLY : qmmm_add_mm_field
!  USE fde,       ONLY : strf_fde

  !
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT)  :: pot(dfftl%nnr)
  complex(dp) , intent(in) :: strf(ngm,ntyp)
  COMPLEX(DP), ALLOCATABLE :: aux (:), v_corr(:)!, strf_embed(:,:)
  ! auxiliary variable
  INTEGER :: nt, ng
  ! counter on atom types
  ! counter on g vectors
  !
  ALLOCATE (aux( dfftl%nnr))
  !ALLOCATE ( strf_embed(ngm, ntyp) )
  aux(:)=(0.d0,0.d0)
  !strf_embed(:,:) = strf_fde(:,:)
  !strf_embed(:,:) = strf_fde(:,:) - strf(:,:)
  !
  IF (do_comp_mt) THEN
      ALLOCATE(v_corr(ngm))
      CALL wg_corr_loc(omega,ntyp,ngm,zv,strf,v_corr)
      aux(dfftp%nl(:)) = v_corr(:)
      DEALLOCATE(v_corr)
  END IF
  !
  do nt = 1, ntyp
     do ng = 1, ngm
        aux (dfftp%nl(ng))=aux(dfftp%nl(ng)) + vloc (igtongl (ng), nt) * strf (ng, nt)
     enddo
  enddo
  IF (gamma_only) THEN
      DO ng = 1, ngm
          aux (dfftp%nlm(ng)) = CONJG(aux (dfftp%nl(ng)))
      END DO
  END IF
  !
!  IF ( do_comp_esm .AND. ( esm_bc .NE. 'pbc' ) ) THEN
!     !
!     ! ... Perform ESM correction to local potential
!     !
!      CALL esm_local ( aux )
!  ENDIF
  !
  ! ... v_of_0 is (Vloc)(G=0)
  !
!  v_of_0=0.0_DP
!  IF (gg(1) < eps8) v_of_0 = DBLE ( aux (nl(1)) )
  !
!  CALL mp_sum( v_of_0, intra_bgrp_comm )
  !
  ! ... aux = potential in G-space . FFT to real space
  !
  CALL invfft ('Rho', aux, dfftl)
  !
  pot (:) =  DBLE (aux (:) )
  !
  ! ... If required add an electric field to the local potential
  !
!  IF ( tefield .AND. ( .NOT. dipfield ) )  &
!      CALL add_efield(vltot,etotefield,rho%of_r,.TRUE.)
  !
  !  ... Add the electrostatic field generated by MM atoms
  !  in a QM/MM calculation to the local potential
  !
!  CALL qmmm_add_mm_field()
  !
  DEALLOCATE(aux)
  !DEALLOCATE(strf_embed)
  !
  RETURN
END SUBROUTINE setlocal_fde_large

SUBROUTINE calc_frag_beta (newmix,origmix,dr20,dr2)
!
! just trying to have a balanced mixing
! among the subsystems
!
USE kinds,     ONLY : DP
!USE fde,       ONLY : nfragments

implicit none

real(DP), intent(in)    :: dr20,dr2,origmix
real(DP), intent(inout) :: newmix

real(DP) :: alpha

     if (dr2 > 1.0d-10 ) then

      alpha = min(1.0d+4*dr2,1.0d0)

      newmix = alpha*(origmix*(1.0d0 - dr20/dr2)) + &
               origmix*(1.0d0 - alpha)

     else
      newmix = origmix
     endif

    newmix=min(abs(newmix),origmix)
    newmix=max(newmix,0.20d0*origmix) ! this is to avoid issued with DIIS
    return

END SUBROUTINE calc_frag_beta


! SUBROUTINE PartitionWeight(RegulRho, RegulRhoG, unitary)
!   !-----------------------------------------------------------------------
!   !
!   ! ADAPTED FROM TDDFT/molecular_operators
!   !
!   !
!   !  This routine will construct a sum of atomic densities
!   !  which properly decays exponentially away from the
!   !  nuclei.
!   !
!   !
!   !
!   !
!   USE kinds,        ONLY : dp
!   USE mp_global,    ONLY : me_pool, intra_bgrp_comm
!   USE mp,           ONLY : mp_sum
!   USE fft_base,     ONLY : dfftp, dffts
!   USE ions_base,    ONLY : ityp, zv, tau, nat
!   USE cell_base,    ONLY : at, bg, alat, omega, tpiba2
!   USE constants,    ONLY : bohr_radius_angs, tpi, tpi2
!   USE fde,          ONLY : tau_fde !, currfrag
!   USE gvect,        ONLY : g, gg, gstart

!   implicit none

!   real(dp), intent(out) :: RegulRho(dfftp%nnr) ! the partition function
!   complex(dp), intent(out) :: RegulRhoG(ngm) ! the partition function
!   logical, intent(in)   :: unitary

!   real(dp) :: wg(dfftp%nnr,nat) ! the partition function
!   real(dp) :: wgsum(nat), wgt
!   real(dp) :: alpha_exp, dist, domega, distx, disty, distz
!   real(dp) :: r(3)
!   real(dp) :: inv_nr1, inv_nr2, inv_nr3
!   real(dp) :: inv_nr1s, inv_nr2s, inv_nr3s
!   complex(dp) :: factor, imag
!   integer :: ia, i, j, k, index, index0, ir, ipol
! !  integer :: iu
! !  iu = 100+currfrag

!   inv_nr1 = 1.d0 / real(dfftp%nr1,dp)
!   inv_nr2 = 1.d0 / real(dfftp%nr2,dp)
!   inv_nr3 = 1.d0 / real(dfftp%nr3,dp)

!   index0 = 0
! #ifdef __MPI
!   do i = 1, me_pool
!     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
!   enddo
! #endif

!   RegulRho(:)   = 0.0d0
!   wg(:,:)  = 0.0d0
!   wgsum(:) = 0.0d0

!   domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

!   ! loop over real space grid

! !$omp parallel do
!   do ir = 1, dfftp%nnr
!     index = index0 + ir - 1
!     k     = index / (dfftp%nr1x*dfftp%nr2x)
!     index = index - (dfftp%nr1x*dfftp%nr2x)*k
!     j     = index / dfftp%nr1x
!     index = index - dfftp%nr1x*j
!     i     = index

!     do ipol = 1, 3
!       r(ipol) = real(i,dp)*inv_nr1*at(ipol,1) + &
!                 real(j,dp)*inv_nr2*at(ipol,2) + &
!                 real(k,dp)*inv_nr3*at(ipol,3)
!     enddo

!     ! minimum image convention
!     call cryst_to_cart( 1, r, bg, -1 )
!     r = r - anint(r)
!     call cryst_to_cart( 1, r, at, 1 )

! !  This is the density of the carbon atom
! !  but is made to integrate to the nuber
! !  of electrons of the given atom.
! !  If "unitary" the density is let integrate to 1


!      do ia = 1, nat
!       distx = r(1)-tau(1,ia)
!       distx = distx - anint(distx)
!       disty = r(2)-tau(2,ia)
!       disty = disty - anint(disty)
!       distz = r(3)-tau(3,ia)
!       distz = distz - anint(distz)
!       dist  = SQRT( distx*distx + disty*disty + distz*distz )

!       dist = dist * alat ! now dist is in Bohrs
!       alpha_exp = 3.5d0 ! effective decay for carbon's density
!       wg(ir,ia) = dist*dist*exp(-alpha_exp*dist)
!       wgsum(ia) = wgsum(ia) + wg(ir,ia) * domega
!      enddo
!   enddo
! !$omp end parallel do

!  ! is this right?
!  ! call mp_sum( wgsum, intra_bgrp_comm )

!   do ia = 1, nat
!    if (unitary) then
!      wg(:,ia) = wg(:,ia) / wgsum(ia)
!    else
!      wg(:,ia) = zv(ityp(ia)) * wg(:,ia) / wgsum(ia)
!    endif
!    RegulRho(:) = RegulRho(:) + wg(:,ia)
!   enddo

! if (.not.unitary) then
!  imag=(0.0d0,1.0d0)
!  do ia = 1, nat
!   do ir = 1,ngm

!    factor = exp( -tpi*imag* &
!             ( g(1,ir)*tau(1,ia) + &
!               g(2,ir)*tau(2,ia) + &
!               g(3,ir)*tau(3,ia) ) )

!    RegulRhoG(ir) = factor * zv(ityp(ia))/wgsum(ia) * 4.0d0*alpha_exp / &
!                    (alpha_exp*alpha_exp-tpi2*gg(ir)*tpiba2)**2.0d0
!   enddo
!  enddo
! endif
! !
! END SUBROUTINE PartitionWeight


! SUBROUTINE GlueRegRhoToRho(rho,RegularDensity,RegularDensityG)
!   USE kinds,        ONLY : dp
!   USE gvect,        ONLY : ngm, gg, gcutm
!   USE fft_base,     ONLY : dfftp
!   USE scf,          ONLY : scf_type
!   USE lsda_mod,     ONLY : nspin
!   USE spin_orb,     ONLY : domag
!   USE fde,          ONLY : currfrag


!   implicit none

!   type(scf_type), intent(inout)  :: rho

!   real(dp),    intent(in) :: RegularDensity(dfftp%nnr) ! regrho in real space
!   complex(dp), intent(in) :: RegularDensityG(ngm)      ! regrho in G space

!   integer  :: nspin0, ipnt, ispin
!   real(dp) :: fac

!   real(dp) :: thr_rho=1.0d-5 , thr_rhog=0.95d0 ! to be set at input?

!   ! TODO: non-collinear
!   nspin0 = nspin
!   if (nspin == 4) nspin0 = 1
!   if (nspin==4 .and. domag) nspin0 = 2

!   fac = 1.d0 / dble(nspin0)

!  ! regularize rho
!   do ispin = 1, nspin0
!    do ipnt=1,dfftp%nnr
!    if (abs(rho%of_r(ipnt,ispin)) < thr_rho) then
!      rho%of_r(ipnt,ispin) = fac * RegularDensity(ipnt)
!     endif
!    enddo
!      do ipnt=1,ngm
!       if (gg(ipnt)>gcutm*thr_rhog) then
!        rho%of_g(ipnt,1) = RegularDensityG(ipnt)*fac
!       endif
!      enddo
!   enddo

!   ! RegularDensity(1:ngm) = 0.0d0
!   ! RegularDensity(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr,1)
!   ! trash = 0.0d0
!   ! trash = SUM(RegularDensity(1:dfftp%nnr))*omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
!   ! write(100,*) 'CHARGE_FROM_REGRHO=',trash

! END SUBROUTINE GlueRegRhoToRho



! SUBROUTINE GradRegularRho(RegularRho, grad, is_fde)
!   !-----------------------------------------------------------------------
!   !
!   !  This routine will construct the derivative of
!   !  a sum of atomic densities
!   !  which properly decays exponentially away from the
!   !  nuclei.
!   !
!   !
!   !
!   !
!   USE kinds,        ONLY : dp
!   USE mp_global,    ONLY : me_pool, intra_bgrp_comm
!   USE mp,           ONLY : mp_sum
!   USE fft_base,     ONLY : dfftp
!   USE ions_base,    ONLY : ityp, zv, tau, nat
!   USE cell_base,    ONLY : at, bg, alat, omega
!   USE constants,    ONLY : bohr_radius_angs, tpi, tpi2
!   USE fde,          only : currfrag, tau_fde, nat_fde

!   implicit none

!   logical , intent(in)  :: is_fde
!   real(dp), intent(in)  :: RegularRho(dfftp%nnr) ! the partition function
!   real(dp), intent(out) :: grad(3,dfftp%nnr) ! the partition function

!   real(dp), allocatable :: tau_tmp(:,:)
!   real(dp), allocatable :: wg(:,:), gwg(:,:,:)
!   real(dp), allocatable :: wgsum(:)
!   real(dp) :: alpha_exp, dist, domega, distx, disty, distz, wgt
!   real(dp) :: r(3)
!   real(dp) :: inv_nr1, inv_nr2, inv_nr3
!   real(dp) :: inv_nr1s, inv_nr2s, inv_nr3s
!   integer :: ia, i, j, k, index, index0, ir, ipol
!   integer :: iu, nat_tmp


!   iu = 100+currfrag

!  if (is_fde) then
!   allocate(tau_tmp(3,nat_fde))
!   tau_tmp(1:3,1:nat_fde) = tau_fde(1:3,1:nat_fde)
!   nat_tmp = nat_fde
!  else
!   allocate(tau_tmp(3,nat))
!   tau_tmp(1:3,1:nat) = tau(1:3,1:nat)
!   nat_tmp = nat
!  endif

!  allocate (wg(dfftp%nnr,nat_tmp),gwg(3,dfftp%nnr,nat_tmp),wgsum(nat_tmp))



!   inv_nr1 = 1.d0 / real(dfftp%nr1,dp)
!   inv_nr2 = 1.d0 / real(dfftp%nr2,dp)
!   inv_nr3 = 1.d0 / real(dfftp%nr3,dp)

!   index0 = 0
! #ifdef __MPI
!   do i = 1, me_pool
!     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
!   enddo
! #endif

!   wg(:,:)  = 0.0d0
!   wgsum(:) = 0.0d0

!   domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

!   ! loop over real space grid

! !$omp parallel do
!   do ir = 1, dfftp%nnr
!     index = index0 + ir - 1
!     k     = index / (dfftp%nr1x*dfftp%nr2x)
!     index = index - (dfftp%nr1x*dfftp%nr2x)*k
!     j     = index / dfftp%nr1x
!     index = index - dfftp%nr1x*j
!     i     = index

!     do ipol = 1, 3
!       r(ipol) = real(i,dp)*inv_nr1*at(ipol,1) + &
!                 real(j,dp)*inv_nr2*at(ipol,2) + &
!                 real(k,dp)*inv_nr3*at(ipol,3)
!     enddo

!     ! minimum image convention
!     call cryst_to_cart( 1, r, bg, -1 )
!     r = r - anint(r)
!     call cryst_to_cart( 1, r, at, 1 )

! !  This is the density of the carbon atom
! !  but is made to integrate to the nuber
! !  of electrons of the given atom.
! !  If "unitary" the density is let integrate to 1


!      do ia = 1, nat_tmp
!       distx = r(1)-tau_tmp(1,ia)
!       distx = distx - anint(distx)
!       disty = r(2)-tau_tmp(2,ia)
!       disty = disty - anint(disty)
!       distz = r(3)-tau_tmp(3,ia)
!       distz = distz - anint(distz)
!       dist  = SQRT( distx*distx + disty*disty + distz*distz )

!       dist = dist * alat ! now dist is in Bohrs
!       alpha_exp = 3.5d0 ! effective decay for carbon's density

!       wg(ir,ia) = dist*dist*exp(-alpha_exp*dist)

!       gwg(1,ir,ia) = alat*distx*(2.0d0-alpha_exp*dist)*exp(-alpha_exp*dist)
!       gwg(2,ir,ia) = alat*disty*(2.0d0-alpha_exp*dist)*exp(-alpha_exp*dist)
!       gwg(3,ir,ia) = alat*distz*(2.0d0-alpha_exp*dist)*exp(-alpha_exp*dist)

!       wgsum(ia) = wgsum(ia) + wg(ir,ia) * domega
!      enddo
!   enddo
! !$omp end parallel do

!   do ia = 1, nat_tmp
!    gwg(:,:,ia) = zv(ityp(ia)) * gwg(:,:,ia) / wgsum(ia)
!   enddo

! !$omp parallel do
!   do ir=1,dfftp%nnr
!    !write(iu,*) 'Rho,Grad,RegGrad', RegularRho(ir), grad(1,ir), sum(gwg(1,ir,:))
!    if (RegularRho(ir) < 1.0d-3) then
!       grad(1,ir) = sum(gwg(1,ir,:))
!       grad(2,ir) = sum(gwg(2,ir,:))
!       grad(3,ir) = sum(gwg(3,ir,:))
!    endif
!   enddo
! !$omp end parallel do


!  deallocate (wg,gwg,wgsum)
!  deallocate(tau_tmp)
! !
! END SUBROUTINE GradRegularRho
!
!
!--------------------------------------------------------------------------------------------------------------
SUBROUTINE get_local_offset(offset,shift,split)
!--------------------------------------------------------------------------------------------------------------
!
USE kinds,              ONLY: DP                 !double-precision kind (selected_real_kind(14,200))
USE ions_base,          ONLY: tau                !positions of atoms/ions
USE ions_base,          ONLY: nat                !number of total atoms/ions (all atomic species)
USE fde,                ONLY: linterlock         !logical determining whether local and global grids are the same
USE fde,                only: tau_large       ! positions of the fragment's ions
                                              ! in the periodicity of the small cell
USE large_cell_base,    ONLY: atl => at !h           !h matrix for converting between r and s coordinates via r = h s (global grid) ***analog of h***
USE large_cell_base,    ONLY: atlinv => bg     !h^-1 matrix for converting between r and s coordinates via s = h^-1 r (global grid) ***analog of ainv***
USE large_cell_base,    ONLY: alat
USE fft_base,           ONLY: dfftl              !FFT derived data type for the large grid
USE fft_base,           ONLY: dfftp              !FFT derived data type for the fragment grid
use io_global, only : stdout
!
IMPLICIT NONE
!
! I/O variables
!
INTEGER, INTENT(INOUT) :: offset(3)
INTEGER, INTENT(INOUT) :: shift(3)
real(dp), INTENT(IN) :: split(3)
!
! Local variables
!
INTEGER :: ia, nr(3), nrl(3), ipol, ja, ref_atm, &
           next_atm, prev_atm, closest_pair(2)
REAL(DP) :: offset_r(3),offset_s(3), min_atxyzs(3), max_atxyzs(3), center_atxyzs(3), &
            size_atxyzs(3), ref_atm_coord(3)
!REAL(DP), DIMENSION(:,:), ALLOCATABLE :: atxyz,atxyzs
REAL(DP) :: atxyz(3,nat)
REAL(DP) :: atxyzs(3,nat)
!REAL(DP) :: dA(3), dAs(3), origin(3), originS(3)
real(dp) :: dr_matrix(nat,nat,3), dd_matrix(nat,nat), dd, dr(3), min_max_dd
logical :: coord_not_set(nat), mask_matrix(nat,nat)

!
! Initialization of offset vector...
!
offset=0
shift=0
atxyz=0.0_DP
atxyzs=0.0_DP
dr_matrix = 0.d0
dd_matrix = 0.d0
mask_matrix = .false.
!
tau = tau_large
!
! Determine the atom which is most representative for the center of the molecule


do ia = 1, nat
  do ja = ia+1, nat
    dr = dr_mic(tau(:,ia), tau(:,ja), atl, atlinv, alat)
    dd = sqrt( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) )
    !write(stdout,*) dr, dd
    dr_matrix(ia, ja, :) =  dr
    dr_matrix(ja, ia, :) = -dr
    dd_matrix(ia, ja) = dd
    dd_matrix(ja, ia) = dd
  enddo
enddo

ref_atm = -1
min_max_dd = 1000.d0
do ia = 1, nat
  !write(stdout,*) tau_large(:,ia)
  !write(stdout,*) d_matrix(ia,:)
  dd = maxval(dd_matrix(ia,:))
  dd_matrix(ia,ia) = 1000.d0
  if (dd - min_max_dd < -1.d-4) then
    min_max_dd = dd
    ref_atm = ia
  endif
enddo

write(stdout,*) "ref_atm",ref_atm
coord_not_set(:) = .true.
ref_atm_coord = dr_mic(tau(:,ref_atm), (/0.d0,0.d0,0.d0/), atl, atlinv, alat)
atxyz(:,ref_atm) = ref_atm_coord
atxyzs(:,ref_atm) = r2s( atxyz(:,ref_atm), atlinv, alat )
coord_not_set(ref_atm) = .false.


! Reconstruct the connectivity of the molecule...
do while (any(coord_not_set(:)))
  call update_mask(coord_not_set, mask_matrix, nat)
  closest_pair = minloc(dd_matrix, mask=mask_matrix)
  !write(*,*) closest_pair
  if (coord_not_set(closest_pair(1))) then
    next_atm = closest_pair(1)
    prev_atm = closest_pair(2)
  else
    next_atm = closest_pair(2)
    prev_atm = closest_pair(1)
  endif

  atxyz(:, next_atm) = atxyz(:,prev_atm) + dr_matrix(next_atm, prev_atm, :)
  atxyzs(:,next_atm) = r2s( atxyz(:,next_atm), atlinv, alat )
  coord_not_set(next_atm) = .false.

enddo




!do ia = 1 , nat
!  dr = dr_mic(tau(:,ia), ref_atm_coord, atl, atlinv, alat)
!  !write(stdout,*) dr
!  atxyz(:,ia) = dr + ref_atm_coord
!  atxyzs(:,ia) = r2s( atxyz(:,ia), atlinv, alat )
!enddo

! Calculate the optimal positioning of the fragment cell...

nrl(1)=dfftl%nr1; nrl(2)=dfftl%nr2; nrl(3)=dfftl%nr3
nr(1)=dfftp%nr1; nr(2)=dfftp%nr2; nr(3)=dfftp%nr3

do ipol = 1, 3
  min_atxyzs(ipol) = minval(atxyzs(ipol,:))
  max_atxyzs(ipol) = maxval(atxyzs(ipol,:))
  center_atxyzs(ipol) = 0.5d0 * ( min_atxyzs(ipol) + max_atxyzs(ipol) )
  size_atxyzs(ipol) = ( max_atxyzs(ipol) - min_atxyzs(ipol) )
  offset(ipol) = nint( center_atxyzs(ipol) * nrl(ipol)) - (nr(ipol) / 2)
  offset_s(ipol) = dble(offset(ipol))/dble(nrl(ipol))
enddo

offset_r(:) = s2r(offset_s, atl, alat)

do ia = 1, nat
  tau(:,ia) = atxyz(:,ia) - offset_r(:)
enddo
!
!sf
if ( split(1) == 1.d0 .and. split(2) == 1.d0 .and. split(3) == 1.d0) then
!if (.false.) then
  do ipol = 1, 3
    shift(ipol) = 0
    offset(ipol) = 0
    tau(ipol,:) = tau_large(ipol,:)
  enddo
  endif

!write(stdout,*) min_atxyzs(:)
!write(stdout,*) max_atxyzs(:)
write(stdout,*) size_atxyzs(:)
!write(stdout,*) center_atxyzs(:)
write(stdout,*) offset(:)

RETURN

CONTAINS

subroutine update_mask(bool_vec, bool_mat, n)
  implicit none
  integer, intent(in) :: n
  logical, intent(in) :: bool_vec(n)
  logical, intent(inout) :: bool_mat(n,n)

  integer :: i, j

  do i = 1, n
    bool_mat(i,i) = .false.
    do j = i+1, n
       if ( (bool_vec(i) .and. .not. bool_vec(j)) .or. &
            (bool_vec(j) .and. .not. bool_vec(i)) ) then
         bool_mat(i,j) = .true.
         bool_mat(j,i) = .true.
       else
         bool_mat(i,j) = .false.
         bool_mat(j,i) = .false.
       endif
    enddo
  enddo

end subroutine update_mask

function r2s( rpos1, bg, alat )

  use kinds, only : dp
  implicit none

  real(dp) :: rpos1(3), bg(3,3), alat
  real(dp) :: r2s(3)
  integer  :: i,j

  r2s(:) = 0.d0

  do i = 1, 3
    do j = 1,3
      ! for some reason bg's columns and rows are swapped...
      r2s(i) = r2s(i) + bg(j,i) * rpos1(j)
    enddo
  enddo

  !r2s = r2s / alat

  return

end function


function s2r( spos1, at, alat )

  use kinds, only : dp
  implicit none

  real(dp) :: spos1(3), at(3,3), alat
  real(dp) :: s2r(3)
  integer  :: i,j

  s2r(:) = 0.d0

  do i = 1, 3
    do j = 1,3
      s2r(i) = s2r(i) + at(i,j) * spos1(j)
    enddo
  enddo

  !s2r = s2r * alat

  return

end function


!--------------------------------------------------------------------------------------------------------------
function dr_mic( rpos1, rpos2, at, bg, alat )
!--------------------------------------------------------------------------------------------------------------
  use kinds, only : dp
  implicit none

  real(dp) :: rpos1(3), rpos2(3), at(3,3), bg(3,3), alat
  real(dp) :: dr_mic(3)
  !
  real(dp) :: spos1(3), spos2(3), ds(3)
  integer :: i,j

  spos1(:) = r2s( rpos1, bg, alat )
  spos2(:) = r2s( rpos2, bg, alat )

  ds(:) = spos1(:) - spos2(:)
  do i = 1, 3
    ds(i)=ds(i)-IDNINT(ds(i))
  enddo

  dr_mic(:) = s2r( ds, at, alat )

  return

end function dr_mic


!
!--------------------------------------------------------------------------------------------------------------
END SUBROUTINE get_local_offset
!--------------------------------------------------------------------------------------------------------------

subroutine potential_wall( dfftp, pot )
!
!
!
  use kinds , only : dp
  use io_global, only : stdout, ionode
  USE fft_types,  ONLY : fft_type_descriptor
  use scatter_mod,          only : gather_grid, scatter_grid
  implicit none

  type(fft_type_descriptor), intent(in) :: dfftp
  real(dp), intent(inout) :: pot(:)


  integer :: i, j, k, dummy, face_points, counter
  integer :: nr(3)
  integer, allocatable :: faces_idx(:)
  real(dp), allocatable :: raux(:)


  nr(1) = dfftp%nr1
  nr(2) = dfftp%nr2
  nr(3) = dfftp%nr3

  if (ionode) then
  face_points = 4*nr(1)*nr(2) + 4*nr(1)*nr(3) + 4*nr(2)*nr(3)
  allocate(faces_idx(face_points))

  counter = 0

  ! xy faces
  do dummy = 1 , 4
     k = dummy
     if ( dummy == 3) k = nr(3)-1
     if ( dummy == 4) k = nr(3)
     do j = 1, nr(2)
        do i = 1, nr(1)
           counter = counter + 1
           faces_idx(counter) = i + (j-1)*nr(1) + (k-1)*nr(1)*nr(2)
        enddo
     enddo
  enddo




  ! xz faces
  do dummy = 1 , 4
     j = dummy
     if ( dummy == 3) j = nr(2)-1
     if ( dummy == 4) j = nr(2)
     do k = 1, nr(3)
        do i = 1, nr(1)
           counter = counter + 1
           faces_idx(counter) = i + (j-1)*nr(1) + (k-1)*nr(1)*nr(2)
        enddo
     enddo
  enddo

  ! yz faces
  do dummy = 1 , 4
     i = dummy
     if ( dummy == 3) i = nr(1)-1
     if ( dummy == 4) i = nr(1)
     do k = 1, nr(3)
        do j = 1, nr(2)
           counter = counter + 1
           faces_idx(counter) = i + (j-1)*nr(1) + (k-1)*nr(1)*nr(2)
        enddo
     enddo
  enddo

  write(stdout,*) "POTENTIAL WALL  ", counter, face_points

  allocate(raux(nr(1)*nr(2)*nr(3)))
  raux=0.d0
  raux(faces_idx(:)) = 1.d3

  deallocate(faces_idx)

  endif

  call scatter_grid( dfftp,raux, pot)

  if (ionode) deallocate(raux)


end subroutine potential_wall







!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grididx_f2g ( fragidx , globidx, dfftp, dfftl, shift_, offset_ )
!
! Given the index of the grid point in the fragment cell, calculate the index in the global cell
!
  USE fft_types,  ONLY : fft_type_descriptor
  use io_global , only : stdout
  implicit none

!  integer :: grididx_f2g
  integer , intent(in) :: fragidx
  integer , intent(out) :: globidx
  integer  :: shift(3)
  integer   :: offset(3)
  integer   :: lolimit(3), hilimit(3)
  integer , intent(in) :: shift_(3)
  integer , intent(in) :: offset_(3)
  type(fft_type_descriptor) , intent(in) :: dfftp, dfftl
  integer :: ijk_f(3) , ijk_g(3)

   ijk_f = 0
   ijk_g = 0

   offset = offset_
   shift = shift_

   if (dfftp%nr1 == dfftl%nr1) then
      !shift(1) = 0
      !offset(1) = 0
   endif

   if (dfftp%nr2 == dfftl%nr2) then
      !shift(2) = 0
      !offset(2) = 0
   endif

   if (dfftp%nr3 == dfftl%nr3) then
      !shift(3) = 0
      !offset(3) = 0
   endif


   ijk_f(1) = mod(fragidx , dfftp%nr1x)
   if ( ijk_f(1) == 0 ) ijk_f(1) = dfftp%nr1x

   ijk_f(2) = mod(((fragidx-1)/dfftp%nr1x + 1) , dfftp%nr2x)
   if ( ijk_f(2) == 0 ) ijk_f(2) = dfftp%nr2x

   ijk_f(3) = mod(((fragidx-1)/(dfftp%nr1x * dfftp%nr2x) + 1) , dfftp%nr3x)
   if ( ijk_f(3) == 0 ) ijk_f(3) = dfftp%nr3x

   ijk_f = ijk_f - 1

!  ijk_g(1) = mod(ijk_f(1) + shift(1), dfftl%nr1x )
!  ijk_g(2) = mod(ijk_f(2) + shift(2), dfftl%nr2x )
!  ijk_g(3) = mod(ijk_f(3) + shift(3), dfftl%nr3x )

!  ijk_g(1) = mod((ijk_f(1) + shift(1)), dfftl%nr1)
!  ijk_g(2) = mod((ijk_f(2) + shift(2)), dfftl%nr2)
!  ijk_g(3) = mod((ijk_f(3) + shift(3)), dfftl%nr3)

!  ijk_g(1) = (ijk_f(1) + shift(1))
!  ijk_g(2) = (ijk_f(2) + shift(2))
!  ijk_g(3) = (ijk_f(3) + shift(3))

  ijk_g(1) = (ijk_f(1) + offset(1))
  ijk_g(2) = (ijk_f(2) + offset(2))
  ijk_g(3) = (ijk_f(3) + offset(3))
!  limit(1) = mod((offset(1) + dfftp%nr1), dfftl%nr1)
!  limit(2) = mod((offset(2) + dfftp%nr2), dfftl%nr2)
!  limit(3) = mod((offset(3) + dfftp%nr3), dfftl%nr3)
  !
!  hilimit(1) = (offset(1) + dfftp%nr1)
!  hilimit(2) = (offset(2) + dfftp%nr2)
!  hilimit(3) = (offset(3) + dfftp%nr3)
  !
!  lolimit(1) = offset(1)
!  lolimit(2) = offset(2)
!  lolimit(3) = offset(3)

!  if (ijk_g(1) >= hilimit(1) ) then
!    ijk_g(1) = ijk_g(1) - dfftp%nr1
!  elseif (ijk_g(1) < lolimit(1) ) then
!    ijk_g(1) = ijk_g(1) + dfftp%nr1
!  endif
!
!  if (ijk_g(2) >= hilimit(2) ) then
!    ijk_g(2) = ijk_g(2) - dfftp%nr2
!  elseif (ijk_g(2) < lolimit(2) ) then
!    ijk_g(2) = ijk_g(2) + dfftp%nr2
!  endif
!
!  if (ijk_g(3) >= hilimit(3) ) then
!    ijk_g(3) = ijk_g(3) - dfftp%nr3
!  elseif (ijk_g(3) < lolimit(3) ) then
!    ijk_g(3) = ijk_g(3) + dfftp%nr3
!  endif

  ijk_g(1) = modulo(ijk_g(1) , dfftl%nr1 )
  ijk_g(2) = modulo(ijk_g(2) , dfftl%nr2 )
  ijk_g(3) = modulo(ijk_g(3) , dfftl%nr3 )

  ijk_g = ijk_g +1

!  if (ijk_g(1) == 0 ) ijk_g(1) = dfftl%nr1
!  if (ijk_g(2) == 0 ) ijk_g(2) = dfftl%nr2
!  if (ijk_g(3) == 0 ) ijk_g(3) = dfftl%nr3

  globidx = ijk_g(1) + &
                ( ijk_g(2) - 1 ) * dfftl%nr1x + &
                ( ijk_g(3) - 1 ) * dfftl%nr1x * dfftl%nr2x

!  write(stdout,*) ' shift: ', shift, ' offs : ' , offset
!  write(stdout,*) ' dense: ', ijk_f, ' large: ' , ijk_g

  return ! grididx_f2g

end subroutine grididx_f2g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grididx_g2f ( globidx , fragidx, dfftp, dfftl, shift_, offset_ )
! Given the index of the grid point in the global cell, calculate the index in the fragment cell
!
! DO NOT USE, use grididx_f2g instead
  USE fft_types,  ONLY : fft_type_descriptor

  implicit none

  !integer :: grididx_g2f
  integer , intent(in) :: globidx
  integer , intent(out) :: fragidx
  integer , intent(in) :: shift_(3)
  integer , intent(in) :: offset_(3)
  integer :: offset(3)
  integer :: shift(3)
  type(fft_type_descriptor) , intent(in) :: dfftp, dfftl
  integer :: ijk_f(3) , ijk_g(3)

   offset = offset_
   shift = shift_

   if (dfftp%nr1 == dfftl%nr1) then
!      shift(1) = 0
!      offset(1) = 0
   endif

   if (dfftp%nr2 == dfftl%nr2) then
!      shift(2) = 0
!      offset(2) = 0
   endif

   if (dfftp%nr3 == dfftl%nr3) then
!      shift(3) = 0
!      offset(3) = 0
   endif

   ijk_f = 0
   ijk_g = 0

   ijk_g(1) = mod(globidx , dfftl%nr1x)
   if ( ijk_g(1) == 0 ) ijk_g(1) = dfftl%nr1x

   ijk_g(2) = mod(((globidx-1)/dfftl%nr1x + 1) , dfftl%nr2x)
   if ( ijk_g(2) == 0 ) ijk_g(2) = dfftl%nr2x

   ijk_g(3) = mod(((globidx-1)/(dfftl%nr1x * dfftl%nr2x) + 1) , dfftl%nr3x)
   if ( ijk_g(3) == 0 ) ijk_g(3) = dfftl%nr3x

   ijk_f(1) = mod(ijk_g(1) - offset(1), dfftl%nr1)
   ijk_f(2) = mod(ijk_g(2) - offset(2), dfftl%nr2)
   ijk_f(3) = mod(ijk_g(3) - offset(3), dfftl%nr3)
   if ( ijk_f(1) == 0 ) ijk_f(1) = dfftl%nr1x
   if ( ijk_f(2) == 0 ) ijk_f(2) = dfftl%nr2x
   if ( ijk_f(3) == 0 ) ijk_f(3) = dfftl%nr3x

   if ( ijk_f(1) < 1 .or. &
        ijk_f(2) < 1 .or. &
        ijk_f(3) < 1 .or. &
        ijk_f(1) > dfftp%nr1x .or. &
        ijk_f(2) > dfftp%nr2x .or. &
        ijk_f(3) > dfftp%nr3x ) &
   then
      fragidx = -1
   else
!      if( ijk_g(1) < shift(1) ) then
!        ijk_f(1) = ijk_g(1) + dfftp%nr1x
!      else
        ijk_f(1) = mod(ijk_g(1) - shift(1), dfftp%nr1)
!      endif
!      if( ijk_g(2) < shift(2) ) then
!        ijk_f(2) = ijk_g(2) + dfftp%nr2x
!      else
        ijk_f(2) = mod(ijk_g(2) - shift(2), dfftp%nr2)
       ijk_f(2) = ijk_g(2)
!      endif
!      if( ijk_g(3) < shift(3) ) then
!        ijk_f(3) = ijk_g(3) + dfftp%nr3x
!      else
        ijk_f(3) = mod(ijk_g(3) - shift(3), dfftp%nr3)
       ijk_f(3) = ijk_g(3)
!      endif
      fragidx = ijk_f(1) + &
                ( ijk_f(2) - 1 ) * dfftp%nr1x + &
                ( ijk_f(3) - 1 ) * dfftp%nr1x * dfftp%nr2x

   endif

  return ! grididx_g2f

end subroutine grididx_g2f


subroutine calc_f2l(f2l, dfftp, dfftl, shift, offset)
   USE fft_types,  ONLY : fft_type_descriptor
   use scatter_mod,          only : gather_grid, scatter_grid

   implicit none

   type(fft_type_descriptor) , intent(in) :: dfftp, dfftl
   integer, intent(in)  :: shift(3)
   integer, intent(in)  :: offset(3)
   integer, intent(out) :: f2l(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
   integer :: idx, idx2, discarded, npoints

   discarded = 0
   npoints = 0

!   do idx = 1, dfftl%nr1 * dfftl%nr2 * dfftl%nr3
!     call grididx_g2f( idx , idx2, dfftp, dfftl, shift, offset )
     ! if idx2 is -1 it means that the large grid point does not belong to the fragment grid
!     if ( idx2 > 0 ) l2f(idx2) = idx
!     if ( idx2 < 1 ) discarded = discarded + 1
!     npoints = npoints + 1
!   enddo

!  write(stdout,*) 'NPOINTS = ' , npoints , discarded

   do idx = 1, dfftp%nr1 * dfftp%nr2 * dfftp%nr3
     call grididx_f2g ( idx, idx2, dfftp, dfftl, shift, offset )
     f2l(idx) = idx2
   enddo

   return

end subroutine calc_f2l

subroutine copy_pot_f2l( vf, vl )
!  This subroutine copies scattered quantities (dens or pot) in frag cell array over to scattered quantities to the supersystem array.
!
  use kinds, only : dp
  use fft_base , only: dfftp, dfftl
  use fde , only : f2l
  use io_global, only : ionode
  use mp, only: mp_sum
  use mp_images, only: inter_fragment_comm
  use scatter_mod,          only : gather_grid, scatter_grid

  real(dp) , intent(in) :: vf(dfftp%nnr)
  real(dp) , intent(inout) :: vl(dfftl%nnr)
  real(dp), allocatable :: raux(:), rauxl(:)
  real(dp), allocatable :: aux(:), auxl(:)


  if (ionode) allocate(raux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  call gather_grid( dfftp,vf, raux)

  if (ionode) then
    allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
    rauxl = 0.d0
    rauxl(f2l(:)) = raux(:)
    deallocate(raux)
  endif

!  if (ionode) call mp_sum(raux1, inter_fragment_comm)
  call scatter_grid( dfftl,rauxl, vl)
  if (ionode) deallocate(rauxl)

!
!  vl = 0.d0
!
!  do ir = 1 , dfftp%nr1*dfftp%nr2*dfftp%nr3
!    vl( f2l(ir) ) = vf(ir)
!  enddo
!
  return
!
end subroutine copy_pot_f2l


subroutine copy_pot_l2f( vl, vf )
!  This subroutine copies scattered quantities (dens or pot) in large cell array over to scattered quantities to the fragment array.
!
  use kinds, only : dp
  use fft_base , only: dfftp, dfftl
  use fde , only : f2l
  use io_global, only : ionode
  use mp, only: mp_bcast
  use mp_images, only: inter_fragment_comm
  use command_line_options,    only : fancy_parallel_
  use scatter_mod,          only : gather_grid, scatter_grid


  real(dp) , intent(inout) :: vf(dfftp%nnr)
  real(dp) , intent(inout) :: vl(dfftl%nnr)
  real(dp), allocatable :: raux(:), rauxl(:)
  real(dp), allocatable :: aux(:), auxl(:)


  if (ionode) then
    allocate(rauxl(dfftl%nr1x*dfftl%nr2x*dfftl%nr3x))
    !rauxl = 0.d0
  endif

  call gather_grid( dfftl,vl, rauxl)

  if (ionode) then
! if we're not using fancy parallelization, each image has its own rauxl
    if (fancy_parallel_) call mp_bcast( rauxl, 0, inter_fragment_comm )
    allocate(raux(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
    raux = 0.d0
    raux(:) = rauxl(f2l(:))
    deallocate(rauxl)
  endif

!  if (ionode) call mp_sum(raux1, inter_fragment_comm)
  call scatter_grid( dfftp,raux, vf)
  if (ionode) deallocate(raux)
!
  return
!
end subroutine copy_pot_l2f

subroutine reg_tau( tau, nat, atl, atlinv )
 use kinds, only : dp

 implicit none
 real(dp) , intent(inout) :: tau(3,nat)
 integer , intent(in) :: nat
 real(dp) , intent(in) :: atl(3,3), atlinv(3,3)

 integer ia
 real(dp) :: atxyzs(3,nat), atxyz(3,nat)

DO ia = 1, nat
  !
  atxyzs(1,ia)=atlinv(1,1)*tau(1,ia)+atlinv(1,2)*tau(2,ia)+atlinv(1,3)*tau(3,ia)! s = h^-1 r
  atxyzs(2,ia)=atlinv(2,1)*tau(1,ia)+atlinv(2,2)*tau(2,ia)+atlinv(2,3)*tau(3,ia)! s = h^-1 r
  atxyzs(3,ia)=atlinv(3,1)*tau(1,ia)+atlinv(3,2)*tau(2,ia)+atlinv(3,3)*tau(3,ia)! s = h^-1 r
  !
  atxyzs(1,ia)=atxyzs(1,ia)-nint(atxyzs(1,ia))   ! impose PBC on s in range:[0,1)
  atxyzs(2,ia)=atxyzs(2,ia)-nint(atxyzs(2,ia))   ! impose PBC on s in range:[0,1)
  atxyzs(3,ia)=atxyzs(3,ia)-nint(atxyzs(3,ia))   ! impose PBC on s in range:[0,1)
  !
  atxyz(1,ia)=atl(1,1)*atxyzs(1,ia)+atl(1,2)*atxyzs(2,ia)+atl(1,3)*atxyzs(3,ia)! r = h s
  atxyz(2,ia)=atl(2,1)*atxyzs(1,ia)+atl(2,2)*atxyzs(2,ia)+atl(2,3)*atxyzs(3,ia)! r = h s
  atxyz(3,ia)=atl(3,1)*atxyzs(1,ia)+atl(3,2)*atxyzs(2,ia)+atl(3,3)*atxyzs(3,ia)! r = h s
  !
END DO

tau = atxyz

return

end subroutine reg_tau


end module fde_routines
