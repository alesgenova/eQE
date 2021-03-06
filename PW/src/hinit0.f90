!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hinit0()
  !-----------------------------------------------------------------------
  !
  ! ... hamiltonian initialization: 
  ! ... atomic position independent initialization for nonlocal PP,
  ! ... structure factors, local potential, core charge
  !
  USE kinds,        ONLY : dp
  USE ions_base,    ONLY : nat, nsp, ityp, tau
  USE basis,        ONLY : startingconfig
  USE cell_base,    ONLY : at, bg, omega, tpiba2
  USE large_cell_base,    ONLY : atl => at, bgl => bg
  USE cellmd,       ONLY : omega_old, at_old, lmovecell
  USE klist,        ONLY : init_igk
  USE wvfct,        ONLY : npwx
  USE fft_base,     ONLY : dfftp
  USE fft_base,     ONLY : dfftl
  USE gvect,        ONLY : ngm, g, eigts1, eigts2, eigts3
  USE gvecl,        ONLY : ngml => ngm, &
                            gl => g, &
                            eigts1l => eigts1, &
                            eigts2l => eigts2, &
                            eigts3l => eigts3
  USE vlocal,       ONLY : strf
  USE gvecw,        ONLY : gcutw
  USE realus,       ONLY : generate_qpointlist,betapointlist,init_realspace_vars,real_space
  use ldaU,         ONLY : lda_plus_U, U_projection
  USE control_flags,ONLY : tqr, tq_smoothing, tbeta_smoothing
  USE io_global,    ONLY : stdout
  USE scf,          ONLY : vltot
  USE scf_large,    ONLY : vltot_large => vltot
  use fde,          only : do_fde, linterlock, fde_cell_offset, &
                            fde_cell_shift, frag_cell_split, &
                            nat_fde, tau_fde, ityp_fde, strf_fde, strf_fde_large, &
                            f2l, tau_large
use fde_routines
  !
  IMPLICIT NONE
  !
  INTEGER :: ik                 ! counter on k points
  REAL(dp), ALLOCATABLE :: gk(:) ! work space
  !
  CALL start_clock( 'hinit0' )
  !
  ! ... calculate the Fourier coefficients of the local part of the PP
  !
  CALL init_vloc()
  !
  if (do_fde) call init_vloc_large()
  !
  ! ... k-point independent parameters of non-local pseudopotentials
  !
  if (tbeta_smoothing) CALL init_us_b0()
  if (tq_smoothing) CALL init_us_0()
  CALL init_us_1()
  IF ( lda_plus_U .AND. ( U_projection == 'pseudo' ) ) CALL init_q_aeps()
  CALL init_at_1()
  !
  CALL init_igk ( npwx, ngm, g, gcutw )
  !
  IF ( lmovecell .AND. startingconfig == 'file' ) THEN
     !
     ! ... If lmovecell and restart are both true the cell shape is read from
     ! ... the restart file and stored. The xxx_old variables are used instead 
     ! ... of the current (read from input) ones.
     ! ... xxx and xxx_old are swapped, the atomic positions rescaled and 
     ! ... the hamiltonian scaled.
     !
     CALL cryst_to_cart( nat, tau, bg, - 1 )
     !
     CALL dswap( 9, at, 1, at_old, 1 )
     CALL dswap( 1, omega, 1, omega_old, 1 )
     !
     CALL cryst_to_cart( nat, tau, at, + 1 )
     !
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     CALL scale_h()
     !
  END IF
  !
  ! ... initialize the structure factor
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  if (do_fde ) then 
    !strf_fde(:,:) = strf(:,:)
    CALL struc_fact( nat_fde, tau_fde, nsp, ityp_fde, ngml, gl, bgl, &
              dfftl%nr1, dfftl%nr2, dfftl%nr3, strf_fde_large, eigts1l, eigts2l, eigts3l )
    call calc_f2l(f2l, dfftp, dfftl, fde_cell_shift, fde_cell_offset)
    call setlocal_fde_large(vltot_large, strf_fde_large)
  endif
  !
  ! these routines can be used to patch quantities that are dependent
  ! on the ions and cell parameters
  !
  CALL plugin_init_ions()
  CALL plugin_init_cell()
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()

  if (do_fde .and. linterlock) call copy_pot_l2f(vltot_large, vltot)
  !
  if (do_fde) call copy_pot_l2f(vltot_large, vltot)
  !
  ! ... calculate the core charge (if any) for the nonlinear core correction
  !
  CALL set_rhoc()
  !
  IF ( tqr ) CALL generate_qpointlist()
  
  IF (real_space ) then
   call betapointlist()
   call init_realspace_vars()
   write(stdout,'(5X,"Real space initialisation completed")')    
  endif
  !
  CALL stop_clock( 'hinit0' )
  !
  RETURN
  !
END SUBROUTINE hinit0

