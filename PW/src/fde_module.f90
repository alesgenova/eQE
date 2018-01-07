!
! Copyright (C) 2007-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This module contains the variables needed for the Frozen Density Embedding
! method [1,2]. The FDE method is implemented thanks to image parallelization,
! and a serial version compilation is not possible. The SCF is performed on
! each fragment in parallel, and the non-additive potentials and energies
! are calculated at every SCF step.
!
! Coded by D. Ceresoli and A. Genova, during the April 2013 PRG hackathon
! in Rutgers-Newark (http://www.michelepavanello.com/pdf/hackathon.pdf)
!
! [1] T. A. Wesolowksi and A. Warshel, J. Phys. Chem. 97, 8050 (1993) 
! [2] see: $QEDIR/FDE-Notes/fde-equations.tex

!----------------------------------------------------------------------------
MODULE fde
  !--------------------------------------------------------------------------
  !
  use kinds,                    only: dp
  use fde_types,                only : simulation_cell
  use scf,                      only: scf_type
  use scf_large,                only: scf_type_large => scf_type
  use kernel_module,            only: KernelType

  implicit none 
  save


  ! .true. is nfragments > 0
  logical :: do_fde = .false.
  !
  logical :: fde_init_rho
  logical :: fde_fat 
  logical :: fde_thaw = .true.

  ! number of fragments
  integer :: nfragments = 1

  ! current fragment (=my_image_id+1)
  integer :: currfrag

  ! System and current fragment nspin
  integer :: fde_nspin, fde_frag_nspin

  ! FDE embedding potential is only updated after each fragment converged
  integer :: fde_fat_iter, fde_fat_maxstep

  ! the structure factor of the full system on the small cell
  complex(dp), allocatable :: strf_fde(:,:)
  
  ! the mixing factor of old and new density when using 'freeze and thaw' FDE
  real(dp) :: fde_fat_mixing

  ! the charge density of the full system
  type(scf_type) :: rho_fde
  type(scf_type) :: rho_fde_old
  ! the charge density of the full system in the large grid
  type(scf_type_large) :: rho_fde_large
  
  ! the valence density of the single fragment at the previous scf step
  type(scf_type) :: rho_old 
  ! the valence density of the single fragment at the previous thaw step
  type(scf_type) :: rho_thaw_old
 
  ! the core charge of the full system, for NLCC
  real(dp), allocatable :: rho_core_fde(:)
  complex(dp), allocatable :: rhog_core_fde(:)

  ! the core charge of the full system, for NLCC on the large grid
  real(dp), allocatable :: rho_core_fde_large(:)
  complex(dp), allocatable :: rhog_core_fde_large(:)

  ! the gaussian core charge of the fragment
  real(dp), allocatable :: rho_gauss(:)
  complex(dp), allocatable :: rhog_gauss(:)

  ! the gaussian core charge of the full system
  real(dp), allocatable :: rho_gauss_fde(:)
  complex(dp), allocatable :: rhog_gauss_fde(:)
  
  ! the gaussian core charge of the full system on the large grid
  real(dp), allocatable :: rho_gauss_fde_large(:)
  complex(dp), allocatable :: rhog_gauss_fde_large(:)

  ! the energies:
  real(dp) :: etot_fde                ! total FDE energy
  real(dp) :: etot_sum                ! fragment energy sum
  real(dp) :: edc_fde                 ! double counting
  real(dp) :: ekin_nadd, ekin0_nadd   ! non-additive kinetic energy (total and single fragment)
  real(dp) :: etxc_nadd, etxc0_nadd   ! non-additive XC energy
  real(dp) :: fde_true_ehart          ! the true total Hartree energy

  ! convergence variables
  real(dp) :: dr20
  real(dp) :: dr2max
  real(dp) :: dr2_fat = 1.0_dp
  real(dp) :: tr2_fat
  integer  :: fde_cycle = 0
  integer  :: current_scf_iteration = 0
  real(dp) :: current_scf_nelec
  real(dp) :: current_scf_nelec_up
  real(dp) :: current_scf_nelec_dw

  ! the atoms of the full system
  integer :: nat_fde, fixatom_fde
  real(dp), allocatable :: tau_fde(:,:)
  real(dp), allocatable :: tau_large(:,:)
  real(dp), allocatable :: force_fde(:,:)
  integer, allocatable :: ityp_fde(:)
  integer, allocatable :: if_pos_fde(:,:)
  integer, allocatable :: atm_offset(:)

  ! Diffusion Coeff.
  real(dp) :: diff_coeff_fde

  ! total number of electrons
  real(dp) :: nelec_fde
  ! nelec at start of calculation (needed only if fde_fractional = .true.)
  real(dp) :: original_nelec
  real(dp) :: original_nelec_up
  real(dp) :: original_nelec_dw
  ! charge (nelec - z_nuc)
  real(dp) :: fde_frag_charge, fde_total_charge

  ! non-local KEDF parameters
!  real(dp), parameter :: alpha_kin = (5.d0 + sqrt(5.d0))/6.d0 
!  real(dp), parameter :: beta_kin = (5.d0 - sqrt(5.d0))/6.d0
  real(dp), parameter :: alpha_kin = (5.d0 )/6.d0 
  real(dp), parameter :: beta_kin = (5.d0 )/6.d0

  ! export the total density and embedding potential flags
  logical :: fde_print_density = .false.
  logical :: fde_print_density_frag = .false.
  logical :: fde_print_density_frag_large = .false.
  logical :: fde_print_embedpot = .false.
  logical :: fde_print_allpot = .false.
  logical :: fde_print_electro = .false.
  integer, allocatable  :: fde_plotemb_vec(:)
  integer, allocatable  :: fde_plotgrad_vec(:)
  integer, allocatable  :: fde_printdensity_vec(:)

  ! use or not compensating gassians
  logical :: use_gaussians = .false.
  
  ! constant for Vandevondele's SIE correction 
  logical  :: fde_si = .false.
  integer, allocatable  :: fde_si_vec(:)
  logical  :: fde_si_all2all = .false.
  real(dp) :: fde_si_alpha = 0.d0

  ! rho0 in nonlocal Ts
  real(dp) :: fde_r0 = 0.d0

  ! amout of Ts_TF in GP functional
  real(dp) :: fde_gp = 0.d0

  ! amout of Ts_TF in GP functional
  real(dp) :: fde_gp_rhot = 0.d0

  ! amout of Ts_TF in GP functional
  real(dp) :: fde_gp_alpha = 0.d0

  ! flag for anisotropic fragment mixing
  logical  :: fde_split_mix = .false.

  ! flag for anisotropic fragment mixing
  logical  :: fde_regrho = .false.

  ! flag for the kin_funct type: gga or nl
  logical :: fde_kin_is_nl = .false.

  ! flag for inclusion of repulsive overlap potential
  logical  :: fde_overlap = .false.
  real(dp) :: fde_overlap_c = 0.0d0
  ! flag for fractional subsystem occupations
  logical :: fde_fractional = .false.
  logical :: fde_fractional_onlyalpha = .false.
  logical :: fde_fractional_onlybeta = .true.
  real(dp):: fde_fractional_mixing = 0.0d0
  real(dp):: fde_fractional_minEtransfer = 0.0d0
  real(dp):: fde_fractional_maxEtransfer = 0.0d0
  integer :: fde_fractional_cycle = 1

  ! Kernel for nonlocal Ts
  type(KernelType) :: NonlocalKernel
  type(KernelType) :: NonlocalKernel_fde

  ! variables for the interlocking fragments cells
  logical :: linterlock = .false.
  real(dp) :: frag_cell_split(3) = 1.d0  
  integer, parameter, dimension(10,2) :: fde_split_types = reshape( &
                                               (/ 1, 1, 1, 2, 1, 3, 2, 3, 4, 1, &
                                                  5, 4, 3, 5, 2, 5, 3, 4, 5, 1 /),&
                                               (/10,2/) )

  integer :: fde_cell_offset(3) = 0
  integer :: fde_cell_shift(3) = 0
  integer :: fde_frag_split_type(3) = 10
  integer :: fde_max_divisor(3) = 1
  integer , allocatable :: f2l(:) , l2f(:) ! masks to map the fragment point to the points in the large grid

  ! strf of total system on large cell 
  complex(dp), allocatable :: strf_fde_large(:,:)
  ! strf of subsystem system on large cell 
  complex(dp), allocatable :: strf_large(:,:)
  !
  real(dp) , allocatable :: pot_wall(:)
  !
  type(simulation_cell) :: native_cell, reduced_cell
  !
  logical :: fde_dotsonlarge = .true.
  !REAL(DP) :: atl(3,3)    = 0.0_DP ! h matrix of the supersystem cell (a.u.) 
  !REAL(DP) :: atlinv(3,3) = 0.0_DP ! inverse of the above matrix

  ! calculate v_of_rho becomes compatible with linterlock if calc_ener and icalc_ener ask to compute the energy
  logical  :: calc_ener_in_v_of_rho=.false. 

  ! SAOP Parameters
  logical :: saop_add ! use saop for the additive part
  logical :: saop_nadd  ! use saop for the nonadditive part
  real(dp) :: saop_hirho = -0.5d0
  real(dp) :: saop_lorho = -7.0d0
  real(dp) :: saop_pow = -1.0d0
  real(dp) :: saop_frac = -1.0d0

  ! An auxiliary array of the large always available
  complex(dp), allocatable :: psic_large(:)


contains
SUBROUTINE fde_fake_nspin ( flag , nspin0 )
  !----------------------------------------------------------------------------
  !
  !  Hacky and dangerous subroutine to temporarely change nspin of eac fragment to that of rho_fde
  !  needed when computing nonadditive potential and energies
  !
  use lsda_mod,                 only : nspin

  implicit none

  logical, intent(in)     :: flag
  integer, optional       :: nspin0

  if ( present(nspin0) ) then
     nspin = nspin0
  else
    if ( flag ) then
       nspin = fde_nspin
    else
       nspin = fde_frag_nspin
    endif

  endif

  return
END SUBROUTINE fde_fake_nspin
END MODULE fde


