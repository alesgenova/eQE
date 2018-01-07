!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE scf_large
  !  
  !  This module contains variables and auxiliary routines needed for
  !  the self-consistent cycle
  !
  !  ROUTINES: allocate_scf_type
  !
  USE kinds,      ONLY : DP
  !
  USE lsda_mod,     ONLY : nspin
  USE ldaU,         ONLY : lda_plus_u, Hubbard_lmax
  USE fde_module,    ONLY : nat_fde
  USE buffers,      ONLY : open_buffer, close_buffer, get_buffer, save_buffer
  USE funct,        ONLY : dft_is_meta
  USE fft_base,     ONLY : dfftp => dfftl
  USE fft_interfaces,ONLY: invfft
  USE gvecl,        ONLY : ngm
  !USE gvecs,        ONLY : ngms
  USE paw_variables,ONLY : okpaw
  USE uspp_param,   ONLY : nhm
  USE extfield,     ONLY : dipfield, emaxpos, eopreg, edir
  USE control_flags,ONLY : lxdm
  !
  SAVE
  INTERFACE charge
   Module Procedure charge_mix
   Module Procedure charge_scf
  END INTERFACE
  !
! Details of PAW implementation:
! NOTE: scf_type is used for two different quantities: density and potential.
!       These correspond, for PAW, to becsum and D coefficients.
!       Due to interference with the ultrasoft routines only the becsum part
!       is stored in the structure (at the moment).
!       This only holds for scf_type; mix_type is not affected.
! NOTE: rho%bec is different from becsum for two reasons:
!       1. rho%bec is mixed, while becsum is not
!       2. for npool > 1 rho%bec is collected, becsum is not
!          ( this is necessary to make the stress work)

  TYPE scf_type
     REAL(DP),   ALLOCATABLE :: of_r(:,:)  ! the charge density in R-space
     COMPLEX(DP),ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     REAL(DP),   ALLOCATABLE :: kin_r(:,:) ! the kinetic energy density in R-space
     COMPLEX(DP),ALLOCATABLE :: kin_g(:,:) ! the kinetic energy density in G-space
     REAL(DP),   ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix
     COMPLEX(DP),ALLOCATABLE :: ns_nc(:,:,:,:)!     ---       noncollinear case
     REAL(DP),   ALLOCATABLE :: bec(:,:,:) ! the PAW hamiltonian elements
     LOGICAL             :: is_fde
  END TYPE scf_type
  !
  TYPE mix_type
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)  ! the charge density in G-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:) ! the charge density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)! the LDA+U occupation matrix 
     COMPLEX(DP), ALLOCATABLE :: ns_nc(:,:,:,:)!     ---     noncollinear case 
     REAL(DP),    ALLOCATABLE :: bec(:,:,:) ! PAW corrections to hamiltonian
     REAL(DP)                 :: el_dipole  ! electrons dipole
  END TYPE mix_type

  type (scf_type) :: rho  ! the charge density and its other components
  type (scf_type) :: v    ! the scf potential
  type (scf_type) :: vnew ! used to correct the forces

  REAL(DP) :: v_of_0    ! vltot(G=0)      
  REAL(DP), ALLOCATABLE :: &
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grid)
       rho_core(:),    &! the core charge in real space
       kedtau(:,:)      ! position dependent kinetic energy enhancement factor
  COMPLEX(DP), ALLOCATABLE :: &
       rhog_core(:)     ! the core charge in reciprocal space

  INTEGER, PRIVATE  :: record_length, &
                       rlen_rho=0,  rlen_kin=0,  rlen_ldaU=0,  rlen_bec=0,&
                       rlen_dip=0, &
                       start_rho=0, start_kin=0, start_ldaU=0, start_bec=0, &
                       start_dipole=0
  ! DFT+U, colinear and noncolinear cases
  LOGICAL, PRIVATE :: lda_plus_u_co, lda_plus_u_nc
  COMPLEX(DP), PRIVATE, ALLOCATABLE:: io_buffer(:)
CONTAINS

 SUBROUTINE create_scf_type ( rho, do_not_allocate_becsum )
   IMPLICIT NONE
   TYPE (scf_type) :: rho
   LOGICAL,INTENT(IN),OPTIONAL :: do_not_allocate_becsum ! PAW hack
   LOGICAL                     :: allocate_becsum        ! PAW hack
   allocate ( rho%of_r( dfftp%nnr, nspin) )
   allocate ( rho%of_g( ngm, nspin ) )
   if (dft_is_meta() .or. lxdm) then
      allocate ( rho%kin_r( dfftp%nnr, nspin) )
      allocate ( rho%kin_g( ngm, nspin ) )
   else
      allocate ( rho%kin_r(1,1) )
      allocate ( rho%kin_g(1,1) )
   endif

   lda_plus_u_co = lda_plus_u .and. .not. (nspin == 4 )
   lda_plus_u_nc = lda_plus_u .and.       (nspin == 4 )
   if (lda_plus_u_co) allocate (rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
   if (lda_plus_u_nc) allocate (rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))

   if (okpaw) then ! See the top of the file for clarification
      if(present(do_not_allocate_becsum)) then
         allocate_becsum = .not. do_not_allocate_becsum
      else
         allocate_becsum = .true.
      endif
      if(allocate_becsum) allocate (rho%bec(nhm*(nhm+1)/2,nat,nspin))
   endif

  ! assume it is not the FDE density
  rho%is_fde=.false.
   
 return
 END SUBROUTINE create_scf_type

 SUBROUTINE destroy_scf_type ( rho )
   IMPLICIT NONE
   TYPE (scf_type) :: rho

   if (ALLOCATED(rho%of_r))  deallocate(rho%of_r)
   if (ALLOCATED(rho%of_g))  deallocate(rho%of_g)
   if (ALLOCATED(rho%kin_r)) deallocate(rho%kin_r)
   if (ALLOCATED(rho%kin_g)) deallocate(rho%kin_g)
   if (ALLOCATED(rho%ns))    deallocate(rho%ns)
   if (ALLOCATED(rho%ns_nc))    deallocate(rho%ns_nc)
   if (ALLOCATED(rho%bec))   deallocate(rho%bec)

   return
 END SUBROUTINE destroy_scf_type
 !
!----------------------------------------
END MODULE scf_large
