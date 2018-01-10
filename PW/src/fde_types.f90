module fde_types
!
! just a module containing a few custom types that point to 
! quantities in the QE global modules.
! makes the code a bit clearer, but it's still not ideal
!
use kinds, only: dp
use fft_types, only : fft_type_descriptor
implicit none

type simulation_cell
! this derived dt contains all the information about the fragment
! in a specific cell (i.e. either the supersystem native cell
! or the fragment electronic cell)
! for backward compatibility we use pointers instead of creating
! a properly private type
  
  real(dp), pointer :: tau(:,:)         ! positions
  integer, pointer :: nat, nsp, ityp(:)
  real(dp), pointer :: at(:,:), bg(:,:) ! lattice vectors 
  complex(dp), pointer :: strf(:,:), strf_fde(:,:) ! structure factors
  real(dp), pointer :: omega, alat, celldm(:), tpiba, tpiba2 ! volume
  type(fft_type_descriptor), pointer :: dfftp
  !
  ! all the gvect stuff
  !
  integer, pointer :: ngm, ngm_g, ngl, ngmx
  real(dp), pointer :: ecutrho, gcutm
  integer, pointer :: nl(:), nlm(:)
  integer, pointer :: gstart
  real(dp), pointer :: gg(:), gl(:), g(:,:)
  integer, pointer :: igtongl(:), mill(:,:), ig_l2g(:), &
                      sortedig_l2g(:), mill_g(:,:)
  complex(dp), pointer :: eigts1(:,:), eigts2(:,:), eigts3(:,:)
  !
  logical :: is_native_cell  ! is it the supersystem cell or the 
                             ! fragment cell?
  logical :: is_allocated=.false. 
end type

interface new
  module procedure new_cell
end interface

contains

subroutine new_cell(this, at, bg, omega, alat, celldm, &
                    tpiba, tpiba2, &
                    nat, tau, nsp, ityp, & 
                    strf, strf_fde,  dfftp, &
                    ngm, ngm_g, ngl, ngmx, &
                    ecutrho, gcutm, & 
                    nl, nlm, gstart, gg, gl, g, &
                    igtongl, mill, ig_l2g, sortedig_l2g, mill_g, &
                    eigts1, eigts2, eigts3, &
                    is_native_cell)
 
  implicit none

  type(simulation_cell) :: this
  real(dp), target, intent(in) :: at(:,:), bg(:,:), omega, alat, celldm(:)
  integer, target, intent(in) :: nat, nsp, ityp(:)
  real(dp), target, intent(in) :: tau(:,:), tpiba, tpiba2
  type(fft_type_descriptor), target, intent(in) :: dfftp
  complex(dp), target, intent(in) :: strf(:,:), strf_fde(:,:)
  integer, target, intent(in) :: ngm, ngm_g, ngl, ngmx
  real(dp), target, intent(in) :: ecutrho, gcutm
  integer, target, intent(in) :: gstart, nl(:), nlm(:)
  integer, target, intent(in) :: igtongl(:), mill(:,:), ig_l2g(:), &
                                 sortedig_l2g(:), mill_g(:,:)
  real(dp), target, intent(in) :: gg(:), gl(:), g(:,:)
  complex(dp), target, intent(in) :: eigts1(:,:), eigts2(:,:), eigts3(:,:)
  logical, intent(in) :: is_native_cell

  this%at => at
  this%bg => bg
  this%omega => omega
  this%alat => alat
  this%tpiba => tpiba
  this%tpiba2 => tpiba2
  this%celldm => celldm
  this%nat => nat
  this%tau => tau
  this%nsp => nsp
  this%ityp => ityp
  this%ngm => ngm
  this%ngm_g => ngm_g
  this%ngl => ngl
  this%ngmx => ngmx
  this%ecutrho => ecutrho
  this%gcutm => gcutm
  this%strf => strf
  this%strf_fde => strf_fde
  this%dfftp => dfftp
  this%gstart => gstart
  this%nl => nl
  this%nlm => nlm
  this%igtongl => igtongl
  this%mill => mill
  this%ig_l2g => ig_l2g
  this%sortedig_l2g => sortedig_l2g
  this%mill_g => mill_g
  this%gg => gg
  this%gl => gl
  this%g => g
  this%eigts1 => eigts1
  this%eigts2 => eigts2
  this%eigts3 => eigts3
  this%is_native_cell = is_native_cell

  this%is_allocated = .true.

end subroutine new_cell

end module fde_types
