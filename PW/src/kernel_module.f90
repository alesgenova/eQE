module kernel_module
!
! Module to build nonlocal Ts kernels
!
! M.P., Sept. 2015 
!

  use kinds,                    only : dp

  implicit none

!================================================================
!           KINETIC LUMP TYPE
!================================================================

  Type KineticLumpType

   ! is kinetic lump?
   logical :: isLump

   ! the type of lump
   character :: TheLump 

   ! fraction of lump to be included (a parameter in input)
   real(dp) :: LumpFraction, LumpAlpha

   ! Lump in G-space
   real(dp), allocatable :: LumpFunction(:)

  End Type KineticLumpType

!================================================================
!               KERNEL TYPE
!================================================================

  Type KernelType

   ! is already generated?
   logical :: isAllocated = .false.

   ! Maximum of functional integration points
   integer :: MaxPoints

   ! number of different LOCAL Fermi wavevectors (GPLDA only)
   integer :: nkf

   ! the charge of the (sub) system
   real(dp) :: charge

   ! the rho_0 parameter
   real(dp) :: rho_star

   ! type of nonlocal Ts
   real(dp) :: ThomasFermiFraction

   ! type of nonlocal Ts
   character(len=80) :: TheFunctional

   ! The Kernel in G-space
   !
   ! Must be used in conjunction with 
   ! c_grid_gather_sum_scatter( w(:,kf) )
   !
   real(dp), allocatable :: KernelFunction(:,:)

   ! if GPLDA need a series of kf values
   real(dp), allocatable :: kfs(:)

   ! this is for handling the Kinetic Lump
   type (KineticLumpType) :: KL

  End Type KernelType

!================================================================
end module kernel_module








