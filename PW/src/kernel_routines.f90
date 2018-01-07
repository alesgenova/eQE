!
!
! Routines for the nonlocal Ts kernel
! MP, Sept. 2015
!

!=========================================================================
!=========================================================================
module kernel_routines
implicit none
contains
SUBROUTINE CreateLump(Lump,charge,cell)
  use kinds,                    only : dp
  use constants,                only : pi
  use kernel_module,            only : KineticLumpType
  use fde_types,                only : simulation_cell
  !use gvect,                    only : ngm, gg

  real(dp), intent(in)                  :: charge
  type(simulation_cell), intent(in)     :: cell
  type (KineticLumpType), intent(inout) :: Lump
  !
  integer :: ngm
  real(dp) :: tpiba2
  real(dp), pointer :: gg(:)

  ngm = cell%ngm
  tpiba2 = cell%tpiba2
  gg => cell%gg

   Lump%isLump = .false.
   if (Lump%LumpAlpha.ne.0.0d0) Lump%isLump = .true.
   
   if (Lump%isLump) then
       ! to be read from input in future versions
       Lump%TheLump      = 'GAUSSIAN'
       
       allocate(Lump%LumpFunction(ngm))
       
       select case (Lump%TheLump)
         case ('GAUSSIAN')
           ! lump = exp(-a g^2)
           Lump%LumpFunction(1:ngm) = (Lump%LumpFraction/charge)* &
                                     exp(-Lump%LumpAlpha*gg(1:ngm)*tpiba2)/4.d0/pi 
       end select 
   endif

END SUBROUTINE CreateLump


!=========================================================================

SUBROUTINE AddLump (Kernel)
  use kernel_module,            only : KernelType 

  type(KernelType), intent(inout) :: Kernel

  integer :: ikf

  if (Kernel%KL%isLump) then
  do ikf = 1, Kernel%nkf
    !Kernel%KernelFunction(1:ngm,ikf) = Kernel%KernelFunction(1:ngm,ikf) + Kernel%KL%LumpFunction(1:ngm)
    Kernel%KernelFunction(:,ikf) = Kernel%KernelFunction(:,ikf) + Kernel%KL%LumpFunction(:)
  enddo
  endif

END SUBROUTINE AddLump

!=========================================================================

SUBROUTINE CreateKernel(Kernel, rho, NAKE, gp, rhot, r0, alpha, cell)
  use kinds,                    only : dp
  use fft_types,                only : fft_dlay_descriptor
  use fde_types,                only : simulation_cell
  use constants,                only : pi
!  use fft_base,                 only : dfftp
  USE fft_interfaces,           ONLY : fwfft, invfft
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum, mp_max
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
!  use cell_base,                only : at, alat, tpiba, tpiba2, omega
!  use gvect,                    only : ngm, nl, g, gg, gstart
  use io_global,                only : ionode, stdout
  use kernel_module,            only : KernelType


  type (KernelType), intent(inout) :: Kernel
  type (scf_type),intent(in)       :: rho
  character(len=80), intent(in)    :: NAKE
  real(dp), intent(inout)          :: gp, rhot, r0, alpha
  type(simulation_cell)            :: cell

  CHARACTER(len=80) :: nonlocal_list(2)
  DATA nonlocal_list / 'GPRHO0', 'GPLDA' /
  
  integer :: i
  logical :: isNonlocal

  type(fft_dlay_descriptor), pointer :: dfftp
  real(dp), pointer :: at(:,:), alat, tpiba, tpiba2, omega
  integer, pointer :: ngm, nl(:), gstart
  real(dp), pointer :: g(:,:), gg(:)

  at => cell%at
  alat => cell%alat
  tpiba => cell%tpiba
  tpiba2 => cell%tpiba2
  omega => cell%omega
  ngm => cell%ngm
  nl => cell%nl
  gstart => cell%gstart
  g => cell%g
  gg => cell%gg
  

  isNonlocal = .false.
  DO i = 1, SIZE( nonlocal_list )
     IF( TRIM( NAKE ) .eq. nonlocal_list(i) ) then
      !this is a nonlocal calculation!
      isNonlocal = .true.
     ENDIF
  END DO

  if (isNonlocal) then
  !   if (ionode) write (stdout,*) " ", trim(NAKE), " is a NONLOCAL NAKE"
     if (ionode) write (stdout,*) ""
  else
  !   if (ionode) write (stdout,*) " ", trim(NAKE), " is a (semi)LOCAL NAKE"
     return
  endif

  Kernel%TheFunctional = NAKE

  ! to be set in input
  Kernel%MaxPoints     = 500

  ! calculate charge of fragment
   Kernel%charge = 0.d0
   if ( gstart == 2 ) then
      Kernel%charge = omega*real( rho%of_g(1,1) )
      if ( nspin == 2 ) Kernel%charge = Kernel%charge + omega*real( rho%of_g(1,2) )
   end if

  call mp_sum(Kernel%charge, intra_bgrp_comm)

  Kernel%ThomasFermiFraction = gp
  Kernel%KL%LumpFraction = rhot
  Kernel%rho_star        = r0
  Kernel%KL%LumpAlpha    = alpha

  call CreateLump (Kernel%KL, Kernel%charge, cell)

  select case (Kernel%TheFunctional)

     case ('GPRHO0')
       Kernel%nkf = 1
       allocate ( Kernel%KernelFunction(ngm,1) )
       call SimpleKernel(Kernel%KernelFunction,Kernel%MaxPoints,Kernel%rho_star, &
                         Kernel%ThomasFermiFraction, cell)

     case ('GPLDA')
       ! will need to be read from input
       Kernel%nkf = 20
       allocate (Kernel%kfs(Kernel%nkf))
       allocate ( Kernel%KernelFunction(ngm,Kernel%nkf) )
       call GPLDAKernel ( Kernel%KernelFunction, rho, Kernel%kfs, Kernel%nkf,  &
                          Kernel%Maxpoints, Kernel%ThomasFermiFraction, cell) 

  end select

  call AddLump(Kernel)

  Kernel%isAllocated = .true.

END SUBROUTINE CreateKernel

!=========================================================================

SUBROUTINE CopyLump (LumpOrig,LumpCopy, cell)
  use kernel_module,            only : KineticLumpType
  use fde_types,                only : simulation_cell
   type (KineticLumpType), intent(in)  :: LumpOrig
   type (KineticLumpType), intent(out) :: LumpCopy
   type (simulation_cell), intent(in) :: cell

   integer :: ngm

   ngm = cell%ngm

   if (.not.LumpOrig%isLump) then 
       LumpCopy%isLump = .false.
       return
   endif

   LumpCopy%isLump       = .true.
   LumpCopy%TheLump      = LumpOrig%TheLump
   LumpCopy%LumpFraction = LumpOrig%LumpFraction 
   LumpCopy%LumpAlpha    = LumpOrig%LumpAlpha

   if (allocated(LumpOrig%LumpFunction)) then
      if (.not.allocated(LumpCopy%LumpFunction)) allocate(LumpCopy%LumpFunction(ngm))
      LumpCopy%LumpFunction = LumpOrig%LumpFunction
   endif
END SUBROUTINE CopyLump

!=========================================================================

SUBROUTINE CopyKernel (KernelOrig,KernelCopy, cell)
  use kernel_module,            only : KernelType
  use fde_types,                only : simulation_cell

   type (KernelType), intent(in)  :: KernelOrig
   type (KernelType), intent(out) :: KernelCopy
   type (simulation_cell), intent(in) :: cell

   integer :: ngm

   ngm = cell%ngm

    KernelCopy%MaxPoints     = KernelOrig%MaxPoints
    KernelCopy%nkf           = KernelOrig%nkf
    KernelCopy%charge        = KernelOrig%charge
    KernelCopy%TheFunctional = KernelOrig%TheFunctional

    if (KernelOrig%nkf .gt. 1) then
     if (allocated(KernelOrig%kfs)) then
       allocate(KernelCopy%kfs(KernelCopy%nkf))
       if (.not. allocated(KernelCopy%kfs)) KernelCopy%kfs = KernelOrig%kfs
     endif
    endif  

    if (allocated(KernelOrig%KernelFunction)) then
         allocate(KernelCopy%KernelFunction(ngm,KernelCopy%nkf))
         if (.not.allocated(KernelCopy%KernelFunction)) &
                            KernelCopy%KernelFunction = KernelOrig%KernelFunction
    endif   

    call CopyLump(KernelOrig%KL,KernelCopy%KL, cell)

END SUBROUTINE CopyKernel

!=========================================================================

SUBROUTINE DeleteKernel (Kernel)
  use kernel_module,            only : KernelType

  type (KernelType) :: Kernel

   if (allocated(Kernel%kfs)) deallocate(Kernel%kfs)
   if (allocated(Kernel%KernelFunction)) deallocate(Kernel%KernelFunction)

   call DeleteLump (Kernel%KL)

END SUBROUTINE DeleteKernel


!=========================================================================

SUBROUTINE DeleteLump (Lump)
  use kernel_module,            only : KineticLumpType

   type (KineticLumpType) :: Lump

    if (Lump%isLump) then
      if (allocated(Lump%LumpFunction)) deallocate (Lump%LumpFunction)
    endif

END SUBROUTINE DeleteLump


!=========================================================================

FUNCTION LindG(eta,lambda,mu)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the q-dependent part of WT kernel, as described
!   in [1] eq. (13) (and (11).)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
! Taken straight from PROFESS 3.0 - Thank you, Emily!
!
!------------------------------------------------------------------------------

  use kinds,       only : dp

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!  

  REAL(KIND=dp), INTENT(IN) :: eta
  ! the point at which the function gets evaluated.
  !
  REAL(KIND=dp), INTENT(IN) :: lambda
  ! the TF multiplier for compensating.
  !
  REAL(KIND=dp), INTENT(IN) :: mu
  ! the vW multiplier
  !

  REAL(KIND=dp) :: LindG
  ! The Lindhard G function as described in [1]
  !

                    !>> INTERNAL VARIABLES <<!  
  REAL(KIND=dp) :: eta2
  ! eta ** 2
  !
  REAL(KIND=dp) :: invEta2
  ! 1/eta**2
  !

                     !>> INITIALIZATION <<!
                     !>> FUNCTION BODY << 

  IF (eta<0._dp) THEN
    LindG = 0._dp

  ! Limit for small eta
  ELSE IF (eta < 1E-10_dp) THEN
    LindG = 1._dp - lambda + eta**2 * (1._dp / 3._dp - 3._dp * mu)

  ! Around the singularity
  ELSE IF (ABS(eta-1._dp) < 1E-10_dp) THEN
    LindG = 2._dp - lambda - 3._dp * mu + 20._dp * (eta-1._dp)

  ! Taylor expansion for high eta
  ELSE IF (eta > 3.65_dp) THEN ! we determined empircally that 3.65 was a 
                               ! good crossover point to the taylor expansion
    eta2 = eta**2
    invEta2 = 1._dp/eta2
    LindG = 3._dp*(1._dp-mu)*eta2 &
            -lambda-0.6_dp &
            + invEta2 * (-0.13714285714285712_dp &
            + invEta2 * (-6.39999999999999875E-2_dp &
            + invEta2 * (-3.77825602968460128E-2_dp &
            + invEta2 * (-2.51824061652633074E-2_dp &
            + invEta2 * (-1.80879839616166146E-2_dp &
            + invEta2 * (-1.36715733124818332E-2_dp &
            + invEta2 * (-1.07236045520990083E-2_dp &
            + invEta2 * (-8.65192783339199453E-3_dp &
            + invEta2 * (-7.1372762502456763E-3_dp &
            + invEta2 * (-5.9945117538835746E-3_dp &
            + invEta2 * (-5.10997527675418131E-3_dp &
            + invEta2 * (-4.41060829979912465E-3_dp &
            + invEta2 * (-3.84763737842981233E-3_dp &
            + invEta2 * (-3.38745061493813488E-3_dp &
            + invEta2 * (-3.00624946457977689E-3_dp)))))))))))))))

  ELSE
    LindG = 1._dp / (0.5_dp + 0.25_dp * (1._dp-eta**2) * LOG((1._dp + eta)/&
             abs(1._dp-eta))/eta) - 3._dp * mu * eta**2 - lambda
  END IF

END FUNCTION LindG


!=========================================================================


!=========================================================================
!
!             Collection of Kernels for nonlocal Ts
!
! MP, Sept. 2015
!
!=========================================================================

SUBROUTINE SimpleKernel (w,MaxPoints,rho_star,ThomasFermiFraction,cell)
  use kinds,                    only : dp
  use constants,                only : pi
  use fft_types,                only : fft_dlay_descriptor
  use fde_types,                only : simulation_cell
  USE fft_interfaces,           ONLY : fwfft, invfft
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum, mp_max
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
  use io_global,                only : ionode, stdout

  type(simulation_cell)            :: cell
  real(dp), intent(inout) :: w(cell%ngm)
  real(dp), intent(inout) :: rho_star,ThomasFermiFraction

  integer, intent(in)     :: MaxPoints

  real(dp)              :: t_tmp, deltat, kf_tmp, kf_star
  real(dp)              :: eta(cell%ngm)
  integer               :: ipnt, i 
 
  integer, pointer :: ngm
  real(dp), pointer :: tpiba2, gg(:)

  ngm => cell%ngm
  tpiba2 => cell%tpiba2
  gg => cell%gg

  kf_star  = (3.d0*pi*pi*rho_star)**(1.d0/3.d0)

  ! Calculate the Kernel in reciprocal space v
  w(:)   = 0.0d0
  deltat = 1.d0/dble(MaxPoints)
  t_tmp  = 1.d0/dble(MaxPoints)

  hypercorrelation_: do ipnt=1,MaxPoints
      kf_tmp = kf_star*(t_tmp**(1.d0/3.d0))
      eta(:) = sqrt(gg(:)*tpiba2)/(2.d0*kf_tmp)
       do i=1,ngm
         w(i) = w(i) + LindG(eta(i),-0.6_dp,1._dp) * (pi*pi/kf_tmp)
       enddo
       t_tmp = t_tmp + deltat
  enddo hypercorrelation_

  w(1:ngm) = w(1:ngm)/dble(MaxPoints)

! re-add the TF component
! \int_0^1 (rho0 t)^1/3 dt = 3/2 rho0^1/3 ==> 3/5*3/2=9/10 this term must be
! removed, otherwise the CZQL is not satisfied
  w(1:ngm) = w(1:ngm) - (9._dp / 10._dp) * (pi*pi/kf_star) 

! regular GP here
  w(1:ngm) = w(1:ngm)*gg(1:ngm)*tpiba2/4.d0/pi

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

END SUBROUTINE SimpleKernel 


!=========================================================================

SUBROUTINE GPLDAKernel ( w, rho, kfs, nkf, Maxpoints,ThomasFermiFraction,cell)
  use kinds,                    only : dp
  use constants,                only : pi
  use fft_types,                only : fft_dlay_descriptor
  use fde_types,                only : simulation_cell
!  use fft_base,                 only : dfftp
  USE fft_interfaces,           ONLY : fwfft, invfft
  use lsda_mod,                 only : nspin
  use scf,                      only : scf_type
  use mp,                       only : mp_sum, mp_max
  use mp_global,                only : intra_pool_comm, intra_bgrp_comm
!  use cell_base,                only : at, alat, tpiba, tpiba2, omega
!  use gvect,                    only : ngm, nl, g, gg, gstart
  use io_global,                only : ionode, stdout

  type(scf_type), intent(in) :: rho
  integer, intent(in)     :: MaxPoints, nkf
  type(simulation_cell)            :: cell
  real(dp), intent(inout) :: w(cell%ngm,nkf)
  real(dp), intent(inout) :: kfs(nkf)
  real(dp), intent(inout) :: ThomasFermiFraction


  real(dp) :: charge, rho_max, kf_max, deltat, kf_tmp
  real(dp) :: eta(cell%ngm), kf_star, t_tmp
  integer  :: ir, ir_ofmax, ikf, itmp, i, ipnt

  type(fft_dlay_descriptor), pointer :: dfftp
  real(dp), pointer :: at(:,:), alat, tpiba, tpiba2, omega
  integer, pointer :: ngm, nl(:), gstart
  real(dp), pointer :: g(:,:), gg(:)

  at => cell%at
  alat => cell%alat
  tpiba => cell%tpiba
  tpiba2 => cell%tpiba2
  omega => cell%omega
  ngm => cell%ngm
  nl => cell%nl
  gstart => cell%gstart
  g => cell%g
  gg => cell%gg

  if ( nkf == 1 ) then
   call errore ('GPLDAKernel','nkf=1 not allowed',1)
  else
  !
  ! find the maximum value of the density
  rho_max = 0.d0
    do ir = 1 , dfftp%nnr
      if (rho%of_r(ir,1) > rho_max ) then
         rho_max = rho%of_r(ir,1)
         ir_ofmax = ir
      endif
      if ( nspin == 2 ) rho_max = max( rho_max , rho%of_r(ir,2) )
    enddo
    call mp_max(rho_max, intra_bgrp_comm )
  endif

  kf_max  = (3.d0*pi*pi*rho_max)**(1.d0/3.d0)


  ! =======================================================
  !
  ! To get rid of the dependency on the rho_star parameter, 
  ! and in order to obtain a position dependent kf(r), 
  ! we calculate the potential several times, according to 
  ! a discre number of values from 
  !
  !  kf = 0     to    kf = (3.d0*pi*pi*rho_max)**(1.d0/3.d0)
  !
  ! =======================================================

  kf_loop : do ikf = 1 , nkf  ! will probably have to include ikf = 0 later on (large q limit?)

    kf_star = kf_max * dble(ikf) / dble(nkf)
    kfs(ikf) = kf_star

    ! Calculate the Kernel in reciprocal space v
    w(:,ikf) = 0.0d0

    deltat = 1.d0/dble(MaxPoints)
    t_tmp  = 1.d0/dble(MaxPoints)


    hypercorrelation_: do ipnt=1,MaxPoints
        kf_tmp = kf_star*(t_tmp**(1.d0/3.d0))
        eta(1:ngm) = sqrt(gg(1:ngm)*tpiba2)/(2.d0*kf_tmp)
         do i=1,ngm
           w(i,ikf) = w(i,ikf) + &
               LindG(eta(i),ThomasFermiFraction,1.0_dp)*(pi*pi/kf_tmp)
         enddo
         t_tmp = t_tmp + deltat
    enddo hypercorrelation_


    ! regular GP here
    w(1:ngm,ikf) = w(1:ngm,ikf)*gg(i)*tpiba2/4.d0/pi/dble(MaxPoints)

  end do kf_loop

END SUBROUTINE GPLDAKernel

!=========================================================================

SUBROUTINE AddTfPartToKernel (w,deltat,k_f,tfp)
  use kinds,                    only : dp
  use constants,                only : pi
  use gvect,                    only : ngm

!
! Here obtain the first part of the integrant between t=0 and t=deltat
!
!
   real(dp), intent(inout) :: w(ngm)
   real(dp), intent(inout) :: deltat, k_f,tfp
   real(dp) :: Coefficient, dTf


  dTf = deltat**(2.0d0/3.0d0)

  Coefficient = (3.0d0/5.0d0 - tfp)

  w=w+Coefficient*dTf*pi*pi/k_f

END SUBROUTINE AddTfPartToKernel


!=========================================================================
SUBROUTINE CreateRhoStar ( rho, rho_fde, rho_star, cell )
! in the making....

  USE kinds,                    only : dp
  use fft_types,                only : fft_dlay_descriptor
  USE fde_types,                only : simulation_cell
  USE fft_interfaces,           only : fwfft, invfft
  USE lsda_mod,                 only : nspin
  USE scf,                      only : scf_type
  USE mp,                       only : mp_sum, mp_max
  USE mp_global,                only : intra_pool_comm, intra_bgrp_comm
  USE io_global,                only : ionode, stdout


! input/output
  type(simulation_cell), intent(in)  :: cell
  type(scf_type), intent(in)         :: rho, rho_fde
  real(dp),      intent(out)         :: rho_star

! local vars
  type(fft_dlay_descriptor), pointer :: dfftp
  real(dp)                           :: rho_tmp, charge, charge_fde
  real(dp), allocatable              :: w(:)


  allocate(w(dfftp%nnr))
     
  ! calculate charge of system
   charge = 0.d0
   if ( cell%gstart == 2 ) then
      charge = cell%omega*real( rho%of_g(1,1) )
      if ( nspin == 2 ) charge = charge + cell%omega*real(rho%of_g(1,2))
   end if
  ! calculate charge of fragment
   charge = 0.d0
   if ( cell%gstart == 2 ) then
      charge = cell%omega*real( rho%of_g(1,1) )
      if ( nspin == 2 ) charge = charge + cell%omega*real(rho%of_g(1,2))
   end if

!  w = 

  deallocate(w)

END SUBROUTINE CreateRhoStar
!=========================================================================


end module kernel_routines







