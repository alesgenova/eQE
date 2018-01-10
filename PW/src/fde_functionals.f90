!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This files contains routines to evaluate the non-additive kinetice
! energy semilocal functionals for the Frozen Density Embedding method [1,2].
!
! Coded by D. Ceresoli and A. Genova, during the April 2013 PRG hackathon
! in Rutgers-Newark (http://www.michelepavanello.com/pdf/hackathon.pdf)
!
! [1] T. A. Wesolowksi and A. Warshel, J. Phys. Chem. 97, 8050 (1993) 
! [2] see: $QEDIR/FDE-Notes/fde-equations.tex

Module fde_functionals
implicit none
contains
! Thomas-Fermi
SUBROUTINE fde_kin_tf (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  Fs = 1.d0
  dF_ds = 0.d0
END SUBROUTINE fde_kin_tf


! von Weizsacker
SUBROUTINE fde_kin_vw (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  Fs = 5.d0/3.d0 * s*s
  dF_ds = 10.d0/3.d0 * s
END SUBROUTINE fde_kin_vw

! Damped von Weizsacker
SUBROUTINE fde_kin_dvw (nr, s, Fs, dF_ds)
  use kinds, only : dp
  use io_global, only : stdout
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter  :: B = 5.d0
  real(dp), parameter  :: C = 3.d0
  write(stdout,*) 'Damped VW!!'
  Fs = (5.d0*s**2*(1.d0 - Tanh(B*(-C + s))))/6.d0 
  dF_ds = (-5.d0*B*s**2*(1.d0/Cosh(B*(-C + s)))**2)/6.d0 - (5.d0*s*(1.d0 - Tanh(B*(-C + s))))/3.d0
END SUBROUTINE fde_kin_dvw


! 2nd order DGE
SUBROUTINE fde_kin_dge2 (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  Fs = 1.d0 + 5.d0/27.d0 * s**2
  dF_ds = 10.d0/27.d0 * s
END SUBROUTINE fde_kin_dge2


! Lee, Lee, Parr
SUBROUTINE fde_kin_llp (nr, s, Fs, dF_ds)
  use kinds, only : dp
  use constants, only : pi
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter :: beta = 2.d0*(6.d0*pi*pi)**(1.d0/3.d0)
  real(dp), parameter :: b1 = 0.0044188d0
  real(dp), parameter :: b2 = 0.0253d0
  Fs = 1.d0 + (b1 * (beta*s)**2) / (1.d0 + b2*beta*s*asinh(beta*s))
  dF_ds = (2*beta*beta*b1*2)/(1.d0+beta*b2*s*asinh(beta*s)) - &
     (b1*(beta*s)**2 * (beta*beta*b2*s/sqrt(1.d0+(beta*s)**2) + beta*b2*asinh(beta*s))) / &
     (1.d0 + b2*beta*s*asinh(beta*s))**2
#ifdef __XLF
CONTAINS
FUNCTION asinh ( x )
  use kinds, only : dp
  implicit none
  real(dp), intent(in) :: x(:)
  real(dp), allocatable :: asinh(:)
  allocate( asinh(size(x)) )
  asinh = log( x + sqrt(x**2 + 1) )
  return
END FUNCTION asinh
#endif
END SUBROUTINE fde_kin_llp


! cojoint Perdew, Wang
SUBROUTINE fde_kin_pw86 (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter :: b1 = 1.296d0
  real(dp), parameter :: b2 = 14.d0
  real(dp), parameter :: b3 = 0.2d0
  Fs = (1.d0 + b1*s**2 + b2*s**4 + b3*s**6)**(1.d0/15.d0)
  dF_ds = (2*b1*s + 4*b2*s**3 + 6*b3*s**5) / (15.*(1 + b1*s**2 + b2*s**4 + &
           b3*s**6)**(14.d0/15.d0))
END SUBROUTINE fde_kin_pw86


! Lembarki, Chermette
SUBROUTINE fde_kin_lc94 (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter :: A0 = 76.32d0
  real(dp), parameter :: A1 = 0.093907d0
  real(dp), parameter :: A2 = 0.26608d0
  real(dp), parameter :: A3 = 0.0809615d0
  real(dp), parameter :: A4 = 100.d0
  real(dp), parameter :: A5 = 0.57767d-4
  Fs = (1.d0 + (A2 - A3 * exp(-A4*s**2))*s**2 + A1*s*asinh(A0*s)) / &
       (1 + A5*s**4 + A1*s*asinh(A0*s))
  dF_ds = -(((4.d0*A5*s**3 + (A0*A1*s)/sqrt(1.d0 + A0**2*s**2) + A1*asinh(A0*s))* &
         (1.d0 + (A2 - A3*exp(-A4*s**2))*s**2 + A1*s*asinh(A0*s)))/(1.d0 + A5*s**4 + &
         A1*s*asinh(A0*s))**2) + &
         (2.d0*(A2 - A3*exp(-A4*s**2))*s + 2.d0*A3*A4*s**3*exp(-A4*s**2) + &
         (A0*A1*s)/sqrt(1.d0 + A0**2*s**2) + A1*asinh(A0*s)) / &
         (1.d0 + A5*s**4 + A1*s*asinh(A0*s))
#ifdef __XLF
CONTAINS
FUNCTION asinh ( x )
  use kinds, only : dp
  implicit none
  real(dp), intent(in) :: x(:)
  real(dp), allocatable :: asinh(:)
  allocate( asinh(size(x)) )
  asinh = log( x + sqrt(x**2 + 1) )
  return
END FUNCTION asinh
#endif
END SUBROUTINE fde_kin_lc94


! asymptotic PBE (APBEK)
SUBROUTINE fde_kin_apbek (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter :: a1 = 0.23899d0
  real(dp), parameter :: a2 = 0.804d0
  Fs = 1.d0 + a1*s*s/(1.d0 + (a1/a2)*s*s)
  dF_ds = (2.d0 * a1*a2*a2*s) / (a1*s*s + a2)**2
END SUBROUTINE fde_kin_apbek

! Rev asymptotic PBE (revAPBEK)
SUBROUTINE fde_kin_rapbek (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter :: a1 = 0.23899d0
  real(dp), parameter :: a2 = 1.245d0
  Fs = 1.d0 + a1*s*s/(1.d0 + (a1/a2)*s*s)
  dF_ds = (2.d0 * a1*a2*a2*s) / (a1*s*s + a2)**2
END SUBROUTINE fde_kin_rapbek


! TFP
SUBROUTINE fde_kin_tfp (nr, s, Fs, dF_ds)
  use kinds, only : dp
  implicit none
  integer, intent(in) :: nr
  real(dp), intent(in) :: s(nr)
  real(dp), intent(out) :: Fs(nr), dF_ds(nr)
  real(dp), parameter :: a1 = 5.d0/27.d0
  !real(dp), parameter :: a1 = 0.23899d0
  real(dp), parameter :: a2 = 5.0d0
  !real(dp), parameter :: a2 = 1.245d0
  Fs = 1.d0 + a1*s*s/(1.d0 + (a1/a2)*s*s)
  dF_ds = (2.d0 * a1*a2*a2*s) / (a1*s*s + a2)**2
END SUBROUTINE fde_kin_tfp

END MODULE fde_functionals
