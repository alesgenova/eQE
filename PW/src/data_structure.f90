!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE data_structure( gamma_only )
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays
  ! (both the smooth and the dense grid)
  ! In the parallel case, it distributes columns to processes, too
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_max
  USE mp_bands,   ONLY : nproc_bgrp, intra_bgrp_comm, nyfft, ntask_groups
  USE mp_large,   ONLY : me_lgrp, nproc_lgrp, root_lgrp, intra_lgrp_comm
  USE mp_pools,   ONLY : inter_pool_comm
  USE fft_base,   ONLY : dfftp, dffts, fft_base_info, smap
  USE fft_base,   ONLY : dfftl, smap_large, fft_base_info_large
  USE fft_types,  ONLY : fft_type_init
  USE cell_base,  ONLY : at, bg, tpiba
  USE large_cell_base,  ONLY : atl => at, bgl => bg, tpibal => tpiba
  USE klist,      ONLY : xk, nks
  USE gvect,      ONLY : gcutm, gvect_init
  USE gvecs,      ONLY : gcutms, gvecs_init, doublegrid
  USE gvecl,      ONLY : gcutml => gcutm, gvecl_init
  use fde,        only : do_fde, linterlock
  use io_global,  only : stdout!, flush_unit
  USE gvecw,      ONLY : gcutw, gkcut
  USE realus,     ONLY : real_space
  USE io_global,  ONLY : stdout, ionode
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in) :: gamma_only
  INTEGER :: ik, ngm_, ngs_
  LOGICAL :: lpara
  !
  lpara =  ( nproc_bgrp > 1 )
  !
  ! ... calculate gkcut = max |k+G|^2, in (2pi/a)^2 units
  !
  IF (nks == 0) THEN
     !
     ! if k-points are automatically generated (which happens later)
     ! use max(bg)/2 as an estimate of the largest k-point
     !
     gkcut = 0.5d0 * max ( &
        sqrt (sum(bg (1:3, 1)**2) ), &
        sqrt (sum(bg (1:3, 2)**2) ), &
        sqrt (sum(bg (1:3, 3)**2) ) )
  ELSE
     gkcut = 0.0d0
     DO ik = 1, nks
        gkcut = max (gkcut, sqrt ( sum(xk (1:3, ik)**2) ) )
     ENDDO
  ENDIF
  gkcut = (sqrt (gcutw) + gkcut)**2
  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max (gkcut, inter_pool_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  ! task group are disabled if real_space calculation of calbec is used
  dffts%has_task_groups = (ntask_groups >1) .and. .not. real_space
  CALL fft_type_init( dffts, smap, "wave", gamma_only, lpara, intra_bgrp_comm, at, bg, gkcut, gcutms/gkcut, nyfft=nyfft )
  CALL fft_type_init( dfftp, smap, "rho" , gamma_only, lpara, intra_bgrp_comm, at, bg, gcutm , 4.d0, nyfft=nyfft )
  ! define the clock labels ( this enables the corresponding fft too ! )
  dffts%rho_clock_label='ffts' ; dffts%wave_clock_label='fftw'
  dfftp%rho_clock_label='fft'
  if (.not.doublegrid) dfftp%grid_id = dffts%grid_id  ! this makes so that interpolation is just a copy.

  CALL fft_base_info( ionode, stdout )
  ngs_ = dffts%ngl( dffts%mype + 1 )
  ngm_ = dfftp%ngl( dfftp%mype + 1 )
  IF( gamma_only ) THEN
     ngs_ = (ngs_ + 1)/2
     ngm_ = (ngm_ + 1)/2
  END IF
  !
  !     on output, ngm_ and ngs_ contain the local number of G-vectors
  !     for the two grids. Initialize local and global number of G-vectors
  !
  call gvect_init ( ngm_ , intra_bgrp_comm )
  call gvecs_init ( ngs_ , intra_bgrp_comm )
  !
  if (do_fde) then
    ! gkcut = 0.d0
    ! CALL pstickset_large( gamma_only, bgl, gcutml, gkcut, gcutml, &
    !              dfftl, dfftl, ngwl_ , ngml_ , ngsl_ , me_lgrp, &
    !              root_lgrp, nproc_lgrp, intra_lgrp_comm, ntask_groups )
    ! call gvecl_init( ngml_ , intra_lgrp_comm )
    lpara =  ( nproc_lgrp > 1 )
    CALL fft_type_init( dfftl, smap_large, "rho" , gamma_only, lpara, intra_lgrp_comm, atl, bgl, gcutml , 4.d0, nyfft=1 )
    dfftl%rho_clock_label='fftl'
    CALL fft_base_info_large( ionode, stdout )
    ngm_ = dfftl%ngl( dfftl%mype + 1 )
    IF( gamma_only ) THEN
       ngm_ = (ngm_ + 1)/2
    END IF
    call gvecl_init ( ngm_ , intra_lgrp_comm )
  endif

END SUBROUTINE data_structure

