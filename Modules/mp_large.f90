!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_large
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_bcast, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Group for the large cell parallelization, it's basically the same as intra_bgrp_comm for now. We use the large cell to calculate non-local potentials functional of the total density.
  !
  INTEGER :: nlrp       = 1  ! number of band groups
  INTEGER :: nproc_lgrp  = 1  ! number of processors within a band group
  INTEGER :: me_lgrp     = 0  ! index of the processor within a band group
  INTEGER :: root_lgrp   = 0  ! index of the root processor within a band group
  INTEGER :: my_lgrp_id  = 0  ! index of my band group
  INTEGER :: intra_lgrp_comm  = 0  ! intra band group communicator  
  !
  ! ... "task" groups (for band parallelization of FFT)
  !
!  INTEGER :: ntask_groups = 1  ! number of proc. in an orbital "task group"
  !
  ! ... The following variables not set during initialization but later
  !
!  INTEGER :: ibnd_start = 0 ! starting band index
!  INTEGER :: ibnd_end = 0   ! ending band index
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_large( parent_comm )
    !---------------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into nband_ pools
    ! ... Requires: nband_, read from command line
    ! ...           parent_comm, typically processors of a k-point pool
    ! ...           (intra_pool_comm)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: parent_comm
    !
    INTEGER :: parent_nproc = 1, parent_mype = 0
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    !
    ! ... Set number of processors per large group
    !
    nproc_lgrp = parent_nproc 
    !
    ! ... set index of processor within the image ( 0 : nproc_image - 1 )
    !
    me_lgrp    = parent_mype
    !
    my_lgrp_id = 0
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the intra_bgrp_comm communicator is created
    !
    CALL mp_comm_split( parent_comm, my_lgrp_id, parent_mype, intra_lgrp_comm )
    !
    CALL mp_barrier( parent_comm )
    !
#endif
    RETURN
    !
  END SUBROUTINE mp_start_large
  !
  !
END MODULE mp_large
