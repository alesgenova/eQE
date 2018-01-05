!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_images
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_bcast, mp_size, mp_rank, mp_comm_split
  USE io_global, ONLY : ionode, ionode_id
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Image groups (processors within an image). Images are used for
  ! ... coarse-grid parallelization of semi-independent calculations,
  ! ... e.g. points along the reaction path (NEB) or phonon irreps 
  !
  INTEGER :: nimage = 1 ! number of images
  INTEGER :: nproc_image=1 ! number of processors within an image
  INTEGER :: me_image  = 0 ! index of the processor within an image
  INTEGER :: root_image= 0 ! index of the root processor within an image
  INTEGER :: my_image_id=0 ! index of my image
  INTEGER :: inter_image_comm = 0  ! inter image communicator
  INTEGER :: intra_image_comm = 0  ! intra image communicator
  INTEGER :: inter_fragment_comm = 0 ! inter fragment communicator (FDE)
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_start_images ( nimage_, parent_comm )
    !-----------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into "images". 
    ! ... Requires: nimage_, read from command line
    ! ...           parent_comm, typically world_comm = group of all processors
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nimage_, parent_comm
    !
#if defined (__MPI)
    INTEGER :: parent_nproc, parent_mype
    !
    ! ... nothing needed to be done in serial calculation
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... nimage_ must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    nimage = nimage_
    !
    IF ( nimage < 1 .OR. nimage > parent_nproc ) &
       CALL errore( 'mp_start_images', 'invalid number of images, out of range', 1 )
    IF ( MOD( parent_nproc, nimage ) /= 0 ) &
       CALL errore( 'mp_start_images', 'n. of images must be divisor of nproc', 1 )
    !
    ! ... set number of cpus per image
    !
    nproc_image = parent_nproc / nimage
    !
    ! ... set index of image for this processor   ( 0 : nimage - 1 )
    !
    my_image_id = parent_mype / nproc_image
    !
    ! ... set index of processor within the image ( 0 : nproc_image - 1 )
    !
    me_image    = MOD( parent_mype, nproc_image )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the intra_image_comm communicator is created
    !
    CALL mp_comm_split ( parent_comm, my_image_id, parent_mype, &
                          intra_image_comm )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the inter_image_comm communicator is created
    !
    CALL mp_comm_split ( parent_comm, me_image, parent_mype, &
                         inter_image_comm )
    !
    ! ... set processor that performs I/O
    !
    ionode = ( me_image == root_image )
    ionode_id = root_image
    !
#endif
    RETURN
    !
  END SUBROUTINE mp_start_images
    !
  SUBROUTINE mp_start_fragments ( nimage_, parent_comm )
    !-----------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into "images". 
    ! ... Requires: nimage_, read from command line
    ! ...           parent_comm, typically world_comm = group of all processors
    !
    use mp,          only: mp_comm_group, mp_group_create, mp_group_free, mp_comm_create
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nimage_, parent_comm
    !
    INTEGER :: parent_nproc, parent_mype, ioerr, i
    INTEGER :: images_procs(0:nimage_ - 1), &
                            images_procs_cum(0:nimage_ - 1), &
                            inter_group_list(0:nimage_ - 1)
    LOGICAL :: exst
    !
    INTEGER :: parent_group, inter_group
    !
    ! ... nothing needed to be done in serial calculation
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... nimage_ must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    nimage = nimage_ 
    !
    IF ( nimage < 1 .OR. nimage > parent_nproc ) &
       CALL errore( 'mp_start_fragments', 'invalid number of fragments, out of range', 1 )
    !
    inquire(file='fragments_procs.in', exist=exst)
    !
    call mp_bcast( exst, 0, parent_comm )
    !
    ! If fragments_procs.in does not exist, the processors are split evenly between the fragments
    !
    if (.not. exst) then
    !
       IF ( MOD( parent_nproc, nimage ) /= 0 ) &
          CALL errore( 'mp_start_fragments', 'n. of images must be divisor of nproc', 1 )
       !
       ! ... set number of cpus per image
       !
       nproc_image = parent_nproc / nimage
       !
       ! ... set index of image for this processor   ( 0 : nimage - 1 )
       !
       my_image_id = parent_mype / nproc_image
       !
       ! ... set index of processor within the image ( 0 : nproc_image - 1 )
       !
       me_image    = MOD( parent_mype, nproc_image )
       !
       inter_group_list(0) = root_image
       !
       do i = 1, nimage - 1
          inter_group_list(i) = inter_group_list(i-1) + nproc_image
       enddo       
       !
    !
    !  If fragments_procs.in exists, the processors are split according to it
    !
    else
       if ( parent_mype == 0 ) then
       images_procs(:) = 1
       images_procs_cum(:) = 0
       open(unit=87, file='fragments_procs.in', status='OLD', action='READ', iostat=ioerr)
       do i=0, nimage - 1
          read(87, *)  images_procs(i)
          if ( i > 0) then
             images_procs_cum(i) = images_procs(i) + images_procs_cum(i-1)
          else
             images_procs_cum(i) = images_procs(i)
          end if
          inter_group_list(i) = images_procs_cum(i) - images_procs(i)
       end do

       close(unit=87, iostat=ioerr)
       endif

       call mp_bcast( images_procs, 0, parent_comm  )
       call mp_bcast( images_procs_cum, 0, parent_comm )
       call mp_bcast( inter_group_list, 0, parent_comm )

       IF ( images_procs_cum(nimage-1) /= parent_nproc ) &
           CALL errore( 'mp_start_fragments', &
           'Processors specified in fragments_procs.dat does not match the total number of processors', 1 )

       if ( parent_mype < images_procs(0) ) then
          my_image_id = 0
          me_image = parent_mype
       else
          do i=1, nimage -1
             if ( parent_mype < images_procs_cum(i) .and. parent_mype >= images_procs_cum(i-1) ) then
                my_image_id = i
                me_image = parent_mype - images_procs_cum(i-1)
             end if
          end do
       end if

       nproc_image = images_procs(my_image_id)

    end if
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the intra_image_comm communicator is created
    !
    CALL mp_comm_split ( parent_comm, my_image_id, parent_mype, &
                          intra_image_comm )
    !
    call mp_comm_group( parent_comm, parent_group )

    call mp_group_create( inter_group_list, nimage, parent_group, inter_group )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the inter_fragment_comm communicator is created
    !
    call mp_comm_create( parent_comm, inter_group, inter_fragment_comm )
    !
    call mp_group_free( parent_group )
    !
    call mp_group_free( inter_group )
    !
    ! ... set processor that performs I/O
    !
    ionode = ( me_image == root_image )
    ionode_id = root_image
    !
#endif
    RETURN
    !
  END SUBROUTINE mp_start_fragments
    !
  SUBROUTINE mp_init_image ( parent_comm )
    !
    use mp,          only: mp_comm_group, mp_group_create, mp_group_free, mp_comm_create
    ! ... There is just one image: set it to the same as parent_comm (world)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: parent_comm
    INTEGER :: inter_group_list(0:1), inter_group, parent_group
    !
    intra_image_comm = parent_comm 
    nproc_image = mp_size( parent_comm )
    me_image    = mp_rank( parent_comm )
    !
    ! ... no need to set inter_image_comm,  my_image_id, root_image
    ! ... set processor that performs I/O
    !
    ionode = ( me_image == root_image )
    ionode_id = root_image
! fake inter_fragment_comm
    inter_group_list(0)=root_image
    call mp_comm_group(parent_comm,parent_group)
    call mp_group_create(inter_group_list, 1, parent_group, inter_group)
    call mp_comm_create(parent_comm,inter_group,inter_fragment_comm)
    call mp_group_free(inter_group)
    call mp_group_free(parent_group)
    !
  END SUBROUTINE mp_init_image
  !
END MODULE mp_images
