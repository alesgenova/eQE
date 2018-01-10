!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM fdepw
  !----------------------------------------------------------------------------
  !
  ! ... FDE pw.x parallel launcher. Usage (for mpirun):
  ! ...    mpirun -np Np fdepw.x -ni Ni -in filename [other options]
  ! ... or whatever is appropriate for your parallel environment
  ! ... Starts Ni pw.x instances each running on Np/Ni processors
  ! ... Each pw.x instance:
  ! ... * reads input data from filename_N.in or pw_N.in if no input
  ! ...   file is specified
  ! ... * saves temporary and final data to "outdir"_N/ directory
  ! ...   (or to tmp_N/ if outdir='./')
  ! ... * writes output to pw_N.out in the current directory if no input
  ! ...   file is specified via the -i option; to "input_file"_N.out
  ! ...   if command-line options -i "input_file" is specified
  !
  USE input_parameters,  ONLY : outdir
  USE environment,       ONLY : environment_start, environment_end
  USE io_global,         ONLY : ionode, ionode_id, stdout
  USE mp_global,         ONLY : mp_startup,  nimage
  USE mp_images,         ONLY : my_image_id
  USE read_input,        ONLY : read_input_file
  USE fde,               ONLY : do_fde, nfragments, currfrag
  USE command_line_options, ONLY: input_file_

  implicit none
  integer :: i, exit_status
  logical :: opnd
  character(len=256) :: filin, filout
  character(len=7) :: image_label
  character(len=6), external :: int_to_char


  ! initialize code
  call mp_startup(start_fragments=.true.)

  !call environment_start('PWSCF')
  call environment_start('eQE')
#ifndef __MPI
  call errore('fdepw', 'FDE cannot be run in serial', 1)
  return
#endif
 
  ! open input files
  image_label = '_' // int_to_char(my_image_id)
  if ( trim (input_file_) == ' ') then
     filin = 'pw' // trim(image_label)  // '.in'
  else
     filin = trim(input_file_) // trim(image_label) // '.in'
  end if

  ! open image-specific output files
  if (ionode) then
     inquire(unit=stdout, opened=opnd)
     if (opnd) close(unit=stdout)
     filout = 'pw' // trim(image_label) // '.out'
     if (trim(input_file_) /= ' ') &
        filout = trim(input_file_) // trim(image_label) // '.out'
     open(unit=stdout, file=trim(filout), status='unknown')
  end if

  ! read input file and check for (nimages == nfragments)
  call start_clock('PWSCF')
  call read_input_file(prog='PW', input_file_=filin)
  
  ! ... set image-specific value for "outdir", starting from input value
  ! ... (read in read_input_file)
  do i=len_trim(outdir),1,-1
     if ( outdir(i:i) /= '/' .and. outdir(i:i) /= '.' ) exit
  end do
  ! ... i = position of last character different from '/' and '.'
  if ( i == 0 ) then
     outdir = 'tmp' // trim(image_label) // '/'
  else
     outdir = outdir(1:i) // trim(image_label) // '/'
  end if

  ! ... perform actual calculation
  do_fde = .false.
  nfragments = nimage
  if (nfragments > 0) do_fde = .true.
  currfrag = my_image_id + 1
  call run_pwscf ( exit_status )
  !
  call stop_run( exit_status )
  call do_stop( exit_status )
  !
  stop
  !
END PROGRAM fdepw
