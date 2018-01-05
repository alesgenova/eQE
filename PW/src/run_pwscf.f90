!
! Copyright (C) 2013-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE run_pwscf ( exit_status ) 
  !----------------------------------------------------------------------------
  !
  !! author: Paolo Giannozzi
  !! license: GNU 
  !! summary: Run an instance of the Plane Wave Self-Consistent Field code
  !!
  !! Run an instance of the Plane Wave Self-Consistent Field code 
  !! MPI initialization and input data reading is performed in the 
  !! calling code - returns in exit_status the exit code for pw.x, 
  !! returned in the shell. Values are:
  !! * 0: completed successfully
  !! * 1: an error has occurred (value returned by the errore() routine)
  !! * 2-127: convergence error
  !!   * 2: scf convergence error
  !!   * 3: ion convergence error
  !! * 128-255: code exited due to specific trigger
  !!   * 255: exit due to user request, or signal trapped,
  !!          or time > max_seconds
  !!     (note: in the future, check_stop_now could also return a value
  !!     to specify the reason of exiting, and the value could be used
  !!     to return a different value for different reasons)
  !! @Note
  !! 10/01/17 Samuel Ponce: Add Ford documentation
  !! @endnote
  !!
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, ethr, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE command_line_options, ONLY : command_line
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE extrapolation,    ONLY : update_file, update_pot
  USE scf,              ONLY : rho
  USE lsda_mod,         ONLY : nspin
  USE fft_base,         ONLY : dfftp
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  USE qexsd_module,     ONLY:   qexsd_set_status
  USE mp,               ONLY : mp_barrier
  USE mp_world,         ONLY : world_comm
  USE fde,              ONLY : do_fde, fde_print_density, fde_print_embedpot,& 
                               fde_print_density_frag, fde_print_allpot, &
                               rho_fde_large, rho_core_fde_large, currfrag, &
                               fde_print_density_frag_large, &
                               fde_printdensity_vec, &
                               fde_plotemb_vec, fde_print_electro, &
                               nonlocalkernel, nonlocalkernel_fde
  use fde_routines
  use kernel_routines,  only : deletekernel
  USE noncollin_module, ONLY : report, i_cons, noncolin
  USE spin_orb,         ONLY : domag
  use lsda_mod ,        only : nspin
  use dynamics_module,  only : trajectory
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_status
  !! Gives the exit status at the end
  LOGICAL, external :: matches
  !! checks if first string is contained in the second
  INTEGER :: idone 
  ! counter of electronic + ionic steps done in this run
  !
  exit_status = 0
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... needs to come before iosys() so some input flags can be
  !     overridden without needing to write PWscf specific code.
  ! 
  CALL qmmm_initialization()
  !
  ! ... convert to internal variables
  !
  CALL iosys()
  !
  ! ... If executable names is "dist.x", compute atomic distances, angles,
  ! ... nearest neighbors, write them to file "dist.out", exit
  !
  IF ( matches('dist.x',command_line) ) THEN
     IF (ionode) CALL run_dist ( exit_status )
     RETURN
  END IF
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  ! call to void routine for user defined / plugin patches initializations
  !
  CALL plugin_initialization()
  !
  CALL check_stop_init()
  !
  CALL setup ()
  !
  CALL qmmm_update_positions()
  !
  CALL init_run()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( check_stop_now() ) THEN
     CALL qexsd_set_status(255)
     CALL punch( 'config' )
     exit_status = 255
     RETURN
  ENDIF
  !
  main_loop: DO idone = 1, nstep
     !
     ! ... electronic self-consistency or band structure calculation
     !
     IF ( .NOT. lscf) THEN
        CALL non_scf ()
     ELSE
        CALL electrons()
     END IF
     ! FDE
     if (do_fde) call mp_barrier(world_comm)
     !
     ! ... code stopped by user or not converged
     !
     IF ( check_stop_now() .OR. .NOT. conv_elec ) THEN
        IF ( check_stop_now() ) exit_status = 255
        IF ( .NOT. conv_elec )  exit_status =  2
        CALL qexsd_set_status(exit_status)
        ! workaround for the case of a single k-point
        twfcollect = .FALSE.
        ! Uncomment if debugging and want to know how the density looks like even when
        ! calculation fails
        ! if (do_fde) then
        ! if (fde_print_density.or.fde_print_density_frag) call fde_plot_density
        ! if (fde_print_embedpot) call fde_plot_embedpot
        ! endif
        CALL punch( 'config' )
        RETURN
     ENDIF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' ); !write(*,*)' start ions' ; FLUSH(6)
     conv_ions = .TRUE.
     !
     ! ... recover from a previous run, if appropriate
     !
     !IF ( restart .AND. lscf ) CALL restart_in_ions()
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) THEN
        CALL pw2casino( istep )
     ELSE
        CALL pw2casino( 0 )
     END IF
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( lstres ) CALL stress ( sigma )
     !
     ! ... send out forces to MM code in QM/MM run
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        if (fix_volume) CALL impose_deviatoric_stress(sigma)
        if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
        !
        ! ... save data needed for potential and wavefunction extrapolation
        !
        CALL update_file ( )
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions ( idone )
        !
        ! write the coordinates of the whole system to file
        if (do_fde) call trajectory()
        !
        ! ... since the ions moved, we need to call make_pointlists again
        !
        IF ( ( (report.ne.0).or.(i_cons.ne.0) ) .and. (noncolin.and.domag) &
                      .or. (i_cons.eq.1) .or. nspin==2 ) THEN
        !
        ! In order to print out local quantities, integrated around the atoms,
        ! we need the following variables
        !
           CALL make_pointlists ( )
        ENDIF
        !
        ! ... then we save restart information for the new configuration
        !
        IF ( idone <= nstep .AND. .NOT. conv_ions ) THEN 
            CALL qexsd_set_status(255)
            CALL punch( 'config' )
        END IF
        !
     END IF
     !
     CALL stop_clock( 'ions' ); !write(*,*)' stop ions' ; FLUSH(6)
     !
     CALL qmmm_update_forces( force, rho%of_r, nspin, dfftp)
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( lmd .OR. lbfgs ) CALL add_qexsd_step(idone)
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... receive new positions from MM code in QM/MM run
     !
     CALL qmmm_update_positions()
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        ! ... update the wavefunctions, charge density, potential
        ! ... update_pot initializes structure factor array as well
        !
        CALL update_pot()
        !
        ! ... re-initialize atomic position-dependent quantities
        !
        CALL hinit1()
        !
     END IF
     ! ... Reset convergence threshold of iterative diagonalization for
     ! ... the first scf iteration of each ionic step (after the first)
     !
     ethr = 1.0D-6
     !
  END DO main_loop
  !
  ! ... save final data file
  !
  CALL qexsd_set_status(exit_status)
  CALL punch('all')
  if (do_fde) then
    call mp_barrier(world_comm)
    !
    !call test_sum_band_selective()
    !
    if ( sum(fde_printdensity_vec(:)) > 0 .or. &
         fde_print_density .or. &
         fde_print_density_frag .or. &
         fde_print_density_frag_large) call fde_plot_density
    if ( sum(fde_plotemb_vec(:)) > 0 ) call fde_plot_embedpot
    !if (fde_print_embedpot) call fde_plot_embedpot
    if (fde_print_allpot) call fde_plot_allpot
    if (fde_print_electro) call fde_plot_electrostatics

    call DeleteKernel(NonlocalKernel)
    call DeleteKernel(NonlocalKernel_fde)

  endif
  !
  CALL qmmm_shutdown()
  !
  IF ( .NOT. conv_ions )  exit_status =  3
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE run_pwscf
