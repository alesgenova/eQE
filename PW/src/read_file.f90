!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__OLDXML)
!----------------------------------------------------------------------------
! TB
! included allocation of the force field of the gate, search for 'TB'
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE read_file()
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for compatibility
  !
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, tmp_dir, wfc_dir
  USE io_global,            ONLY : stdout, ionode
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE uspp,                 ONLY : becsum
  USE scf,                  ONLY : rho
  USE realus,               ONLY : betapointlist, &
                                   init_realspace_vars,real_space
  USE dfunct,               ONLY : newd
  USE ldaU,                 ONLY : lda_plus_u, U_projection
  USE pw_restart,           ONLY : pw_readfile
  USE control_flags,        ONLY : io_level
  USE input_parameters,     ONLY : fde_xc_funct
! for debugging
!  USE mp_images,       ONLY : nimage, my_image_id
!  USE mp_bands,        ONLY : inter_bgrp_comm, nbgrp
!  USE mp_pools,        ONLY : nproc_pool, me_pool
  USE klist,                ONLY : init_igk
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
#if defined (__HDF5)
  USE hdf5_qe
#endif
  !
  IMPLICIT NONE 
  INTEGER :: ierr
  LOGICAL :: exst
  CHARACTER( 256 )  :: dirname
  !
  !
  ierr = 0 
  !
  ! ... Read the contents of the xml data file
  !
  IF ( ionode ) WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:', TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
#if defined(__HDF5)
  CALL initialize_hdf5()
#endif
  !
  CALL read_xml_file ( )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  wfc_dir = tmp_dir
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  ! ... Allocate and compute k+G indices and number of plane waves
  ! ... FIXME: should be read from file, not re-computed
  !
  CALL init_igk ( npwx, ngm, g, gcutw ) 
  !
  CALL pw_readfile( 'wave', ierr )
  !
  ! ... Assorted initialization: pseudopotentials, PAW
  ! ... Not sure which ones (if any) should be done here
  !
  CALL init_us_1()
  !
  IF (lda_plus_u .AND. (U_projection == 'pseudo')) CALL init_q_aeps()
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  !
  IF ( real_space ) THEN
    CALL betapointlist()
    CALL init_realspace_vars()
    IF( ionode ) WRITE(stdout,'(5x,"Real space initialisation completed")')
  ENDIF
  CALL newd()
  !
  CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE read_file
!
SUBROUTINE read_xml_file()
  ! wrapper routine to call the default behavior
  call read_xml_file_internal(.true.)
END SUBROUTINE read_xml_file

SUBROUTINE read_xml_file_nobs()
  ! wrapper routine to load everything except for the band structure
  call read_xml_file_internal(.false.)
END SUBROUTINE read_xml_file_nobs

!----------------------------------------------------------------------------
SUBROUTINE read_xml_file_internal(withbs)
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  ! ... All quantities that are initialized in subroutine "setup" when
  ! ... starting from scratch should be initialized here when restarting
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, nsp, ityp, tau, if_pos, extfor
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE large_cell_base,      ONLY : alatl => alat, omegal => omega, atl => at, bgl => bg, ibravl => ibrav!, tpiba2 => tpiba2
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE symm_base,            ONLY : irt, d1, d2, d3, checkallsym, nsym
  USE extfield,             ONLY : forcefield, tefield, gate, forcegate
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE fft_base,             ONLY : dfftp
  USE fft_base,             ONLY : dfftl, dffts
  USE fft_interfaces,       ONLY : fwfft
  USE fft_types,            ONLY : fft_type_allocate
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : gg, ngm, g, gcutm, mill, ngm_g, ig_l2g, &
                                   eigts1, eigts2, eigts3, gstart
  USE gvecl,                ONLY : gg_l => gg, ngm_l => ngm, g_l => g, gcutm_l => gcutm, &
                                   eigts1_l => eigts1, eigts2_l => eigts2, eigts3_l => eigts3, &
                                   nl_l => nl, gstart_l => gstart, ig_l2g_l=>ig_l2g, &
                                   mill_l => mill
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact
  USE fft_base,             ONLY : dfftp, dffts
  USE gvecs,                ONLY : ngms, gcutms 
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE scf,                  ONLY : vltot
  use scf_large,            only : vltot_large => vltot
  USE wavefunctions_module, ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE pw_restart,           ONLY : pw_readfile, pp_check_file
  USE io_rho_xml,           ONLY : read_scf
  USE read_pseudo_mod,      ONLY : readpp
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : upf
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE ldaU,                 ONLY : lda_plus_u, eth, init_lda_plus_u
  USE control_flags,        ONLY : gamma_only
  USE funct,                ONLY : get_inlc, get_dft_name
  USE funct,                ONLY : dft_is_nonlocc
  USE kernel_table,         ONLY : initialize_kernel_table
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE mp_bands,             ONLY : intra_bgrp_comm, nyfft
  USE fde
  USE fde_routines
  USE input_parameters,     ONLY : fde_xc_funct
! for debugging
  USE mp_images,            ONLY : nimage, my_image_id, intra_image_comm, inter_fragment_comm
!  USE mp_bands,        ONLY : inter_bgrp_comm, nbgrp
!  USE mp_pools,        ONLY : nproc_pool, me_pool
  use mp,                   only : mp_sum, mp_barrier, mp_bcast
!  use mp_world,                 only : world_comm
  use io_global,            only : stdout, ionode, ionode_id
  !
  IMPLICIT NONE

  ! Used to specify whether to read the band structure (files 
  ! K??????/eigenval.xml), so one can skip it if not needed by
  ! the post-processing tool. 
  ! Set to True for the 'default' behavior of reading these files.
  LOGICAL :: withbs

  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr, inlc
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  CHARACTER(LEN=20) dft_name
  !
  !
  ! ... first we get the version of the qexml file
  !     if not already read
  CALL pw_readfile( 'header', ierr )
  CALL errore( 'read_xml_file ', 'unable to determine qexml version', ABS(ierr) )
  !
  ! ... then we check if the file can be used for post-processing
  !
  IF ( .NOT. pp_check_file() ) CALL infomsg( 'read_xml_file', &
               & 'file ' // TRIM( tmp_dir ) // TRIM( prefix ) &
               & // '.save not guaranteed to be safe for post-processing' )
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
  CALL pw_readfile( 'reset', ierr )
  CALL pw_readfile( 'dim',   ierr )
  CALL errore( 'read_xml_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  flush(stdout)
  !
  ! ... allocate space for atomic positions, symmetries, forces
  !
  IF ( nat < 0 ) CALL errore( 'read_xml_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau(    3, nat ) )
  if (do_fde)  ALLOCATE(tau_large(3, nat))
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( extfor(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  IF ( gate ) ALLOCATE( forcegate( 3, nat ) ) ! TB
  !
  ALLOCATE( irt( 48, nat ) )
  !
  CALL set_dimensions()
  CALL fft_type_allocate ( dfftp, at, bg, gcutm, intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft )
  if (do_fde) CALL fft_type_allocate ( dfftl, atl, bgl, gcutml, intra_lgrp_comm, nyfft=1 )
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  if (cell_factor == 0.d0) cell_factor = 1.D0
  !
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  ! ... here we read all the variables defining the system
  !
  IF  ( withbs .EQV. .TRUE. ) THEN  
     CALL pw_readfile( 'nowave', ierr )
  ELSE
     CALL pw_readfile( 'nowavenobs', ierr )
  END IF
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... check on symmetry
  !
  IF (nat > 0) CALL checkallsym( nat, tau, ityp )
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF
  !
  ! ... read pseudopotentials
  !
  CALL pw_readfile( 'pseudo', ierr )

  dft_name = get_dft_name () ! already set, should not be set again
  CALL readpp ( dft_name )
  !
  ! ... read the vdw kernel table if needed
  !
  if (dft_is_nonlocc()) then
  if (do_fde) then
        if (trim(fde_xc_funct) == 'SAME' ) then
           inlc = get_inlc()
        elseif (trim(fde_xc_funct) == 'RVV10' ) then
           inlc = 3
        elseif (trim(fde_xc_funct) == 'REV-VDW-DF2' .or. &
                trim(fde_xc_funct) == 'VDW-DF2-C09' .or. &
                trim(fde_xc_funct) == 'VDW-DF2' ) then
           inlc = 2
        elseif (trim(fde_xc_funct) == 'VDW-DF4'    .or. &
                trim(fde_xc_funct) == 'VDW-DF3'    .or. &
                trim(fde_xc_funct) == 'VDW-DF-C09' .or. &
                trim(fde_xc_funct) == 'VDW-DF' ) then
           inlc = 1
        else
           inlc = 0
        endif
        if (inlc > 0) call initialize_kernel_table(inlc)
  else 
        inlc = get_inlc()
        if (inlc > 0) call initialize_kernel_table(inlc)
  endif

  !if (trim(fde_xc_funct) /= 'SAME' .and. inlc .ne. 0) inlc = 0

  endif
  !
  okpaw = ANY ( upf(1:nsp)%tpawp )
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL pre_init()
  CALL data_structure ( gamma_only )
  CALL allocate_fft()
  if ( do_fde ) call allocate_fft_large()
  ! moved allocate_fde earlier, because set_rhoc needs some "large". cannot be earlier than allocate_fft because allocate_fde needs ngm.
  if (do_fde)  then
       call allocate_fde
       !
       ! fde_cell_offset , tau , tau_large are all read from xml, no need to call local_offset again
       !if ( linterlock) call get_local_offset(fde_cell_offset,fde_cell_shift,frag_cell_split)
  endif
  allocate( fde_si_vec(nfragments) )
  fde_si_vec=0
  if (fde_si) fde_si_vec(currfrag) = 1
  if (ionode) call mp_sum( fde_si_vec, inter_fragment_comm )
  call mp_bcast( fde_si_vec, ionode_id, intra_image_comm )
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
              g, gg, mill, ig_l2g, gstart ) 
  if ( do_fde ) CALL ggen ( dfftl, gamma_only, atl, bgl, gcutm_l, ngm_g_l, ngm_l, &
                            g_l, gg_l, mill_l, ig_l2g_l, gstart_l ) 
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms ) 
  IF (do_comp_esm) THEN
    CALL pw_readfile( 'esm', ierr )
    CALL esm_init()
  END IF
  CALL gshells ( lmovecell ) 
  if ( do_fde ) call gshells_large( lmovecell )
  ! Not sure here is the best place to initialize the datatype, but let's se...
  if (do_fde) call gen_simulation_cells()
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  CALL allocate_locpot()
  if (do_fde ) call allocate_locpot_large()
  !
  IF ( lda_plus_u ) THEN
     CALL init_lda_plus_u ( upf(1:nsp)%psd, noncolin )
  ENDIF
  !
  CALL allocate_wfc()
  !
  ! ... read the charge density
  !
  CALL read_scf( rho, nspin )
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  IF (do_fde) CALL fde_summary
  ! 2D calculations: re-initialize cutoff fact before calculating potentials
  IF(do_cutoff_2D) CALL cutoff_fact()
  !
  CALL init_vloc()
  if (do_fde) then
!     tau_large = tau
     !call get_local_offset(fde_cell_offset,fde_cell_shift,frag_cell_split)
     call init_vloc_large()
  endif
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2, &
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
! this part copies from hinit0, where it's between struc_fact and setlocal
  if (do_fde ) then 
    !strf_fde(:,:) = strf(:,:)
    CALL struc_fact( nat_fde, tau_fde, nsp, ityp_fde, ngml, gl, bgl, &
          dfftl%nr1, dfftl%nr2, dfftl%nr3, strf_fde_large, eigts1l, eigts2l, eigts3l )
    call calc_f2l(f2l, dfftp, dfftl, fde_cell_shift, fde_cell_offset)
    call setlocal_fde_large(vltot_large, strf_fde_large)
  endif
! END: this part copies from hinit0, where it's between struc_fact and setlocal
  CALL setlocal()
  if (do_fde) call copy_pot_l2f(vltot_large, vltot)
  CALL set_rhoc()
  !
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     CALL fwfft ('Rho', psic, dfftp)
     rho%of_g(:,is) = psic(dfftp%nl(:))
     !
  END DO
  if (do_fde) call update_rho_fde(rho, .true.)
  !
  ! ... read info needed for hybrid functionals
  flush(stdout)
  !
  CALL pw_readfile('exx', ierr)
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_dimensions()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi, eps8
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutrho, gcutm
      USE gvecs,     ONLY : gcutms, dual, doublegrid
      use gvecl,     ONLY : gcutml => gcutm
      USE gvecw,     ONLY : gcutw, ecutwfc
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutw =        ecutwfc / tpiba2
      gcutm = dual * ecutwfc / tpiba2
      ecutrho=dual * ecutwfc
      !
      doublegrid = ( dual > 4.0_dp + eps8 )
      IF ( doublegrid ) THEN
         gcutms = 4.D0 * ecutwfc / tpiba2
      ELSE
         gcutms = gcutm
      END IF
      if (do_fde .and. linterlock) gcutml = gcutm
      !
    END SUBROUTINE set_dimensions
    !
  END SUBROUTINE read_xml_file_internal
#else
SUBROUTINE read_file_dummy()
END SUBROUTINE read_file_dummy
#endif
