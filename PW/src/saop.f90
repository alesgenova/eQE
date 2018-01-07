#ifdef __SAOP
module saop

  use kinds,                only: dp
  use fde_types, only: simulation_cell
  USE mp,            ONLY : mp_sum
  !USE gvect,                ONLY : nl, nlm, ngm, g, gg
  !USE lsda_mod,             ONLY : nspin
  !USE fft_base,             ONLY : dfftp
  use fft_types,            only : fft_type_descriptor
  USE control_flags, ONLY : gamma_only
  use fft_interfaces,       only: invfft, fwfft
  USE scf,                  ONLY : scf_type, create_scf_type, destroy_scf_type
  use io_global,            only : stdout
  !use cell_base,            only: omega, alat, tpiba
  USE wvfct,                ONLY : nbnd, wg, et
  USE klist,                ONLY : wk, nelec
  !use fde_routines, only: plot_dense
  use fde, only : saop_hirho, saop_lorho, saop_pow, saop_frac

  implicit none

    REAL(DP), PARAMETER :: epsr = 1.D-6, vanishing_charge = 1.d-10, epsg = 1.D-10
    REAL(DP), PARAMETER :: epsmin = 1.D3, epsmax = 1.D5, epsr2 = 1.d-4

  public v_rho_saop

contains
  ! Initial driver to calculate sigma grad rho and whatnot
  subroutine v_rho_saop(rho, rho_core, rhog_core, etxc, vtxc, v, cell, nspin, v_pbe, v_lb) !dfftp, nspin, nl, nlm, ngm, g, gg, omega, alat, tpiba)

    ! REFERENCE PAPER:
    ! Schipper, Gritsenko, van Gisbergen, Baerends, JCP 112, 1344, (2000)

    use xc_f90_types_m
    use xc_f90_lib_m

    implicit none
    TYPE(scf_type), INTENT(INout) :: rho
  !  TYPE(fft_type_descriptor), INTENT(IN) :: dfftp
  !  integer, intent(in) :: ngm
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP),    INTENT(IN)    :: rho_core(cell%dfftp%nnr)
    COMPLEX(DP), INTENT(IN)    :: rhog_core(cell%ngm)
    REAL(DP),    INTENT(INOUT) :: v(cell%dfftp%nnr,nspin)
    REAL(DP),    INTENT(INOUT) :: vtxc, etxc
    REAL(DP), optional,    INTENT(INOUT) :: v_pbe(cell%dfftp%nnr,nspin)
    REAL(DP), optional,    INTENT(INOUT) :: v_lb(cell%dfftp%nnr,nspin)
    

    real(dp), allocatable :: rhosum(:,:)  ! sum of rho and rhocore
    real(dp), allocatable :: auxrho(:,:)  ! sum of rho and rhocore
    complex(dp), allocatable :: rhogsum(:,:) ! same, in gspace
    complex(dp), allocatable :: rhogsum_lb(:,:) ! same, in gspace
    real(dp), allocatable :: regrhosum(:,:) ! regularized rho (interpolation between rho, and rho_smooth)
    real(dp), allocatable :: grho(:,:,:) ! the gradient of rhosum
    real(dp), allocatable :: logrho(:,:) ! the gradient of rhosum
    real(dp), allocatable :: glogrho(:,:,:) ! the gradient of rhosum
    real(dp), allocatable :: gauxrho(:,:,:) ! the gradient of rhosum
    !real(dp), allocatable :: prod(:,:,:) ! the gradient of rhosum
    !real(dp), allocatable :: gdot_prod(:) ! the gradient of rhosum
    real(dp), allocatable :: sigma(:,:)   ! dot product of rhosum gradient.
    real(dp), allocatable :: sigma_alt(:,:)   ! dot product of rhosum gradient.
                                          ! if spin polarized:
                                          ! sigma(1,:) = gradrhoup dot gradrhoup
                                          ! sigma(2,:) = gradrhoup dot gradrhodown
                                          ! sigma(3,:) = gradrhodown dot gradrhodown
    real(dp), allocatable :: v_outer(:,:), v_inner(:,:), v_xc(:,:)
    complex(dp), allocatable :: gaux(:)

    type(scf_type) :: rho_iorb

    real(dp), allocatable, dimension(:,:) :: rho_smooth
    real(dp), allocatable, dimension(:,:) :: rho_iorb_smooth

    real(dp) :: fac, domega, arhox, rhox, vmod, occupation, lnrho, e_inner, de
    real(dp) :: ehomo, eigenval, exp_coeff, weight, hirho, lorho, pow
    integer :: ik, iorb, ir, ispin, nocc, i
    integer :: band_vec(1)
    real(dp), parameter :: hart_to_ryd = 2.d0
    REAL(DP), PARAMETER :: occ_thr = 1.d-4
    character(64) :: orb_label

    if (nspin>1) then
      call errore('v_rho_saop', 'spin polarized case not implemented yet', 1)
    endif


    associate( &
      dfftp => cell%dfftp, &
      ngm => cell%ngm, &
      nl => cell%nl, &
      nlm => cell%nlm, &
      g => cell%g, &
      gg => cell%gg, &
      omega => cell%omega, &
      tpiba => cell%tpiba, &
      alat => cell%alat, &
      llarge => cell%is_native_cell &
    )

    allocate( rhosum(nspin,dfftp%nnr), &
              regrhosum(nspin,dfftp%nnr), &
              rhogsum(nspin,ngm), &
              rhogsum_lb(nspin,ngm), &
              auxrho(nspin,dfftp%nnr), &
              logrho(nspin,dfftp%nnr), &
              sigma(nspin+nspin-1, dfftp%nnr), &
              sigma_alt(nspin+nspin-1, dfftp%nnr), &
              grho(3,dfftp%nnr, nspin), &
              glogrho(3,dfftp%nnr, nspin), &
              gauxrho(3,dfftp%nnr, nspin), &
              v_xc(nspin, dfftp%nnr), &
              v_outer(nspin, dfftp%nnr), &
              v_inner(nspin, dfftp%nnr), &
              gaux(dfftp%nnr), & 
              rho_smooth(dfftp%nnr,nspin), &
              rho_iorb_smooth(dfftp%nnr,1) &
            )

    domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

    sigma = 0.d0
    v_xc = 0.d0

    fac = 1.d0/dble(nspin)

    ! first calculate rho_smooth,
    ! useful to use as a real space mask, and to calculate regularized gradients

    call smoothen(rho%of_r, rho_smooth, cell, nspin)

    ! Calculate the gradient and sigma to be used later in v_lb and v_gllb
    ! loop over spins
    do ispin = 1, nspin
      ! add core charge
      rhosum(ispin,:)  = fac * rho_core(:)  + rho%of_r(:,ispin)
      rhogsum(ispin,:) = fac * rhog_core(:) + rho%of_g(:,ispin)

      !rhogsum_lb(ispin,:) = rhogsum(ispin,:) * exp(-0.5*gg(:)*(tpiba/5.0d0)**2)

      !auxrho(ispin,:) = sqrt(abs(rhosum(ispin,:)))
      auxrho(ispin,:) = sqrt(abs(rho_smooth(:,ispin)))
      

      ! compute \nabla\rho
      if (llarge) then
        call gradrho_large(dfftp%nnr, rhogsum(ispin,:), ngm, g, nl, grho(:,:,ispin))
        !call gradrho(dfftp%nnr, rhogsum_lb(ispin,:), ngm, g, nl, glogrho(:,:,ispin))
        call gradient_large(dfftp%nnr, auxrho(ispin,:), ngm, g, nl, gauxrho(:,:,ispin))
      else
        call gradrho(dfftp%nnr, rhogsum(ispin,:), ngm, g, nl, grho(:,:,ispin))
        !call gradrho(dfftp%nnr, rhogsum_lb(ispin,:), ngm, g, nl, glogrho(:,:,ispin))
        call gradient(dfftp%nnr, auxrho(ispin,:), ngm, g, nl, gauxrho(:,:,ispin))
      endif
      
      do i=1,3; gauxrho(i,:,1) = 2.0d0*gauxrho(i,:,1) * auxrho(1,:); enddo

      hirho = -3
      lorho = -5
      do ir = 1, dfftp%nnr
        lnrho = log(abs(rho_smooth(ir,1)))
        if ( lnrho >= hirho) then
          grho(:,ir,ispin) = grho(:,ir,ispin)
          regrhosum(ispin,ir) = rhosum(ispin,ir)
        elseif (hirho > lnrho .and. lnrho >= lorho) then
          fac = (lnrho-lorho)/(hirho-lorho)
          grho(:,ir,ispin) = fac*grho(:,ir,ispin) + (1.d0-fac)*gauxrho(:,ir,ispin)
          regrhosum(ispin,ir) = fac*rhosum(ispin,ir) + (1.d0-fac)*rho_smooth(ir,ispin)
        elseif (lorho > lnrho) then
          grho(:,ir,ispin) = gauxrho(:,ir,ispin)
          regrhosum(ispin,ir) = rho_smooth(ir,ispin)
        endif
      enddo


      !gaux(nl(1:ngm)) = rhogsum_lb(ispin,1:ngm)
      !IF ( gamma_only ) gaux(nlm(1:ngm)) = CONJG(rhogsum_lb(ispin,1:ngm))
      !CALL invfft ('Dense', gaux, dfftp)
      !logrho(ispin,:) = abs(dble(gaux))


      !do ir = 1, dfftp%nnr
      !  if ((abs(rhosum(ispin,ir)) < epsmax*epsr).and.&
      !      (abs(rhosum(ispin,ir)).ge.epsmin*epsr)) then
      !      fac = ((abs(rhosum(ispin,ir))-epsmin*epsr)/(epsmax*epsr-epsmin*epsr))**2
      !      grho(:,ir,ispin)=((1.0d0-fac)*gauxrho(:,ir,ispin)+fac*grho(:,ir,ispin))
      !  elseif (abs(rhosum(ispin,ir)) < epsmin*epsr) then
      !      fac = ((abs(rhosum(ispin,ir))-epsr)/(epsmin*epsr-epsr))**2
      !      grho(:,ir,ispin)=(1.0d0-fac)*glogrho(:,ir,ispin)+fac*gauxrho(:,ir,ispin)
      !  elseif (abs(rhosum(ispin,ir)) < epsr) then
      !      grho(:,ir,ispin)=glogrho(:,ir,ispin)
      !  endif
      !enddo
  
      if (ispin == 1) then
        do ir = 1, dfftp%nnr
                 sigma(1,ir) = grho(1,ir,1)*grho(1,ir,1) &
                             + grho(2,ir,1)*grho(2,ir,1) &
                             + grho(3,ir,1)*grho(3,ir,1)
        end do
      else if (ispin == 2) then
        sigma(2,:) = grho(1,:,1)*grho(1,:,2) &
                  + grho(2,:,1)*grho(2,:,2) &
                  + grho(3,:,1)*grho(3,:,2)
        !
        sigma(3,:) = grho(1,:,2)*grho(1,:,2) &
                  + grho(2,:,2)*grho(2,:,2) &
                  + grho(3,:,2)*grho(3,:,2)
      endif
    enddo

    open(unit=668, file='sigma.dat', status="REPLACE")
    write(668,"(E15.7)") (sigma(1,ir) , ir = 1, dfftp%nr1)
    close(unit=668)

    open(unit=668, file='density.dat', status="REPLACE")
    write(668,"(E15.7)") (rho%of_r(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    !open(unit=668, file='sigma_alt.dat', status="REPLACE")
    !write(668,"(E15.7)") (sigma_alt(1,ir) , ir = 1, dfftp%nr1)
    !close(unit=668)

!    open(unit=668, file='density_alt.dat', status="REPLACE")
!    write(668,"(E15.7)") (exp(logrho(1,ir)) , ir = 1, dfftp%nr1)
!    close(unit=668)


    ! recalculate vxc through libxc, just for debugging and to make sure my implementation is right
    call v_rho_xc(regrhosum, grho, sigma, v_inner, rhosum, e_inner, cell, nspin)
    ! calculate vlb
    call v_rho_lb(regrhosum, grho, sigma, v_outer, cell, nspin)
    ! calculate vgllb
    !call v_rho_gllb(rhosum, grho, sigma, v_gllb)

    !open(unit=668, file='vxc2.dat', status="REPLACE")
    !write(668,"(E15.7)") (v_xc(1,ir) , ir = 1, dfftp%nr1)
    !close(unit=668)

    open(unit=668, file='vouter.dat', status="REPLACE")
    write(668,"(E15.7)") (v_outer(1,ir) , ir = 1, dfftp%nr1)
    close(unit=668)

    open(unit=668, file='vinner.dat', status="REPLACE")
    write(668,"(E15.7)") (v_inner(1,ir) , ir = 1, dfftp%nr1)
    close(unit=668)

    open(unit=668, file='density_smooth.dat', status="REPLACE")
    write(668,"(E15.7)") (rho_smooth(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    orbital_saop: if (.false.) then
    ! loop over the ks orbitals to build the averaged potential
    ! EQ. 2.8
      write(stdout,*) "POST FILEWRITE"

      !nocc = int(nelec/2)
      ik = 1
      !ehomo = et(nocc,ik)
      ehomo = max_occ_eigen(ik)
      call create_scf_type ( rho_iorb )
      v = 0.d0
      band_loop: do iorb = 1, nbnd
        occupation = wg(iorb, ik)/wk(ik)
        if (occupation < occ_thr) cycle band_loop
        ! et is in Ry, convert it to Ha first
        eigenval = et(iorb,ik)*0.5d0
        band_vec(1) = iorb
        rho_iorb%of_r = 0.d0
        rho_iorb%of_g = (0.d0, 0.d0)
        !write(stdout, *) ehomo, nocc, eigenval, band_vec
        !call flush_unit(stdout)
        call sum_band_selective(1, band_vec, rho_iorb)
        call regrho(rho_iorb%of_r, rho_iorb_smooth, 1)
        ! parameter to smoothly switch the summation from vinner (vgllb) to vouter (vlb)
        ! EQ. 2.9
        exp_coeff = exp(-2.d0*((ehomo-eigenval )**2))
        !write (orb_label, "('orb_',I3.3)") iorb
        !call plot_dense(rho_iorb%of_r(:,1), orb_label)
        write(stdout,*) "V_SAOP"
        write(stdout,*) iorb, sum(rho_iorb%of_r(:,1))*domega, eigenval, occupation, ehomo, exp_coeff, 1.d0-exp_coeff
        do ir = 1, dfftp%nnr
          ! EQ. 2.9
          !
          !let's try something...
          !if (rho%of_r(ir,1) > epsr) then
          !  vmod = (exp_coeff)*v_outer(1,ir) + (1.d0-exp_coeff)*v_inner(1,ir)
          !else
          !  vmod = v_outer(1,ir)
          !endif

          vmod = (exp_coeff)*v_outer(1,ir) + (1.d0-exp_coeff)*v_inner(1,ir)
          
          ! hard confine between 0 and 1
          !weight = min(max((rho_iorb%of_r(ir,1) / rhosum(1,ir))*cutoff(abs(rhosum(1,ir)),1.d-4, 3.d-5), 0.d0), 1.d0)
          !!weight = (rho_iorb%of_r(ir,1) / rhosum(1,ir))*cutoff(abs(rho_iorb%of_r(ir,1)),1.d-4, 3.d-5)
          
          !weight = regratio(rho_iorb%of_r(ir,1), rho%of_r(ir,1),1.d-10, 3.d-11, exp_coeff)
          !weight = rho_iorb_smooth(ir,1)/rho_smooth(ir,1)
          weight = regratio2(rho_iorb_smooth(ir,1),rho_smooth(ir,1))
          v(ir,1) = v(ir,1) + (occupation*vmod*weight)
          !if ( (rho_iorb%of_r(ir,1) > vanishing_charge) .and. (rhosum(1,ir) > vanishing_charge) ) then
          !  vmod = (exp_coeff)*v_lb(1,ir) + (1.d0-exp_coeff)*v_gllb(1,ir)
          !  v(ir,1) = v(ir,1) + (occupation*vmod*rho_iorb%of_r(ir,1) / rhosum(1,ir))*cutoff(abs(rhosum(1,ir)),1.d-4, 3.d-5)
          !endif
        enddo
      enddo band_loop
      call destroy_scf_type ( rho_iorb )
    endif orbital_saop

    !call saop_mixing_real(v_inner, v_outer, rho_smooth, v, cell, nspin)
    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 10.d0)
    open(unit=668, file='vsaop10.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 1.d0)
    open(unit=668, file='vsaop1.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 0.5d0)
    open(unit=668, file='vsaop05.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 0.1d0)
    open(unit=668, file='vsaop01.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 0.01d0)
    open(unit=668, file='vsaop001.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 0.001d0)
    open(unit=668, file='vsaop0001.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 0.0001d0)
    open(unit=668, file='vsaop00001.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    call saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, 0.d0)
    open(unit=668, file='vsaop.dat', status="REPLACE")
    write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    close(unit=668)

    !open(unit=668, file='vsaop2.dat', status="REPLACE")
    !write(668,"(E15.7)") (v(ir,1) , ir = 1, dfftp%nr1)
    !close(unit=668)
    de = 0.d0
    do ir = 1, dfftp%nnr
      if (rhosum(1,ir) > epsr) then
        de = de + (v(ir,1)-v_inner(1,ir))*rhosum(1,ir)
      endif
    enddo

    CALL mp_sum(  de , dfftp%comm )
    de = de*hart_to_ryd*domega

    write(stdout,*) "E_INNER: ", e_inner
    write(stdout,*) "dE: ", de
    etxc = e_inner + de

    !Manual overrides, just for debug
    !etxc = e_inner
    do ir = 1, dfftp%nnr
      !v(ir,1) = v_outer(1,ir)
      !v(ir,1) = v_inner(1,ir)
      !if (rho%of_r(ir,1) < epsr) then
        !v(ir,1) = v_outer(1,ir)
      !endif
    enddo

    if (present(v_pbe)) then
      do ir = 1, dfftp%nnr
       v_pbe(ir,1) = v_inner(1,ir)
      enddo
    endif

    if (present(v_lb)) then
      do ir = 1, dfftp%nnr
       v_lb(ir,1) = v_outer(1,ir)
      enddo
    endif

    !de = sum(rho%of_r(:,1))
    !CALL mp_sum(  de , dfftp%comm )
    !de = de*domega
    !write(*,*) "intrho: ", de

    !de = sum(rho%of_r(:,1)*rho%of_r(:,1))
    !CALL mp_sum(  de , dfftp%comm )
    !de = de*domega
    !write(*,*) "intrho2: ", de

    !de = sum(rho_smooth(:,1))
    !CALL mp_sum(  de , dfftp%comm )
    !de = de*domega
    !write(*,*) "intrhosmooth: ", de

    !de = sum(rho_smooth(:,1)*rho_smooth(:,1))
    !CALL mp_sum(  de , dfftp%comm )
    !de = de*domega
    !write(*,*) "intrhosmooth2: ", de

    !de = sum(rho_core)
    !CALL mp_sum(  de , dfftp%comm )
    !de = de*domega
    !write(*,*) "intrhocore: ", de

    !de = sum(rho_core*rho_core)
    !CALL mp_sum(  de , dfftp%comm )
    !de = de*domega
    !write(*,*) "intrhocore2: ", de


  end associate

  end subroutine v_rho_saop


  subroutine saop_mixing_real(v_inner, v_outer, rho_smooth, v, cell, nspin)
    implicit none

    ! parameters
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP),    INTENT(IN) :: v_inner(nspin, cell%dfftp%nnr)
    REAL(DP),    INTENT(IN) :: v_outer(nspin, cell%dfftp%nnr)
    REAL(DP),    INTENT(IN) :: rho_smooth(nspin, cell%dfftp%nnr)
    REAL(DP),    INTENT(INOUT) :: v(cell%dfftp%nnr, nspin)

    ! local vars
    integer :: is, ir
    real(dp) :: fac, lnrho


    associate( &
      dfftp => cell%dfftp, &
      ngm => cell%ngm, &
      nl => cell%nl, &
      nlm => cell%nlm, &
      g => cell%g, &
      gg => cell%gg, &
      omega => cell%omega, &
      tpiba => cell%tpiba, &
      alat => cell%alat, &
      llarge => cell%is_native_cell &
    )

    do is = 1, nspin
      do ir = 1, dfftp%nnr
        lnrho = log(abs(rho_smooth(ir,is)))
        if ( lnrho >= saop_hirho) then
          v(ir,is) = v_inner(is,ir)
        elseif (saop_hirho > lnrho .and. lnrho >= saop_lorho) then
          fac = ((saop_hirho-lnrho)/(saop_hirho-saop_lorho))**saop_pow
          fac = fac * saop_frac
          v(ir,is) = (1.d0-fac)*v_inner(is,ir) + fac*v_outer(is,ir)
        elseif (saop_lorho > lnrho) then
          v(ir,is) = saop_frac * v_outer(is,ir) + (1.d0-saop_frac)*v_inner(is,ir)
        endif
      enddo
    enddo

    end associate

  end subroutine saop_mixing_real

  subroutine saop_mixing_reciprocal(v_inner, v_outer, v, cell, nspin, a)
    USE constants, ONLY : fpi, e2

    implicit none

    ! parameters
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP),    INTENT(IN) :: v_inner(nspin, cell%dfftp%nnr)
    REAL(DP),    INTENT(IN) :: v_outer(nspin, cell%dfftp%nnr)
    REAL(DP),    INTENT(INOUT) :: v(cell%dfftp%nnr, nspin)
    real(dp), intent(in) :: a

    ! local vars
    integer :: is, ir, ig
    real(dp) :: fac, invfac, lnrho
    complex(dp), dimension(:), allocatable :: aux
    complex(dp), dimension(:), allocatable :: v_inner_g, v_outer_g
    complex(dp), dimension(:), allocatable :: r12_sr, r12_lr !, r12_tot
    !real(dp), parameter :: a = 0.0d0
    real(dp) :: aa


    ! only numerically unstable case should be a=0,
    ! but this is simply v = vinner
    if (a<1.d-6) then
      do is=1, nspin
        v(:,is) = v_inner(is,:)
      enddo
      return
    endif


    associate( &
      dfftp => cell%dfftp, &
      ngm => cell%ngm, &
      nl => cell%nl, &
      nlm => cell%nlm, &
      g => cell%g, &
      gg => cell%gg, &
      gstart => cell%gstart, &
      omega => cell%omega, &
      tpiba => cell%tpiba, &
      tpiba2 => cell%tpiba2, &
      alat => cell%alat, &
      llarge => cell%is_native_cell &
    )

    aa = a*a

    allocate(aux(dfftp%nnr), v_inner_g(ngm), v_outer_g(ngm))
    allocate(r12_sr(ngm), r12_lr(ngm))!, r12_tot(ngm))

    do is=1, nspin
      aux(:) = CMPLX(v_inner(is,:),0.d0,kind=dp)
      call fwfft ('Custom', aux, dfftp)
      v_inner_g(:) = aux(nl(:))

      aux(:) = CMPLX(v_outer(is,:),0.d0,kind=dp)
      call fwfft ('Custom', aux, dfftp)
      v_outer_g(:) = aux(nl(:))

      ! check for gstart, just in case
      if (gstart==1)then
        r12_sr(1) = 0.d0
        r12_lr(1) = v_outer_g(ig)
      endif
      do ig = gstart, ngm
        !r12_tot(ig) = 1.D0 / gg(ig)
        !r12_lr(ig) = r12_tot(ig) - r12_sr(ig)
        ! calculate rho_{PBE} * 1/r_{SR}
        r12_sr(ig) = v_inner_g(ig) * gg(ig) / (gg(ig)+aa)
        r12_lr(ig) = v_outer_g(ig) * (1.d0 - gg(ig) / (gg(ig)+aa))
      enddo

      !fac = e2 * fpi / tpiba2
      !invfac = tpiba2/(e2 * fpi)
      !invfac = 1.d0/fac
      !r12_sr = r12_sr * fac
      !r12_lr = r12_lr * fac
      !r12_tot = r12_tot * fac

      !aux(nl(1:ngm)) = CMPLX ( r12_sr(1:ngm)*real(v_inner_g(1:ngm,is)), r12_sr(1:ngm)*aimag(v_inner_g(1:ngm,is)), KIND=dp )
      ! shouldn't this be equivalent?
      !aux(nl(1:ngm)) = r12_sr(1:ngm)*v_inner_g(1:ngm,is)*gg(1:ngm)*invfac
      !aux(nl(1:ngm)) = r12_lr(1:ngm)*v_outer_g(1:ngm,is)*gg(1:ngm)*invfac + aux(nl(1:ngm))
      aux(nl(1:ngm)) = r12_sr(1:ngm)
      aux(nl(1:ngm)) = r12_lr(1:ngm) + aux(nl(1:ngm))
      IF ( gamma_only ) then
        aux(nlm(1:ngm)) = CONJG(aux(nl(1:ngm)))
      endif
      
      CALL invfft ('Custom', aux, dfftp)
      !
      v(:,is) = dble(aux(:))

    enddo

    deallocate(aux, v_inner_g, v_outer_g)
    deallocate(r12_sr, r12_lr)!, r12_tot)

    end associate

  end subroutine saop_mixing_reciprocal


! Debug Subroutine to calculate the XC Potential and compare it to the reference
  subroutine v_rho_xc(rho, grho, sigma, v, rho_orig, e, cell, nspin)

    use xc_f90_types_m
    use xc_f90_lib_m


    implicit none
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP), INTENT(INout) :: rho(nspin,cell%dfftp%nnr)
    REAL(DP), INTENT(INout) :: rho_orig(nspin,cell%dfftp%nnr)
    REAL(DP), INTENT(INout) :: grho(3,cell%dfftp%nnr,nspin)
    REAL(DP), INTENT(INout) :: sigma(nspin+nspin-1,cell%dfftp%nnr)
    REAL(DP),    INTENT(INOUT) :: v(nspin, cell%dfftp%nnr)
    real(dp), intent(inout) :: e

    real(dp), allocatable :: prod(:,:,:) ! the gradient of rhosum
    real(dp), allocatable :: gdot_prod(:) ! the gradient of rhosum
    real(dp), allocatable :: aux(:,:),exc(:,:)     ! auxiliary array for the xc energy density
    real(dp), allocatable :: aux2(:,:),vrho(:,:)   ! auxiliary array for the xc energy density
    real(dp), allocatable :: aux3(:,:),vsigma(:,:) ! auxiliary array for the xc energy density

    integer :: ifunc, ispin, funcs_id(2), ir

    real(dp) :: fac, domega, arhox, rhox
    real(dp), parameter :: hart_to_ryd = 2.d0
    real(dp), parameter :: occ_thr = 1.d-4

    !-- libxc stuff
    TYPE(xc_f90_pointer_t) :: xc_func
    TYPE(xc_f90_pointer_t) :: xc_info

    associate( &
      dfftp => cell%dfftp, &
      ngm => cell%ngm, &
      nl => cell%nl, &
      nlm => cell%nlm, &
      g => cell%g, &
      gg => cell%gg, &
      omega => cell%omega, &
      tpiba => cell%tpiba, &
      alat => cell%alat, &
      llarge => cell%is_native_cell &
    )

    allocate( &
              aux(nspin,dfftp%nnr), &
              aux2(nspin,dfftp%nnr), &
              aux3(nspin,dfftp%nnr), &
              exc(nspin,dfftp%nnr), &
              vrho(nspin,dfftp%nnr), &
              vsigma(nspin,dfftp%nnr), &
              prod(3,dfftp%nnr,nspin), &
              gdot_prod(dfftp%nnr) )

    funcs_id = (/XC_GGA_X_PBE, XC_GGA_C_PBE/)
    exc = 0.d0
    vrho = 0.d0
    vsigma = 0.d0

    do ifunc = 1, 2
      aux = 0.d0
      aux2 = 0.d0
      aux3= 0.d0
      call xc_f90_func_init(xc_func, xc_info, funcs_id(ifunc), nspin)
      select case (xc_f90_info_family(xc_info))
       case(XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(xc_func, dfftp%nnr, rho(1,1), aux(1,1), aux2(1,1))
       case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
         call xc_f90_gga_exc_vxc(xc_func, dfftp%nnr, rho(1,1), sigma(1,1), aux(1,1), aux2(1,1), aux3(1,1))
         !call xc_f90_gga_exc(xc_func, dfftp%nnr, rhosum(1,1), sigma(1,1), aux(1,1))
      end select
      exc = exc + aux

      vrho = vrho + aux2
      vsigma = vsigma + aux3
    enddo

    v = 0.d0
    prod = 0.d0
    gdot_prod = 0.d0

    !effin numerics
    do ir = 1, dfftp%nnr
      arhox = abs(rho(1,ir))
      if ( arhox > epsr) then
        if ( sigma(1,ir) > epsg) then
          prod(1,ir,1) = hart_to_ryd*hart_to_ryd*vsigma(1,ir) * grho(1,ir,1)!*cutoff(arhox,1.d-4, 3.d-5)! / sqrt(sigma(1,ir))
          prod(2,ir,1) = hart_to_ryd*hart_to_ryd*vsigma(1,ir) * grho(2,ir,1)!*cutoff(arhox,1.d-4, 3.d-5)! / sqrt(sigma(1,ir))
          prod(3,ir,1) = hart_to_ryd*hart_to_ryd*vsigma(1,ir) * grho(3,ir,1)!*cutoff(arhox,1.d-4, 3.d-5)! / sqrt(sigma(1,ir))
        else
          prod(:,ir,1) = 0.d0
        endif
      else
        prod(:,ir,1) = 0.d0
      endif

      if ( arhox < vanishing_charge ) then
      !if (.true.) then
        ! apparently vrho can be unstable too
        vrho(1,ir) = 0.0d0 !vrho(1,ir)*cutoff(arhox,1.d-4, 3.d-5)!0.d0
      endif
    enddo

    if (llarge) then
      call grad_dot_large(dfftp%nnr, prod(1,1,1), ngm, g, nl, alat, gdot_prod)
    else
      call grad_dot(dfftp%nnr, prod(1,1,1), ngm, g, nl, alat, gdot_prod)
    endif

    v(1,:) = hart_to_ryd*vrho(1,:) - gdot_prod(:)
    e = sum(exc*rho_orig)
    CALL mp_sum(  e , dfftp%comm )
    domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
    e = e*hart_to_ryd*domega

    end associate

  end subroutine v_rho_xc

  ! Subroutine to calculate the LB Potential in the short range
  subroutine v_rho_lb(rho, grho, sigma, v, cell, nspin)

    use xc_f90_types_m
    use xc_f90_lib_m


    implicit none
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP), INTENT(INout) :: rho(nspin,cell%dfftp%nnr)
    REAL(DP), INTENT(INout) :: grho(3,cell%dfftp%nnr,nspin)
    REAL(DP), INTENT(INout) :: sigma(nspin+nspin-1,cell%dfftp%nnr)
    REAL(DP),    INTENT(INOUT) :: v(nspin, cell%dfftp%nnr)

    real(dp), allocatable :: prod(:,:,:) ! the gradient of rhosum
    real(dp), allocatable :: gdot_prod(:) ! the gradient of rhosum
    real(dp), allocatable :: aux(:,:),exc(:,:)     ! auxiliary array for the xc energy density
    real(dp), allocatable :: aux2(:,:),vrho(:,:)   ! auxiliary array for the xc energy density
    real(dp), allocatable :: aux3(:,:),vsigma(:,:) ! auxiliary array for the xc energy density

    integer :: ifunc, ispin, funcs_id(2), ir

    real(dp) :: fac, domega, arhox, rhox
    real(dp) :: auxrho(nspin,cell%dfftp%nnr)
    real(dp), parameter :: hart_to_ryd = 2.d0
    real(dp), parameter :: occ_thr = 1.d-4
    REAL(DP), PARAMETER :: epsr = 1.D-10
    REAL(DP), PARAMETER :: epsg = 1.D-10

    !-- libxc stuff
    TYPE(xc_f90_pointer_t) :: xc_func
    TYPE(xc_f90_pointer_t) :: xc_info

    associate( &
      dfftp => cell%dfftp, &
      ngm => cell%ngm, &
      nl => cell%nl, &
      nlm => cell%nlm, &
      g => cell%g, &
      gg => cell%gg, &
      omega => cell%omega, &
      tpiba => cell%tpiba, &
      alat => cell%alat, &
      llarge => cell%is_native_cell &
    )

    allocate( &
              aux(nspin,dfftp%nnr), &
              aux2(nspin,dfftp%nnr), &
              aux3(nspin,dfftp%nnr), &
              exc(nspin,dfftp%nnr), &
              vrho(nspin,dfftp%nnr), &
              vsigma(nspin,dfftp%nnr), &
              prod(3,dfftp%nnr,nspin), &
              gdot_prod(dfftp%nnr) )

    ! the LB (mod) implementation in libxc only includes the exchange part, 
    ! an LDA correlation functional needs to be picked as well.
    funcs_id = (/XC_GGA_X_LBM, XC_LDA_C_PW/)
    exc = 0.d0
    vrho = 0.d0
    vsigma = 0.d0

    !auxrho = abs(rho)
    auxrho = rho
    do ifunc = 1, 2
      aux = 0.d0
      aux2 = 0.d0
      aux3= 0.d0
      call xc_f90_func_init(xc_func, xc_info, funcs_id(ifunc), nspin)
      select case (xc_f90_info_family(xc_info))
       case(XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(xc_func, dfftp%nnr, auxrho(1,1), aux(1,1), aux2(1,1))
       case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
         call xc_f90_gga_vxc(xc_func, dfftp%nnr, auxrho(1,1), sigma(1,1), aux2(1,1), aux3(1,1))
         !call xc_f90_gga_exc(xc_func, dfftp%nnr, rhosum(1,1), sigma(1,1), aux(1,1))
      end select
      do ir=1,dfftp%nnr
      if (auxrho(1,ir)>=epsr) then
         exc(1,ir) = exc(1,ir) + aux(1,ir)
         vrho(1,ir) = vrho(1,ir) + aux2(1,ir)
         if (sigma(1,ir)>=epsg) vsigma(1,ir) = vsigma(1,ir) + aux3(1,ir)
      endif                                                         
      enddo
    enddo

    v = 0.d0
    prod = 0.d0
    gdot_prod = 0.d0

    !effin numerics
    do ir = 1, dfftp%nnr
      arhox = abs(rho(1,ir))
      if ( arhox > epsr) then
        if ( sigma(1,ir) >= epsg) then
          prod(1,ir,1) = hart_to_ryd*hart_to_ryd*vsigma(1,ir) * grho(1,ir,1)!*cutoff(arhox,1.d-4, 3.d-5)! / sqrt(sigma(1,ir))
          prod(2,ir,1) = hart_to_ryd*hart_to_ryd*vsigma(1,ir) * grho(2,ir,1)!*cutoff(arhox,1.d-4, 3.d-5)! / sqrt(sigma(1,ir))
          prod(3,ir,1) = hart_to_ryd*hart_to_ryd*vsigma(1,ir) * grho(3,ir,1)!*cutoff(arhox,1.d-4, 3.d-5)! / sqrt(sigma(1,ir))
        else
          prod(:,ir,1) = 0.d0
        endif
      else
        prod(:,ir,1) = 0.d0
      endif
    enddo

    if (llarge) then
      call grad_dot_large(dfftp%nnr, prod(1,1,1), ngm, g, nl, alat, gdot_prod)
    else
      call grad_dot(dfftp%nnr, prod(1,1,1), ngm, g, nl, alat, gdot_prod)
    endif

    v(1,:) = hart_to_ryd*vrho(1,:) - gdot_prod(:)

  end associate
    
  end subroutine v_rho_lb

  ! Subroutine to calculate the GLLB potential in the long range
  subroutine v_rho_gllb(rho, grho, sigma, v, cell, nspin)

    use xc_f90_types_m
    use xc_f90_lib_m


    implicit none
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP), INTENT(INout) :: rho(nspin,cell%dfftp%nnr)
    REAL(DP), INTENT(INout) :: grho(3,cell%dfftp%nnr,nspin)
    REAL(DP), INTENT(INout) :: sigma(nspin+nspin-1,cell%dfftp%nnr)
    REAL(DP),    INTENT(INOUT) :: v(nspin,cell%dfftp%nnr)

    real(dp), allocatable :: aux(:,:),exc(:,:)   ! auxiliary array for the xc energy density

    type(scf_type) :: rho_iorb
    real(dp), allocatable :: v_hole(:,:), v_resp(:,:)

    integer :: ifunc, ispin, funcs_id(2), ir
    integer :: iocc, ik, iorb, nocc
    real(dp) :: ehomo, eigenval
    integer :: band_vec(1)

    real(dp) :: fac, domega, arhox, rhox, occupation
    real(dp), parameter :: hart_to_ryd = 2.d0
    real(dp), parameter :: occ_thr = 1.d-4



    !-- libxc stuf
    TYPE(xc_f90_pointer_t) :: xc_func
    TYPE(xc_f90_pointer_t) :: xc_info

    associate( &
      dfftp => cell%dfftp, &
      ngm => cell%ngm, &
      nl => cell%nl, &
      nlm => cell%nlm, &
      g => cell%g, &
      gg => cell%gg, &
      omega => cell%omega, &
      tpiba => cell%tpiba, &
      alat => cell%alat &
    )

    allocate( &
              aux(nspin,dfftp%nnr), &
              exc(nspin,dfftp%nnr), &
              v_hole(nspin,dfftp%nnr), &
              v_resp(nspin,dfftp%nnr) &
            )


    domega = omega / dble(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

    ! V_gllb = V_hole + V_resp
    ! EQ. 2.2

    ! V Hole = 2\epsilon_x(B) + 2\epsilon_c(PW) (I think GPAW implementation uses PBEsol instead of PW91)
    funcs_id = (/XC_GGA_X_B88, XC_GGA_C_PW91/)
    exc = 0.d0
    ! aux stores exc
    ! aux2 stores vrho := dexc/drho
    ! aux3 stores vsigma := dexc/dsigma
    do ifunc = 1, 2
      aux = 0.d0
      call xc_f90_func_init(xc_func, xc_info, funcs_id(ifunc), nspin)
      select case (xc_f90_info_family(xc_info))
       case(XC_FAMILY_LDA)
         call xc_f90_lda_exc(xc_func, dfftp%nnr, rho(1,1), aux(1,1))
       case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
         call xc_f90_gga_exc(xc_func, dfftp%nnr, rho(1,1), sigma(1,1), aux(1,1))
         !call xc_f90_gga_exc(xc_func, dfftp%nnr, rhosum(1,1), sigma(1,1), aux(1,1))
      end select
      exc = exc + aux
    enddo

    !effin numerics
    do ir = 1, dfftp%nnr
      arhox = abs(rho(1,ir))
      exc(1,ir) = exc(1,ir)*cutoff(arhox,1.d-4, 3.d-5)
    enddo

    v_hole = 2.d0*exc*hart_to_ryd



    ! V resp 
    ! EQ. 2.6

    nocc = int(nelec/2)
    ik = 1
    !ehomo = et(nocc,ik)
    ehomo = max_occ_eigen(ik)
    v_resp = 0.d0

    call create_scf_type ( rho_iorb )

    band_loop: do iorb = 1, nbnd
      occupation = wg(iorb, ik)/wk(ik)
      if (occupation < occ_thr) cycle band_loop
      ! et is in Ry, convert it to Ha first
      eigenval = et(iorb,ik)*0.5d0
      band_vec(1) = iorb
      call sum_band_selective(1, band_vec, rho_iorb)
      aux = 0.d0
      write(stdout,*) "V_GLLB"
      write(stdout,*) iorb, sum(rho_iorb%of_r(:,1))*domega, eigenval, occupation, ehomo, sqrt(ehomo - eigenval)
      do ir = 1, dfftp%nnr
        !if (rho_iorb%of_r(ir,1) > vanishing_charge .and. rho(1,ir) > vanishing_charge) then
        if (.true.) then
          aux(1,ir) = regratio(rho_iorb%of_r(ir,1), rho(1,ir),1.d-4, 3.d-5, 1.d0)
        else
          aux(1,ir) = 0.d0
        endif
      enddo
      aux = occupation*sqrt(ehomo - eigenval) * aux
      v_resp = v_resp + aux
    enddo band_loop

    call destroy_scf_type ( rho_iorb )

    ! K\sigma = 0.42
    v_resp = v_resp * 0.42d0 * hart_to_ryd

    v = v_hole + v_resp

    end associate

  end subroutine v_rho_gllb

  subroutine atom_mask(tau, strf, nat, ntyp, ngm, dfftp, msh, rgrid, &
                       nl, nlm, ngl, gl, igtongl, omega, tpiba2, maskr)

    use fft_interfaces, only: invfft
    use fft_types, only: fft_type_descriptor
    use radial_grids, only: radial_grid_type
    USE control_flags, ONLY : gamma_only
    USE uspp_param,ONLY : upf

    implicit none
    ! inputs
    integer, intent(in) :: nat, ntyp, ngm, msh(ntyp)
    real(dp), intent(in) :: tau(3,nat)
    complex(dp), intent(in) :: strf(ntyp, ngm)
    type(fft_type_descriptor), intent(in) :: dfftp
    type(radial_grid_type), intent(in) :: rgrid(:)
    integer, intent(in) :: ngl
    integer, intent(in) :: nl(ngm), nlm(ngm), igtongl(ngm)
    real(dp), intent(in) :: gl(ngl)
    real(dp), intent(in) :: omega, tpiba2

    ! outputs
    real(dp), intent(inout) :: maskr(dfftp%nnr)

    ! local vars
    complex(dp), allocatable :: gaux(:), maskg(:)
    real(dp), allocatable :: aux(:,:), gauss(:)
    integer :: nt, ir, ng, it
    real(dp) :: r, rcut, sigma

    allocate (maskg( ngl))
    allocate (gaux( dfftp%nnr))
    gaux (:) = (0.0_DP, 0.0_DP)
    do nt = 1, ntyp
      allocate (gauss(1:msh(nt)))
      sigma = 0.1d0
      rcut = maxval(upf(nt)%rcut(:))
      !rcut = 1.d0
      do it = 1, msh(nt)
        r = rgrid(nt)%r(it)
        if (r <= rcut) then
          gauss(it) = 1.d0
        else
          gauss(it) = exp(-((r-rcut)/(sigma))**2)
        endif
      enddo

      ! a little comment here could have been nice...
      call drhoc (ngl, gl, omega, tpiba2, msh (nt), rgrid(nt)%r, &
           rgrid(nt)%rab, gauss, maskg)

      !if (nt == 1) then
      !do ir = 1 , msh(nt) - 1
      !   write(666,*) rgrid(nt)%r(ir), gauss(ir), rgrid(nt)%r(ir+1) - rgrid(nt)%r(ir)
      !enddo
      !endif

      deallocate (gauss)

      !     multiply by the structure factor and sum
      do ng = 1, ngm
         gaux(nl(ng)) = gaux(nl(ng)) + strf(ng,nt) * maskg(igtongl(ng))
      enddo
    enddo

    if (gamma_only) then
      do ng = 1, ngm
         gaux(nlm(ng)) = CONJG(gaux(nl(ng)))
      end do
    end if

    !mask_g(:) = aux(nl(:)) # no need for the mask in g space

    !   the core charge in real space
    CALL invfft ('Custom', gaux, dfftp)

    do ir = 1, dfftp%nnr
       !rhoneg = rhoneg + min (0.d0,  DBLE (aux (ir) ) )
       !rhoima = rhoima + abs (AIMAG (aux (ir) ) )
       maskr(ir) =  DBLE (gaux(ir))
       !corecharge = corecharge + rho_gauss(ir)
       !
    enddo

  end subroutine atom_mask
  ! Mix LB and GLLB potentials to obtain the statistically averaged SAOP.

  real(dp) function max_occ_eigen(ik)

    use kinds, only: DP
    USE wvfct, ONLY : nbnd, wg, et
    USE klist, ONLY : wk

    implicit none

    integer, intent(in) :: ik
    integer :: ibnd
    real(dp), parameter :: occ_thr = 1.d-4
    real(dp) :: occupation, eigenval

    max_occ_eigen = -1.d20

    ! et is in Ry, convert it to Ha first

    band_loop: do ibnd=1, nbnd
      occupation = wg(ibnd, ik)/wk(ik)
      if (occupation < occ_thr) cycle band_loop
      eigenval = et(ibnd,ik)*0.5d0
      if (eigenval > max_occ_eigen) max_occ_eigen = eigenval
    end do band_loop

  end function max_occ_eigen

  real(dp) function cutoff(rho, m, w)
    implicit none

    real(dp), intent(in) :: rho, m, w

    !cutoff = 0.5*(erf((rho-m)/w) + 1.d0)
    cutoff = 1.0d0

  end function cutoff


  real(dp) function newrho(rho)
    implicit none

    real(dp), intent(in) :: rho

    newrho = log(rho)

  end function newrho



  real(dp) function regratio(nom, denom, m, w, exp_coeff)
    implicit none

    real(dp), intent(in) :: nom, denom, m, w, exp_coeff
    
    if (abs(nom) <= abs(denom) .and. abs(denom) > epsr) then
      !regratio = abs(nom*cutoff(nom,m,w)/denom)
      regratio = abs(nom*cutoff(nom,m,w)/denom)
    else
      !regratio = 0.d0!cutoff(nom, m,w)
      ! this might be good
      regratio = exp_coeff
    endif
    !regratio = nom/denom
    !regratio = min(1.0d0,regratio)
    !regratio = max(0.0d0,regratio)

  end function regratio


  real(dp) function regratio2(nom, denom)
    implicit none

    real(dp), intent(in) :: nom, denom
    
    if (nom >= 1.d-7) then
      !regratio = abs(nom*cutoff(nom,m,w)/denom)
      regratio2 = abs(nom/denom)
    else
      !regratio = 0.d0!cutoff(nom, m,w)
      ! this might be good
      regratio2 = 0.d0
    endif
    !regratio = nom/denom
    !regratio2 = min(1.0d0,regratio)
    !regratio2 = max(0.0d0,regratio)

  end function regratio2


  subroutine smoothen(rhoin, rhoout, cell, nspin)

  USE control_flags, ONLY : gamma_only
  !USE wavefunctions_module, only: psic

  implicit none

  type(simulation_cell), intent(in) :: cell
  integer, intent(in) :: nspin
  real(dp), dimension(:,:), intent(in) :: rhoin
  real(dp), dimension(:,:), intent(inout) :: rhoout

  real(dp), parameter :: epsr = 1.0e-4
  real(dp), dimension(:,:), allocatable :: rhosmooth
  complex(dp), dimension(:), allocatable :: aux, gaux
  integer :: is, ir

  associate( &
    dfftp => cell%dfftp, &
    ngm => cell%ngm, &
    nl => cell%nl, &
    nlm => cell%nlm, &
    g => cell%g, &
    gg => cell%gg, &
    omega => cell%omega, &
    tpiba => cell%tpiba, &
    alat => cell%alat &
  )

  allocate(aux(dfftp%nnr), gaux(ngm), rhosmooth(dfftp%nnr, nspin))

  do is=1, nspin
    aux(:) = rhoin(:,is)
    call fwfft ('Custom', aux, dfftp)
    gaux(:) = aux(nl(:))

    ! smoothen
    !gaux(1:ngm) = gaux(1:ngm) * exp(-0.5*gg(1:ngm)*(tpiba/5.0d0)**2)
    gaux(1:ngm) = gaux(1:ngm) * exp(-0.5*gg(1:ngm)*(tpiba/3.5d0)**2)
    ! invfft
    aux(nl(1:ngm)) = gaux(1:ngm)

    IF ( gamma_only ) then
      !aux(nlm(1:ngm)) = CONJG(gaux(1:ngm))
    endif

    CALL invfft ('Custom', aux, dfftp)
    !
    rhosmooth(:,is) = dble(aux(:))
  enddo
  
  do is=1, nspin
    do ir = 1, dfftp%nnr
      if (rhoin(ir,is) < epsr) then
        rhoout(ir,is) = rhosmooth(ir,is)
      else
        rhoout(ir,is) = rhoin(ir,is)
      endif
    enddo
  enddo
  rhoout = rhosmooth

  end associate

  end subroutine smoothen

end module saop
#else
module saop
  use kinds,                only: dp
  use fde_types, only: simulation_cell
  USE scf,                  ONLY : scf_type
  
  implicit none

  public v_rho_saop

contains
  ! Initial driver to calculate sigma grad rho and whatnot
  subroutine v_rho_saop(rho, rho_core, rhog_core, etxc, vtxc, v, cell, nspin, v_pbe, v_lb) !dfftp, nspin, nl, nlm, ngm, g, gg, omega, alat, tpiba)

    implicit none
    TYPE(scf_type), INTENT(INout) :: rho
  !  TYPE(fft_type_descriptor), INTENT(IN) :: dfftp
  !  integer, intent(in) :: ngm
    type(simulation_cell), intent(in) :: cell
    integer, intent(in) :: nspin
    REAL(DP),    INTENT(IN)    :: rho_core(cell%dfftp%nnr)
    COMPLEX(DP), INTENT(IN)    :: rhog_core(cell%ngm)
    REAL(DP),    INTENT(INOUT) :: v(cell%dfftp%nnr,nspin)
    REAL(DP),    INTENT(INOUT) :: vtxc, etxc
    REAL(DP), optional,    INTENT(INOUT) :: v_pbe(cell%dfftp%nnr,nspin)
    REAL(DP), optional,    INTENT(INOUT) :: v_lb(cell%dfftp%nnr,nspin)

    return
  end subroutine v_rho_saop
end module saop
#endif
