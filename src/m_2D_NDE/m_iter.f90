module m_iter
  use m_diffusion_matrix
  use m_PCG_solver
  use m_outpVTK
  use m_constants
  implicit none
  contains
  subroutine iter(PCG_mode, max_iter)
      implicit none
      integer, intent(in) :: PCG_mode, max_iter
      integer, allocatable :: ia(:), ja(:)
      real(8), allocatable :: aa(:), phi(:), phi1(:), phi_ptr(:), lambda(:), src(:)
      integer :: nx = 100, ny = 200, N, ii = 1
      real(8) :: x_domain = 20.0, y_domain = 20.0, D = 1/0.6, Sigma_a = 0.2, nuSigma_f = 1.0, norm, dx, dy
      !real(8) :: x_domain = 2.0, y_domain = 2.0, D = 1/0.06, Sigma_a = 0.02, nuSigma_f = 0.12, norm, dx, dy
      call system("pkill gnuplot")

      N = nx*ny
      allocate(ia(ny+1), ja(N), aa(N), phi(N), phi1(N), phi_ptr(N), lambda(100), src(N))
      
      call assemble_2D_diffusion_matrix(nx, ny, x_domain, y_domain, D, Sigma_a, ia, ja, aa, dx, dy)

      phi(1:N) = 1; lambda(1) = 1

      do
          ii = ii + 1
          src = (nuSigma_f * phi) / lambda(ii-1)

          phi_ptr = phi(1:N)
          call PCG_algorithm(AA, ja, IA, phi_ptr, src, PCG_mode, max_iter)
          phi1(1:N) = phi_ptr

          lambda(ii) = lambda(ii-1) * sum(nuSigma_f*phi1*dx*dy)/sum(nuSigma_f*phi*dx*dy)

          norm = (sum(phi1*nuSigma_f*dx*dy))
          phi = phi1/norm

          !print*,phi1(10),phi(10),lambda(ii), lambda(ii-1)
          print*,"Iteration =", ii, "Keff =", lambda(ii), "norm_check =",sum(phi*nuSigma_f*dx*dy)/(lambda(ii))

          if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-13) exit
          if (maxval(abs(phi1 - phi)) < 1.0d-12) exit
      end do

      !print'(20F6.2)',phi1
      call outpVTK(phi1, nx, ny, dx, dy)

  end subroutine iter

  subroutine iter_rz(PCG_mode, max_iter)
    implicit none
    integer, intent(in) :: PCG_mode, max_iter
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: aa(:), phi(:), phi1(:), phi_ptr(:), lambda(:), src(:), r(:)
    real(8) :: dr, dz
    integer :: nr = 200, nz = 200, N, ii = 1
    !real(8) :: r_domain = 2.0, z_domain = 2.0, D = 0.04, Sigma_a = 0.02, nuSigma_f = 3.0, norm
    real(8) :: r_domain = 20.0, z_domain = 20.0, D = 1/0.6, Sigma_a = 0.2, nuSigma_f = 1.0, norm
    call system("pkill gnuplot")

    N = nr*nz
    allocate(phi(N), phi1(N), phi_ptr(N), lambda(100000), src(N))
    allocate(r(nr))

    call assemble_rz_diffusion_matrix(nr, nz, r_domain, z_domain, D, Sigma_a, ia, ja, aa, dr, dz)
    
    phi(1:N) = 1; lambda(1) = 1

    do
        ii = ii + 1
        src = (nuSigma_f * phi) / lambda(ii-1)

        phi_ptr = phi(1:N)
        call PCG_algorithm(aa, ja, ia, phi_ptr, src, PCG_mode, max_iter)
        phi1(1:N) = phi_ptr

        lambda(ii) = lambda(ii-1) * sum(nuSigma_f*phi1*dr*dz)/sum(nuSigma_f*phi*dr*dz)

        norm = sum(phi1*nuSigma_f*dr*dz)
        phi = phi1/norm

        !print*,phi1(10),phi(10),lambda(ii), lambda(ii-1)
        print*,"Iteration =", ii, "Keff =", lambda(ii), "norm_check =",sum(phi*nuSigma_f*dr*dz)/(lambda(ii))

        if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-8) exit
    end do

    !print'(20F6.2)',phi1
    call outpVTK_rz(phi, nr, nz, dr, dz)

  end subroutine iter_rz

  subroutine iter_rth(PCG_mode, max_iter)
    implicit none
    integer, intent(in) :: PCG_mode, max_iter
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: aa(:), phi(:), phi1(:), x(:), lambda(:), src(:), ri(:), r_node(:)
    integer :: nr, nth, N, iter_idx, max_outer, ii, jj
    real(8) :: R_domain, D, Sigma_a, nuSigma_f, norm, dr, dth
      
    nr = 100
    nth = 100
    R_domain = 10
    D = 3
    Sigma_a = 0.12
    nuSigma_f = 0.35

    N = nr * nth
    max_outer = 200

    allocate(phi(N), phi1(N), src(N), ri(nr), x(N), lambda(max_outer+2), r_node(N))

    dr = R_domain/real(nr)
    dth = 2.0d0*PI/real(nth)

    do ii = 1, nr
      ri(ii) = (ii-0.5d0) * dr
      do jj = 1, nth
          r_node(jj+(ii-1)*nth) = ri(ii)
      end do
    end do
  
    call assemble_rth_diffusion_matrix(nr, nth, ri, dr, dth, D, Sigma_a, ia, ja, aa)
    print*,"Matrix Assembly Successful"

    phi(1:N) = 1
    lambda(1) = 1

    iter_idx = 1
    do while (iter_idx < max_outer)
      src = (nuSigma_f * phi) / lambda(iter_idx)

      x = phi

      call PCG_algorithm(aa, ja, ia, x, src, PCG_mode, max_iter)
      
      phi1 = x

      lambda(iter_idx+1) = lambda(iter_idx) * sum(nuSigma_f*phi1*dr*r_node*dth)/sum(nuSigma_f*phi*dr*r_node*dth)

      norm = sqrt(sum((phi1**2) * r_node)*dr*dth)
      !print'(4F6.2)',phi1(10),phi(10),lambda(ii), lambda(ii-1)
      print*,"Iteration =", iter_idx, "Keff =", lambda(iter_idx), "norm_check =",sum(phi*nuSigma_f*dr*r_node*dth)/(lambda(iter_idx))
      
      if (abs((lambda(iter_idx+1) - lambda(iter_idx)) / lambda(iter_idx)) < 1.0d-13) exit
      if (maxval(abs(phi1 - phi)) < 1.0d-8) exit

      iter_idx = iter_idx + 1
      phi = phi1
  
    end do

    call outpVTK_rth(phi1, nr, nth, dr, dth)

  end subroutine iter_rth
end module