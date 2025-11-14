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
        real(8), allocatable :: aa(:), phi(:), phi1(:), phi_ptr(:), lambda(:), src(:), dx, dy
        integer :: nx = 50, ny = 50, N, ii = 1
        real(8) :: x_domain = 1.0, y_domain = 1.0, D = 0.001, Sigma_a = 0.3, nuSigma_f = 0.2, norm
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

            norm = (sum(phi1*nuSigma_f*dx*dy)/(lambda(ii)))
            phi = phi1/norm

            !print*,phi1(10),phi(10),lambda(ii), lambda(ii-1)
            print*,"Iteration =", ii, "Keff =", lambda(ii), "norm_check =",sum(phi*nuSigma_f*dx*dy)/(lambda(ii))

            if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-13) exit
            if (maxval(abs(phi1 - phi)) < 1.0d-8) exit
        end do

        !print'(20F6.2)',phi1
        call outpVTK(phi1, nx, ny, dx, dy)

    end subroutine iter

    subroutine iter_cyl(PCG_mode, max_iter)
      implicit none
      integer, intent(in) :: PCG_mode, max_iter
      integer, allocatable :: ia(:), ja(:)
      real(8), allocatable :: aa(:), phi(:), phi1(:), x(:), lambda(:), src(:), ri(:), r_node(:)
      integer :: nr, nth, N, iter_idx, max_outer, ii, jj
      real(8) :: R_domain, D, Sigma_a, nuSigma_f, norm, dr, dth
        
      nr = 10
      nth = 10
      R_domain = 1.0d0
      D = 0.001d0
      Sigma_a = 0.3d0
      nuSigma_f = 0.2d0

      N = nr * nth
      max_outer = 200

      allocate(phi(N), phi1(N), src(N), ri(nr))

      dr = R_domain/real(nr)
      dth = 2.0d0*PI/real(nth)

      do ii = 1, nr
        ri(ii) = (ii-0.5d0) * dr
        do jj = 1, nth
            r_node(jj+(ii-1)*nr) = ri(ii)
        end do
      end do
    
      call assemble_rth_diffusion_matrix(nr, nth, ri, dr, dth, D, Sigma_a, ia, ja, aa)

      allocate(x(N), lambda(max_outer+2), r_node(N))

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
        !print*,phi1(10),phi(10),lambda(ii), lambda(ii-1)
        print*,"Iteration =", ii, "Keff =", lambda(ii), "norm_check =",sum(phi*nuSigma_f*dr*r_node*dth)/(lambda(ii))
        
        if (abs((lambda(ii) - lambda(ii-1)) / lambda(ii-1)) < 1.0d-13) exit
        if (maxval(abs(phi1 - phi)) < 1.0d-8) exit

        iter_idx = iter_idx + 1
    
      end do

      !print'(20F6.2)',phi1
      call outpVTK(phi1, nr, nth, dr, dth)

    end subroutine iter_cyl
end module