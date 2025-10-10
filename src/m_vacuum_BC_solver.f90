module m_Vacuum_BC_solver
!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Contains the subroutines required to solve the tridiagonal matrix    -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 08/10/2025    T. Charlton      Original code                          -!
!------------------------------------------------------------------------! 

implicit none
contains
    subroutine build_matrix_A_vacuum(a, b, c, N, alpha, phi, phi_analytical, L2Err)
        use m_constants
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: a(N), b(N), c(N)
        real(8), intent(out) :: phi(N), L2Err, phi_analytical(N)

        real(8), parameter :: R_domain = 1.0
        real(8), parameter :: Dif = 1.0/3.0
        real(8), parameter :: Sigma_a = 0.01
        real(8), parameter :: S0 = 1.0
        real(8) :: dx, L
        real(8) :: AA
        real(8), dimension(N) :: d
        integer :: ii

        dx = R_domain/(N-1)
        L=sqrt(Dif/Sigma_a)

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0
                b(ii) = Sigma_a + 2*Dif/dx**2 + (1/dx)*(1-alpha)/(1+alpha)
                c(ii) = -2*Dif/dx**2
                d(ii) = S0
            elseif (ii == N) then
                a(ii) = -2*Dif/dx**2
                b(ii) = Sigma_a + 2.0*Dif/dx**2 + (1/dx)*(1-alpha)/(1+alpha)
                c(ii) = 0.0
                d(ii) = S0
            else
                a(ii) = -Dif/dx**2
                b(ii) = Sigma_a + 2.0*Dif/dx**2
                c(ii) = -Dif/dx**2
                d(ii) = S0
            end if
        end do

        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")

        do ii = 1, N
            AA = (-((S0/(2*Sigma_a))*((L/(2*Dif))*sinh(R_domain/L)+1+cosh(R_domain/L)))/(cosh(R_domain/L)+((L/(4*Dif)+(Dif/L)))*sinh(R_domain/L)))
            phi_analytical(ii) = AA*cosh(dx*(ii-1)/L) + (L/(2*Dif))*(AA + S0/Sigma_a)*sinh(dx*(ii-1)/L) + S0/Sigma_a
            write(991,'(F10.5,2(1X,E15.8))') (ii-1)*dx, phi(ii), phi_analytical(ii)
        end do
        close(991)

        L2Err = sqrt(sum((phi - phi_analytical)**2)/N)

         !call inverse_power_iteration(a, b, c, N, max_iter, tol, lambda, phi)

    end subroutine build_matrix_A_vacuum

    !------------------------------------------------------------!
    ! Standard Thomas algorithm for tridiagonal systems         !
    !------------------------------------------------------------!

    subroutine thomas_algorithm(a, b, c, d, x, n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: a(n-1), b(n), c(n-1), d(n)
        real(8), intent(out) :: x(n)
        real(8), allocatable :: c_star(:), d_star(:)
        integer :: ii

        allocate(c_star(n-1))
        allocate(d_star(n))

        ! Forward sweep
        c_star(1) = c(1)/b(1)
        d_star(1) = d(1)/b(1)

        do ii = 2, n-1
            c_star(ii) = c(ii)/(b(ii)-a(ii-1)*c_star(ii-1))
        end do

        do ii = 2, n
            d_star(ii) = (d(ii)-a(ii-1)*d_star(ii-1))/(b(ii)-a(ii-1)*c_star(ii-1))
        end do

        ! Back substitution
        x(n) = d_star(n)
        do ii = n-1, 1, -1
            x(ii) = d_star(ii) - c_star(ii)*x(ii+1)
        end do

        deallocate(c_star, d_star)
    end subroutine thomas_algorithm
end module m_Vacuum_BC_solver
