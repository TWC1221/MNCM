module m_iterative_power_solver
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
    subroutine build_matrix_A_iterative_power(a, b, c, N, alpha, phi, dx, Dif, Sigma_a, S0)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: a(N), b(N), c(N), phi(N), dx(N-1), Dif(N), Sigma_a(N), S0(N)

        real(8), dimension(N) :: L, d, x
        integer :: ii

        L = sqrt(Dif(1:N)/Sigma_a(1:N))

        x(1) = 0.0d0
        do ii = 2, N
            x(ii) = x(ii-1) + dx(ii-1)
        end do

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0d0
                b(ii) = ((Dif(ii)+Dif(ii+1))/(dx(ii)*dx(ii)))+Sigma_a(ii) + 1/dx(ii)*(1-alpha)/(1+alpha)
                c(ii) = -(Dif(ii)+Dif(ii+1))/(dx(ii)*(dx(ii)))
                d(ii) = S0(ii)
            elseif (ii == N) then
                a(ii) = -(Dif(ii-1)+Dif(ii))/(dx(ii-1)*dx(ii-1))
                b(ii) = ((Dif(ii)+Dif(ii-1))/dx(ii-1))/(dx(ii-1))+Sigma_a(ii) + 1/dx(ii-1)*(1-alpha)/(1+alpha)
                c(ii) = 0.0d0
                d(ii) = S0(ii)
            else
                a(ii) = -(Dif(ii-1)+Dif(ii))/(dx(ii-1)*(dx(ii)+dx(ii-1)))
                b(ii) = ((Dif(ii)+Dif(ii-1))/dx(ii-1)+(Dif(ii)+Dif(ii+1))/dx(ii))/(dx(ii)+dx(ii-1))+Sigma_a(ii) 
                c(ii) = -(Dif(ii)+Dif(ii+1))/(dx(ii)*(dx(ii)+dx(ii-1)))
                d(ii) = S0(ii)
        end if

        end do

        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

        ! Write all columns at once after the loop
        open(unit=991, file="flux.dat", status="replace", action="write")
        do ii = 1, N
            write(991,'(F10.5,2(1X,E15.8))') x(ii), phi(ii)
        end do
        close(991)

    end subroutine build_matrix_A_iterative_power

end module m_iterative_power_solver
