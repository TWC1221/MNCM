module m_multigroup_1D_diffusion
implicit none
contains
    subroutine build_matrix_multigroup(N, alpha, G, phi, dx, Dif, A_scatter, S0)
        use m_constants
        use m_thomas_algorithm
        implicit none
        integer, intent(in) :: N, G
        real(8), intent(in) :: alpha
        real(8), intent(inout) :: phi(N), dx(N-1), Dif(N), A_scatter(G), S0(N)

        real(8), dimension(N) :: a, b, c, d, x
        integer :: ii

        x(1) = 0.0d0
        do ii = 2, N
            x(ii) =  sum(dx(1:ii-1))
        end do

        do ii = 1, N
            if (ii == 1) then
                a(ii) = 0.0d0
                b(ii) = ((Dif(ii)+Dif(ii+1))/(dx(ii)*dx(ii))) + 1/dx(ii)*(1-alpha)/(1+alpha) + sum(A_scatter) 
                c(ii) = -(Dif(ii)+Dif(ii+1))/(dx(ii)*(dx(ii)))
                d(ii) = S0(ii)
            elseif (ii == N) then
                a(ii) = -(Dif(ii-1)+Dif(ii))/(dx(ii-1)*dx(ii-1))
                b(ii) = ((Dif(ii)+Dif(ii-1))/dx(ii-1))/(dx(ii-1)) + 1/dx(ii-1)*(1-alpha)/(1+alpha) + sum(A_scatter)
                c(ii) = 0.0d0
                d(ii) = S0(ii)
            else
                a(ii) = -(Dif(ii-1)+Dif(ii))/(dx(ii-1)*(dx(ii)+dx(ii-1)))
                b(ii) = ((Dif(ii)+Dif(ii-1))/dx(ii-1)+(Dif(ii)+Dif(ii+1))/dx(ii))/(dx(ii)+dx(ii-1)) + sum(A_scatter)
                c(ii) = -(Dif(ii)+Dif(ii+1))/(dx(ii)*(dx(ii)+dx(ii-1)))
                d(ii) = S0(ii)
            end if
        end do
        ! Solve system
        call thomas_algorithm(a(2:N), b, c(1:N-1), d, phi, N)

    end subroutine

end module m_multigroup_1D_diffusion