module m_driver
!------------------------------------------------------------------------!
!  Purpose:                                                             -!
!  Contains the subroutines required to drive solvers                   -!
!------------------------------------------------------------------------!
!  Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 09/10/2025    T. Charlton      Original code                          -!
!------------------------------------------------------------------------! 
    use m_constants
    use m_MMS_solver
    use m_vacuum_BC_solver
    contains
    subroutine MMS_driver()
        implicit none
        integer :: N, ios, ii
        integer :: index = 0
        real(8), parameter :: alpha = 0.0d0  ! 0.0 = vacuum, 1.0 = reflective
        real(8), allocatable :: a(:), b(:), c(:)    ! tridiagonal matrix
        real(8), allocatable :: phi(:)             ! eigenvector
        real, allocatable :: x(:), y(:)
        real(8) :: L2Err, sum_x, sum_y, sum_xy, sum_xx, slope
        
        !-----------------------------
        ! Discretization Unit / Nodal Iteration Unit
        !-----------------------------
        ii=0
        open(unit=993, file="error.dat", status="replace", action="write")
        do N = 5, 5010, 1
            allocate(a(N), b(N), c(N), phi(N))
            call build_matrix_A_MMS(a, b, c, N, phi, L2Err)
            write(993,'(F10.5,1(1X,E15.8))') log10(real(N)), log10(L2Err)
            deallocate(a, b, c, phi)
            ii = ii + 1
        end do
        close(993)

        allocate(x(ii),y(ii))
        open(unit=993, file="error.dat", status='old', action='read', iostat=ios)
        do index = 1, ii 
            read(993, *, iostat=ios) x(index), y(index)
            sum_x = sum_x + x(index)
            sum_y = sum_y + y(index)
            sum_xy = sum_xy + x(index)*y(index)
            sum_xx = sum_xx + x(index)*x(index)
        end do
        close(993)

        slope = (ii*sum_xy - sum_x*sum_y) / (ii*sum_xx - sum_x*sum_x)

        print *, 'Log-Log Linear Fit: dy/dx = ', slope

        !-----------------------------
        ! Plot Error vs N
        !-----------------------------
        open(unit=994, file="plot_error.gp", status="replace", action="write")
        write(994,*) "set term wxt 0 title 'Error vs N'"
        write(994,*) "set title 'Error vs N'"
        write(994,*) "set xlabel 'log10(N)'"
        write(994,*) "set ylabel 'log10(Error)'"
        write(994,*) "plot 'error.dat' using 1:2 with lines title 'Error'"
        close(994)
        call system("gnuplot -persist plot_error.gp")

        !-----------------------------
        ! Plot Flux vs x
        !-----------------------------
        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x for different alpha'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical', \"
        write(995,*) "     'flux.dat' using 1:3 with lines title 'exact'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")

    end subroutine MMS_driver

    subroutine vacuum_BC_driver()
        implicit none
        integer :: N, ios, ii, index
        real(8), parameter :: alpha = 0.0d0  ! 0.0 = vacuum, 1.0 = reflective
        real(8), allocatable :: a(:), b(:), c(:)    ! tridiagonal matrix
        real(8), allocatable :: phi(:), phi_analytical(:)             ! eigenvector
        real, allocatable :: x(:), y(:)
        real(8) :: L2Err, sum_x, sum_y, sum_xy, sum_xx, slope
        
        !-----------------------------
        ! Discretization Unit / Nodal Iteration Unit
        !-----------------------------
        open(unit=993, file="error.dat", status="replace", action="write")
        ii=0
        do N = 5, 40, 1
            allocate(a(N), b(N), c(N), phi(N), phi_analytical(N))
            call build_matrix_A_vacuum(a, b, c, N, alpha, phi, phi_analytical, L2Err)
            write(993,'(F10.5,1(1X,E15.8))') log10(real(N)), log10(L2Err)
            deallocate(a, b, c, phi, phi_analytical)
            ii = ii + 1
        end do
        close(993)

        allocate(x(ii),y(ii))
        open(unit=993, file="error.dat", status='old', action='read', iostat=ios)
        do index = 1, ii 
            read(993, *, iostat=ios) x(index), y(index)
            sum_x = sum_x + x(index)
            sum_y = sum_y + y(index)
            sum_xy = sum_xy + x(index)*y(index)
            sum_xx = sum_xx + x(index)*x(index)
        end do
        close(993)

        slope = (ii*sum_xy - sum_x*sum_y) / (ii*sum_xx - sum_x*sum_x)

        print *, 'Log-Log Linear Fit: dy/dx = ', slope

        !-----------------------------
        ! Plot Error vs N
        !-----------------------------
        open(unit=994, file="plot_error.gp", status="replace", action="write")
        write(994,*) "set term wxt 0 title 'Error vs N'"
        write(994,*) "set title 'Error vs N'"
        write(994,*) "set xlabel 'log10(N)'"
        write(994,*) "set ylabel 'log10(Error)'"
        write(994,*) "plot 'error.dat' using 1:2 with lines title 'Error'"
        close(994)
        call system("gnuplot -persist plot_error.gp")

        !-----------------------------
        ! Plot Flux vs x
        !-----------------------------
        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical', \"
        write(995,*) "     'flux.dat' using 1:3 with lines title 'exact'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")

    end subroutine Vacuum_BC_driver

    subroutine heterogeneous_driver_m1()
        implicit none
        integer :: N = 1000, N_interface = 500
        real(8) :: alpha = 0.0
        real(8), allocatable :: a(:), b(:), c(:), phi(:)

        allocate(a(N), b(N), c(N), phi(N))

        call build_matrix_A_heterogeneous_m1(a, b, c, N, N_interface, alpha, phi)

        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x for different alpha'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")
    end subroutine heterogeneous_driver_m1

     subroutine heterogeneous_driver_m2()
        implicit none
        integer :: N = 1000, N_interface = 500
        real(8) :: alpha = 0.0, x_interface = 1
        real(8), allocatable :: a(:), b(:), c(:), phi(:)

        allocate(a(N), b(N), c(N), phi(N))

        call build_matrix_A_heterogeneous_m2(a, b, c, N, N_interface, x_interface, alpha, phi)

        open(unit=995, file="plot_flux.gp", status="replace", action="write")
        write(995,*) "set term wxt 1 title 'Flux vs x'"
        write(995,*) "set title 'Flux vs x for different alpha'"
        write(995,*) "set xlabel 'x'"
        write(995,*) "set ylabel 'phi(x)'"
        write(995,*) "plot 'flux.dat' using 1:2 with lines title 'numerical', 'flux.dat' using 1:3 with lines title 'analytical'"
        close(995)
        call system("gnuplot -persist plot_flux.gp")
    end subroutine heterogeneous_driver_m2

end module m_driver

