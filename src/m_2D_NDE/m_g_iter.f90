module CSR_types
implicit none
    private
    public :: CSRMatrix

    type :: CSRMatrix
        real(8), allocatable :: AA(:)
        integer, allocatable :: JA(:)
        integer, allocatable :: IA(:)
        integer :: N
    end type CSRMatrix
end module CSR_types

module m_g_iter
    implicit none
    contains

    subroutine multigroup_diffusion_iter(mode, PCG_mode)
        use CSR_types, only: CSRMatrix
        implicit none
        type(CSRMatrix) :: A_all(3)
        integer, intent(in) :: mode, PCG_mode ! 1 = Fixed Source No Upscatter, 2 = Fixed Source Upscatter Iteration, 3 = Multigroup Iterative Fission No Upscatter, 4 = Multigroup Iterative Fission with Upscatter Iteration
        integer :: N = 1000, G = 3, gg, jj, ii
        real(8) :: alpha = 0.0, R_domain = 10
        real(8), allocatable :: phi(:,:), phi_prime(:,:), phi_ptr(:), Dif(:), Sigma_t(:), Sigma_a(:), nu_Sigma_f(:), Sigma_s(:,:), Sigma_sr(:,:), S_f(:), S_f_iter(:), K_eff(:), x(:), dx(:), dx_V(:), chi(:)
        real(8), allocatable :: Sigma_s_upscatter(:,:), Sigma_s_L(:,:), Sigma_s_U(:,:), Sigma_s_L_sum(:,:), Sigma_s_U_sum(:,:) ! Experimental

        allocate(phi(G,N), phi_prime(G,N), phi_ptr(N), Dif(G), Sigma_t(G), Sigma_a(G), nu_Sigma_f(G), Sigma_s(G,G), Sigma_sr(G,G), S_f(N), S_f_iter(N), K_eff(1000), x(N), dx(N-1), dx_V(N), chi(G))
        allocate(Sigma_s_upscatter(1:G,1:G), Sigma_s_L(1:G,1:G), Sigma_s_L_sum(1:G,1:N), Sigma_s_U(1:G,1:G), Sigma_s_U_sum(1:G,1:N))
        call system("pkill gnuplot")

        dx(1:N) = R_domain / real(N, 8) 

        chi = [0.9, 0.1, 0.0]
        Sigma_a = [0.015, 0.04, 0.12]
        nu_Sigma_f = [0.02, 0.1, 0.35] !0.35]

        S_f = 1 ; K_eff(1) = 1; 
        ! MULTIGROUP ITERATIVE FISSION SOLVER with Up-Scatter Iteration
        print *, "MULTIGROUP ITERATIVE FISSION SOLVER WITH UP-SCATTER ITERATION MODE"
        Sigma_s_upscatter = reshape([0.20d0, 0.10d0, 0.17d0, 0.05d0, 0.25d0, 0.10d0, 0.10d0, 0.07d0, 0.3d0], shape=[3,3]) ! TEST DATA (UNREALISTIC UPSCATTER)
        !Sigma_s_upscatter = reshape([0.30d0, 0.20d0, 0.085d0,0.02d0, 0.25d0, 0.09d0, 0.001d0, 0.009d0, 0.17d0], shape=[3,3]) ! REALISTIC DATA

        if (mode == 5) Sigma_s_upscatter = transpose(Sigma_s_upscatter)

        Sigma_s_L(:,:) = 0 ; Sigma_s_L(2:3,1:2) = Sigma_s_upscatter(2:3,1:2) ; Sigma_s_L(2,2) = 0 ;
        Sigma_s_U(:,:) = 0 ; Sigma_s_U(1:2,2:3) = Sigma_s_upscatter(1:2,2:3) ; Sigma_s_U(2,2) = 0 ;

        do gg = 1,G
            Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G))
            Dif(gg) = 1.0d0/(3.0d0*Sigma_t(gg))
            call build_CSR_matrix_multigroup(N, alpha, G, gg, dx, Dif(gg), Sigma_s_upscatter(1:G,1:G), Sigma_a(gg), A_all(gg))
        end do

        do ii = 2,1000
            S_f_iter(1:N) = S_f(1:N)
            do 
                do gg = 1,3
                    do jj = 1,3
                        Sigma_s_L_sum(jj,1:N) = Sigma_s_L(jj,gg) * phi_prime(jj,1:N)
                        Sigma_s_U_sum(jj,1:N) = Sigma_s_U(jj,gg) * phi(jj,1:N)
                    end do

                    if (mode == 4) S_f(1:N) = chi(gg)/K_eff(ii-1) * S_f_iter(1:N) + sum(Sigma_s_L_sum,dim=1) + sum(Sigma_s_U_sum,dim=1)
                    if (mode == 5) S_f(1:N) = nu_sigma_f(gg)/K_eff(ii-1) * S_f_iter(1:N) + sum(Sigma_s_L_sum,dim=1) + sum(Sigma_s_U_sum,dim=1)

                    phi_ptr = phi(gg,1:N)
                    call PCG_algorithm(A_all(gg)%AA, A_all(gg)%JA, A_all(gg)%IA, phi_ptr, S_f, PCG_mode, N, dx)
                    phi(gg,1:N) = phi_ptr
                end do

                !print'(10F10.5)',phi(1,N/2),phi(2,N/2),phi(3,N/2)
                if (maxval(abs(phi(3,1:N) - phi_prime(3,1:N))) < 1.0d-7 .and. maxval(abs(phi(2,1:N) - phi_prime(2,1:N))) < 1.0d-7) exit
                phi_prime(3,1:N) = phi(3,1:N)
                phi_prime(2,1:N) = phi(2,1:N) 
                exit ! Accelerates convergence by performing many more outer iterations and negleting inner iterations
            end do
            
            S_f(1:N) = matmul(transpose(phi),nu_Sigma_f)
            K_eff(ii) = K_eff(ii-1)*sum(S_f*dx_V)/sum(S_f_iter*dx_V)
            print*,"Iteration =", ii-1, "Keff =", K_eff(ii-1)
            if (abs((K_eff(ii) - K_eff(ii-1)) / K_eff(ii-1)) < 1.0d-6 .and. maxval(abs(S_f - S_f_iter)) < 1.0d-6) exit
        end do

        print*,"Iteration =", ii-1, "Keff =", K_eff(ii-1)

        x(1) = 0.0d0
        do jj = 2, N
            x(jj) =  sum(dx(1:jj-1))
        end do

        open(unit=991, file="flux.dat", status="replace", action="write")
        do jj = 1, N
            write(991,'(F10.5,3(1X,E15.8))') x(jj), phi(1,jj), phi(2,jj), phi(3,jj)
        end do
        close(991)

        open(unit=991, file="plot_flux.gp", status="replace", action="write")
        write(991,*) "set term qt 1 noraise title 'Multigroup Diffusion Flux Profile'"
        if (mode == 5) then
            write(991,*) "set title 'Adjoint Flux (Importance Function) vs x for three energy groups (g_{3} < g_{1})'"
        else
            write(991,*) "set title 'Flux vs x for three energy groups'"
        end if
        write(991,*) "set xlabel 'x'"
        if (mode == 5) then
            write(991,*) "set ylabel 'Adjoint Flux (x)'"
        else
            write(991,*) "set ylabel 'phi(x)'"
        end if
        write(991,*) "plot 'flux.dat' using 1:2 with lines title 'g=1',\"
        write(991,*) "     'flux.dat' using 1:3 with lines title 'g=2',\"
        write(991,*) "     'flux.dat' using 1:4 with lines title 'g=3'"
        close(991)
        call system("env QT_QPA_PLATFORM=xcb /usr/bin/gnuplot -persist plot_flux.gp")

    end subroutine multigroup_diffusion_iter
end module