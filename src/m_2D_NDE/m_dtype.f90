!--------------------------------------------------------------------------------------------!
!  Purpose:                                                                                 -!
!  Contains inner/outer iteration on multigroup diffusion eigenvalue problem                -!
!           no-upscatter fixed fission source solver                                        -! 
!--------------------------------------------------------------------------------------------!
!  Record of revisions:                                                                     -!
!   Date       Programmer     Description of change                                         -!
!   ====       ==========     =====================                                         -!
! 27/11/2025    T. W Charlton     Original code                                             -!
!--------------------------------------------------------------------------------------------! 
module multigroup_multidimensional
    USE m_PCG_solver
    use m_diffusion_matrix
    implicit none
    contains

 subroutine multigroup_diffusion_iter_3d(PCG_mode, max_iter)
    use CSR_types, only: CSRMatrix
    implicit none

    type(CSRMatrix) :: A_all(3)

    integer, intent(in) :: PCG_mode, max_iter
    integer :: N, G, nx, ny, nz
    integer :: ii, jj, gg

    real(8) :: X_domain, Y_domain, Z_domain
    real(8) :: dx, dy, dz, alpha

    real(8), allocatable :: phi(:,:), phi_prime(:,:), phi_ptr(:)
    real(8), allocatable :: D(:), Sigma_t(:), Sigma_a(:), Sigma_r(:), nu_Sigma_f(:)
    real(8), allocatable :: Sigma_s(:,:), Sigma_s_upscatter(:,:), Sigma_s_L(:,:), Sigma_s_U(:,:)
    real(8), allocatable :: Sigma_s_L_sum(:,:), Sigma_s_U_sum(:,:)
    real(8), allocatable :: S_f(:), S_f_iter(:), K_eff(:), x(:), chi(:)

    !------------------------------
    ! Problem size and groups
    !------------------------------
    G = 3
    nx = 4
    ny = 4
    nz = 4
    N = nx * ny * nz

    X_domain = 10.0d0
    Y_domain = 10.0d0
    Z_domain = 10.0d0
    dx = X_domain / nx
    dy = Y_domain / ny
    dz = Z_domain / nz
    alpha = 0.0d0

    !------------------------------
    ! Allocate arrays
    !------------------------------
    allocate(phi(G,N), phi_prime(G,N), phi_ptr(N))
    allocate(D(G), Sigma_t(G), Sigma_a(G), Sigma_r(G), nu_Sigma_f(G))
    allocate(Sigma_s(G,G), Sigma_s_upscatter(G,G), Sigma_s_L(G,G), Sigma_s_U(G,G))
    allocate(Sigma_s_L_sum(G,N), Sigma_s_U_sum(G,N))
    allocate(S_f(N), S_f_iter(N), K_eff(1000), x(N), chi(G))

    !------------------------------
    ! Initialize physics
    !------------------------------
    chi = [0.9d0, 0.1d0, 0.0d0]
    Sigma_a = [0.015d0, 0.04d0, 0.12d0]
    nu_Sigma_f = [0.02d0, 0.1d0, 0.35d0]
    S_f = 1.0d0
    K_eff(1) = 1.0d0
    phi_prime = 0.0d0

    ! Test scattering matrix
    Sigma_s_upscatter = reshape([0.20d0, 0.10d0, 0.17d0, &
                                 0.05d0, 0.25d0, 0.10d0, &
                                 0.10d0, 0.07d0, 0.30d0], shape=[3,3])
    Sigma_s_L(:,:) = 0.0d0
    Sigma_s_L(2:3,1:2) = Sigma_s_upscatter(2:3,1:2)
    Sigma_s_L(2,2) = 0.0d0
    Sigma_s_U(:,:) = 0.0d0
    Sigma_s_U(1:2,2:3) = Sigma_s_upscatter(1:2,2:3)
    Sigma_s_U(2,2) = 0.0d0

    !------------------------------
    ! Assemble diffusion matrices
    !------------------------------
    do gg = 1,G
        Sigma_t(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G))
        Sigma_r(gg) = Sigma_a(gg) + sum(Sigma_s_upscatter(gg,1:G)) - Sigma_s_upscatter(gg,gg)
        D(gg) = 1.0d0 / (3.0d0 * Sigma_t(gg))
        call assemble_3D_diffusion_matrix_g(nx, ny, nz, N, G, gg, X_domain, Y_domain, Z_domain, D(gg), Sigma_r(gg), alpha, A_all(gg), dx, dy, dz)
    end do

    !------------------------------
    ! Iterative solver
    !------------------------------
    do ii = 2,max_iter
        S_f_iter(1:N) = S_f(1:N)
        do
            do gg = 1,G
                do jj = 1,G
                    Sigma_s_L_sum(jj,1:N) = Sigma_s_L(jj,gg) * phi_prime(jj,1:N)
                    Sigma_s_U_sum(jj,1:N) = Sigma_s_U(jj,gg) * phi(jj,1:N)
                end do
                S_f(1:N) = chi(gg)/K_eff(ii-1) * S_f_iter(1:N) + sum(Sigma_s_L_sum, dim=1) + sum(Sigma_s_U_sum, dim=1)

                phi_ptr = phi(gg,:)
                call PCG_algorithm(A_all(gg)%AA, A_all(gg)%JA, A_all(gg)%IA, phi_ptr, S_f, PCG_mode, max_iter)
                phi(gg,:) = phi_ptr

            end do

            if (maxval(abs(phi(3,:) - phi_prime(3,:))) < 1.0d-7 .and. &
                maxval(abs(phi(2,:) - phi_prime(2,:))) < 1.0d-7) exit

            phi_prime(3,:) = phi(3,:)
            phi_prime(2,:) = phi(2,:)
        end do

        S_f(1:N) = matmul(transpose(phi), nu_Sigma_f)
        K_eff(ii) = K_eff(ii-1) * sum(S_f*dx*dy*dz) / sum(S_f_iter*dx*dy*dz)
        print *, "Iteration =", ii-1, "Keff =", K_eff(ii-1)

        if (abs((K_eff(ii) - K_eff(ii-1)) / K_eff(ii-1)) < 1.0d-6 .and. &
            maxval(abs(S_f - S_f_iter)) < 1.0d-6) exit
    end do

end subroutine multigroup_diffusion_iter_3d


end module multigroup_multidimensional

! module vector_module
!     implicit none
!     type :: vector
!         real :: x
!         real :: y
!     end type vector
!     contains
!     type (vector) function vector_add ( v1, v2 )
!         implicit none
!         type (vector), intent(in) :: v1
!         type (vector), intent(in) :: v2

!         vector_add%x = v1%x + v2%x
!         vector_add%y = v1%y + v2%y
!     end function vector_add
! end module vector_module

! PROGRAM test_vectors
!     use vector_module
!     implicit none

!     type(vector) :: v1
!     type(vector) :: v2
!     write (*,*) 'Enter the first vector (x1,y1):'
!     read (*,*) v1%x, v1%y

!     write (*,*) 'Enter the second vector (x2,y2):'
!     read (*,*) v2%x, v2%y

!     write (*,1000) vector_add(v1,v2)
!     1000 FORMAT('The sum of the points (v1, v2) is (',F8.2,','F8.2,')')
! END PROGRAM



