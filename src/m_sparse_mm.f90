!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Carry out sparse matrix multiplication                               -!  
!  Carry out Conjugate Gradient method for solving linear systems       -!  
!  Generate Preconditioning matrices M for use in Preconditioned CG     -!  
!  using Incomplete Cholesky factorization (IC(0)) in CSR format        -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 24/10/2025    T. Charlton      Original code                          -!      
! 28/10/2025    T. Charlton      Added true CSR Cholesky preconditioner -!
! 29/10/2025    T. Charlton      Integrated CG with 1D multigroup solver-!
!------------------------------------------------------------------------! 

module m_Sparse_mm
implicit none
contains 
    function CSR_dot_product(AA,JA,IA,x) result(Y)
        implicit none
        real(8), intent(in) :: AA(:), x(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: Y(SIZE(x))
        integer :: ii, k1, k2

        do ii = 1,SIZE(x)
            k1 = IA(ii)
            k2 = IA(ii+1)-1
            Y(ii) = dot_product(AA(k1:k2),x(JA(k1:k2)))
        end do
        return
    end function

  !-------------------------
  ! Incomplete Cholesky (lower, CSR) Decomposition of CSR matrix A
  ! AA, JA, IA define original A in CSR. 
  ! L_AA, L_JA, L_IA define the lower triangular L. L transformation into Cholesky preconditioning matrix M done in place
  ! Returns L_AA, L_JA, L_IA. The CSR format of ICD(0) matrix M
  !-------------------------
    subroutine Cholesky_CSR(AA, JA, IA, L_AA, L_JA, L_IA)
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)

        real(8), allocatable, intent(out) :: L_AA(:)
        integer, allocatable, intent(out) :: L_JA(:), L_IA(:)

        real(8), allocatable :: L_D0(:)
        real(8) :: s
        integer :: ii, jj, kk, k1, k2, nnz, jcol, kcol, kk1, kk2, n
        
        n = size(IA) - 1
        nnz = 0

        do ii = 1, n
            k1 = IA(ii)
            k2 = IA(ii+1) - 1
            do jj = k1,k2
                if (JA(jj) <= ii) nnz = nnz + 1
            end do
        end do

        allocate(L_AA(nnz))
        allocate(L_JA(nnz))
        allocate(L_IA(n+1))
        allocate(L_D0(n))

        kk = 0
        do ii = 1, n
            L_IA(ii) = kk + 1
            k1 = IA(ii)
            k2 = IA(ii+1) - 1
            do jj = k1, k2
                if (JA(jj) <= ii) then
                    kk = kk + 1
                    L_AA(kk) = AA(jj)
                    L_JA(kk) = JA(jj)
                end if
            end do
        end do
        L_IA(n+1) = nnz + 1  ! last element points past the end

        do ii = 1, n
            k1 = L_IA(ii)
            k2 = L_IA(ii+1)-1
            do jj = k1, k2
                jcol = L_JA(jj)
                s = L_AA(jj)
                do kk1 = k1, jj-1
                    kcol = L_JA(kk1)
                    do kk2 = L_IA(kcol), L_IA(kcol+1)-1
                        if (L_JA(kk2) == jcol) then
                            s = s - L_AA(kk1)*L_AA(kk2)
                            exit
                        end if
                    end do
                end do

                if (ii == jcol) then
                    ! diagonal
                    s = max(s, 1.0d-20)
                    L_D0(ii) = sqrt(s)
                    L_AA(jj) = L_D0(ii)
                else
                    ! off-diagonal
                    L_AA(jj) = s / L_D0(jcol)
                end if
            end do
        end do
    end subroutine

    subroutine Jacobi_CSR(AA, JA, IA, diag)
        implicit none
        ! Inputs
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)
        ! Output
        real(8), intent(out), allocatable :: diag(:)

        integer :: n, ii, jj

        n = size(IA) - 1
        allocate(diag(n))

        do ii = 1, n
            ! Loop over row ii
            do jj = IA(ii), IA(ii+1)-1
                if (JA(jj) == ii) then
                    diag(ii) = AA(jj)
                    exit
                end if
            end do
        end do
    end subroutine

  !-------------------------
  ! Function to compute d = M^-1 r using preconditioning matrix M from IC(0) factorization
  ! L_AA, L_JA, L_IA (CSR format) define preconditioning matrix M computed by Cholesky Decomposition
  ! Returns search vector d
  !-------------------------
    function PCG_CSR(L_AA, L_JA, L_IA, r) result(d)
        real(8), intent(in) :: L_AA(:), r(:)
        integer, intent(in) :: L_JA(:), L_IA(:)
        real(8), allocatable :: d(:), y(:)

        allocate(y(size(L_IA)-1))
        call forward_solver_CSR(L_AA,L_JA,L_IA,r,y)

        allocate(d(size(L_IA)-1))
        call backward_solver_CSR(L_AA,L_JA,L_IA,y,d)

    end function

    !-------------------------
    ! Forward solve: L y = r  (L stored in CSR; diagonal is last element per row)
    !-------------------------
    subroutine forward_solver_CSR(L_AA, L_JA, L_IA, r0, y0)
        implicit none
        real(8), intent(in) :: L_AA(:), r0(:)
        integer, intent(in) :: L_JA(:), L_IA(:)
        real(8), intent(out) :: y0(:)
        integer :: ii, k1, k2

        do ii = 1, size(L_IA) - 1
            k1 = L_IA(ii)
            k2 = L_IA(ii+1) - 1
            y0(ii) = (r0(ii) - dot_product(L_AA(k1:k2-1), y0(L_JA(k1:k2-1))))/L_AA(k2)
        end do
        return
    end subroutine

    !-------------------------
    ! Backward solve: L^T d = y  (L stored in CSR; diagonal is last element per row)
    !-------------------------
    subroutine backward_solver_CSR(L_AA,L_JA,L_IA,y0,d0)
        implicit none
        real(8), intent(in) :: L_AA(:), y0(:)
        integer, intent(in) :: L_JA(:), L_IA(:)

        real(8),intent(out) :: d0(size(L_IA)-1)

        integer :: ii, jj, k1, k2, n
        real(8) :: L_D0

        n = SIZE((L_IA))-1
        d0 = y0

        do ii = n, 1, -1
            k1 = L_IA(ii)
            k2 = L_IA(ii+1) - 1
            L_D0 = 0.0d0
            do jj = k1, k2
                if (L_JA(jj) == ii) then
                    L_D0 = L_AA(jj)
                    exit
                end if
            end do

            d0(ii) = d0(ii)/L_D0

            do jj = k1, k2
                if (L_JA(jj) < ii) then
                    d0(L_JA(jj)) = d0(L_JA(jj)) - L_AA(jj)*d0(ii)
                end if
            end do
        end do
        return
    end subroutine

    subroutine build_CSR_matrix_multigroup(N, alpha, G, phi_ptr, dx, Dif, A_scatter, S0, gg, Sigma_a)
        integer, intent(in) :: N, G, gg
        real(8), intent(in) :: Dif, dx(:), alpha, A_scatter(:,:), Sigma_a, S0(:)
        
        real(8), intent(out) :: phi_ptr(:)
        real(8), allocatable :: AA(:), b(:)
        integer, allocatable :: JA(:), IA(:)
        integer :: ii, nnz

        allocate(AA(3*N-2), JA(3*N-2), IA(N+1))

        do ii = 1,N
            if (ii == 1) then
                AA(1:2) = [((2*Dif)/(dx(ii)**2) + 1/dx(ii)*(1-alpha)/(1+alpha) + sum(A_scatter(gg,1:G)) - A_scatter(gg,gg) + Sigma_a), (-(2*Dif)/(dx(ii)**2))]
                JA(1:2) = [1, 2]
                IA(1) = 1
                nnz = 2
            elseif (ii == N) then
                AA(nnz+1:nnz+2) = [(-(2*Dif)/(dx(ii-1)**2)), ((2*Dif)/(dx(ii-1)**2) + 1/dx(ii-1)*(1-alpha)/(1+alpha) + sum(A_scatter(gg,1:G)) - A_scatter(gg,gg) + Sigma_a)]
                JA(nnz+1:nnz+2) = [N-1, N]
                IA(N) = nnz + 1
                nnz = nnz + 2
            else
                AA(nnz+1:nnz+3) = [(-(2*Dif)/(dx(ii-1)*(dx(ii)+dx(ii-1)))), ((2*Dif)/dx(ii-1)+(2*Dif)/dx(ii))/(dx(ii)+dx(ii-1)) + sum(A_scatter(gg,1:G)) - A_scatter(gg,gg) + Sigma_a, (-(2*Dif)/(dx(ii)*(dx(ii)+dx(ii-1))))]
                JA(nnz+1:nnz+3) = [ii-1, ii, ii+1]
                IA(ii) = nnz + 1
                nnz = nnz + 3
            end if
        end do
        IA(N+1) = nnz + 1 

        b = S0

        call PCG_algorithm(AA, JA, IA, phi_ptr, b)

    end subroutine build_CSR_matrix_multigroup

    subroutine PCG_algorithm(AA, JA, IA, x, b)
        implicit none
        real(8), intent(in) :: AA(:), b(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8), intent(out) :: x(:)
        real(8) :: alpha, beta
        integer :: ii, max_iter
        real(8), allocatable :: L_AA(:), diag(:), r(:), d(:), z(:), q(:)
        integer, allocatable :: L_IA(:), L_JA(:)

        call Cholesky_CSR(AA, JA, IA, L_AA, L_JA, L_IA)
        !call Jacobi_CSR(AA, JA, IA, diag)
        
        r = b - CSR_dot_product(AA, JA, IA, x)
        z = PCG_CSR(L_AA, L_JA, L_IA, r)! r/diag !
        d = z

        max_iter = 100
        do ii = 1, max_iter
            q = CSR_dot_product(AA, JA, IA, d)
            alpha = dot_product(r, z)/dot_product(d, q)
            x = x + alpha*d
            r = r - alpha*q
            z = PCG_CSR(L_AA, L_JA, L_IA, r) !r/diag 
            beta = dot_product(r, z)/dot_product(r - alpha*q, z - alpha*d) ! optional check
            d = z + beta*d
        end do

        !print*, "Solution:", x
        !print*, "aX = ", CSR_dot_product(AA, JA, IA, x)

    end subroutine

    subroutine PCG_solver_test
        implicit none

        integer :: n, ii
        real(8), allocatable :: AA(:), b(:), x(:)
        integer, allocatable :: JA(:), IA(:)

        n = 4
        allocate(AA(10), JA(10), IA(n+1), b(n), x(n))

        ! CSR representation of 4x4 SPD matrix
        AA = [4.0d0, -1.0d0, -1.0d0, 4.0d0, -1.0d0, -1.0d0, 4.0d0, -1.0d0, -1.0d0, 3.0d0]
        JA = [1,2,1,2,3,2,3,4,3,4]
        IA = [1,3,6,9,11]

        ! Right-hand side and initial guess
        b = [1.0d0, 2.0d0, 2.0d0, 1.0d0]
        x = 0.0d0

        ! Solve
        call PCG_algorithm(AA, JA, IA, x, b)

        ! Print solution
        print*, "PCG Solution:"
        do ii = 1, n
            print*, x(ii)
        end do
    end subroutine

end module m_sparse_mm
