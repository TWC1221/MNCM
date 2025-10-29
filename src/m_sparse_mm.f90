!------------------------------------------------------------------------!
!! Purpose:                                                             -!
!  Carry out sparse matrix multiplicatiom                               -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 22/10/2025    T. Charlton      Original code                          -!
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

    subroutine CSR_rebuild(AA,JA,IA)
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: A(SIZE(IA)-1,SIZE(IA)-1), L(SIZE(IA)-1,SIZE(IA)-1)
        integer :: ii, jj, k1, k2

        do ii = 1, SIZE(IA)-1
            k1 = IA(ii) !Start column index of ith row
            k2 = IA(ii+1) - 1 !End column index of ith row
            do jj = k1,k2
                A(ii,JA(jj)) = (AA(jj))
            end do
        end do

        L = transpose(A)*A
        do ii = 1, 5
            print '(5F10.5)', (L(ii, jj), jj = 1, 5)
        end do
    end subroutine CSR_rebuild

  !-------------------------
  ! Incomplete Cholesky (lower, CSR) where each row's diagonal stored last
  ! AA,JA,IA define original A in CSR.
  ! Returns d0 = (L^T \ (L \ r0)) i.e. apply preconditioner to r0 (M^{-1} r0).
  !-------------------------
    function Incomplete_Cholesky_CSR_2(AA,JA,IA,r0) result(d0)
        real(8), intent(in) :: AA(:), r0(:)
        integer, intent(in) :: JA(:), IA(:)

        real(8), allocatable :: L_AA(:), d0(:), y0(:)
        integer, allocatable :: L_JA(:), L_IA(:), L_D0(:)
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
                    L_D0(ii) = sqrt(s)
                    L_AA(jj) = L_d0(ii)
                else
                    ! off-diagonal
                    L_AA(jj) = s / L_D0(jcol)
                end if
            end do
        end do

        ! ! Print real array L_AA
        ! print '(13F5.2)', L_AA

        ! ! Print integer arrays L_JA and L_IA
        ! print '(13I5)', L_JA
        ! print '(6I5)', L_IA

        !call CSR_rebuild(L_AA,L_JA,L_IA)
        allocate(y0(n))
        call forward_solver_CSR(L_AA,L_JA,L_IA,r0,y0)
        
        allocate(d0(n))
        call backward_solver_CSR(L_AA,L_JA,L_IA,y0,d0)
        print*,d0
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

    subroutine backward_solver_CSR(L_AA,L_JA,L_IA,y0,d0)
        implicit none
        real(8), intent(in) :: L_AA(:), y0(:)
        integer, intent(in) :: L_JA(:), L_IA(:)

        real(8) :: d0(size(L_IA)-1)

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

    subroutine PCG_solver
        real(8) :: AA(13), b(5), alpha, beta, x(5), r(5), d(5), z(5), q(5)
        integer :: JA(13), IA(6), ii, max_iter

        AA = (/4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0/)
        JA = (/1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5/)
        IA = (/1, 3, 6, 9, 12, 14/)

            ! RHS and initial guess
        b = (/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)
        x = 0.0d0

        ! Initial residual and preconditioned residual
        r = b - CSR_dot_product(AA, JA, IA, x)
        z = Incomplete_Cholesky_CSR_2(AA, JA, IA, r)
        d = z

        max_iter = 5
        do ii = 1, max_iter
            q = CSR_dot_product(AA, JA, IA, d)
            alpha = dot_product(r, z)/dot_product(d, q)
            x = x + alpha*d
            r = r - alpha*q
            z = Incomplete_Cholesky_CSR_2(AA, JA, IA, r)
            beta = dot_product(r, z)/dot_product(r - alpha*q, z - alpha*d) ! optional check
            d = z + beta*d
        end do

        print*, "Solution:", x
        print*, "aX = ", CSR_dot_product(AA, JA, IA, x)

    end subroutine PCG_solver

    ! function ILU(AA,JA,IA, r0) result(d0)
    !     real(8), intent(in) :: AA(:), r0(:)
    !     integer, intent(in) :: JA(:), IA(:)
    !     integer :: n, ii, jj, kk, k1, k2, jj_pos
    !     real(8), allocatable :: L(:,:), U(:,:), y0(:), s, d0(:)

    !     n = SIZE(IA)-1
    !     allocate(L(n,n))
    !     allocate(U(n,n))
    !     allocate(y0(n))

    !     ! Initialize L and U
    !     L = 0.0d0
    !     U = 0.0d0

    !     ! Construct ILU(0)
    !     do ii = 1,n
    !         k1 = IA(ii)
    !         k2 = IA(ii+1) - 1

    !         do jj = k1,k2
    !             jj_pos = JA(jj)
    !             s = AA(jj)

    !             ! Subtract previous contributions
    !             do kk = 1, jj_pos-1
    !                 s = s - L(ii,kk) * U(kk,jj_pos)
    !             end do

    !             if (jj_pos < ii) then
    !                 ! Lower triangle
    !                 if (U(jj_pos,jj_pos) == 0.0d0) then
    !                     print *, "Zero pivot at row ", jj_pos
    !                     stop
    !                 end if
    !                 L(ii,jj_pos) = s / U(jj_pos,jj_pos)
    !             else
    !                 ! Upper triangle (including diagonal)
    !                 U(ii,jj_pos) = s
    !             end if
    !         end do

    !         ! Set unit diagonal for L
    !         L(ii,ii) = 1.0d0
    !     end do

    !     ! Solve Ly = r0 (forward substitution)
    !     y0 = forward_solver(L, r0)

    !     ! Solve Ux = y0 (backward substitution)
    !     d0 = backward_solver(U, y0)

    !     ! Deallocate temporary arrays
    !     deallocate(L,U,y0)
    ! end function

    ! function Jacobi_CRS(AA,JA,IA, r0) result(d0)
    !     real(8), intent(in) :: AA(:), r0(:)
    !     integer, intent(in) :: JA(:), IA(:)
    !     real(8) :: J(SIZE(IA)-1), d0(SIZE(IA)-1),M_minus(SIZE(IA)-1)
    !     integer :: ii, jj, k1, k2

    !     do ii = 1, SIZE(IA)-1
    !         k1 = IA(ii) !Start column index of ith row
    !         k2 = IA(ii+1) - 1 !End column index of ith row
    !         do jj = k1,k2
    !             if (JA(jj) == ii) then !Enforce Lower Triangular Matrix Formulation
    !                 J(ii) = (AA(jj))
    !             end if
    !         end do
    !     end do

    !     do jj = 1,SIZE(IA)-1
    !         M_minus(jj) = 1/J(jj)
    !     end do
    !     d0 = M_minus*r0
    ! end function

    ! subroutine CG_solver
    !     real(8) :: AA(13), b(5), alpha, beta, p(5), r(5), x(5), p0(5), r0(5), x0(5)
    !     integer :: JA(13), IA(6), ii

    !     AA = (/4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0/)
    !     JA = (/1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5/)
    !     IA = (/1, 3, 6, 9, 12, 14/)

    !     b = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0/)
    !     x0 = ([3,6,9,12,15])

    !     r0 = b - CSR_dot_product(AA,JA,IA,x0)
    !     p0 = r0
        
    !     do ii = 1,5

    !         alpha = dot_product(r0,r0)/dot_product(p0, CSR_dot_product(AA,JA,IA,p0))
        
    !         x = x0 + alpha*p0

    !         r = r0 - alpha * CSR_dot_product(AA,JA,IA,p0)

    !         beta = dot_product(r,r)/dot_product(r0,r0)

    !         p = r + beta*p0

    !         r0 = r      
    !         p0 = p
    !         x0 = x
    !     end do

    !     print '(A, 5F10.5)', "Ax =", CSR_dot_product(AA, JA, IA, x)


    ! end subroutine CG_solver

end module m_sparse_mm
