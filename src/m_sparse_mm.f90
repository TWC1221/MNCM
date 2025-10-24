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

    function forward_solver(L,r0) result(y0)
        implicit none
        real(8), intent(in) :: L(:,:), r0(:)
        real(8) :: y0(SIZE(r0))
        integer :: ii

        y0(1) = r0(1) / L(1,1)
        do ii = 2, size(r0)
            y0(ii) = (r0(ii) - sum(L(ii,1:ii-1) * y0(1:ii-1))) / L(ii,ii)
        end do
        return
    end function

    function backward_solver(L,y0) result(d0)
        implicit none
        real(8), intent(in) :: L(:,:), y0(:)
        real(8) :: d0(SIZE(y0))
        integer :: ii, n

        n = SIZE(y0)
        d0(n) = y0(n)/L(n,n)
        do ii = n-1, 1, -1
            d0(ii) = (y0(ii) - sum(L(ii+1:n,ii) * d0(ii+1:n))) / L(ii,ii)
        end do
        return
    end function

    subroutine CSR_rebuild(AA,JA,IA)
        real(8), intent(in) :: AA(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: A(SIZE(IA)-1,SIZE(IA)-1)
        integer :: ii, jj, k1, k2

        do ii = 1, SIZE(IA)-1
            k1 = IA(ii) !Start column index of ith row
            k2 = IA(ii+1) - 1 !End column index of ith row
            do jj = k1,k2
                A(ii,JA(jj)) = (AA(jj))
            end do
        end do
        do ii = 1, 5
            print '(5F10.5)', (A(ii, jj), jj = 1, 5)
        end do
    end subroutine CSR_rebuild

    function Incomplete_Cholesky_CRS(AA,JA,IA, r0) result(d0)
        real(8), intent(in) :: AA(:), r0(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: L(SIZE(IA)-1,SIZE(IA)-1), s, y0(SIZE(IA)-1), d0(SIZE(IA)-1)
        integer :: ii, jj, k, k1, k2

        L = 0
        do ii = 1, SIZE(IA)-1
            k1 = IA(ii) !Start column index of ith row
            k2 = IA(ii+1) - 1 !End column index of ith row
            do jj = k1,k2
                if (JA(jj) <= ii) then !Enforce Lower Triangular Matrix Formulation
                    L(ii,JA(jj)) = (AA(jj))
                end if
            end do
        end do

        do ii = 1, SIZE(IA)-1
            do jj = 1,ii
                if (L(ii,jj) /= 0.0d0) then
                    s = L(ii,jj)
                    do k =1, jj-1
                        s = s - L(ii,k)*L(jj,k)
                    end do
                    if (ii == jj) then
                        L(ii,jj) = sqrt(s)
                    else
                        L(ii,jj) = s/L(jj,jj)
                    end if
                end if
            end do
        end do
        y0 = forward_solver(L,r0)
        d0 = backward_solver(L,y0)
    end function

    function ILU(AA,JA,IA, r0) result(d0)
        real(8), intent(in) :: AA(:), r0(:)
        integer, intent(in) :: JA(:), IA(:)
        integer :: n, ii, jj, kk, k1, k2, jj_pos
        real(8), allocatable :: L(:,:), U(:,:), y0(:), s, d0(:)

        n = SIZE(IA)-1
        allocate(L(n,n))
        allocate(U(n,n))
        allocate(y0(n))

        ! Initialize L and U
        L = 0.0d0
        U = 0.0d0

        ! Construct ILU(0)
        do ii = 1,n
            k1 = IA(ii)
            k2 = IA(ii+1) - 1

            do jj = k1,k2
                jj_pos = JA(jj)
                s = AA(jj)

                ! Subtract previous contributions
                do kk = 1, jj_pos-1
                    s = s - L(ii,kk) * U(kk,jj_pos)
                end do

                if (jj_pos < ii) then
                    ! Lower triangle
                    if (U(jj_pos,jj_pos) == 0.0d0) then
                        print *, "Zero pivot at row ", jj_pos
                        stop
                    end if
                    L(ii,jj_pos) = s / U(jj_pos,jj_pos)
                else
                    ! Upper triangle (including diagonal)
                    U(ii,jj_pos) = s
                end if
            end do

            ! Set unit diagonal for L
            L(ii,ii) = 1.0d0
        end do

        ! Solve Ly = r0 (forward substitution)
        y0 = forward_solver(L, r0)

        ! Solve Ux = y0 (backward substitution)
        d0 = backward_solver(U, y0)

        ! Deallocate temporary arrays
        deallocate(L,U,y0)
    end function


    function Jacobi_CRS(AA,JA,IA, r0) result(d0)
        real(8), intent(in) :: AA(:), r0(:)
        integer, intent(in) :: JA(:), IA(:)
        real(8) :: J(SIZE(IA)-1), d0(SIZE(IA)-1),M_minus(SIZE(IA)-1)
        integer :: ii, jj, k1, k2

        do ii = 1, SIZE(IA)-1
            k1 = IA(ii) !Start column index of ith row
            k2 = IA(ii+1) - 1 !End column index of ith row
            do jj = k1,k2
                if (JA(jj) == ii) then !Enforce Lower Triangular Matrix Formulation
                    J(ii) = (AA(jj))
                end if
            end do
        end do

        do jj = 1,SIZE(IA)-1
            M_minus(jj) = 1/J(jj)
        end do
        d0 = M_minus*r0
    end function

    subroutine CG_solver
        real(8) :: AA(13), b(5), alpha, beta, p(5), r(5), x(5), p0(5), r0(5), x0(5)
        integer :: JA(13), IA(6), ii

        AA = (/4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0/)
        JA = (/1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5/)
        IA = (/1, 3, 6, 9, 12, 14/)

        b = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0/)
        x0 = ([3,6,9,12,15])

        r0 = b - CSR_dot_product(AA,JA,IA,x0)
        p0 = r0
        
        do ii = 1,5

            alpha = dot_product(r0,r0)/dot_product(p0, CSR_dot_product(AA,JA,IA,p0))
        
            x = x0 + alpha*p0

            r = r0 - alpha * CSR_dot_product(AA,JA,IA,p0)

            beta = dot_product(r,r)/dot_product(r0,r0)

            p = r + beta*p0

            r0 = r      
            p0 = p
            x0 = x
        end do

        print '(5F10.5)', "Ax =", CSR_dot_product(AA, JA, IA, x)
        print '(5F10.5)', "x =", x


    end subroutine CG_solver

    subroutine PCG_solver
        real(8) :: AA(13), b(5), alpha, beta, d(5), r(5), x(5), d0(5), r0(5), x0(5)
        integer :: JA(13), IA(6), ii

        AA = (/4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0/)
        JA = (/1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5/)
        IA = (/1, 3, 6, 9, 12, 14/)

        b = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0/)
        x0 = ([3,6,9,12,15])

        r0 = b - CSR_dot_product(AA,JA,IA,x0)
        d0 = Jacobi_CRS(AA,JA,IA,r0) !Incomplete_Cholesky_CRS(AA,JA,IA, r0) !ILU(AA,JA,IA, r0)
        
        do ii = 1,5

            alpha = dot_product(r0,d0)/dot_product(d0, CSR_dot_product(AA,JA,IA,d0))
        
            x = x0 + alpha*d0

            r = r0 - alpha * CSR_dot_product(AA,JA,IA,d0)

            beta = dot_product(r,d)/dot_product(r0,d0)

            !d = ILU(AA,JA,IA, r) + beta*d0
            
            !d = Jacobi_CRS(AA,JA,IA, r) + beta*d0

            !d = Incomplete_Cholesky_CRS(AA,JA,IA, r) + beta*d0          

            r0 = r      
            d0 = d
            x0 = x
        end do

        print*,CSR_dot_product(AA,JA,IA,x)

    end subroutine PCG_solver

end module m_sparse_mm
