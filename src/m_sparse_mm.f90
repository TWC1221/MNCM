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

    subroutine CG_solver
        real(8) :: AA(13), b(5), alpha, beta, p(5), r(5), x(5), p0(5), r0(5), x0(5)
        integer :: JA(13), IA(6), ii

        AA = (/4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0, 1.0d0, 1.0d0, 4.0d0/)
        JA = (/1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5/)
        IA = (/1, 3, 6, 9, 12, 14/)

        b = (/1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0/)
        x0 = ([3,6,9,12,15])

        call CSR_rebuild(AA,JA,IA)

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

        print*,CSR_dot_product(AA,JA,IA,x)

    end subroutine CG_solver

end module m_sparse_mm
