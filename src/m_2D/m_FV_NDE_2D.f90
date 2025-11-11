module diffusion_matrix
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    contains

    subroutine assemble_diffusion_matrix(nx, ny, dx, dy, D, Sigma_a, ia, ja, a)
    !---------------------------------------------------------------
    ! Assembles CSR matrix for 2D diffusion operator:
    !   -D*(d2/dx2 + d2/dy2) + Sigma_a
    ! Dirichlet boundary conditions (phi=0) on all edges.
    !
    ! Input:
    !   nx, ny     - grid points (including boundaries)
    !   dx, dy     - grid spacing
    !   D          - diffusion coefficient
    !   Sigma_a    - absorption coefficient
    !
    ! Output:
    !   ia(:)      - row pointer (CSR)
    !   ja(:)      - column indices (CSR)
    !   a(:)       - nonzero values (CSR)
    !
    ! Notes:
    !   Lexicographic ordering: k = i + (j-1)*nx
    !---------------------------------------------------------------
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: dx, dy, D, Sigma_a
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(dp), allocatable, intent(out) :: a(:)

    integer :: N, i, j, row, ptr, est_nnz
    real(dp) :: ax, ay, diag, offx, offy

    N = nx * ny
    est_nnz = 5 * N
    allocate(ia(N+1), ja(est_nnz), a(est_nnz))

    ax = D / dx**2
    ay = D / dy**2
    diag = 2.0_dp * (ax + ay) + Sigma_a
    offx = -ax
    offy = -ay

    ptr = 1
    ia(1) = 1
    do j = 1, ny
      do i = 1, nx
        row = i + (j-1)*nx

        if (i == 1 .or. i == nx .or. j == 1 .or. j == ny) then
          ! Dirichlet boundary node
          ja(ptr) = row
          a(ptr)  = 1.0_dp
          ptr = ptr + 1
          ia(row+1) = ptr
          cycle
        end if

        ! Center
        ja(ptr) = row;     a(ptr) = diag; ptr = ptr + 1
        ! Left
        ja(ptr) = row-1;   a(ptr) = offx; ptr = ptr + 1
        ! Right
        ja(ptr) = row+1;   a(ptr) = offx; ptr = ptr + 1
        ! Down
        ja(ptr) = row-nx;  a(ptr) = offy; ptr = ptr + 1
        ! Up
        ja(ptr) = row+nx;  a(ptr) = offy; ptr = ptr + 1

        ia(row+1) = ptr
      end do
    end do

    ! Finalize CSR structure
    ia(N+1) = ptr
    if (ptr-1 < size(a)) then
      a = a(1:ptr-1)
      ja = ja(1:ptr-1)
    end if
  end subroutine assemble_diffusion_matrix

end module diffusion_matrix
