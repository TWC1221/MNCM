module FV_NDE_2D
    contains
    subroutine mat()
        real(8), allocatable :: AA(:,:), dx(:), dy(:), D(:,:)
        real(8) :: x_domain=1.0, y_domain=1.0
        integer :: ii, jj, Nx=10, Ny=10
        
        allocate(AA(Nx,Ny),dx(Nx),dy(Ny),D(Nx,Ny))

        D(:,:) = 1

        dx(1:Nx) = x_domain/real(Nx)
        dy(1:Ny) = y_domain/real(Ny)

        do ii = 1,Nx
            do jj = 1, Ny
                if (ii == jj) then
                    AA(ii,jj) = 2*(1/(dx(ii)*(dx(ii-1)/D(ii-1,jj)+dx(ii)/D(ii,jj)))+1/(dx(ii)*(dx(ii+1)/D(ii+1,jj)+dx(ii)/D(ii,jj)))+1/(dy(ii)*(dy(ii-1)/D(ii-1,jj)+dy(ii)/D(ii,jj)))+1/(dy(ii)*(dy(ii+1)/D(ii+1,jj)+dy(ii)/D(ii,jj))))
                elseif (ii+1 == jj) then

                elseif (ii+3 == jj) then
                
                elseif (ii-1 == jj) then
                
                elseif (ii-3 == jj) then

                end if
            end do
        end do
        print'(10F10.5)',
    end subroutine
end module