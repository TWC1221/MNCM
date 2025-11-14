module m_outpVTK
    implicit none
    contains
    subroutine outpVTK(phi, nx, ny, dx, dy)
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: phi(:), dx, dy 
        integer :: ii, jj, kk, n0, n1, n2, n3
        real(8) :: x, y
        
        open(unit=11, file="flux_points.vtk", status="replace", action="write")

        ! Write header
        write(11,'(A)') '# vtk DataFile Version 2.0'
        write(11,'(A)') 'cell-centred output'
        write(11,'(A)') 'ASCII'
        write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'
        write(11,'(A,I0,A)') 'POINTS ', (nx+1)*(ny+1), ' double'

        ! Loop over 2D domain nodes
        do jj = 0, ny
            y = jj*dy        ! y-direction
            do ii = 0, nx    ! x-direction
                x = ii*dx
                write(11,'(F8.5 F8.5 ES16.8)') x, y, 0.0
            end do
        end do

        write(11,'(A)') ''
        write(11,'(A,I0,A,I0)') 'CELLS ', nx*ny,' ', 5*nx*ny

        do jj = 0, (ny-1)
            do ii = 0, (nx-1)
                ! point indices for quad cell
                ! node numbering: row-major
                n0 = jj*(nx+1) + ii
                n1 = n0 + 1
                n2 = n1 + (nx+1)
                n3 = n0 + (nx+1)
                write(11,'(I5,4I5)') 4, n0, n1, n2, n3
            end do
        end do

        write(11,'(A)') ''
        write(11,'(A,I0,A,I0)') 'CELL_TYPES ', nx*ny

        do jj = 1, nx*ny
            write(11,'(I5)') 9    
        end do
        
        write(11,'(A)') ''
        write(11,'(A,I0,A,I0)') 'CELL_DATA ', nx*ny
        write(11,'(A)') 'SCALARS phi double'
        write(11,'(A)') 'LOOKUP_TABLE default'
        do jj = 1, ny
            do ii = 1, nx
                kk = (jj-1)*nx + ii
                write(11,'(ES16.8)') phi(kk)
            end do
        end do

        close(11)
    end subroutine outpVTK
end module m_outpVTK