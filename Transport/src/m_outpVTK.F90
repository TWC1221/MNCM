
!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  Produce VTK files for 2D discretized NDE, xy, rz and rth geometries  -!  
!  Writes : Points, Cells, Cell_Types, Cell_Data                        -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 14/11/25     T. Charlton    Implemented 2D matrix                     -!
! 15/11/25     T. Charlton    Implemented rth matrix                    -!
! 17/11/25     T. Charlton    Implemented rz matrix                     -! 
!------------------------------------------------------------------------! 

module m_outpVTK
   implicit none
   contains
   subroutine outpVTK_xyz_transport(filename, phi, I, J, K, dx, dy, dz, use_adjoint)
      logical, optional, intent(in) :: use_adjoint
      character(len=*), intent(in) :: filename
      real(8), intent(in) :: phi(:,:,:)
      integer, intent(in) :: I, J, K

      real(8), allocatable :: dx(:), dy(:), dz(:)
      character(len=256) :: filename_final
      integer :: ii, jj, kk, n0, n1, n2, n3, n4, n5, n6, n7

      filename_final = trim(filename)
      if (present(use_adjoint)) then
         filename_final = trim(filename) // "_adjoint"
      endif

      open(unit=11, file="../" // trim(filename_final) // ".vtk", status="replace", action="write")

      ! Write header
      write(11,'(A)') '# vtk DataFile Version 2.0'
      write(11,'(A)') 'cell-centred output'
      write(11,'(A)') 'ASCII'
      write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'
      write(11,'(A,I0,A)') 'POINTS ', (I+1)*(J+1)*(K+1), ' double'

      ! Loop over 2D domain nodes
      do kk = 0, K
         do jj = 0, J 
            do ii = 0, I
               write(11,'(F9.5, F9.5, ES16.8)') ii*dx(1), jj*dy(1), kk*dz(1)
            end do
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELLS ', I*J*K,' ', 9*I*J*K

      do kk = 0, (K-1)
         do jj = 0, (J-1)
            do ii = 0, (I-1)
               n0 = kk*(I+1)*(J+1) + jj*(I+1) + ii
               n1 = n0 + 1
               n2 = n1 + (I+1)
               n3 = n0 + (I+1)
               n4 = n0 + (I+1)*(J+1)
               n5 = n4 + 1
               n7 = n4 + (I+1)
               n6 = n7 + 1

               write(11,'(I4,8I8)') 8, n0, n1, n2, n3, n4, n5, n6, n7
            end do
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_TYPES ', I*J*K
      do ii = 1, I*J*K
         write(11,'(I5)') 12    
      end do
        
      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_DATA ', I*J*K
      write(11,'(A)')   'SCALARS phi double'
      write(11,'(A)') 'LOOKUP_TABLE default'
      do kk = 1, K
         do jj = 1, J
            do ii = 1, I
               write(11,'(ES16.8)') phi(ii,jj,kk)
            end do
         end do
      end do

      close(11)
   end subroutine outpVTK_xyz_transport

  subroutine outpVTK_xyz_vector(filename, psi, I, J, K, dx, dy, dz, Nang, mu, eta, zeta)
    implicit none

    character(len=*), intent(in) :: filename
    real(8), intent(in) :: psi(:,:,:,:), mu(:), eta(:), zeta(:)
    integer, intent(in) :: I, J, K, Nang
    real(8), intent(in) :: dx(:), dy(:), dz(:)

    integer :: nn, jj, kk, ii, ios
    real(8), allocatable :: Jx(:,:,:), Jy(:,:,:), Jz(:,:,:)
    character(len=256) :: fname

    ! Build output filename in parent directory
    fname = trim(filename)
    if (len_trim(fname) >= 4) then
      if (fname(len_trim(fname)-3:len_trim(fname)) == '.vtk') then
        fname = fname(1:len_trim(fname)-4)
      end if
    end if

    ! Allocate vector fields
    allocate(Jx(I,J,K), Jy(I,J,K), Jz(I,J,K))
    Jx = 0.0d0
    Jy = 0.0d0
    Jz = 0.0d0

    ! Compute Jx, Jy, Jz
    do nn = 1, Nang
      do kk = 1, K
        do jj = 1, J
          do ii = 1, I
            Jx(ii,jj,kk) = Jx(ii,jj,kk) + mu(nn)  * psi(nn,ii,jj,kk)
            Jy(ii,jj,kk) = Jy(ii,jj,kk) + eta(nn) * psi(nn,ii,jj,kk)
            Jz(ii,jj,kk) = Jz(ii,jj,kk) + zeta(nn)* psi(nn,ii,jj,kk)
          end do
        end do
      end do
    end do

    ! Write VTK file in parent directory
    open(unit=11, file='../' // trim(fname) // '.vtk', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      print*, 'ERROR: Unable to open VTK file: ', '../' // trim(fname) // '.vtk'
      stop
    end if

   write(11,'(A)')'# vtk DataFile Version 3.0'
   write(11,'(A)')'Current vector field'
   write(11,'(A)')'ASCII'
   write(11,'(A)')'DATASET STRUCTURED_POINTS'
   write(11,'(A,3I8)')        'DIMENSIONS', I, J, K
   write(11,'(A,3F16.8)')     'ORIGIN', 0.0d0, 0.0d0, 0.0d0
   write(11,'(A,3F16.8)')     'SPACING', dx(1), dy(1), dz(1)
   write(11,'(A,I8)')         'POINT_DATA', I*J*K

   ! Treat each cell as a "point" for visualization.
   write(11,*) 'VECTORS J double'

    do kk = 1, K
      do jj = 1, J
        do ii = 1, I
          write(11,'(3(1X,F12.6))') Jx(ii,jj,kk), Jy(ii,jj,kk), Jz(ii,jj,kk)
        end do
      end do
    end do

    close(11)

    deallocate(Jx, Jy, Jz)

  end subroutine outpVTK_xyz_vector

  subroutine outpVTK_xyz_vector_allangles(filename, psi, I, J, K, dx, dy, dz, Nang, mu, eta, zeta)
  implicit none

  character(len=*), intent(in) :: filename
  real(8), intent(in) :: psi(:,:,:,:), mu(:), eta(:), zeta(:)
  integer, intent(in) :: I, J, K, Nang
  real(8), intent(in) :: dx(:), dy(:), dz(:)

  integer :: nn, jj, kk, ii, ios
  integer :: npoints
  character(len=256) :: fname

  ! Build output filename
  fname = trim(filename)
  if (len_trim(fname) >= 4) then
    if (fname(len_trim(fname)-3:len_trim(fname)) == '.vtk') then
      fname = fname(1:len_trim(fname)-4)
    end if
  end if

  ! Write VTK file
  open(unit=11, file='../' // trim(fname) // '.vtk', status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    print*, 'ERROR: Unable to open VTK file: ', '../' // trim(fname) // '.vtk'
    stop
  end if

  write(11,'(A)') '# vtk DataFile Version 3.0'
  write(11,'(A)') 'Multiple angles vector field'
  write(11,'(A)') 'ASCII'
  write(11,'(A)') 'DATASET STRUCTURED_POINTS'
  write(11,'(A,3I8)') 'DIMENSIONS', I, J, K
  write(11,'(A,3F16.8)') 'ORIGIN', dx(1)/2, dy(1)/2, dz(1)/2
  write(11,'(A,3F16.8)') 'SPACING', dx(1), dy(1), dz(1)

  npoints = I * J * K
  write(11,'(A,I8)') 'POINT_DATA', npoints

  ! ------------------------------------------------------------------
  ! Write one VECTOR field per angle
  ! ------------------------------------------------------------------
  do nn = 1, Nang, 1
    write(11,'(A)') 'VECTORS J_angle_'//trim(adjustl(itoa(nn)))//' double'

    do kk = 1, K
      do jj = 1, J
        do ii = 1, I
          write(11,'(3(1X,F12.6))') mu(nn)*psi(nn,ii,jj,kk), &
                                    eta(nn)*psi(nn,ii,jj,kk), &
                                    zeta(nn)*psi(nn,ii,jj,kk)
        end do
      end do
    end do
  end do

  close(11)

end subroutine outpVTK_xyz_vector_allangles

pure function itoa(i) result(str)
  integer, intent(in) :: i
  character(len=12) :: str
  write(str,'(I0)') i
end function itoa


end module m_outpVTK