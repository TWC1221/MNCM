module m_boundaries
  implicit none
  TYPE :: BoundaryCondition
    !-----------------------------
    ! Enumerated boundary types
    !-----------------------------
    integer :: BOUNDARY_NONE = 0
    integer :: BOUNDARY_ZERO = 1
    integer :: BOUNDARY_REFLECTIVE = 2
    integer :: BOUNDARY_VACUUM = 3
    integer :: BOUNDARY_ALBEDO = 4
    integer :: BOUNDARY_EXTRAPOLATED = 5
  end TYPE
end module