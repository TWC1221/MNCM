module m_system

    use m_constants
    use m_global
    implicit none
    private

    public :: t_system, CalcNodeBoundaryOrder
    public :: XABORT, XWARNING

    type :: t_system

        real(dp)                                    :: xmin, xmax, ymin, ymax, zmin, zmax                 
        
        integer, dimension(:,:), allocatable        :: initial_setup                   
        integer, dimension(:,:), allocatable        :: MatID      
        integer, dimension(:,:), allocatable        :: source_setup             
        integer                                     :: nCellx, nCelly               ! number of cells in x and y
        integer                                     :: nCellTotal                   ! total number of cells
        integer                                     :: No1DElements                 ! number of 1D elements in the system
        integer                                     :: No2DElements                 ! number of 2D elements in the system
        integer                                     :: No3DElements                 ! number of 3D elements in the system
        integer                                     :: NoNodes                      ! number of nodes in the system
        integer                                     :: NoCtrlPnts                   ! number of control points in the system

        integer, dimension(:,:), allocatable        :: NE_List
        integer, dimension(:,:), allocatable        :: Sweep                        ! Sweep order of elements
        integer, dimension(:,:), allocatable        :: Sweep_half                   ! Sweep order for half angles in 2D cylindrical geometry
        integer, dimension(:), allocatable          :: ElementKinds                 ! ID of elements in the system
        integer, dimension(:), allocatable          :: PolyOrders                   ! Polynomial orders of elements
        integer, dimension(:,:), allocatable        :: NBOrder                      ! Node boundary order for FE/SE
        integer, dimension(:,:,:), allocatable      :: NBOrder3D
        integer, dimension(:,:), allocatable        :: NBORder3Dfinv, NBORder3Df
        integer, dimension(:,:,:,:), allocatable    :: NBOrder3Dfexp
        integer, dimension(:,:), allocatable        :: NBOrderNURB                      ! Node boundary order for NURBS or NURBS-Enhanced
        

        real(dp), dimension(:), allocatable       :: Chi
        real(dp), dimension(4)                      :: Albedo, SurfaceSource
        real(dp)                                    :: albedo_surf

        integer                                     :: BC_left, BC_right, BC_bottom, BC_top, BC_front, BC_back
        integer, dimension(:,:), allocatable        :: REFLMAP     ! Reflected angles mapping
        real(dp), dimension(:,:,:), allocatable     :: omega_dot_n  ! Omega dot n mapping
        integer, dimension(:), allocatable          :: face_indices ! Reflected angles mapping for 1D
        real(dp), dimension(:,:,:,:), allocatable   :: face_nodes   ! Omega dot n mapping for 3D

        integer, dimension(:,:), allocatable        :: NOUT     ! number of outgoing directions
        integer, dimension(:,:), allocatable        :: NIN      ! number of incoming directions
        integer, dimension(:,:,:), allocatable      :: OUTIND   ! outgoing indices
        integer, dimension(:,:,:), allocatable      :: ININD    ! incoming indices

        integer, dimension(:,:), allocatable        :: NBOUND   ! number of boundary faces for each element
        integer, dimension(:,:,:), allocatable      :: BOUND    ! boundary face indices for each element

        integer, dimension(:,:), allocatable        :: BNDINDX      ! boundary face global index for each element

        integer, dimension(:,:), allocatable        :: BNDELEM      ! element index for each boundary 
        integer, dimension(:,:), allocatable        :: BNDELEMBC    ! number of faces for each element

        integer, dimension(:,:,:), allocatable      :: neighbours   ! neighbouring elements for each element
        integer, dimension(:,:,:), allocatable      :: neigh_faces  ! face indices of neighbouring elements for each element
        integer, dimension(:,:), allocatable        :: bound_faces  ! boundary faces for each element

        integer, dimension(:,:), allocatable        :: EN_LIST_1D   ! element to node connectivity for 1D elements
        integer, dimension(:,:), allocatable        :: EN_LIST_2D   ! element to node connectivity for 2D elements
        integer, dimension(:,:), allocatable        :: EN_LIST_3D   ! element to node connectivity for 3D elements

        real(dp), dimension(:,:,:), allocatable     :: FLUX         ! Scalar flux (NELEM,NGRP,NNODE)

        contains

            procedure :: write_system_info
            procedure :: get_sweep_order
            procedure :: get_omega_dot_n
            procedure :: get_face_indices
            procedure :: get_face_nodes
            procedure :: get_nb_order
            procedure :: get_NOUT
            procedure :: get_NIN
            procedure :: get_ININD
            procedure :: get_OUTIND
            procedure :: get_neighbours
            procedure :: get_neigh_faces
            procedure :: get_NBOUND
            procedure :: get_BOUND
            procedure :: get_BNDINDX
            procedure :: get_BNDELEM
            procedure :: get_BNDELEMBC
            procedure :: get_REFLMAP
            procedure :: get_chi_spectrum
            procedure :: get_sweep
            procedure :: get_NBORDER3DF
            procedure :: get_EN_LIST_1D
            procedure :: get_EN_LIST_2D
            procedure :: get_EN_LIST_3D
            procedure :: collect_flux_transport
            procedure :: collect_flux_diffusion

    end type t_system
contains 

    function get_EN_LIST_1D(this) result (EN_LIST)

        class(t_system), intent(in)             :: this
        integer, dimension(:,:), allocatable    :: EN_LIST

        if (allocated(this%EN_LIST_1D)) then
            EN_LIST = this%EN_LIST_1D
        else
            ! Raise an error 
            write(*,*) 'Error: Element to node connectivity is not allocated in the system.'
            stop
        end if

    end function get_EN_LIST_1D

    function get_EN_LIST_2D(this) result (EN_LIST)

        class(t_system), intent(in)             :: this
        integer, dimension(:,:), allocatable    :: EN_LIST

        if (allocated(this%EN_LIST_2D)) then
            EN_LIST = this%EN_LIST_2D
        else
            ! Raise an error 
            write(*,*) 'Error: Element to node connectivity is not allocated in the system.'
            stop
        end if

    end function get_EN_LIST_2D

    function get_EN_LIST_3D(this) result (EN_LIST)

        class(t_system), intent(in)             :: this
        integer, dimension(:,:), allocatable    :: EN_LIST

        if (allocated(this%EN_LIST_3D)) then
            EN_LIST = this%EN_LIST_3D
        else
            ! Raise an error 
            write(*,*) 'Error: Element to node connectivity is not allocated in the system.'
            stop
        end if

    end function get_EN_LIST_3D

    function get_sweep_order(this) result (sweep_order)

        class(t_system), intent(in)             :: this
        integer, dimension(:,:), allocatable    :: sweep_order

        if (allocated(this%Sweep)) then
            sweep_order = this%Sweep
        else
            ! Raise an error 
            write(*,*) 'Error: Sweep order is not allocated in the system.'
            stop
        end if

    end function get_sweep_order

    function get_omega_dot_n(this) result (omega_dot_n)

        class(t_system), intent(in)             :: this
        real(dp), dimension(:,:,:), allocatable :: omega_dot_n

        if (allocated(this%omega_dot_n)) then
            omega_dot_n = this%omega_dot_n
        else
            call XABORT('SYS: Omega_dot_n is not allocated.')
        end if

    end function get_omega_dot_n

    function get_face_indices(this) result (face_indices)

        class(t_system), intent(in)             :: this
        integer, dimension(:), allocatable       :: face_indices

        if (allocated(this%face_indices)) then
            face_indices = this%face_indices
        else
            call XABORT('SYS: Face indices are not allocated.')
        end if

    end function get_face_indices

    function get_face_nodes(this) result (face_nodes)

        class(t_system), intent(in) :: this
        real(dp), allocatable       :: face_nodes(:,:,:,:)

        if (allocated(this%face_nodes)) then
            face_nodes = this%face_nodes
        else
            call XABORT('SYS: Face nodes are not allocated.')
        end if

    end function get_face_nodes

    function get_nb_order(this) result (nb_order)

        class(t_system), intent(in) :: this
        integer, dimension(:,:), allocatable :: nb_order

        if (allocated(this%NBOrder)) then
            nb_order = this%NBOrder
        else
            call XABORT('SYS: Node boundary order is not allocated.')
        end if

    end function get_nb_order

    function get_NOUT(this) result (NOUT)

        class(t_system), intent(in) :: this
        integer, dimension(:,:), allocatable :: NOUT

        if (allocated(this%NOUT)) then
            NOUT = this%NOUT
        else
            call XABORT('SYS: NOUT is not allocated.')
        end if

    end function get_NOUT

    function get_NIN(this) result (NIN)

        class(t_system), intent(in) :: this
        integer, dimension(:,:), allocatable :: NIN

        if (allocated(this%NIN)) then
            NIN = this%NIN
        else
            call XABORT('SYS: NIN is not allocated.')
        end if

    end function get_NIN

    function get_ININD(this) result (ININD)

        class(t_system), intent(in) :: this
        integer, allocatable :: ININD(:,:,:)

        if (allocated(this%ININD)) then
            ININD = this%ININD
        else
            call XABORT('SYS: ININD is not allocated.')
        end if

    end function get_ININD

    function get_OUTIND(this) result (OUTIND)

        class(t_system), intent(in) :: this
        integer, allocatable :: OUTIND(:,:,:)

        if (allocated(this%OUTIND)) then
            OUTIND = this%OUTIND
        else
            call XABORT('SYS: OUTIND is not allocated.')
        end if

    end function get_OUTIND

    function get_neighbours(this) result (neighbours)

        class(t_system), intent(in) :: this
        integer, allocatable :: neighbours(:,:,:)

        if (allocated(this%neighbours)) then
            neighbours = this%neighbours
        else
            call XABORT('SYS: Neighbours is not allocated.')
        end if

    end function get_neighbours

    function get_neigh_faces(this) result (neigh_faces)

        class(t_system), intent(in) :: this
        integer, allocatable :: neigh_faces(:,:,:)

        if (allocated(this%neigh_faces)) then
            neigh_faces = this%neigh_faces
        else
            call XABORT('SYS: Neighbour faces are not allocated.')
        end if

    end function get_neigh_faces

    function get_NBOUND(this) result (NBOUND)

        class(t_system), intent(in) :: this
        integer, dimension(:,:), allocatable :: NBOUND

        if (allocated(this%NBOUND)) then
            NBOUND = this%NBOUND
        else
            call XABORT('SYS: NBOUND is not allocated.')
        end if

    end function get_NBOUND

    function get_BOUND(this) result (BOUND)

        class(t_system), intent(in) :: this
        integer, dimension(:,:,:), allocatable :: BOUND

        if (allocated(this%BOUND)) then
            BOUND = this%BOUND
        else
            call XABORT('SYS: BOUND is not allocated.')
        end if

    end function get_BOUND

    function get_BNDINDX(this) result (BNDINDX)

        class(t_system), intent(in) :: this
        integer, allocatable :: BNDINDX(:,:)

        if (allocated(this%BNDINDX)) then
            BNDINDX = this%BNDINDX
        else
            call XABORT('SYS: BNDINDX is not allocated.')
        end if

    end function get_BNDINDX

    function get_BNDELEM(this) result (BNDELEM)

        class(t_system), intent(in) :: this
        integer, allocatable :: BNDELEM(:,:)

        if (allocated(this%BNDELEM)) then
            BNDELEM = this%BNDELEM
        else
            call XABORT('SYS: BNDELEM is not allocated.')
        end if

    end function get_BNDELEM

    function get_chi_spectrum(this) result (CHI)

        class(t_system), intent(in)             :: this
        real(dp), dimension(:), allocatable   :: CHI

        if (allocated(this%CHI)) then
            CHI = this%CHI
        else
            call XABORT('SYS: Chi spectrum is not allocated.')
        end if

    end function get_chi_spectrum

    function get_BNDELEMBC(this) result (BNDELEMBC)

        class(t_system), intent(in) :: this
        integer, allocatable :: BNDELEMBC(:,:)

        if (allocated(this%BNDELEMBC)) then
            BNDELEMBC = this%BNDELEMBC
        else
            call XABORT('SYS: BNDELEMBC is not allocated.')
        end if

    end function get_BNDELEMBC

    function get_REFLMAP(this) result (REFLMAP)

        class(t_system), intent(in) :: this
        integer, allocatable :: REFLMAP(:,:)

        if (allocated(this%REFLMAP)) then
            REFLMAP = this%REFLMAP
        else
            call XABORT('SYS: Reflective map is not allocated.')
        end if

    end function get_REFLMAP

    function get_sweep(this) result (SWEEP)

        class(t_system), intent(in) :: this
        integer, allocatable        :: SWEEP(:,:)

        if (allocated(this%SWEEP)) then
            SWEEP = this%SWEEP
        else
            call XABORT('SYS: Sweep is not allocated.')
        end if

    end function get_sweep

    function get_NBORDER3DF(this) result (NBORDER3DF)

        class(t_system), intent(in) :: this
        integer, allocatable        :: NBORDER3DF(:,:,:,:)

        if (allocated(this%NBOrder3Dfexp)) then
            NBORDER3DF = this%NBOrder3Dfexp
        else
            call XABORT('SYS: NBORDER3DF is not allocated.')
        end if

    end function get_NBORDER3DF

    subroutine collect_flux_transport(this, FLUX, NELEM, NGRP, NNODE)

        class(t_system), intent(inout)  :: this
        real(dp)                        :: FLUX(:,:,:)
        integer                         :: NELEM, NGRP, NNODE

        allocate(this%FLUX(NELEM, NGRP, NNODE))
        this%FLUX = FLUX

    end subroutine collect_flux_transport
    
    subroutine collect_flux_diffusion(this, FLUX, EN_LIST, NELEM, NGRP, NNODE)

        ! Input variables
        class (t_system), intent(inout) :: this
        real(dp)                        :: FLUX(:,:)
        integer, intent(in)             :: EN_LIST(:,:)
        integer, intent(in)             :: NELEM, NGRP, NNODE

        ! Local variables
        integer :: EINDX, GINDX

        allocate(this%FLUX(NELEM, NGRP, NNODE))
        this%FLUX = 0.0_dp

        do EINDX = 1, NELEM
            do GINDX = 1, NGRP
                this%FLUX(EINDX, GINDX, :) = FLUX(EN_LIST(EINDX, :), GINDX)
            end do
        end do

    end subroutine collect_flux_diffusion

    subroutine write_system_info(this)
            
            class(t_system), intent(inout) :: this

            write(*,*) '------------------Geometry info------------------'
            ! write(*,*) 'Number of regions in x: ', this%no_x_regions
            ! write(*,*) 'Number of regions in y: ', this%no_y_regions
            ! write(*,*) 'dx: ', this%dx(1)
            ! write(*,*) 'dy: ', this%dy(1)

            write(*,*) 'Boundary conditions'
            if (this%BC_left == 1) then
                write(*,*) 'Left boundary:   Zero flux'
            else if (this%BC_left == 2) then
                write(*,*) 'Left boundary:   Reflective'
            else if (this%BC_left == 3) then
                write(*,*) 'Left boundary:   Vacuum'
            else if (this%BC_left == 4) then
                write(*,*) 'Left boundary:   Albedo'
            else if (this%BC_left == 5) then
                write(*,*) 'Left boundary:   Periodic'
            end if

            if (this%BC_right == 1) then
                write(*,*) 'Right boundary:  Zero flux'
            else if (this%BC_right == 2) then
                write(*,*) 'Right boundary:  Reflective'
            else if (this%BC_right == 3) then
                write(*,*) 'Right boundary:  Vacuum'
            else if (this%BC_right == 4) then
                write(*,*) 'Right boundary:  Albedo'
            else if (this%BC_right == 5) then
                write(*,*) 'Right boundary:  Periodic'
            end if

            if (this%BC_bottom == 1) then
                write(*,*) 'Bottom boundary:  Zero flux'
            else if (this%BC_bottom == 2) then
                write(*,*) 'Bottom boundary: Reflective'
            else if (this%BC_bottom == 3) then
                write(*,*) 'Bottom boundary: Vacuum'
            else if (this%BC_bottom == 4) then
                write(*,*) 'Bottom boundary: Albedo'
            else if (this%BC_bottom == 5) then
                write(*,*) 'Bottom boundary: Periodic'
            end if

            if (this%BC_top == 1) then
                write(*,*) 'Top boundary:    Zero flux'
            else if (this%BC_top == 2) then
                write(*,*) 'Top boundary:    Reflective'
            else if (this%BC_top == 3) then
                write(*,*) 'Top boundary:    Vacuum'
            else if (this%BC_top == 4) then
                write(*,*) 'Top boundary:    Albedo'
            else if (this%BC_top == 5) then
                write(*,*) 'Top boundary:    Periodic'
            end if

            if (this%BC_front == 1) then
                write(*,*) 'Front boundary:  Zero flux'
            else if (this%BC_front == 2) then
                write(*,*) 'Front boundary:  Reflective'
            else if (this%BC_front == 3) then
                write(*,*) 'Front boundary:  Vacuum'
            else if (this%BC_front == 4) then
                write(*,*) 'Front boundary:  Albedo'
            else if (this%BC_front == 5) then
                write(*,*) 'Front boundary:  Periodic'
            end if

            if (this%BC_back == 1) then
                write(*,*) 'Back boundary:   Zero flux'
            else if (this%BC_back == 2) then
                write(*,*) 'Back boundary:   Reflective'
            else if (this%BC_back == 3) then
                write(*,*) 'Back boundary:   Vacuum'
            else if (this%BC_back == 4) then
                write(*,*) 'Back boundary:   Albedo'
            else if (this%BC_back == 5) then
                write(*,*) 'Back boundary:   Periodic'
            end if

            write(*,*) '-------------------------------------------------'
            

    end subroutine write_system_info

    subroutine CalcNodeBoundaryOrder(sys, NDIM, NBASIS)
        
        type(t_system), intent(inout)   :: Sys
        integer, intent(in)             :: NDIM, NBASIS
        integer                         :: i, j, k, n
        
        if (NDIM==2) then
            allocate(Sys%NBOrder(4, NBASIS))
            ! First face of NBOrder
            do i = 1, NBASIS
                Sys%NBOrder(1,i) = i
            end do

            ! Second face of NBOrder
            Sys%NBOrder(2,1) = NBASIS
            do i = 1, NBASIS - 1
                Sys%NBOrder(2,i+1) = Sys%NBOrder(2,i) + NBASIS 
            end do

            ! Third face of NBOrder
            Sys%NBOrder(3,1) = NBASIS**2
            do i = NBASIS - 1, 1, -1
                Sys%NBOrder(3,i+1) = NBASIS**2 - i
            end do

            ! Fourth face of NBOrder
            Sys%NBOrder(4,1) = NBASIS**2 - NBASIS + 1
            do i = 1, NBASIS - 1
                Sys%NBOrder(4,i+1) = Sys%NBOrder(4,i) - NBASIS
            end do
            
        elseif (NDIM==3) then
            allocate(Sys%NBOrder(6, NBASIS * NBASIS))
            allocate(Sys%NBORder3Dfinv(6, NBASIS*NBASIS))

            allocate(Sys%NBOrder3D(NBASIS, NBASIS, NBASIS))

            n = 1
            do i = 1, NBASIS
                do j = 1, NBASIS
                    do k = 1, NBASIS
                        sys%NBOrder3D(i,j,k) = n
                        n = n + 1
                    end do
                end do
            end do

            sys%NBORder3Dfinv = 0

            n = 1
            do i = 1, NBASIS
                do j = 1, NBASIS

                    sys%NBOrder(1,n) = Sys%NBOrder3D(i,j,1)
                    sys%NBOrder(2,n) = Sys%NBOrder3D(i,j,NBASIS)

                    sys%NBOrder(3,n) = Sys%NBOrder3D(i,1,j)
                    sys%NBOrder(4,n) = Sys%NBOrder3D(i,NBASIS,j)

                    sys%NBOrder(5,n) = Sys%NBOrder3D(1,i,j)
                    sys%NBOrder(6,n) = Sys%NBOrder3D(NBASIS,i,j)

                    sys%NBORder3Dfinv(1,n) = Sys%NBOrder3D(i,j,1)
                    sys%NBORder3Dfinv(2,n) = Sys%NBOrder3D(i,j,NBASIS)

                    sys%NBORder3Dfinv(3,n) = Sys%NBOrder3D(i,1,j)
                    sys%NBORder3Dfinv(4,n) = Sys%NBOrder3D(i,NBASIS,j)
                    
                    sys%NBORder3Dfinv(5,n) = Sys%NBOrder3D(1,i,j)
                    sys%NBORder3Dfinv(6,n) = Sys%NBOrder3D(NBASIS,i,j)
                    n = n + 1
                end do
            end do
        end if

    end subroutine CalcNodeBoundaryOrder


    subroutine XABORT(msg)
        character(len=*), intent(in) :: msg
        write(*,*) msg
        stop
    end subroutine XABORT

    subroutine XWARNING(msg)
        character(len=*), intent(in) :: msg
        write(*,*) msg
    end subroutine XWARNING

end module m_system