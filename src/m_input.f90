module m_ReadInput
    use m_system
    use m_global
    use m_constants
    implicit none

    contains

    subroutine ReadInput(GLOB, InputFile)

        type(t_global)  :: GLOB
        character(60), intent(in)       :: InputFile
        character(15)                   :: line
        integer                         :: iunit, ierr

        iunit = 1
        ! Open the file
        open(newunit=iunit, file=InputFile, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, 'Error opening file ', InputFile
            stop
        end if

        GLOB%OutputVTK = .false.
        GLOB%OutputData = .false.
        GLOB%OutMaterial = .false.
        GLOB%CalcDiffCoeff = .false.
        GLOB%CalcCondition = .false.
        GLOB%TerminalOutput = .false.

        ! ====================================================================================================
        ! Load GLOB parameters
        read(iunit,'(A16,A)', iostat=ierr) line, GLOB%MeshPath
        read(iunit,'(A16,A)', iostat=ierr) line, GLOB%MatPath
        read(iunit,'(A16,A)', iostat=ierr) line, GLOB%OutVTKPath
        read(iunit,'(A16,A)', iostat=ierr) line, GLOB%OutDataPath
        read(iunit,*, iostat=ierr) line, GLOB%OutputVTK
        read(iunit,*, iostat=ierr) line, GLOB%RUN_ANALYSIS
        read(iunit,*, iostat=ierr) line, GLOB%OutputData
        read(iunit,*, iostat=ierr) line, GLOB%OutMaterial
        read(iunit,*, iostat=ierr) line, GLOB%OutFluxInt
        read(iunit,*, iostat=ierr) line, GLOB%Outkeff
        read(iunit,*, iostat=ierr) line, GLOB%CalcDiffCoeff
        read(iunit,*, iostat=ierr) line, GLOB%Method
        read(iunit,*, iostat=ierr) line, GLOB%Coordsys
        read(iunit,*, iostat=ierr) line, GLOB%Preconditioner
        read(iunit,*, iostat=ierr) line, GLOB%Source
        read(iunit,*, iostat=ierr) line, GLOB%NoGroups
        read(iunit,*, iostat=ierr) line, GLOB%SN
        read(iunit,*, iostat=ierr) line, GLOB%Anisotropy
        read(iunit,*, iostat=ierr) line, GLOB%adjoint
        read(iunit,*, iostat=ierr) line, GLOB%Tolerance
        read(iunit,*, iostat=ierr) line, GLOB%MaxIter
        read(iunit,*, iostat=ierr) line, GLOB%CalcCondition
        read(iunit,*, iostat=ierr) line, GLOB%TerminalOutput
                
        if (GLOB%CalcCondition) then
            allocate(GLOB%lambda_min(GLOB%NoGroups))
            allocate(GLOB%lambda_max(GLOB%NoGroups))
            allocate(GLOB%cond_number(GLOB%NoGroups))
            GLOB%lambda_min = 0.0_dp
            GLOB%lambda_max = 0.0_dp
            GLOB%cond_number = 0.0_dp
        end if
        
        write(*,'(A)') '-------------------------------------------------------------------------------'
        write(*,'(A)') '| File paths                                                                  |'
        write(*,'(A)') '|                                                                             |'
        write(*,'(A18,A,A)') '| Input file:       ', InputFile, '|'
        write(*,'(A18,A,A)') '| Material file:    ', GLOB%MatPath, '|'
        write(*,'(A18,A,A)') '| Mesh file:        ', GLOB%MeshPath, '|'
        write(*,'(A18,A,A)') '| VTK file:         ', GLOB%OutputVTK, '|'
        write(*,'(A18,A,A)') '| Data file:        ', GLOB%OutDataPath, '|'
        write(*,'(A)') '|                                                                             |'
        write(*,'(A)') '-------------------------------------------------------------------------------'

    end subroutine ReadInput


    subroutine initialise_boundary_conditions(sys, iunit, ierr)
        type(t_system)       :: sys
        character(15)           :: line
        integer                 :: iunit, ierr
        
        read(iunit,*, iostat=ierr) line, sys%BC_left
        read(iunit,*, iostat=ierr) line, sys%BC_right
        read(iunit,*, iostat=ierr) line, sys%BC_top
        read(iunit,*, iostat=ierr) line, sys%BC_bottom
        read(iunit,*, iostat=ierr) line, sys%BC_front
        read(iunit,*, iostat=ierr) line, sys%BC_back

        
        ! ### Check if boundary conditions are valid for periodic ###

        if (sys%BC_left == 6 .and. sys%BC_right /= 6 .or.  &
            sys%BC_left /= 6 .and. sys%BC_right == 6) then
            print *, 'Read input file error: Left and right boundary conditions must be the same for periodic.'
            stop
        end if

        if (sys%BC_top == 6 .and. sys%BC_bottom /= 6 .or.  &
            sys%BC_top /= 6 .and. sys%BC_bottom == 6) then
            print *, 'Read input file error: Top and bottom boundary conditions must be the same for periodic.'
            stop
        end if

    end subroutine initialise_boundary_conditions

    subroutine RunInputChecks(GLOB)
        type(t_GLOBal)      :: GLOB

        ! write(*,*) 'Running input checks...'

        if (GLOB%Method == 'Transport') then

            ! Check SN 1 or even
            if (mod(GLOB%SN, 2) /= 0) then
                print *, 'Input file error: Number of SN must be even for transport method.'
                stop
            end if

        end if

        if (GLOB%source /= 'Fixed' .and. GLOB%source /= 'Eigenvalue') then
            write(*,*) '!---------------------------------------------------------!'
            write(*,*) '!- Input file error: Source must be Fixed or Eigenvalue. -!'
            write(*,*) '!---------------------------------------------------------!'
            print*
            stop
        end if

        if (GLOB%Coordsys == 2 .and. GLOB%Adjoint .eqv. .true.) then
            print *, 'Input file error: Adjoint not yet supported 1D spherical coordinates.'
            stop
        end if

        ! write(*,*) 'Input checks passed.'


    end subroutine RunInputChecks


end module m_ReadInput