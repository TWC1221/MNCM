Module C_timer
    !
    ! PURPOSE:
    !   To time functions so that the speed of the program can be calculated accuratly
    !
    ! Record of Revsions:
    !   Date        Programer       Change
    !   04/11/22    J.I.M           Original code
    !   23/06/23    J.I.M           Updated the timer code to use date
    !
    use m_constants
    use m_global
    implicit none
    private
    public t_timer, OutputTimesToTerminal

    TYPE t_timer
        REAL(Kind=dp)   :: time
        real(dp)        :: Secs,mins,milli,hours
        integer         :: newsecs,newmins,newmilli,newhours

        real(kind=dp)   :: Total_seconds_start, Total_seconds_end, Total_seconds
        real(kind=dp)   :: Total_seconds_start_temp, Total_seconds_end_temp, Total_seconds_temp

        integer         :: temp_sec,temp_min,temp_hour,temp_milli
        integer         :: temp_newsec,temp_newmin,temp_newhour,temp_newmilli
    CONTAINS
        PROCEDURE,pass :: InitTimer     
        PROCEDURE,pass :: Startdate
        PROCEDURE,pass :: elapseddate
        PROCEDURE,PASS :: start_timer
        PROCEDURE,PASS :: elapsed_time
        PROCEDURE,PASS :: start_temp_date
        PROCEDURE,PASS :: Add_Temp_Time
    end type t_timer
    
contains

    subroutine InitTimer(this)
        class(t_timer)  :: this
        this%Total_seconds = 0.0_dp
        this%Total_seconds_temp = 0.0_dp

    end subroutine

    subroutine Startdate(this)
        class(t_timer)  :: this
        integer         :: values(8)
        call date_and_time(values=values)
        this%hours = values(5)
        this%mins = values(6)
        this%Secs = values(7)
        this%milli = values(8)

        this%Total_seconds_start = real(this%Secs) + real(this%mins * 60,dp) + real(this%hours * 3600,dp) + real(this%milli * 0.001,dp)

    end subroutine

    subroutine elapseddate(this)
        class(t_timer)  :: this
        integer         :: values(8)
        call date_and_time(values=values)
        this%newhours = values(5)
        this%newmins = values(6)
        this%newsecs = values(7)
        this%newmilli = values(8)

        this%Total_seconds_end = real(this%newsecs, dp) + real(this%newmins * 60,dp) + real(this%newhours * 3600, dp) + real(this%newmilli * 0.001,dp)
        this%Total_seconds = this%Total_seconds_end - this%Total_seconds_start
        if (abs(this%Total_seconds) < 1e-5) then
            this%Total_seconds = 0.0_dp
        end if

    end subroutine

    subroutine start_temp_date(this)
        class(t_timer)  :: this
        integer         :: values(8)
        call date_and_time(values=values)
        this%temp_hour = values(5)
        this%temp_min = values(6)
        this%temp_sec = values(7)
        this%temp_milli = values(8)

        this%Total_seconds_start_temp = this%temp_sec + (this%temp_min * 60) + (this%temp_hour * 3600) + (this%temp_milli * 0.001)
        this%Total_seconds_end_temp = 0.0_dp
        this%Total_seconds_temp = 0.0_dp


    end subroutine

    subroutine Add_Temp_Time(this)
        class(t_timer)  :: this
        integer         :: values(8)
        call date_and_time(values=values)
        this%temp_newhour = values(5)
        this%temp_newmin = values(6)
        this%temp_newsec = values(7)
        this%temp_newmilli = values(8)

        this%Total_seconds_end_temp = this%temp_newsec + (this%temp_newmin * 60) + (this%temp_newhour * 3600) + (this%temp_newmilli * 0.001)
        this%Total_seconds_temp = this%Total_seconds_temp + this%Total_seconds_end_temp - this%Total_seconds_start_temp
        this%Total_seconds = this%Total_seconds + this%Total_seconds_temp

    end subroutine


    subroutine start_timer(this) !starts the timer
        CLASS(t_timer),INTENT(INOUT) :: this
        CALL CPU_TIME(this%time)
    end subroutine start_timer

    subroutine elapsed_time(this) !calls the time taken
        CLASS(t_timer),INTENT(INOUT) :: this
        REAL(KIND = dp) :: stopTime
        CALL CPU_TIME(stopTime)
        print *, "THE ELAPSED SERIAL TIME IS: ", (stopTime - this%time)
        this%time = stopTime
    end subroutine elapsed_time

    subroutine OutputTimesToTerminal(glob, Global, Mesh, ElemTimer, MatCreateTimer, Solve, Output)

        type(t_global), intent(in) :: glob
        type(t_timer), intent(in) :: Global, Mesh, ElemTimer, MatCreateTimer, Solve, Output

        write(*,'(A)') ''
        write(*,'(A)') '---------------------------------------------------------------------'
        write(*,'(A)') '- Time Breakdown:                                                   -'
        write(*,'(A)') '- Stage                Total Time (s)             Percentage (%)    -'
        write(*,'(A)') '---------------------------------------------------------------------'
        write(*,'(A11, F17.3, F28.2,A13)') '- Total Run', Global%Total_seconds, (Global%Total_seconds / Global%Total_seconds) * 100, '-'
        write(*,'(A11, F17.3, F28.2,A13)') '- Mesh       ', Mesh%Total_seconds, (Mesh%Total_seconds / Global%Total_seconds) * 100, '-'
        write(*,'(A14, F14.3, F28.2,A13)') '- ElemCreate   ', ElemTimer%Total_seconds, (ElemTimer%Total_seconds / Global%Total_seconds) * 100, '-'

        if (glob%Method=='Diffusion') then
            write(*,'(A14, F14.3, F28.2,A13)') '- MatrixCreate   ', MatCreateTimer%Total_seconds, (MatCreateTimer%Total_seconds / Global%Total_seconds) * 100, '-'
        end if
        
        write(*,'(A11, F17.3, F28.2,A13)') '- Solve      ', Solve%Total_seconds, (Solve%Total_seconds / Global%Total_seconds) * 100, '-'    
        write(*,'(A17, F11.3, F28.2,A13)') '- Post-processing', Output%Total_seconds, (Output%Total_seconds / Global%Total_seconds) * 100, '-'
        write(*,'(A)') '---------------------------------------------------------------------'

    end subroutine OutputTimesToTerminal

end module