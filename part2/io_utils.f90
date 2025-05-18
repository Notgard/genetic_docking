module io_utils
    implicit none
    
    private
    public :: create_csv, add_csv_entry
    
contains
    
    ! Subroutine to create a new CSV file with headers
    subroutine create_csv(filename)
        character(len=*), intent(in) :: filename
        integer :: unit_num, io_status
        
        ! Open the file for writing
        open(newunit=unit_num, file=trim(filename), status='replace', iostat=io_status)
        
        if (io_status /= 0) then
            write(*, '(A,A)') "Error: Could not create file ", trim(filename)
            return
        end if
        
        ! Write the header
        write(unit_num, '(A)') "Generation,Min_Fitness,Max_Fitness,Avg_Fitness"
        
        ! Close the file
        close(unit_num)
    end subroutine create_csv
    
    ! Subroutine to add an entry to an existing CSV file
    subroutine add_csv_entry(filename, generation, min_fitness, max_fitness, avg_fitness)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: generation
        real, intent(in) :: min_fitness, max_fitness, avg_fitness
        integer :: unit_num, io_status
        logical :: file_exists
        
        ! Check if file exists
        inquire(file=trim(filename), exist=file_exists)
        if (.not. file_exists) then
            write(*, '(A,A)') "Error: File does not exist: ", trim(filename)
            return
        end if
        
        ! Open the file for appending
        open(newunit=unit_num, file=trim(filename), status='old', position='append', iostat=io_status)
        
        if (io_status /= 0) then
            write(*, '(A,A)') "Error: Could not open file ", trim(filename)
            return
        end if
        
        ! Write the data
        write(unit_num, '(I0,",",F0.6,",",F0.6,",",F0.6)') &
            generation, min_fitness, max_fitness, avg_fitness
        
        ! Close the file
        close(unit_num)
    end subroutine add_csv_entry
    
end module io_utils