program GameOfLife
!---------------------------------------------------------------------
!
!  This program runs a Collisional N-body simulation
!
!  Uses:  
!
!---------------------------------------------------------------------
        use nbody_common 
        implicit none
        interface 
        ! prototypes
        subroutine accel_update(opt, parts, nparts)
                use nbody_common 
                implicit none
                type(Options), intent(in) :: opt
                integer, intent(in) :: nparts 
                type(Particle), dimension(:), intent(inout) :: parts
        end subroutine
        subroutine velocity_update(opt, parts, nparts)
                use nbody_common 
                implicit none
                type(Options), intent(in) :: opt
                integer, intent(in) :: nparts 
                type(Particle), dimension(:), intent(inout) :: parts
        end subroutine
        subroutine position_update(opt, parts, nparts)
                use nbody_common 
                implicit none
                type(Options), intent(in) :: opt
                integer, intent(in) :: nparts 
                type(Particle), dimension(:), intent(inout) :: parts
        end subroutine
        ! GOL stats protoype
        subroutine nbody_stats(opt, step, parts, nparts)
                use nbody_common 
                implicit none 
                type(Options), intent(in) :: opt
                integer, intent(in) :: nparts, step
                type(Particle), dimension(:), intent(in) :: parts
        end subroutine 
        end interface 
        type(Options) :: opt
        integer :: nparts, nsteps, current_step
        type(Particle), dimension(:), target, allocatable :: parts
        real*8 :: time1, time2

        call getinput(opt)
        call generate_IC(opt, parts, nparts)
        time1 = init_time()
        current_step = 0
        print *, "Running simulation ... "
        do while (current_step .ne. nsteps)
                time2 = init_time()
                call visualise(opt%ivisualisetype, current_step, parts, nparts, opt%period)
                call nbody_stats(opt, current_step, parts, nparts)
                call accel_update(opt, parts, nparts)
                call velocity_update(opt, parts, nparts)
                call position_update(opt, parts, nparts)
                current_step = current_step + 1
                call get_elapsed_time(time2)
                time2 = init_time()
        end do 
        write(*,*) "Finnished NBody simulation"
        call get_elapsed_time(time1);
        deallocate(parts)
end program

! calculate the accelerations on particles 
subroutine accel_update(opt, parts, nparts)
        use nbody_common 
        implicit none
        type(Options), intent(in) :: opt
        integer, intent(in) :: nparts 
        type(Particle), dimension(:), intent(inout) :: parts
        real(8) :: rad2
        integer :: i, j, k

        do i = 1, nparts
                do k = 1,3
                        parts(i)%accel(k) = 0
                end do
        end do
        ! calculate gravitational force 
        do i = 1, nparts
                do j = 1, nparts
                        if (j .ne. i) then 
                                do k = 1,3
                                        ! write gravitational force
                                        ! parts(i)%accel(k) = 
                                end do
                        end if 
                end do 
        end do 
        ! calculate collisional force, other particles must be in the window
        do i = 1, nparts
                do j = 1, nparts
                        if (j .ne. i) then 
                                ! rad = get_dist(parts(i)%position, parts(j)%position)
                                if (rad2 .lt. parts(i)%radius) then 
                                        do k = 1,3
                                                ! write collisional force 
                                                ! parts(i)%accel(k) = 
                                        end do
                                end if 
                        end if 
                end do 
        end do 
end subroutine

! update velocities 
subroutine velocity_update(opt, parts, nparts)
        use nbody_common 
        implicit none
        type(Options), intent(in) :: opt
        integer, intent(in) :: nparts 
        type(Particle), dimension(:), intent(inout) :: parts
        integer :: i, k
        do i = 1, nparts
                do k = 1,3
                        parts(i)%velocity(k) = parts(i)%velocity(k) + parts(i)%accel(k) * opt%time_step
                end do
        end do
end subroutine
! update positions 
subroutine position_update(opt, parts, nparts)
        use nbody_common 
        implicit none
        type(Options), intent(in) :: opt
        integer, intent(in) :: nparts 
        type(Particle), dimension(:), intent(inout) :: parts
        integer :: i, k
        do i = 1, nparts
                do k = 1,3
                        parts(i)%position(k) = parts(i)%position(k) + parts(i)%velocity(k) * opt%time_step
                end do
        end do
end subroutine
! stats protoype
subroutine nbody_stats(opt, step, parts, nparts)
        use nbody_common 
        implicit none 
        type(Options), intent(in) :: opt
        integer, intent(in) :: nparts, step
        type(Particle), dimension(:), intent(in) :: parts
        integer :: i, k

        do i = 1, nparts
                do k = 1,3
                end do 
        end do 
end subroutine 
