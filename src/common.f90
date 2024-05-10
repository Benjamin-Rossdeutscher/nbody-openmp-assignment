module nbody_common
!---------------------------------------------------------------------
!
!  Common routines and functions for Conway's Game of Life
!
!---------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
#ifdef _OPENMP
    use omp_lib
#endif
#ifdef _MPI 
    include "mpif.h"
#endif
   
    integer, parameter :: NUMVISUAL = 4
    integer, parameter :: VisualiseType_VISUAL_NONE = 0
    integer, parameter :: VisualiseType_VISUAL_ASCII = 1
    integer, parameter :: VisualiseType_VISUAL_ASCII_MESH = 2

    integer, parameter :: NUMICS = 3
    integer, parameter :: ICType_IC_RAND = 0
    integer, parameter :: ICType_IC_ORBIT = 1
    integer, parameter :: ICType_IC_FILE = 2

    integer, parameter :: NUMBOUNDARYCHOICES = 2
    integer, parameter :: BoundaryType_BOUNDARY_NONE = 0
    integer, parameter :: BoundaryType_BOUNDARY_PERIODIC = 1

    integer, parameter :: NUMTIMESTEPSTYPES = 2
    integer, parameter :: TimeStepCrit_Adaptive = 1
    integer, parameter :: TimeStepCrit_Static = 0 
    
    type Options 
        ! assuming position in solar masses, velocity in km/s 
        ! and masses in solar masses
        integer :: nparts, nsteps 
        integer :: iictype
        integer :: ivisualisetype
        integer :: iboundarytype
        integer :: itimestepcrit 
        real(8) :: initial_size
        real(8) :: period
        real(8) :: radiusfac
        real(8) :: time, time_step, time_step_fac 
        real(8) :: munit 
        real(8) :: vunit 
        real(8) :: tunit
        real(8) :: lunit
        real(8) :: vlunittolunit  
        real(8) :: grav_unit
        real(8) :: collision_unit
        integer :: vis_res 
        character(len=2000) :: outfile
        character(len=2000) :: asciifile
    end type Options 

    type Particle
        integer(8) :: ID
        real(8) :: mass
        real(8) :: radius
        real(8) :: position(3)
        real(8) :: velocity(3)
        real(8) :: accel(3)
        integer(8) :: PID
    end type Particle

    contains 

    ! get the memory used as reported by system for process
    subroutine report_memory_usage()
        implicit none
        character(len=2000) :: stat_file
        integer :: ios 
        integer, parameter :: read_unit = 99
        character(len=200) :: line
        character(len=200), dimension(4) :: searchstr
        integer :: i, ii
        integer, dimension(4) :: mem_usage
    ! apple does not provide /proc/self/status file so 
    ! just ignore mem usage 
#ifndef __APPLE__
        write(*,*) "Mem reporting still in dev"
        stat_file = "/proc/self/status";
        searchstr(1) = "VmSize:"
        searchstr(2) = "VmPeak:"
        searchstr(3) = "VmRSS:"
        searchstr(4) = "VmHWM:"
        open (unit = read_unit, file = stat_file, iostat = ios, status = "old")
        if ( ios /= 0 ) then 
            write(*,*) "========================================================= "
            write(*,*) "Couldn't open", stat_file, " for memory usage reading"
            write(*,*) "========================================================= "
            return 
        end if
        do
            read(read_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            do i=1,4
                ii = index(line, searchstr(i))
                ! would need to add parsing of line to extract value of memory 
                if (ii > 0) then
                    mem_usage(ii) = 0
                end if 
            end do 
        end do
        close(read_unit)
        write(*,*) "========================================================= "
        write(*,*) "Memory Usage (VM, Peak VM, RSS, Peak RSS) = \n"
        write(*,*) mem_usage(1), mem_usage(2), mem_usage(3), mem_usage(4), "kilobytes"
        write(*,*) "========================================================= "
#endif
    end subroutine

    ! get some basic timing info
    real(8) function init_time()
        integer, dimension(8) :: value
        call date_and_time(VALUES=value)
        init_time = value(5)*3600.0+value(6)*60.0+value(7)+value(8)/1000.0
        return 
    end function
    ! get the elapsed time relative to start
    subroutine get_elapsed_time(start)
        real(8), intent(in) :: start
        real(8) :: finish, delta
        integer, dimension(8) :: value
        call date_and_time(VALUES=value)
        finish = value(5)*3600.0+value(6)*60.0+value(7)+value(8)/1000.0
        delta = finish - start
        write(*,*) "Elapsed time is ", delta, "s"
    end subroutine

    ! get distances 
    function get_dist(pos1, pos2) result(rad)
        implicit none
        real(8), dimension(:), intent(in) :: pos1, pos2 
        real(8) :: rad
        integer :: i 
        rad = 0
        do i = 1,3
            rad = rad + (pos1(i) - pos2(i))**2.0
        end do
        rad = sqrt(rad)
        return 
    end function get_dist

    ! generate orthonal vector
    ! where x*v = 0 and (x cross v) = rad*vel
    ! randomly choose one of the vectors to be zero so 
    ! long as its related position is  
    ! and then calculate the other two given this constraint
    function get_tagential_velocity(opt, mass, xin) result(v)
        implicit none
        type(Options), intent(in) :: opt
        real(8), intent(in) :: mass 
        real(8), dimension(3), intent(in) :: xin
        real(8), dimension(3) :: x, v
        integer, dimension(3) :: izero
        real(8) :: rad, vel, fac, cross 
        integer :: i, j, k, numzeros, seed

        rad = sqrt(xin(1)*xin(1) + xin(2)*xin(2) + xin(3)*xin(3))
        x = xin / rad
        vel = sqrt(opt%grav_unit/opt%vlunittolunit * mass/rad)

        seed = 4224
        call random_seed(seed)
        numzeros = 0
        do i = 1, 3
            if (x(i) .eq. 0) then 
                numzeros = numzeros + 1
                izero(numzeros) = i
            end if 
        end do

        if (numzeros .eq. 0) then
            ! if all positions are non-zero, then calculate
            ! tangential velocity, randomly setting one of the velocity components as 0
            ! get permutations 
            k = floor(rand()*3.0)
            i = MOD(k+1,3)
            j = MOD(i+1,3)
            i = i+1
            j = j+1 
            k = k+1 
            v(k) = 0
            fac = -x(i)/x(j)
            v(i) = 1.0 / sqrt((-fac*x(k))**2.0 + x(k)**2.0 + (x(i)*fac-x(j))**2.0)
            v(j) = fac*v(i)
        else if (numzeros .gt. 0) then
            v = 0 
            v(izero(1)) = 1.0
        end if 
        v = v * vel 
        return 
    end function get_tagential_velocity
        
    ! periodically wrap positions if necessary 
    subroutine period_wrap(opt, pos1)
        implicit none
        type(Options), intent(in) :: opt
        real(8), dimension(:), intent(inout) :: pos1 
        real(8), dimension(3) :: fac 
        if (opt%period .gt. 0) then
            fac = -floor(pos1 / opt%period)*opt%period
            pos1 = pos1 + fac
        end if 
    end subroutine

    ! periodically wrap difference between two positions if necessary 
    subroutine period_wrap_delta(opt, pos1)
        implicit none
        type(Options), intent(in) :: opt
        real(8), dimension(:), intent(inout) :: pos1 
        integer :: i
        if (opt%period .gt. 0) then
            do i = 1,3
                if (pos1(i) .gt. 0.5*opt%period) then 
                    pos1(i) = pos1(i) - opt%period
                else if (pos1(i) .lt. -0.5*opt%period) then
                    pos1(i) = pos1(i) + opt%period
                end if 
            end do
        end if 
    end subroutine

    ! function that determines the limits of the particle distribution
    function get_particle_limits(parts, nparts) result(limits)
        use, intrinsic :: iso_c_binding
        implicit none
        integer, intent(in) :: nparts
        type(Particle), dimension(:), intent(in) :: parts
        real(8), dimension(:,:), allocatable :: limits
        integer :: i,j

        allocate(limits(3,2))
        do j = 1, 3
            limits(j,1) = 0
            limits(j,2) = 0
        end do 
        do i = 1, nparts
            do j = 1, 3
                if (parts(i)%position(j) .lt. limits(j,1)) then
                    limits(j,1) = parts(i)%position(j)
                end if 
                if (parts(i)%position(j) .gt. limits(j,2)) then
                    limits(j,2) = parts(i)%position(j)
                end if 
            end do 
        end do 
        return
    end function get_particle_limits

    ! function that projects particles on to a density mesh
    function get_mesh(opt, parts, limits) result(mesh)
        use, intrinsic :: iso_c_binding
        implicit none
        type(Options), intent(in) :: opt
        type(Particle), dimension(:), intent(in) :: parts
        real(8), dimension(3,2), intent(in) :: limits
        real(8), dimension(:,:), allocatable :: mesh
        real(8), dimension(3) :: delta
        integer :: ix, iy
        integer :: i,j

        ! get delta to map coordinates to mesh
        do i = 1,3
            delta(i) = (limits(i,2) - limits(i,1))/opt%vis_res
        end do 
        ! allocate and init mesh
        allocate(mesh(opt%vis_res,opt%vis_res))
        do i = 1, opt%vis_res
            do j = 1, opt%vis_res 
                mesh(i,j) = 0 
            end do
        end do 
        ! iterate over particles and construct mesh 
        do i = 1, opt%nparts
            ix = floor((parts(i)%position(1)-limits(1,1))/delta(1)) + 1
            iy = floor((parts(i)%position(2)-limits(2,1))/delta(2)) + 1
            mesh(ix,iy) = mesh(ix,iy) + parts(i)%mass
        end do 
        return
    end function get_mesh

    function normalize_mesh(opt, mesh) result(meshlim)
        implicit none
        type(Options), intent(in) :: opt
        real(8), dimension(:,:), intent(inout) :: mesh
        real(8), dimension(2) :: meshlim
        integer :: i,j
        meshlim(1) = mesh(1,1)
        meshlim(2) = mesh(1,1)
        do i = 1, opt%vis_res
            do j = 1, opt%vis_res
                if (meshlim(1) .gt. mesh(j,i)) then 
                    meshlim(1) = mesh(j,i)
                end if 
                if (meshlim(2) .lt. mesh(j,i)) then 
                    meshlim(2) = mesh(j,i)
                end if
            end do
        end do 
        do i = 1, opt%vis_res
            do j = 1, opt%vis_res
                mesh(j,i) = (mesh(j,i)-meshlim(1))/(meshlim(2)-meshlim(1))
            end do
        end do 
    end function normalize_mesh

        !   ascii visualisation
    subroutine visualise_ascii_mesh(opt, step, parts)
        implicit none
        integer, intent(in) :: step
        type(Particle), dimension(:), intent(in) :: parts
        type(Options), intent(in) :: opt
        real(8), dimension(:,:), allocatable :: limits
        real(8), dimension(:,:), allocatable :: mesh
        real(8), dimension(2) :: meshlim
        real(8) :: delta
        integer :: i,j

        ! construct mesh. First find maximum distance or use period
        if (opt%period .gt. 0) then
            allocate(limits(3,2))
            do i = 1,3
                limits(i,1) = 0
                limits(i,2) = opt%period
            end do
        else 
            limits = get_particle_limits(parts, opt%nparts)
            ! minor correction for mesh drawing as want to deal with
            ! particles at the edges
            do j = 1, 3
                delta = 0.01*(limits(j,2)-limits(j,1))
                limits(j,1) = limits(j,1) - delta
                limits(j,2) = limits(j,2) + delta 
            end do 
    
        end if 
        mesh = get_mesh(opt, parts, limits)
        meshlim = normalize_mesh(opt,mesh)

        write(*,*) "Collisional NBody"
        write(*,*) "Step ", step
        write(*,*) "Mesh limits", limits
        write(*,*) "Mesh value limits", meshlim
        write(*,*) "Normlised mesh"
        write(*,*) "-----------------------"
        do i = 1, opt%vis_res
            do j = 1, opt%vis_res
                write(*,"(1F5.2)", advance="no") mesh(j,i)
            end do 
            write(*,*) ""
        end do 
        write(*,*) "-----------------------"
        deallocate(limits)
        deallocate(mesh)
    end subroutine 

    ! ascii visualisation to a file 
    subroutine visualise_ascii(opt, step, parts)
        implicit none 
        integer, intent(in) :: step
        type(Particle), dimension(:), intent(in) :: parts
        type(Options), intent(in) :: opt

        integer :: i 

        if (step .eq. 0) then 
            open(10, file=opt%asciifile, access="sequential")
        else 
            open(10, file=opt%asciifile, access="append")
        end if
        do i = 1, opt%nparts
            write(10,*) step, opt%time, parts(i)
        end do 
        close(10)
    end subroutine 

    ! no visualisation
    subroutine visualise_none(step)
        implicit none 
        integer, intent(in) :: step
        write(*,*) "Collisional NBody, Step ", step
    end subroutine 

    ! visualisation routine
    subroutine visualise(opt, step, parts)
        implicit none 
        type(Options), intent(in) :: opt
        integer, intent(in) :: step
        type(Particle), dimension(:), intent(in) :: parts
        if (opt%ivisualisetype .eq. VisualiseType_VISUAL_ASCII) then 
            call visualise_ascii(opt, step, parts)
        else if (opt%ivisualisetype .eq. VisualiseType_VISUAL_ASCII_MESH) then 
            call visualise_ascii_mesh(opt, step, parts)
        else  
            call visualise_none(step)
        end if 

    end subroutine

    ! generate random IC
    function generate_rand_IC(opt) result(parts)
        use, intrinsic :: iso_c_binding
        use :: omp_lib
        implicit none
        type(Particle), dimension(:), allocatable :: parts
        type(Options), intent(in) :: opt 
        integer :: i, j
        real(8), dimension(3) :: x
        real(8) :: mval, time1
        integer :: seed
       
        ! Initialize variables
        allocate(parts(opt%nparts))
        
        ! Printing information
        print *, "Particle class requires ", sizeof(parts(1)), " bytes"
        print *, "Total number of particles ", opt%nparts
        print *, "Total memory required is ", opt%nparts*sizeof(parts(1))/1024./1024./1024., " GB"
        print *, "Generating random particle positions and velocities"
        time1 = init_time()
        
        seed = 4322
        call random_seed(seed)
        do i = 1, opt%nparts
            do j = 1,3
                x(j)=rand()
            end do
            parts(i)%position = x * opt%initial_size
        end do

        ! Set mass, ID, and PID for all particles
        do i = 1, opt%nparts
            mval=rand()
            parts(i)%mass = mval * opt%munit
            parts(i)%ID = i
            parts(i)%PID = i
            mval=rand()
            parts(i)%radius = opt%initial_size * opt%radiusfac*mval 
        end do
        print *, "Done generating positions in "
        call get_elapsed_time(time1)
    end function generate_rand_IC

    ! uses the randomly generated positions and then 
    ! calculates velocities as if particles are orbiting 
    ! centre-of-mass 
    ! using largest particle to define orbital plane 
    function generate_orbit_IC(opt) result(parts)
        use, intrinsic :: iso_c_binding
        use :: omp_lib
        implicit none
        type(Particle), dimension(:), allocatable :: parts
        type(Options), intent(in) :: opt 
        integer :: i 
        real(8), dimension(3) :: cmx, x
        real(8) :: mtot, time1

        time1 = init_time()

        parts = generate_rand_IC(opt)
        do i = 1, opt%nparts
            mtot = mtot + parts(i)%mass
            cmx = cmx + parts(i)%position * parts(i)%mass
        end do
        cmx = cmx / mtot

        do i = 1, opt%nparts
            x = cmx - parts(i)%position 
           parts(i)%velocity = get_tagential_velocity(opt, mtot - parts(i)%mass, x) 
        end do
        
        call get_elapsed_time(time1)
    end function generate_orbit_IC

    ! generate IC
    subroutine generate_IC(opt, parts)
        implicit none 
        type(Options), intent(in) :: opt
        type(Particle), dimension(:), intent(out), allocatable :: parts
        if (opt%iictype .eq. ICType_IC_RAND) then 
            parts = generate_rand_IC(opt)
        end if 
        if (opt%iictype .eq. ICType_IC_ORBIT) then 
            parts = generate_orbit_IC(opt)
        end if 
        !call write_data(parts)
    end subroutine 


    ! Report run state 
    subroutine run_state(opt)
        implicit none 
        type(Options), intent(inout) :: opt
        integer :: count
        integer*8 :: nbytes
        real*8 :: memfootprint
        type(Particle) :: p
        integer :: mpi_size, omp_size, ierror

        nbytes = sizeof(p) * opt%nparts
        memfootprint = real(nbytes)/1024.0/1024.0/1024.0
        write(*,*) "========================================================= "
        write(*,*) " NBody Running with following "
        write(*,*) "========================================================= "
        write(*,*) "Requesting number of particles ", opt%nparts
        write(*,*) "which requires", memfootprint, " GB "
        write(*,*) "Number of steps ", opt%nsteps
        write(*,*) "Physical parameters : "
        write(*,*) "munit in solar masses", opt%munit
        write(*,*) "lunit in kpc", opt%lunit
        write(*,*) "vunit in km/s", opt%vunit
        write(*,*) "tunit in s", opt%tunit
        write(*,*) "time step", opt%time_step
        write(*,*) "grav in solar masses kpc^2 km/s^2", opt%grav_unit
        write(*,*) "collisional force unit in grav units", opt%collision_unit
        write(*,*) "Boundary rule type ",opt%iboundarytype
        write(*,*) "IC type ",opt%iictype
        write(*,*) "Time-step criterion ",opt%itimestepcrit
        write(*,*) "Visualization type ",opt%ivisualisetype
        write(*,*) "========================================================= "
    
#ifdef _MPI
        call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierror)
        write(*,*)a "========================================================"
        write(*,*) "Running with MPI"
        write(*,*) "Comm World contains ", mpi_size, " processes "
        write(*,*) "========================================================="
#endif
#ifdef _OPENMP
        omp_size = omp_get_max_threads()
        write(*, *) "========================================================"
        write(*, *) "Running with OpenMP, version ", _OPENMP
        write(*, *) "Maximum number of threads ", omp_size
        write(*, *) "========================================================"
#endif
    end subroutine

    ! UI
    subroutine getinput(opt)
        implicit none 
        type(Options), intent(inout) :: opt
        character(len=2000) :: cmd
        character(len=32) :: arg 
        character(len=2000) :: outfilename, asciifilename
        character(len=32) :: value
        integer :: count
        ! get the commands passed and the number of args passed 
        call get_command(cmd)
        count = command_argument_count()
        if (count .lt. 2) then 
            write(*,*) "Usage: <number of particles> <nsteps> "
            write(*,*) "[<Boundary type> <IC type> <Time Step Criterion> <Visualisation type> <Vis res>]"
            call exit();
        end if 
        
        outfilename = "nbody-data.txt"
        asciifilename = "nbody-asciivis.txt"
        call get_command_argument(1,arg)
        read(arg,*) opt%nparts
        call get_command_argument(2,arg)
        read(arg,*) opt%nsteps
        opt%vis_res=100
        opt%ivisualisetype = VisualiseType_VISUAL_NONE
        opt%iictype = ICType_IC_RAND
        opt%iboundarytype = BoundaryType_BOUNDARY_NONE
        opt%itimestepcrit = TimeStepCrit_Adaptive
        opt%initial_size = 1.0 
        opt%period = 0.0
        ! and radius is 1/10th of the initial size of the box 
        opt%radiusfac = 0.1/opt%nparts**(1.0/3.0)
        opt%munit = 1.0
        opt%vunit = 1.0 
        opt%lunit = 1.0 
        ! conversion from km to pc 
        opt%vlunittolunit = 3.24078e-14
        ! grav in km/s^2 pc^2 / solar masses
        opt%grav_unit = 4.30241002e-3 * opt%vlunittolunit
        ! make time unit related to graviational unit in seconds
        opt%tunit = 1.0/sqrt((opt%grav_unit *opt%vlunittolunit * opt%munit/opt%initial_size**3.0)) 
        ! time step should be some fraction of the dynamical time of a close encounter the system which is related 
        ! to tunit by radiusfac **1.5
        opt%time = 0
        opt%time_step_fac = 0.1
        ! opt%time_step = opt%tunit * opt%radiusfac**1.5 * opt%time_step_fac
        opt%time_step = opt%tunit * opt%time_step_fac
        ! make collisional (repulsive) force 10 times stronger than gravity 
        ! and is in units of gravitational forces 
        opt%collision_unit = 2.0
        if (count .ge. 3) then
            call get_command_argument(3,arg)
            read(arg,*) opt%iboundarytype
        end if 
        if (opt%iboundarytype .eq. BoundaryType_BOUNDARY_PERIODIC) then 
            opt%period = opt%initial_size
        end if 
        if (count .ge. 4) then 
            call get_command_argument(4, arg)
            read(arg,*) opt%iictype
        end if 
        if (count .ge. 5) then 
            call get_command_argument(5, arg)
            read(arg,*) opt%itimestepcrit
        end if 
        if (count .ge. 6) then 
            call get_command_argument(6, arg)
            read(arg,*) opt%ivisualisetype
        end if 
        if (count .ge. 7) then 
            call get_command_argument(7, arg)
            read(arg,*) opt%vis_res
        end if 
        if (opt%nparts .le. 0) then 
            value = CHAR(opt%nparts)
            write(*,*) "Requested ", value, " particles"
            write(*,*) "Invalid number of particles."
            call exit(1)
        end if 
        if (opt%nsteps .le. 0) then
            value = CHAR(opt%nsteps)
            write(*,*) "Requested ", value, " steps"
            write(*,*) "Invalid number of steps."
            call exit(1)
       end if 
       opt%outfile = outfilename
       opt%asciifile = asciifilename 
       call run_state(opt)
    end subroutine 

    
end module
