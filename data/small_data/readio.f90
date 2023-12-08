module particle_module
    implicit none

    type :: particle_type
        integer :: id
        real(4) :: pos(3)
        real(4) :: vel(3)
        real(4) :: fld(3)
        real(4) :: uf(3)
    end type particle_type

contains

     function get_num_particles(particle_file) result(npart)
        character(len=*), intent(in) :: particle_file
        integer :: npart, io_stat
        integer(4) :: num_particles

        open(unit=10, file=particle_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "Error opening file: ", particle_file
            stop
        end if

        ! Skip the first 240 bytes
        read(10, pos=241, iostat=io_stat) num_particles

        if (io_stat /= 0) then
            print *, "Error reading from file: ", particle_file
            stop
        end if

        npart = num_particles
        close(10)
    end function get_num_particles
    subroutine read_particle_data(npart, fld_file, id_file, particle_file, uf_file, vel_file, particles)
        implicit none
        integer, intent(in) :: npart
        character(len=*), intent(in) :: fld_file, id_file, particle_file, uf_file, vel_file
        character(len=80) :: tmp
        character(len=244+(4*npart)) :: tmp1
        real(4) :: idtmp
        type(particle_type), intent(out) :: particles(npart)
        integer :: i, io_stat

        open(unit=11, file=id_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=12, file=particle_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=13, file=vel_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=14, file=fld_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=15, file=uf_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        print*,npart
           
        read(11) tmp 
        read(12) tmp1
        read(13) tmp 
        read(14) tmp 
        read(15) tmp
 
        do i = 1, npart
            read(11) idtmp
            particles(i)%id=int(idtmp)
            read(12) particles(i)%pos
            read(13) particles(i)%vel
            read(14) particles(i)%fld
            read(15) particles(i)%uf
        end do

        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
    end subroutine read_particle_data

    subroutine write_particle_data(particles, np, output_file)
        implicit none
        type(particle_type), intent(in) :: particles(:)
        integer, intent(in) :: np
        character(len=*), intent(in) :: output_file
        integer :: i

        open(unit=20, file=output_file, status='replace')

        do i = 1, np
            write(20, '(I5, 3F15.7, 3F15.7, 3F15.7, 3F15.7)') particles(i)%id, &
                particles(i)%pos, particles(i)%vel, particles(i)%fld, particles(i)%uf
        end do

        close(20)
    end subroutine write_particle_data

end module particle_module

program main
    use particle_module
    implicit none

    type(particle_type), allocatable :: particles(:)
    integer :: npart

    ! Determine the number of particles based on the ID file
    npart = get_num_particles('particle.000148')

    ! Allocate particles array
    allocate(particles(npart))

    ! Read particle data from files
    call read_particle_data(npart, 'fld.000148', 'id.000148', 'particle.000148', 'uf.000148', 'vel.000148', particles)

    ! Write particle data to a formatted text file
    call write_particle_data(particles, npart, 'formatted_particles_output.txt')

end program main

