program test
    use particle_stats  
    implicit none
    type(particle_type), allocatable :: particles(:)
    integer :: npart, nb, i
    real(4) :: L, dr, s
    real(4), allocatable :: uul(:), uut(:)

    L = 6.2832  ! Replace with actual system size
    nb = 20  ! Replace with actual number of bins

    npart = count_particles('formatted_particles_output.txt')

    allocate(particles(npart))
    call read_particle_data('formatted_particles_output.txt', particles, npart)

    allocate(uul(nb), uut(nb))
    call uu_lt_cond_r(particles, nb, L, dr, uul, uut, s)

    open(unit=20, file='output_stats.txt', status='replace')
    write(20, '(A, F10.5)') 'dr: ', dr
    write(20, '(A, F10.5)') 's: ', s
    do i = 1, nb
        write(20, '(I5, 2F15.8)') i, uul(i), uut(i)
    end do
    close(20)

 
    call sf_lt_cond_r(particles, nb, L, dr, uul, uut, s)

    open(unit=20, file='output_stats1.txt', status='replace')
    write(20, '(A, F10.5)') 'dr: ', dr
    write(20, '(A, F10.5)') 's: ', s
    do i = 1, nb
        write(20, '(I5, 2F15.8)') i, uul(i), uut(i)
    end do

    call rdf_lt_cond_r(particles, nb, L, dr, uut)

    open(unit=20, file='output_stats2.txt', status='replace')
    write(20, '(A, F10.5)') 'dr: ', dr
    write(20, '(A, F10.5)') 's: ', s
    do i = 1, nb
        write(20, '(I5, F15.8)') i, uut(i)
    end do

contains 

function count_particles(filename) result(npart)
    implicit none
    character(len=*), intent(in) :: filename
    integer :: npart, io_stat
    integer :: temp_id
    real(4) :: temp_data(12)

    open(unit=10, file=filename, status='old', iostat=io_stat)
    npart = 0
    if (io_stat /= 0) then
        print *, "Error opening file: ", filename
        stop
    end if

    do
        read(10, *, iostat=io_stat) temp_id, temp_data
        if (io_stat /= 0) exit
        npart = npart + 1
    end do
    close(10)
end function count_particles

subroutine read_particle_data(filename, particles, npart)
    implicit none
    character(len=*), intent(in) :: filename
    type(particle_type), intent(out) :: particles(:)
    integer, intent(in) :: npart
    integer :: i, io_stat
    real(4) :: temp_data(12)

    open(unit=11, file=filename, status='old', iostat=io_stat)
    if (io_stat /= 0) then
        print *, "Error opening file: ", filename
        stop
    end if

    do i = 1, npart
        read(11, *) particles(i)%id, temp_data
        particles(i)%pos = temp_data(1:3)
        particles(i)%fld = temp_data(10:12)
    end do
    close(11)
end subroutine read_particle_data

end program test
