module particle_class
    implicit none
    private

    public :: part
    public :: particles

    !> Used by set_vec to set the computed two-point quantity
    integer, parameter, public :: fluid_velocity=1
    integer, parameter, public :: particle_velocity=2
    integer, parameter, public :: stochastic_velocity=3
    integer, parameter, public :: total_fluid_velocity=4

    type :: part
        integer :: id         !> Particle id number
        integer :: ind(3)     !> Particle cartesian index
        real(4) :: pos(3)     !> Particle position
        real(4) :: vel(3)     !> Particle velocity
        real(4) :: fld(3)     !> Fluid velocity 'seen'
        real(4) :: uf(3)      !> Modeled velocity
        real(4) :: vec(3)      !> Modeled velocity
    end type part

    type :: particles
        integer :: stat_to_compute=stochastic_velocity   !> What two-point quantity to compute
        integer :: npart                                 !> Number of particles
        type(part), allocatable, dimension(:) :: p       !> Vector containing particles
        character(len=80) :: name="UNNAMED"              !> Name of particle group
        character(len=80) :: dir="UNDEFINED"             !> Location of particle data
        character(len=6)  :: suf="000001"                !> Which time step to load in
      contains
        procedure :: get_npart             !> Routine to get the number of particles
        procedure :: set_vec               !> Set the quantity to compute
        procedure :: read_particle_data    !> Routine to load in particle ensight data
        procedure :: write_particle_data   !> Routine to write particle information to text
    end type particles
 
    interface particles
        procedure :: constructor
    end interface

contains

     function constructor(directory, suffix, name) result(self)
        implicit none
        type(particles) :: self
        character(len=*), intent(in) :: directory
        character(len=*), intent(in) :: suffix
        character(len=*), optional :: name

        self%dir=trim(adjustl(directory))
        self%suf=trim(adjustl(suffix))

        if (present(name)) self%name=trim(adjustl(name))

        call self%get_npart()
      
        allocate(self%p(1:self%npart))

        call self%read_particle_data()

     end function constructor

     !> Read particle.* ensight file for number of particles
     subroutine get_npart(this)
        class(particles), intent(inout) :: this
        character(len=240) :: particle_file
        integer :: io_stat

        particle_file = trim(adjustl(this%dir)) // "particle." // this%suf

        open(unit=10, file=particle_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        if (io_stat /= 0) then
            print *, "[get_npart] Error opening file: ", particle_file
            stop
        end if

        ! Skip the first 240 bytes
        read(10, pos=241, iostat=io_stat) this%npart

        if (io_stat /= 0) then
            print *, "[get_npart] Error reading from file: ", particle_file
            stop
        end if

        close(10)
      end subroutine get_npart


     subroutine set_vec(this)
      implicit none
      class(particles), intent(inout) :: this
      integer :: i
      select case(this%stat_to_compute)
      case (fluid_velocity)
         do i = 1,this%npart
            this%p(i)%vec=this%p(i)%fld
         end do
      case (particle_velocity)
         do i = 1,this%npart
            this%p(i)%vec=this%p(i)%vel
         end do
      case (stochastic_velocity)
         do i = 1,this%npart
            this%p(i)%vec=this%p(i)%uf
         end do
      case (total_fluid_velocity)
         do i = 1,this%npart
            this%p(i)%vec=this%p(i)%uf+this%p(i)%fld
         end do
      end select
     end subroutine set_vec


    subroutine read_particle_data(this)
        implicit none
        class(particles), intent(inout) :: this
        character(len=240) :: fld_file, id_file, particle_file, uf_file, vel_file
        character(len=80) :: buf
        character(len=244+(4*this%npart)) :: larger_buf
        real(4) :: idtmp
        integer :: i, io_stat

        particle_file = trim(adjustl(this%dir)) // "particle." // this%suf
        fld_file      = trim(adjustl(this%dir)) // "fld." // this%suf
        id_file       = trim(adjustl(this%dir)) // "id." // this%suf
        uf_file       = trim(adjustl(this%dir)) // "uf." // this%suf
        vel_file      = trim(adjustl(this%dir)) // "vel." // this%suf

        open(unit=11, file=id_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=12, file=particle_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=13, file=vel_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=14, file=fld_file, form='unformatted', access='stream', status='old', iostat=io_stat)
        open(unit=15, file=uf_file, form='unformatted', access='stream', status='old', iostat=io_stat)
           
        !> Skip header (position file has larger header)
        read(11) buf
        read(12) larger_buf
        read(13) buf 
        read(14) buf 
        read(15) buf
 
        do i = 1, this%npart
            read(11) idtmp; this%p(i)%id=int(idtmp)
            read(12) this%p(i)%pos
            read(13) this%p(i)%vel
            read(14) this%p(i)%fld
            read(15) this%p(i)%uf

            !> Set default stat to uf
            this%p(i)%vec = this%p(i)%uf
        end do

        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
    end subroutine read_particle_data

    subroutine write_particle_data(this, output_file)
        implicit none
        class(particles), intent(inout) :: this
        character(len=*), intent(in) :: output_file
        integer :: i

        open(unit=20, file=output_file, status='replace')

        do i = 1, this%npart
            write(20, '(I5, 3F15.7, 3F15.7, 3F15.7, 3F15.7)') this%p(i)%id, &
                this%p(i)%pos, this%p(i)%vel, this%p(i)%fld, this%p(i)%uf
        end do

        close(20)
    end subroutine write_particle_data

end module particle_class