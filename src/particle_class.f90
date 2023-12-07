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

    type :: grid
       integer :: nx, ny, nz  !> Number of cells in each grid direction
       real(4) :: dx, dy, dz  !> Size of cells in each direction
       real(4) :: Lx, Ly, Lz  !> Length of grid in ecah direction 

       integer :: no=1 !> Number of neighbor cells to search over

       integer, allocatable, dimension(:,:,:) :: npic
       integer, allocatable, dimension(:,:,:,:) :: ipic
    end type grid

    type :: part
        integer :: id         !> Particle id number
        integer :: ind(3)     !> Particle cartesian index
        real(4) :: pos(3)     !> Particle position
        real(4) :: vel(3)     !> Particle velocity
        real(4) :: fld(3)     !> Fluid velocity 'seen'
        real(4) :: uf(3)      !> Modeled velocity
        real(4) :: vec(3)     !> Modeled velocity
    end type part

    type :: particles
        integer :: stat_to_compute=stochastic_velocity   !> What two-point quantity to compute
        integer :: npart                                 !> Number of particles
        type(part), allocatable, dimension(:) :: p       !> Vector containing particles
        character(len=80) :: name="UNNAMED"              !> Name of particle group
        character(len=80) :: dir="UNDEFINED"             !> Location of particle data
        character(len=6)  :: suf="000001"                !> Which time step to load in

        type(grid) :: grid                               !> Cartesian-based particle storage
      contains
        procedure :: locate_particles
        procedure :: grid_count
        procedure :: build_ipic
        procedure :: get_npart             !> Routine to get the number of particles
        procedure :: set_vec               !> Set the quantity to compute
        procedure :: read_particle_data    !> Routine to load in particle ensight data
        procedure :: write_particle_data   !> Routine to write particle information to text
    end type particles
 
    interface particles
        procedure :: constructor
    end interface

contains

     function constructor(directory, suffix, name, Lx, Ly, Lz, nx, ny, nz, nover) result(self)
        implicit none
        type(particles) :: self
        integer, intent(in) :: nx
        integer, optional, intent(in) :: ny, nz, nover
        real(4), intent(in) :: Lx
        real(4), optional, intent(in) :: Ly, Lz
        character(len=*), intent(in) :: directory
        character(len=*), intent(in) :: suffix
        character(len=*), optional :: name

        self%dir=trim(adjustl(directory))
        self%suf=trim(adjustl(suffix))

        if (present(name)) self%name=trim(adjustl(name))

        call self%get_npart()
      
        allocate(self%p(1:self%npart))

        call self%read_particle_data()

        construct_grid: block

         self%grid%nx = nx
         if (present(ny)) then
            self%grid%ny = ny
         else
            self%grid%ny = nx
         end if
         if (present(nz)) then
            self%grid%nz = nz
         else
            self%grid%nz = nx
         end if

         self%grid%Lx = Lx
         if (present(Ly)) then
            self%grid%Ly = Ly
         else
            self%grid%Ly = Lx
         end if
         if (present(Lz)) then
            self%grid%Lz = Lz
         else
            self%grid%Lz = Lx
         end if

         self%grid%no = nover

         self%grid%dx = self%grid%Lx / self%grid%nx
         self%grid%dy = self%grid%Ly / self%grid%ny
         self%grid%dz = self%grid%Lz / self%grid%nz

         !> Allocate npic
         allocate(self%grid%npic(self%grid%nx,self%grid%ny,self%grid%nz)); self%grid%npic=0
         !> Locate particles on the grid (get index)
         call self%locate_particles()
         !> Count up particles per grid cell
         call self%grid_count()
         !> Allocate according to the max number of particles that will appear in one grid cell
         call self%build_ipic()
        end block construct_grid
     end function constructor

     subroutine build_ipic(this)
      implicit none
      class(particles), intent(inout) :: this
      integer :: i, ip, jp, kp

      allocate(this%grid%ipic(maxval(this%grid%npic),this%grid%nx,this%grid%ny,this%grid%nz))
      this%grid%npic = 0
      do i = 1, this%npart
         ip = this%p(i)%ind(1)
         jp = this%p(i)%ind(2)
         kp = this%p(i)%ind(3)
         this%grid%npic(ip,jp,kp) = this%grid%npic(ip,jp,kp) + 1
         this%grid%ipic(this%grid%npic(ip,jp,kp),ip,jp,kp) = i
      end do
     end subroutine build_ipic

     subroutine grid_count(this)
      implicit none
      class(particles), intent(inout) :: this
      integer :: i, ip, jp, kp

      do i = 1, this%npart
         ip = this%p(i)%ind(1)
         jp = this%p(i)%ind(2)
         kp = this%p(i)%ind(3)
         this%grid%npic(ip,jp,kp) = this%grid%npic(ip,jp,kp) + 1
      end do
     end subroutine grid_count

     subroutine locate_particles(this)
         implicit none
         class(particles), intent(inout) :: this
         type(part) :: myp
         integer :: i
         do i = 1, this%npart
            myp = this%p(i)
            myp%ind(1) = INT(FLOOR(myp%pos(1) / this%grid%dx)) + 1
            myp%ind(2) = INT(FLOOR(myp%pos(2) / this%grid%dy)) + 1
            myp%ind(3) = INT(FLOOR(myp%pos(3) / this%grid%dz)) + 1
            if (myp%ind(1).gt.this%grid%nx) myp%ind(1) = this%grid.nx
            if (myp%ind(2).gt.this%grid%ny) myp%ind(2) = this%grid.ny
            if (myp%ind(3).gt.this%grid%nz) myp%ind(3) = this%grid.nz
            if (myp%ind(1).lt.1) myp%ind(1) = 1
            if (myp%ind(2).lt.1) myp%ind(2) = 1
            if (myp%ind(3).lt.1) myp%ind(3) = 1
            this%p(i) = myp
         end do
     end subroutine locate_particles

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
