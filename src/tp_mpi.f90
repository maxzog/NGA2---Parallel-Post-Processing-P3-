module twopoint_mpi
   use particle_class
   implicit none
   include "mpif.h"

   private

   real(4) :: PI = 3.14159265

   public :: mpi_stats

   type mpi_stats
      integer :: nb = 32
      real(4) :: L = 6.2832
      real(4) :: dr
      real(4) :: rmax
      real(4) :: mean = 0.0
      real(4) :: var = 0.0
      real(4) :: dt = 0.0
      real(4) :: bw = 0.5
      real(4), allocatable, dimension(:) :: uul, uut, sfl, sft, rdf, ac
      integer, allocatable, dimension(:) :: c
      real(4), allocatable, dimension(:) :: drift, covar
      integer, allocatable, dimension(:) :: driftc
      real(4), allocatable, dimension(:,:) :: reldrift
      integer, allocatable, dimension(:,:) :: reldriftc
      integer :: imin, imax, jmin, jmax
      integer :: nstep = 1, step = 1
      integer :: comm, nproc, rank
      integer :: ubins
   contains
      procedure :: decomp
      procedure :: compute_var
      procedure :: compute_uu
      procedure :: compute_sf
      procedure :: compute_rdf
      procedure :: compute_ac
      procedure :: compute_dispersion
      procedure :: compute_relative_dispersion

      procedure :: infer_drift
      procedure :: infer_covar
      procedure :: infer_relative_drift

      procedure :: write_uu
      procedure :: write_sf
      procedure :: write_ac
      procedure :: write_rdf
      procedure :: write_covar
      procedure :: write_drift
      procedure :: write_relative_drift
   end type mpi_stats

   interface mpi_stats
      procedure :: constructor
   end interface mpi_stats

contains

   function constructor(numbins, length, nstep, dt, ubins) result(self)
      implicit none
      type(mpi_stats) :: self
      integer, optional :: numbins
      integer, optional :: ubins
      integer, optional :: nstep
      real(4), optional :: length
      real(4), optional :: dt

      if (present(numbins)) self%nb = numbins
      if (present(ubins)) self%ubins = ubins
      if (present(length)) self%L = length
      if (present(nstep)) self%nstep = nstep
      if (present(dt)) self%dt = dt

      allocate (self%uul(1:self%nb)); self%uul   = 0.0
      allocate (self%uut(1:self%nb)); self%uut   = 0.0
      allocate (self%sfl(1:self%nb)); self%sfl   = 0.0
      allocate (self%sft(1:self%nb)); self%sft   = 0.0
      allocate (self%rdf(1:self%nb)); self%rdf   = 0.0
      allocate (self%covar(1:self%nb)); self%covar = 0.0
      allocate (self%ac(1:self%nstep));  self%ac = 0.0
      allocate (self%drift( 1:self%ubins)); self%drift  = 0.0
      allocate (self%driftc(1:self%ubins)); self%driftc = 0
      allocate (self%reldrift(1:self%ubins,1:self%nb)); self%reldrift=0.0
      allocate (self%reldriftc(1:self%ubins,1:self%nb)); self%reldriftc=0

      allocate (self%c(1:self%nb)); self%c = 0

      self%rmax = norm2([0.5*self%L, 0.5*self%L, 0.5*self%L])
      self%dr = self%rmax/self%nb
   end function constructor

   subroutine decomp(this, length, nproc, rank, s, e)
      implicit none
      class(mpi_stats), intent(inout) :: this
      integer, intent(in) :: nproc, rank, length
      integer, intent(inout) :: s, e
      integer :: nlocal, deficit

      nlocal = length/nproc
      s = rank*nlocal + 1
      deficit = mod(length, nproc)

      if (rank .lt. deficit) then
         s = s + rank
      else
         s = s + deficit
      end if

      if (rank .lt. deficit) nlocal = nlocal + 1
      e = s + nlocal - 1
      if (e .gt. length .or. rank .eq. nproc - 1) e = length
   end subroutine decomp

   subroutine compute_var(this, parts)
      type(particles), intent(in) :: parts
      class(mpi_stats), intent(inout) :: this
      real(4) :: mean
      integer :: i, ierr

      !> Reset
      this%mean = 0.0; this%var = 0.0

      !> Compute mean (isotropic)
      do i = this%imin, this%imax
         this%mean = this%mean + sum(parts%p(i)%vec)/3
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%mean, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%mean = this%mean/parts%npart

      !> Compute variance
      do i = this%imin, this%imax
         this%var = this%var + sum((parts%p(i)%vec - this%mean)**2)/3
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%var, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%var = this%var/parts%npart
   end subroutine compute_var

   subroutine compute_uu(this, parts)
      type(particles), intent(in) :: parts
      class(mpi_stats), intent(inout) :: this
      integer :: i, j, ir, ierr
      real(4) :: r(3), rll(3), rt2(3)

      !> Reset
      this%c = 0; this%uul = 0.0; this%uut = 0.0

      do i = this%imin, this%imax
         do j = i + 1, this%imax
            call par_perp_u(this, parts%p(i), parts%p(j), rll, rt2)
            r = parts%p(j)%pos - parts%p(i)%pos
            ir = MIN(floor(norm2(r)/this%dr) + 1, this%nb)
            !! this%uul(ir) = this%uul(ir) + dot_product(parts%p(i)%vec, rll)*dot_product(parts%p(j)%vec, rll)
            this%uul(ir) = this%uul(ir) + dot_product(parts%p(i)%vec, parts%p(j)%vec) / 3.0
            this%uut(ir) = this%uut(ir) + dot_product(parts%p(i)%vec, rt2)*dot_product(parts%p(j)%vec, rt2)
            this%c(ir) = this%c(ir) + 1
         end do
      end do

      do i = 1, this%nb
         !> Prevent division by zero
         if (this%c(i) == 0) then
            this%c(i) = 1
         end if
         !> Normalize bins
         this%uul(i) = this%uul(i)/this%c(i)
         this%uut(i) = this%uut(i)/this%c(i)
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, this%uul(:), this%nb, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%uut(:), this%nb, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%c(:), this%nb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%uul = this%uul/this%nproc
      this%uut = this%uut/this%nproc
   end subroutine compute_uu

   subroutine compute_sf(this, parts)
      type(particles), intent(in) :: parts
      class(mpi_stats), intent(inout) :: this
      integer :: i, j, ir, ierr
      real(4) :: r(3), rll(3), rt2(3)

      !> Reset
      this%c = 0; this%sfl = 0.0; this%sft = 0.0

      do i = this%imin, this%imax
         do j = i + 1, this%imax
            call par_perp_u(this, parts%p(i), parts%p(j), rll, rt2)
            r = parts%p(j)%pos - parts%p(i)%pos
            ir = MIN(floor(norm2(r)/this%dr) + 1, this%nb)
            this%sfl(ir) = this%sfl(ir) + dot_product(parts%p(j)%vec - parts%p(i)%vec, rll) &
                          & *dot_product(parts%p(j)%vec - parts%p(i)%vec, rll)
            this%sft(ir) = this%sft(ir) + dot_product(parts%p(j)%vec - parts%p(i)%vec, rt2) &
                          & *dot_product(parts%p(j)%vec - parts%p(i)%vec, rt2)
            this%c(ir) = this%c(ir) + 1
         end do
      end do

      do i = 1, this%nb
         !> Prevent division by zero
         if (this%c(i) == 0) then
            this%c(i) = 1
         end if
         !> Normalize bins
         this%sfl(i) = this%sfl(i)/this%c(i)
         this%sft(i) = this%sft(i)/this%c(i)
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, this%sfl(:), this%nb, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%sft(:), this%nb, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, this%c(:), this%nb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%sfl = this%sfl/this%nproc
      this%sft = this%sft/this%nproc
   end subroutine compute_sf

   subroutine compute_rdf(this, parts)
      type(particles), intent(in) :: parts
      class(mpi_stats), intent(inout) :: this
      integer :: i, j, ir, ierr
      real(4) :: V, Vshell, myN, N, rho, r(3)

      !> Reset
      this%rdf = 0.0
      myN = 0.0; N = 0.0

      do i = this%imin, this%imax
         do j = i + 1, this%imax
            r = get_minr(this, parts%p(i), parts%p(j))
            ir = MIN(floor(norm2(r)/this%dr) + 1, this%nb)
            this%rdf(ir) = this%rdf(ir) + 1.0
            myN = myN + 1.0
         end do
      end do

      !call MPI_ALLREDUCE(myN, N, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

      !> Prepare RDF normalization
      V = this%L**3
      !N   = (this%imax - this%imin) * (this%imax - this%imin) / 2.0
      rho = myN/V

      do i = 1, this%nb
         Vshell = 4.0/3.0*PI*((this%dr*i)**3 - (this%dr*(i - 1))**3)
         this%rdf(i) = this%rdf(i)/(Vshell*rho)
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, this%rdf(:), this%nb, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%rdf = this%rdf/this%nproc
   end subroutine compute_rdf

   subroutine compute_dispersion(this, partsn, partsm)
      type(particles), intent(in) :: partsn,partsm
      class(mpi_stats), intent(inout) :: this
      integer :: i, ierr
      real(4) :: tmpdis=0.0
      real(4) :: npi

      do i = 1, partsn%npart
         tmpdis=tmpdis+sum((partsm%p(i)%pos - partsn%p(i)%pos)**2)             
      end do

      !call MPI_ALLREDUCE(MPI_IN_PLACE, tmpdis, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%ac(this%step)=tmpdis/partsn%npart
   end subroutine compute_dispersion

   subroutine compute_relative_dispersion(this, partsn, partsm)
      type(particles), intent(in) :: partsn,partsm
      class(mpi_stats), intent(inout) :: this
      integer :: i, j, ierr
      real(4) :: tmpdis=0.0

      do i = this%imin, this%imax
         do j = i + 1, this%imax
            tmpdis=tmpdis+sum(((partsm%p(j)%pos-partsm%p(i)%pos)-(partsn%p(j)%pos-partsn%p(i)%pos))**2)  
         end do
      end do

      tmpdis=tmpdis/((this%imax-this%imin)*(this%imax-this%imin+1.0))*2.0
      call MPI_ALLREDUCE(MPI_IN_PLACE, tmpdis, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%ac(this%step)=tmpdis/this%nproc
   end subroutine compute_relative_dispersion

   subroutine infer_drift(this, partsn, partsm)
      type(particles), intent(in) :: partsn,partsm
      class(mpi_stats), intent(inout) :: this
      real(4) :: umin, umax
      real(4), dimension(3) :: du
      integer :: i, j, ierr, ind

      umax=0.0
      do i = 1, partsn%npart
         umax=maxval(abs(partsm%p(i)%vec)) 
      end do
      this%bw=2.0*umax/this%ubins
      do i = 1, partsn%npart
         du = partsn%p(i)%vec - partsm%p(i)%vec
         do j = 1, 3
            ind = ceiling((partsm%p(i)%vec(j)+umax) / this%bw)
            if (ind.gt.this%ubins.or.ind.lt.0) cycle
            this%drift(ind) = this%drift(ind) + du(j)/this%dt
            this%driftc(ind) = this%driftc(ind) + 1
         end do
      end do

      do i = 1, this%ubins
         if (this%driftc(i).gt.0) then
            this%drift(i) = this%drift(i) / this%driftc(i) 
         end if
      end do
   end subroutine infer_drift

   subroutine infer_covar(this, partsn, partsm)
      type(particles), intent(in) :: partsn,partsm
      class(mpi_stats), intent(inout) :: this
      real(4) :: umin, umax
      real(4), dimension(3) :: dui,duj,r
      integer :: i, j, ir, ierr

      this%c=0
      do i = this%imin, this%imax
         do j = i + 1, this%imax
            dui = partsn%p(i)%vec - partsm%p(i)%vec
            duj = partsn%p(j)%vec - partsm%p(j)%vec
            r   = partsm%p(j)%pos - partsm%p(i)%pos
            ir  = MIN(floor(norm2(r)/this%dr) + 1, this%nb)
            this%covar(ir) = this%covar(ir) + sum(dui*duj)/3.0
            this%c(ir) = this%c(ir) + 1
         end do
      end do

      do i = 1, this%nb
         if (this%c(i).gt.0) then
            this%covar(i) = this%covar(i) / this%c(i) / this%dt 
         end if
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, this%covar(:), this%nb, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%covar = this%covar/this%nproc
   end subroutine infer_covar

   subroutine infer_relative_drift(this, partsn, partsm)
      type(particles), intent(in) :: partsn,partsm
      class(mpi_stats), intent(inout) :: this
      real(4) :: umax, dur, rn, durn, durm
      real(4), dimension(3) :: du, r, rh
      integer :: i, j, ir, iu, ierr
   
      !> Get maximum rel vel for data storage
      umax=0.0
      do i = this%imin, this%imax
         do j = i + 1, this%imax
            du  = partsm%p(j)%vec - partsm%p(i)%vec
            r   = partsm%p(j)%pos - partsm%p(i)%pos
            rn  = norm2(r)
            rh  = r / rn
            dur = dot_product(du, rh)
            if (abs(dur).gt.umax) umax = abs(dur)
         end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, umax, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
      this%bw=2.0*umax/this%ubins

      this%c=0
      do i = this%imin, this%imax
         do j = i + 1, this%imax
            du = partsn%p(j)%vec - partsn%p(i)%vec
            r  = partsn%p(j)%pos - partsn%p(i)%pos
            rn = norm2(r)
            rh = r / rn
            durn = dot_product(du, rh)
            
            du = partsm%p(j)%vec - partsm%p(i)%vec
            r  = partsm%p(j)%pos - partsm%p(i)%pos
            rn = norm2(r)
            rh = r / rn
            durm = dot_product(du, rh)

            r  = partsm%p(j)%pos - partsm%p(i)%pos
            ir = MIN(floor(norm2(r)/this%dr) + 1, this%nb)
            iu = ceiling((durm + umax)/this%bw)
            this%reldrift(iu,ir) = this%reldrift(iu,ir) + durn-durm
            this%reldriftc(iu,ir) = this%reldriftc(iu,ir) + 1
         end do
      end do

      do i = 1, this%ubins
         do j = 1, this%nb
            if (this%reldriftc(i,j).gt.0) then
               this%reldrift(i,j) = this%reldrift(i,j) / this%reldriftc(i,j) / this%dt 
            end if
         end do
      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, this%reldrift(:,:), this%nb*this%ubins, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%reldrift = this%reldrift/this%nproc
   end subroutine infer_relative_drift

   subroutine compute_ac(this, partsn, partsm)
      implicit none
      class(mpi_stats), intent(inout) :: this
      type(particles), intent(in) :: partsn
      type(particles), intent(in) :: partsm
      integer :: i, ierr
      real(4) :: meanm, tmpac

      meanm = 0.0
      !> Compute mean (isotropic)
      do i = this%imin, this%imax
         meanm = meanm + sum(partsm%p(i)%vec)/3
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, meanm, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      meanm = meanm/partsm%npart

      tmpac=0.0
      ! Compute two-time variance
      do i = this%imin, this%imax
         tmpac = tmpac + dot_product(partsn%p(i)%vec - this%mean, partsm%p(i)%vec - meanm)/3
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, tmpac, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      this%ac(this%step) = tmpac/partsm%npart
   end subroutine compute_ac

   function get_minr(this, p, q) result(r)
      type(mpi_stats) :: this
      type(part), intent(in) :: p, q
      real(4) :: r(3)
      r = q%pos - p%pos
      r(1) = r(1) - this%L*nint(r(1)/this%L)
      r(2) = r(2) - this%L*nint(r(2)/this%L)
      r(3) = r(3) - this%L*nint(r(3)/this%L)
   end function get_minr

   function cross_product(a, b) result(cross_val)
      real(4), intent(in) :: a(3), b(3)
      real(4) :: cross_val(3)
      cross_val(1) = a(2)*b(3) - a(3)*b(2)
      cross_val(2) = a(3)*b(1) - a(1)*b(3)
      cross_val(3) = a(1)*b(2) - a(2)*b(1)
   end function cross_product

   subroutine par_perp_u(this, p, q, rll, rt2)
      implicit none
      type(mpi_stats) :: this
      type(part), intent(in) :: p, q
      real(4), intent(out) :: rll(3), rt2(3)
      real(4) :: r(3), rt1(3)
      integer :: i

      r = get_minr(this, p, q)
      rll = r/norm2(r)
      rt1 = cross_product(rll, p%vec)/norm2(cross_product(rll, p%vec))
      rt2 = cross_product(rt1, rll)/norm2(cross_product(rt1, rll))
   end subroutine par_perp_u

   !> I/O routines
   subroutine write_uu(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
         write (20, '(I5, 2F15.8)') i, this%uul(i), this%uut(i)
      end do
      close (20)
   end subroutine write_uu

   subroutine write_sf(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
         write (20, '(I5, 2F15.8)') i, this%sfl(i), this%sft(i)
      end do
      close (20)
   end subroutine write_sf

   subroutine write_rdf(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i
      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
         write (20, '(I5, F15.8)') i, this%rdf(i)
      end do
      close (20)
   end subroutine write_rdf

   subroutine write_covar(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i
      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
         write (20, '(I5, F15.8)') i, this%covar(i)
      end do
      close (20)
   end subroutine write_covar

   subroutine write_ac(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'dt: ', this%dt
      do i = 1, this%nstep
         write (20, '(I5, 2F15.8)') i, this%ac(i)
      end do
      close (20)
   end subroutine write_ac

   subroutine write_drift(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'du: ', this%bw
      do i = 1, this%ubins
         write (20, '(I5, 2F15.8)') i, this%drift(i)
      end do
      close (20)
   end subroutine write_drift
   
   subroutine write_relative_drift(this, outfile)
      implicit none
      class(mpi_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i, j

      open (unit=20, file=outfile, status='replace')
      write (20, '(A, F10.5)') 'du: ', this%bw
      write (20, '(A, F10.5)') 'dr: ', this%dr
      ! Iterate through the 2D array and write to file
      do i = 1, this%ubins
          ! Begin writing a new line for the new row
          write(20, '(F10.6)', advance='no') this%reldrift(i, 1)
          do j = 2, this%nb
              write(20, '(A)', advance='no') ' '
              write(20, '(F10.6)', advance='no') this%reldrift(i, j)
          end do
          ! End the line after the row has been written
          write(20, *)
      end do
      close (20)
   end subroutine write_relative_drift

end module twopoint_mpi

