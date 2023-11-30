module twopoint_serial
    use particle_class
    implicit none
    private
    
    real(4) :: PI=3.14159265

    public :: serial_stats

    type serial_stats
      integer :: nb=32
      real(4) :: L=6.2832
      real(4) :: dr
      real(4) :: rmax
      real(4) :: Ui2
      real(4) :: DeltaU2
      real(4), allocatable, dimension(:) :: uul, uut, sfl, sft, rdf
      integer, allocatable, dimension(:) :: c
    contains
      procedure :: compute_uu
      procedure :: compute_sf
      procedure :: compute_rdf

      procedure :: write_uu
      procedure :: write_sf
      procedure :: write_rdf
    end type serial_stats
    
    interface serial_stats
      procedure :: constructor
    end interface serial_stats

contains

    function constructor(numbins, length) result(self)
      implicit none
      type(serial_stats) :: self
      integer, optional :: numbins
      real(4), optional :: length

      if (present(numbins)) self%nb=numbins
      if (present(length)) self%L=length

      allocate(self%uul(1:self%nb)); self%uul=0.0
      allocate(self%uut(1:self%nb)); self%uut=0.0
      allocate(self%sfl(1:self%nb)); self%sfl=0.0
      allocate(self%sft(1:self%nb)); self%sft=0.0
      allocate(self%rdf(1:self%nb)); self%rdf=0.0

      allocate(self%c(1:self%nb)); self%c=0

      self%rmax=norm2([0.5*self%L, 0.5*self%L, 0.5*self%L])
      self%dr=self%rmax / self%nb
    end function constructor

    subroutine compute_uu(this, parts)
        type(particles), intent(in) :: parts
        class(serial_stats), intent(inout) :: this
        integer :: i, j, ir
        real(4) :: r(3), rll(3), rt1(3), rt2(3)

        !> Reset
        this%c = 0; this%uul=0.0; this%uut=0.0

        do i = 1, parts%npart
            do j = i, parts%npart
                call par_perp_u(this, parts%p(i), parts%p(j), rll, rt2)
                r = parts%p(j)%pos - parts%p(i)%pos
                ir = floor(norm2(r) / this%dr) + 1
                if (ir <= this%nb .and. parts%p(i)%id /= parts%p(j)%id) then
                    this%uul(ir) = this%uul(ir) + dot_product(parts%p(i)%vec, rll) * dot_product(parts%p(j)%vec, rll)
                    this%uut(ir) = this%uut(ir) + dot_product(parts%p(i)%vec, rt2) * dot_product(parts%p(j)%vec, rt2)
                    this%c(ir) = this%c(ir) + 1
                end if
            end do
        end do

        do i = 1, this%nb
            !> Prevent division by zero
            if (this%c(i) == 0) then
                this%c(i) = 1
            end if
            !> Normalize bins
            this%uul(i) = this%uul(i) / this%c(i)
            this%uut(i) = this%uut(i) / this%c(i)
        end do
    end subroutine compute_uu

    subroutine compute_sf(this, parts)
      type(particles), intent(in) :: parts
      class(serial_stats), intent(inout) :: this
      integer :: i, j, ir
      real(4) :: r(3), rll(3), rt1(3), rt2(3)

      !> Reset 
      this%c = 0; this%sfl=0.0; this%sft=0.0

      do i = 1, parts%npart
          do j = i, parts%npart
              call par_perp_u(this, parts%p(i), parts%p(j), rll, rt2)
              r = parts%p(j)%pos - parts%p(i)%pos
              ir = floor(norm2(r) / this%dr) + 1
              if (ir <= this%nb .and. parts%p(i)%id /= parts%p(j)%id) then
                  this%sfl(ir) = this%sfl(ir) + dot_product(parts%p(j)%vec - parts%p(i)%vec, rll) &
                               & * dot_product(parts%p(j)%vec - parts%p(i)%vec, rll)
                  this%sft(ir) = this%sft(ir) + dot_product(parts%p(j)%vec - parts%p(i)%vec, rt2) & 
                               & * dot_product(parts%p(j)%vec - parts%p(i)%vec, rt2)
                  this%c(ir) = this%c(ir) + 1
              end if
          end do
      end do

      do i = 1, this%nb
          !> Prevent division by zero
          if (this%c(i) == 0) then
              this%c(i) = 1
          end if
          !> Normalize bins
          this%sfl(i) = this%sfl(i) / this%c(i)
          this%sft(i) = this%sft(i) / this%c(i)
      end do
    end subroutine compute_sf

    subroutine compute_rdf(this, parts)
      type(particles), intent(in) :: parts
      class(serial_stats), intent(inout) :: this
      integer :: i, j, ir
      real(4) :: V, Vshell, N, rho, r(3)

      !> Reset
      this%rdf = 0.0

      do i = 1, parts%npart
         do j = i, parts%npart
               r = get_minr(this, parts%p(i), parts%p(j))
               ir = floor(norm2(r) / this%dr) + 1
               if (ir <= this%nb .and. parts%p(i)%id /= parts%p(j)%id) then 
                  this%rdf(ir) = this%rdf(ir) + 1.0
               end if
         end do
      end do

      !> Prepare RDF normalization
      V=this%L**3 
      N=parts%npart * (parts%npart - 1) / 2
      rho=N/V

      do i = 1, this%nb
         Vshell=1.33333 * PI * ( (this%dr * i)**3 - (this%dr * (i - 1))**3 )
         this%rdf(i) = this%rdf(i)/(Vshell*rho)
      end do
    end subroutine compute_rdf

   function get_minr(this, p, q) result(r)
      type(serial_stats) :: this
      type(part), intent(in) :: p, q
      real(4) :: r(3)
      r = q%pos - p%pos
      r(1)=r(1)-this%L*nint(r(1)/this%L)
      r(2)=r(2)-this%L*nint(r(2)/this%L)
      r(3)=r(3)-this%L*nint(r(3)/this%L)
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
      type(serial_stats) :: this
      type(part), intent(in) :: p, q
      real(4), intent(out) :: rll(3), rt2(3)
      real(4) :: r(3), rt1(3)
      integer :: i

      r = q%pos - p%pos
      do i = 1, 3
         if (r(i) > 0.5 * this%L) then
            r(i) = r(i) - this%L
         end if
         if (r(i) < -0.5 * this%L) then
            r(i) = r(i) + this%L
         end if
      end do

      rll = r / norm2(r)
      rt1 = cross_product(rll, p%vec) / norm2(cross_product(rll, p%vec))
      rt2 = cross_product(rt1, rll) / norm2(cross_product(rt1, rll))
   end subroutine par_perp_u

   !> I/O routines
   subroutine write_uu(this, outfile)
      implicit none
      class(serial_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open(unit=20, file=outfile, status='replace')
      write(20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
          write(20, '(I5, 2F15.8)') i, this%uul(i), this%uut(i)
      end do
      close(20)
   end subroutine write_uu

   subroutine write_sf(this, outfile)
      implicit none
      class(serial_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open(unit=20, file=outfile, status='replace')
      write(20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
          write(20, '(I5, 2F15.8)') i, this%sfl(i), this%sft(i)
      end do
      close(20)
   end subroutine write_sf

   subroutine write_rdf(this, outfile)
      implicit none
      class(serial_stats), intent(in) :: this
      character(len=*), intent(in) :: outfile
      integer :: i

      open(unit=20, file=outfile, status='replace')
      write(20, '(A, F10.5)') 'dr: ', this%dr
      do i = 1, this%nb
          write(20, '(I5, F15.8)') i, this%rdf(i)
      end do
      close(20)
   end subroutine write_rdf

end module twopoint_serial

