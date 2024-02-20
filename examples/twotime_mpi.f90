program testing
   use particle_class
   use twopoint_mpi
   implicit none
   include "mpif.h"

   integer :: nproc, rank, ierr, i

   type(particles) :: partsn
   type(mpi_stats) :: stats
   real(4) :: length, dt
   integer :: numbins, nstep, start, finish, count_rate
   integer :: stepi, stepf
   real(8) :: elapsed_time
   character(len=6) :: strn

   !> Initialize MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,stats%nproc,ierr)

   !> These are the default values for the serial stats class
   numbins = 512
   length = 6.2832
   nstep = 46
   dt = 0.01
   stepi = 47
   stepf = stepi - nstep

   call intToPaddedString(stepi, strn)
   partsn = particles(directory="/Users/maxzog/DATA/stk1dns/ensight/HIT/particles/",suffix=strn,name="TEST")

   print *, "NPART :: ", partsn%npart

   partsn%stat_to_compute=stochastic_velocity
   call partsn%set_vec()
   call partsn%sort_particles(1,partsn%npart)

   stats = mpi_stats(numbins, length, nstep, dt)
   call MPI_Cart_create(MPI_COMM_WORLD,1,stats%nproc,.false.,.true.,stats%comm,ierr)
   call stats%decomp(partsn%npart, stats%nproc, rank, stats%imin, stats%imax)
   call stats%decomp(partsn%npart, stats%nproc, rank, stats%jmin, stats%jmax)

   call stats%compute_var(partsn)
   
   CALL SYSTEM_CLOCK(start,count_rate) !get start time
   
   !> Loop through time
   do i = 1, stats%nstep
      ! Increment step counter
      stats%step = i

      ! Load data and compute two-time correlation
      load_and_compute: block
         type(particles) :: partsm
         character(len=6) :: strm

         call intToPaddedString(stepi-stats%step, strm)
         partsm = particles(directory="/Users/maxzog/DATA/stk1dns/ensight/HIT/particles/",suffix=strm,name="TEST")
         partsm%stat_to_compute=stochastic_velocity
         call partsm%set_vec()
         call partsm%sort_particles(1,partsm%npart)

         call stats%compute_ac(partsn, partsm)

         call partsm%deallocate_particles()
      end block load_and_compute

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
   end do
 
   CALL SYSTEM_CLOCK(finish) !get finish time

   !Convert time to seconds and print
   elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   if (rank.eq.0) then
      WRITE(*,'(a,f9.3,a)') "    compute took", elapsed_time, " seconds"
      call stats%write_ac("./outs/ac.txt")
   end if

   call partsn%deallocate_particles()

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   print *, "Compute done on rank ", rank

   call MPI_Finalize(ierr)
   contains

   subroutine intToPaddedString(num, paddedString)
       implicit none
       integer, intent(in) :: num
       character(len=6), intent(out) :: paddedString
   
       write(paddedString, '(I6.6)') num
   end subroutine intToPaddedString

end program testing
