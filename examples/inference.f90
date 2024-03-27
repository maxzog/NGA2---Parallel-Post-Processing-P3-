program testing
   use particle_class
   use twopoint_mpi
   implicit none
   include "mpif.h"

   integer :: nproc, rank, ierr, i

   type(particles) :: partsn
   type(mpi_stats) :: stats
   real(4) :: length, dt
   integer :: numbins, N, nstep, start, finish, count_rate
   integer :: stepi, stepf, ubins
   real(8) :: elapsed_time
   character(len=6) :: strn


   !> These are the default values for the serial stats class
   numbins = 512
   N       = 64
   length  = 6.2832
   dt      = 0.1254
   stepi   = 11
   stepf   = 10
   ubins   = 32

   call intToPaddedString(stepi, strn)
   partsn = particles(directory="/home/maxzog/NGA2/examples/hit/ensight/HIT/particles/"&
                     &,suffix=strn,name="TEST",Lx=length,nx=N,has_uf=.false.)

   partsn%stat_to_compute=fluid_velocity
   call partsn%set_vec()
   call partsn%sort_particles(1,partsn%npart)

   stats = mpi_stats(numbins, length, nstep, dt, ubins)
   !> Initialize MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,stats%rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,stats%nproc,ierr)
   call MPI_Cart_create(MPI_COMM_WORLD,1,stats%nproc,.false.,.true.,stats%comm,ierr)
   call stats%decomp(partsn%npart, stats%nproc, stats%rank, stats%imin, stats%imax)
   call stats%decomp(partsn%npart, stats%nproc, stats%rank, stats%jmin, stats%jmax)

   call stats%compute_var(partsn)
   
   CALL SYSTEM_CLOCK(start,count_rate) !get start time
   
   ! Load data and compute two-time correlation
   load_and_compute: block
      type(particles) :: partsm
      character(len=6) :: strm

      call intToPaddedString(stepf, strm)
      partsm = particles(directory="/home/maxzog/NGA2/examples/hit/ensight/HIT/particles/"&
                        &,suffix=strm,name="TEST",Lx=length,nx=N,has_uf=.false.)
      partsm%stat_to_compute=fluid_velocity
      call partsm%set_vec()
      call partsm%sort_particles(1,partsm%npart)

      call stats%infer_drift(partsn, partsm)
      call partsm%deallocate_particles()
   end block load_and_compute

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
 
   CALL SYSTEM_CLOCK(finish) !get finish time

   !Convert time to seconds and print
   elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   if (stats%rank.eq.0) then
      WRITE(*,'(a,f9.3,a)') "    compute took", elapsed_time, " seconds"
      call stats%write_drift("./outs/drift_128.txt")
   end if

   call partsn%deallocate_particles()

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   print *, "Compute done on rank ", stats%rank

   call MPI_Finalize(ierr)
   contains

   subroutine intToPaddedString(num, paddedString)
       implicit none
       integer, intent(in) :: num
       character(len=6), intent(out) :: paddedString
   
       write(paddedString, '(I6.6)') num
   end subroutine intToPaddedString

end program testing
