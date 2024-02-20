program testing
   use particle_class
   use twopoint_mpi
   implicit none
   include "mpif.h"

   integer :: nproc, rank, ierr

   type(particles) :: parts
   type(mpi_stats) :: stats
   real(4) :: length
   integer :: numbins, start, finish, count_rate
   real(8) :: elapsed_time

   !> Initialize MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,stats%nproc,ierr)

   !> These are the default values for the serial stats class
   numbins = 512
   length = 6.2832

   parts = particles(directory="./data/small_data/",suffix="000148",name="TEST")
   !parts = particles(directory="/home/maxzog/NGA2/examples/SDE_tester/ensight/SDE/particles/",suffix="000151",name="TEST")
   !parts = particles(directory="./data/ners_data/",suffix="000049",name="TEST")

   parts%stat_to_compute=stochastic_velocity
   call parts%set_vec()
   call parts%sort_particles(1,parts%npart)

   stats = mpi_stats(numbins, length, nstep, dt)
   call MPI_Cart_create(MPI_COMM_WORLD,1,stats%nproc,.false.,.true.,stats%comm,ierr)
   call stats%decomp(parts%npart, stats%nproc, rank, stats%imin, stats%imax)
   call stats%decomp(parts%npart, stats%nproc, rank, stats%jmin, stats%jmax)
   
   CALL SYSTEM_CLOCK(start,count_rate) !get start time
 
   call stats%compute_rdf(parts)
   call stats%compute_uu(parts)
   call stats%compute_sf(parts)
 
   CALL SYSTEM_CLOCK(finish) !get finish time

   !Convert time to seconds and print
   elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   if (rank.eq.0) then
      WRITE(*,'(a,f9.3,a)') "    compute took", elapsed_time, " seconds"
      call stats%write_rdf("./outs/rdf_inertial_drift.txt")
      call stats%write_sf("./outs/sf_inertial_drift.txt")
      call stats%write_uu("./outs/uu_inertial_drift.txt")
   end if

   call parts%deallocate_particles()

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   print *, "Compute done on rank ", rank

   call MPI_Finalize(ierr)
end program testing
