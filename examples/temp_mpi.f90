program testing
   use particle_class
   use twopoint_mpi
   implicit none
   include "mpif.h"

   integer :: nproc, rank, ierr

   type(particles) :: parts
   type(mpi_stats) :: stats
   integer :: numbins, start, finish, count_rate
   real(4) :: length
   real(8) :: elapsed_time

   !> These are the default values for the serial stats class
   numbins = 512
   length = 6.2832

   !parts = particles(directory="./data/100000_particles_data/SDE/particles/",suffix="000002",name="TEST")
   parts = particles(directory="./data/medium_data/",suffix="000051",name="TEST")
   !parts = particles(directory="./data/ners_data/",suffix="000049",name="TEST")

   parts%stat_to_compute=stochastic_velocity
   call parts%set_vec()

   stats = mpi_stats(numbins, length)

   !> Initialize MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,stats%nproc,ierr)

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
   
   if (rank.eq.0) then
      WRITE(*,'(a,f9.3,a)') "    compute took", elapsed_time, " seconds"
      call stats%write_rdf("./outs/rdf.txt")
      call stats%write_sf("./outs/sf.txt")
      call stats%write_uu("./outs/uu.txt")
   end if

   print *, "YOU MADE IT OUT"
end program testing
