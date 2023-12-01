program testing
   use particle_class
   use twopoint_mpi
   implicit none
   include "mpif.h"

   integer :: nproc, rank, ierr

   type(particles) :: parts
   type(mpi_stats) :: stats
   integer :: numbins
   real(4) :: length

   !> These are the default values for the serial stats class
   numbins = 128
   length = 6.2832

   parts = particles(directory="./data/medium_data/",suffix="000051",name="TEST")

   call parts%write_particle_data("./outs/test_particle_data.txt")

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
  
   call stats%compute_rdf(parts)
   call stats%compute_uu(parts)
   call stats%compute_sf(parts)
   
   if (rank.eq.0) then
      call stats%write_rdf("./outs/rdf.txt")
      call stats%write_sf("./outs/sf.txt")
      call stats%write_uu("./outs/uu.txt")
   end if

end program testing
