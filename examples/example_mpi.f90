program testing
   use particle_class
   use twopoint_mpi
   implicit none
   include "mpif.h"

   integer :: nproc, rank, ierr

   type(particles) :: parts
   type(mpi_stats) :: stats
   real(4) :: length, D
   integer :: numbins, start, finish, count_rate
   real(8) :: elapsed_time


   !> These are the default values for the serial stats class
   numbins = 512
   length = 6.2832

   parts = particles(directory="/home/maxzog/NGA2/examples/SDE_tester/ensight/SDE/particles/"&
                              &,suffix="000023",name="TEST",Lx=length, nx=64, has_uf=.true.)
   ! parts = particles(directory="/home/maxzog/iso1024data/jhtdb_particles/dns/"&
   !                            &,suffix="000001",name="TEST",Lx=length, nx=64, has_uf=.true.)
   parts%stat_to_compute=stochastic_velocity
   call parts%set_vec()
   call parts%sort_particles(1,parts%npart)

   stats = mpi_stats(numbins, length)

   !> Initialize MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
   call MPI_Comm_size(MPI_COMM_WORLD,stats%nproc,ierr)
   call MPI_Cart_create(MPI_COMM_WORLD,1,stats%nproc,.false.,.true.,stats%comm,ierr)
   call stats%decomp(parts%npart, stats%nproc, rank, stats%imin, stats%imax)
   call stats%decomp(parts%npart, stats%nproc, rank, stats%jmin, stats%jmax)
 
   CALL SYSTEM_CLOCK(start,count_rate) !get start time
 
   call stats%compute_uu(parts)
   call stats%compute_sf(parts)
 
   CALL SYSTEM_CLOCK(finish) !get finish time

   !Convert time to seconds and print
   elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   
!!   D_statistic: block
!!      integer :: i, j, k
!!      real(4) :: lambda, s, sp
!!
!!      lambda = real(parts%npart) / (parts%grid%nx*parts%grid%ny*parts%grid%nz)
!!      s = sqrt(lambda)
!!      sp = 0.0
!!      do k = 1, parts%grid%nz
!!         do j = 1, parts%grid%ny
!!            do i = 1, parts%grid%nx
!!               sp = sp + (parts%grid%npic(i,j,k) - lambda)**2
!!            end do
!!         end do
!!      end do
!!      sp = sp / (parts%grid%nx*parts%grid%ny*parts%grid%nz)  
!!      sp = sqrt(sp)
!!      D = (sp - s) / lambda
!!   end block D_statistic

   if (rank.eq.0) then
      print *, "D = ", D
      WRITE(*,'(a,f9.3,a)') "    compute took", elapsed_time, " seconds"
      call stats%write_sf("./outs/sf_rho.txt")
      call stats%write_uu("./outs/uu_rho.txt")
   end if

   call parts%deallocate_particles()

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   print *, "Compute done on rank ", rank

   call MPI_Finalize(ierr)
end program testing
