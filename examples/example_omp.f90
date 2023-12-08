program testing
   use particle_class
   use twopoint_omp
   use omp_lib
   implicit none

   type(particles) :: parts
   type(omp_stats) :: stats
   integer :: numbins, start, finish, count_rate
   real(4) :: length
   real(8) :: elapsed_time
   integer :: i, nthreads, thread_id

   !> These are the default values for the serial stats class
   numbins = 32
   length = 6.2832

   parts = particles(directory="./data/medium_data/",suffix="000051",name="TEST")
   call omp_set_num_threads(4)
   call parts%write_particle_data("./outs/test_particle_data.txt")

   parts%stat_to_compute=stochastic_velocity
   call parts%set_vec()

   stats = omp_stats(numbins, length)

   ! Test parallel execution and print thread IDs
   nthreads = omp_get_max_threads()
   !$omp parallel private(thread_id)
   thread_id = omp_get_thread_num()
   !$omp end parallel
   
   CALL SYSTEM_CLOCK(start,count_rate) !get start time
 
   call stats%compute_rdf(parts)
   call stats%compute_uu(parts)
   call stats%compute_sf(parts)
 
   CALL SYSTEM_CLOCK(finish) !get finish time

   !Convert time to seconds and print
   elapsed_time=REAL(finish-start,8)/REAL(count_rate,8)
   
   WRITE(*,'(a,f9.3,a)') "    compute took", elapsed_time, " seconds"
   
   call stats%write_uu("./outs/uu.txt")
   call stats%write_rdf("./outs/rdf.txt")
   call stats%write_sf("./outs/sf.txt")

end program testing
