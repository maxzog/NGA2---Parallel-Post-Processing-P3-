program testing
   use particle_class
   use twopoint_omp
   use omp_lib
   implicit none

   type(particles) :: parts
   type(omp_stats) :: stats
   integer :: numbins
   real(4) :: length
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
    print *, 'Running with', nthreads, 'threads'

    !$omp parallel private(thread_id)
    thread_id = omp_get_thread_num()
    print *, 'Hello from thread', thread_id
    !$omp end parallel
   call stats%compute_rdf(parts)
   call stats%write_rdf("./outs/rdf.txt")

   call stats%compute_uu(parts)
   call stats%write_uu("./outs/uu.txt")

   call stats%compute_sf(parts)
   call stats%write_sf("./outs/sf.txt")

end program testing
