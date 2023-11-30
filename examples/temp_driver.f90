program testing
   use particle_class
   use twopoint_serial
   implicit none

   type(particles) :: parts
   type(serial_stats) :: stats
   integer :: numbins
   real(4) :: length

   !> These are the default values for the serial stats class
   numbins = 32
   length = 6.2832

   parts = particles(directory="./data/small_data/",suffix="000148",name="TEST")

   call parts%write_particle_data("./outs/test_particle_data.txt")

   parts%stat_to_compute=stochastic_velocity
   call parts%set_vec()

   stats = serial_stats(numbins, length)

   call stats%compute_rdf(parts)
   call stats%write_rdf("./outs/rdf.txt")

   call stats%compute_uu(parts)
   call stats%write_uu("./outs/uu.txt")

   call stats%compute_sf(parts)
   call stats%write_sf("./outs/sf.txt")

end program testing