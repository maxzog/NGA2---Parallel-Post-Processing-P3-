program testing
   use particle_class
   use twopoint_serial_grid
   implicit none

   type(particles) :: parts
   type(grid_stats) :: stats
   integer :: numbins
   real(4) :: length

   !> These are the default values for the serial stats class
   numbins = 512
   length = 6.2832

   parts = particles(directory="./data/ners_data/",suffix="000049",name="TEST",Lx=length,nx=32,nover=1)

   ! call parts%write_particle_data("./outs/test_particle_data.txt")
   parts%stat_to_compute=stochastic_velocity
   call parts%set_vec()

   stats = grid_stats(numbins, length)

   call stats%compute_rdf(parts)
   call stats%compute_uu(parts)
   call stats%compute_sf(parts)
   
   call stats%write_rdf("./outs/rdf.txt")
   call stats%write_sf("./outs/sf.txt")
   call stats%write_uu("./outs/uu.txt")

end program testing
