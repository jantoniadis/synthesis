def evolve_galaxy(end_time, time_step,gal_fraction=1e-4,binfraction=0.5):
    
    stellar_galaxy_mass = 0.0
    time = 0.0 * end_time
    code = SeBa()

    while time < end_time:
        time += time_step

        mass_in_timestep =  star_formation(12,7e9,time.value_in(units.yr)) - star_formation(12,7e9,(time-time_step).value_in(units.yr))
        mass_in_timestep = mass_in_timestep*gal_fraction
        
        N = stars_for_given_Mtot(mass_in_timestep,0.15,125.)
        

        if (time == time_step):
            
            binaries, stars = generate_initial_population_grid(0.15, 125.,N, 1e-5, 1, 1, 1e6, 0.0 , 1.0)
 

        else:
             binaries2, stars2 = generate_initial_population_grid(0.15, 125.,N, 1e-5, 1, 1, 1e6, 0.0 , 1.0)
             binaries.add_particles(stars2)
             binaries.add_particles(binaries2)
             #code.particles.add_particles(stars)
             #code.binaries.add_particles(binaries)
            
        print "Stellar population:", len(binaries), "binaries"
    

            
        print "Stellar population:", len(binaries), "binaries"
    

    
        channel_from_code_to_model_for_binaries = code.binaries.new_channel_to(binaries)
        channel_from_code_to_model_for_stars = code.particles.new_channel_to(stars)
    

        stellar_galaxy_mass += mass_in_timestep
        print "Total stellar mass of the Galaxy [in 1e10 Msun]:", stellar_galaxy_mass/1e10
        code.evolve_model(time)
        print "evolved to time: ", time.as_quantity_in(units.Myr)
        
    channel_from_code_to_model_for_stars.copy()
    channel_from_code_to_model_for_binaries.copy()
    make_hr_diagram(binaries)
