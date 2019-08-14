"""
Generates a realistic population of single and binary stars with 
different, primary mass, mass ratio and separation and evolves these over time.
"""

#import all amuse-related modules
from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import quantities
from amuse import datamodel
from amuse.community.seba.interface import SeBa
from amuse.community.bse.interface import BSE
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.lab import *

from matplotlib import pyplot

import numpy as np
import time
import datetime

from alchemy_funcs import *
from plotting import *


def stellar_population(number_of_stars = 1000,min_mass=0.1, max_mass=125, alpha=-2.35):

    """
    Generates a stellar population that follows the Salpeter mass function
    """
    masses = new_salpeter_mass_distribution(number_of_stars, min_mass | units.MSun, max_mass | units.MSun, alpha = alpha)
    return masses


def salpeter(mass,xi=1,alpha=-2.35):
    return xi*mass**alpha 


def stars_per_unit_mass2(m1,m2,ma=0.125,mb=125,a=2.35):
    frac = (2-a)/(1-a)
    return frac*(m1**(1-a) - m2**(1-a))/(ma**(2-a) - mb**(2-a)) 


def chabrier(mass,xi=1,alpha=-2.35):
    result = np.zeros_like(mass)
    filt = np.where(mass < 1.)
    result[filt] = 0.086*(1/(np.log(10)*mass[filt]))*np.exp(-(np.log10(mass[filt]) -np.log10(0.22))**2. / (2*0.57**2))
    filt2 = np.where(mass >= 1.)
    result[filt2] = 0.019199*mass[filt2]**alpha 
    return result


def chabrier2(mass,xi=1,alpha=-2.35):
    result = np.zeros_like(mass)
    filt = np.where(mass < 1.)
    result[filt] = 0.086*(1/(np.log(10)*mass[filt]))*np.exp(-(np.log10(mass[filt]) -np.log10(0.22))**2. / (2*0.57**2))
    filt2 = np.where(mass >= 1.)
    result[filt2] = 0.019199*mass[filt2]**alpha 
    return result


def stars_per_solar_mass(min_mass,max_mass,ma=0.1,mb=125,xi=1.,alpha=-2.35,imf=salpeter):
    mass = np.arange(min_mass,max_mass+1e-9,0.001)
    mass_tot = np.arange(ma,mb+1e-9,0.001)
    N = np.trapz(imf(mass,xi=xi,alpha=alpha),mass)
    Mtot = np.trapz(imf(mass_tot,xi=xi,alpha=alpha)*mass_tot,mass_tot)
    return N/Mtot
    

def random_power_law(a, b, p, size=1):
    g = p+1 
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)




def stars_for_given_Mtot(Mtot,min_mass=0.1,max_mass=125.,k=2.35,_return_int=True):

    frac = Mtot*(k-2)/(k-1)
    a1 = min_mass**(-(k-1)) - max_mass**(-(k-1))
    a2 = min_mass**(-(k-2)) - max_mass**(-(k-2))

    N = frac*a1/a2
    if _return_int:
        N = np.int(N)
    return N


def roche_lobe(q):
    """
    Calculates the Roche lobe radius in units of the orbital separation
    using Eggleton's formula

    input
    q: mass ratio
   
    output
    (r/a): Roche lobe radius as a fraction or the orbital separation
    """

    ra = 0.49*q**(2./3.) / (0.6*q**(2./3.) + np.log(1 + q**(1./3.)))
    return ra



def adjust_separation(separations,rl1):
    """
    checks whether the orbital separation is at least twice as 
    large as the roche lobe radius of the primary and adjust accordingly
    """
    filt                   = np.where(separations <= (2*rl1*separations))
    separations[filt]      = (2*rl1[filt]*separations[filt])

    return separations 
    

def generate_initial_population_grid(
    min_mass, max_mass, number_of_stars, 
    min_ratio, max_ratio,
    min_separation, max_separation,
    min_eccentricity, max_eccentricity,binary_fraction=0.7
):
    """
    creates a set of binaries and a set of stars, the grid
    will be divided in equal parts amongst all dimensions (primary mass,
    ratio and separation), the total number of binaries will
    be the product of all bins.
    """
    

    number_of_single_stars = np.int((1-binary_fraction)*number_of_stars)

    number_of_binaries     = np.int(binary_fraction*number_of_stars/2)

    number_of_systems      = number_of_binaries + number_of_single_stars

    primary_masses         = stellar_population(number_of_systems,min_mass,max_mass) 

    mass_ratios            = np.random.uniform(0.0,1+1e-5,number_of_systems)

    roche_lobe_1           = roche_lobe(mass_ratios)

    separations            = np.random.uniform(np.log10(min_separation),np.log10(max_separation), number_of_systems)

    separations            = (10**separations) 
    separations            = adjust_separation(separations,roche_lobe_1) | units.RSun

    eccentricities         = eccentricity_distribution(number_of_systems)

    #We have twice as many stars as needed and all in binaries, so we need to take care of that
    mass_ratios[number_of_binaries:] = 0
    

        
    primary_stars   = datamodel.Particles(mass=primary_masses)
    secondary_stars = datamodel.Particles(mass=primary_masses * mass_ratios)
        
    stars = datamodel.Particles()
    primary_stars   = stars.add_particles(primary_stars)
    secondary_stars = stars.add_particles(secondary_stars)
        
    binaries = datamodel.Particles(
            semi_major_axis = separations,
            eccentricity    = eccentricities)
    binaries.child1 = list(primary_stars)
    binaries.child2 = list(secondary_stars)

    return binaries, stars


def star_formation(a,tau,time):
    """
    returns the integral of the star formation rate
    """
    return -a*tau*np.exp(-time/tau)


def eccentricity_distribution(N=1):
    return np.random.random(N)**2.0


def separation_distribution(N,a_min,a_max):
    s = np.random.random(N)


def evolve_galaxy(
        end_time, 
        time_step,
        gal_fraction=1e-6,
        binfraction=1.0, 
        Mmin=0.125,
        Mmax=125.,
        a_min=0.01,
        a_max=1e6,
        e_min=0,
        e_max=1.,
        ce=0.5,
        Mmax_NS=2.1,
        plot = True
):



    
    #initialize some variables for monitoring
    date = str(datetime.datetime.now())
    stellar_galaxy_mass = 0.0
    number_of_binaries = 0
    number_of_systems = 0
    number_of_stars = 0
    time = 0.0 * end_time



    try:
        session_id = session.query(Simulation.id).all()[-1][0] + 1
    except:
        session_id = 0


    
    code = SeBa() #link code
    #code = BSE()
    #configure parameters

    code.parameters.common_envelope_efficiency = ce
    code.parameters.maximum_neutron_star_mass = Mmax_NS | units.MSun

    pyplot.ion()
    pyplot.show()

    while time < end_time:
        time += time_step
        mass_in_timestep =  star_formation(15,7e9,time.value_in(units.yr)) - star_formation(15,7e9,(time-time_step).value_in(units.yr))

        mass_in_timestep = mass_in_timestep*gal_fraction
        
        N = stars_for_given_Mtot(mass_in_timestep,Mmin,Mmax)
        #N = 4 #for de-bugging 

        number_of_stars += N
        number_of_binaries = np.int(number_of_stars*binfraction/2)
        number_of_systems = np.int(number_of_stars*(1-binfraction) + number_of_binaries)

        stellar_galaxy_mass += mass_in_timestep

        print
        print
        print
        print "Total stellar mass of the Galaxy [in 1e10 Msun]:", stellar_galaxy_mass/gal_fraction/1e10
        print "Age of the disk:",time.value_in(units.Gyr), "Gyr"
        print "Number of stellar systems:", number_of_systems, "of which", number_of_binaries ,"are binaries"
        print "Injected", N, "stars in the last timestep. Total number of stars:", number_of_stars
        print
        print
        print

        if (time == time_step):
            
            binaries, stars = generate_initial_population_grid(Mmin, Mmax, N, gal_fraction, binfraction, a_min, a_max, e_min , e_max,binfraction)

            code.particles.add_particles(stars)
            code.binaries.add_particles(binaries)

        else:
             binaries2, stars2 = generate_initial_population_grid(0.125, 125,N, 1e-5, 1, 0.1, 1e6, 0.0 , 1.0, binfraction)
             stars.add_particles(stars2)
             binaries.add_particles(binaries2)
             code.particles.add_particles(stars2)
             code.binaries.add_particles(binaries2)



   

        channel_from_code_to_model_for_binaries = code.binaries.new_channel_to(binaries)
        channel_from_code_to_model_for_stars = code.particles.new_channel_to(stars)

        code.evolve_model(time)
        channel_from_code_to_model_for_stars.copy()
        channel_from_code_to_model_for_binaries.copy()
       
        if plot:
            if (time > time_step):
                make_hr_diagram3(stars)
                pyplot.clf()

    #add_sim_info_to_database(end_time=end_time,time_step=time_step,
    #                         gal_fraction=gal_fraction,binfraction=binfraction,Mmin=Mmin,
    #                         Mmax=Mmax,a_min=a_min,a_max=a_max,date=date,code='BSE',ce=ce,MNSmax = Mmax_NS)

    #add_stars_to_database(binaries,sim_id = session_id)
    write_set_to_file(code.particles,'stellar_properties.bse.hdf5',format='hdf5')
    write_set_to_file(binaries,'binary_properties.bse.hdf5',format='hdf5')

   





if __name__ == '__main__':
    evolve_galaxy(100 | units.Myr, 0.05 | units.Myr,a_min=0.1, ce=1, gal_fraction=1e-5,plot=False, Mmin=5, Mmax=30)

