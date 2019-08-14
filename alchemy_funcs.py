from alchemy import *
from amuse.units import units
import numpy as np

def add_sim_info_to_database(
        end_time,
        time_step,
        gal_fraction,
        binfraction,
        Mmin,
        Mmax,
        a_min,
        a_max,
        date,
        code='BSE',
        ce=0.5,
        MNSmax=2.1,
        session = session,
        table=Simulation):
    
    end_time  = end_time.value_in(units.Gyr)
    time_step = time_step.value_in(units.Gyr)

    session.add(table(
        date = date,
        code = code,
        MaxTime = end_time,
        TimeStep = time_step,
        CEefficiency = ce,
        BinaryFraction = binfraction,
        GalacticFraction = gal_fraction,
        MinimumMass = Mmin,
        MaximumMass = Mmax
    ))

    session.commit()




def add_stars_to_database(binaries,sim_id,session=session,table=Star):
    n = len(binaries)
    for i in range(n):
        binary = binaries[i]
        if (binary.child2.mass.value_in(units.MSun) != 0.0):
            is_binary = True
        else:
            is_binary = False


        star = table(
            key                = int(binary.key),
            is_binary          = is_binary,
            eccentricity       = np.float(binary.eccentricity),
            m1                 = np.float(binary.child1.mass.value_in(units.MSun)),
            temperature_1      = np.float(binary.child1.temperature.value_in(units.K)),
            luminosity_1       = np.float(binary.child1.luminosity.value_in(units.LSun)), 
            R_1                = np.float(binary.child1.radius.value_in(units.RSun)),
            age                = np.float(binary.age.value_in(units.Gyr)),
            stellar_type_1     = np.int(binary.child1.stellar_type.value_in(units.stellar_type)),
            simulation_id      = sim_id)


        if is_binary:
            star.semimajoraxis  = np.float(binary.semi_major_axis.value_in(units.RSun))
            star.m2             = np.float(binary.child2.mass.value_in(units.MSun))
            star.temperature_2  = np.float(binary.child2.temperature.value_in(units.K))
            star.luminosity_2   = np.float(binary.child2.luminosity.value_in(units.LSun)) 
            star.R_2            = np.float(binary.child2.radius.value_in(units.RSun))
            try:
                star.stellar_type_2 = np.int(binary.child2.stellar_type.value_in(units.stellar_type))
            except:
                star.stellar_type_2 = -1
        
        session.add(star)
        session.commit()

    
