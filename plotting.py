from matplotlib import pyplot
from amuse.units import units
import numpy as np

def make_hr_diagram(binaries):
    separation = binaries.semi_major_axis.value_in(units.RSun)
    eccentricity = binaries.eccentricity
    pyplot.hexbin(
        separation,
        eccentricity,
        gridsize = 40,
        bins = 'log',
        extent = (0, 1e6, 0.0, 1.0)
    )
    pyplot.xlabel('semi major axis (RSun)')
    pyplot.ylabel('eccentricity')
    pyplot.pause(2)
    pyplot.draw()


def make_hr_diagram3(stars):
    mass = stars.mass.value_in(units.MSun)
    luminosity = stars.luminosity.value_in(units.LSun)
    pyplot.hexbin(
        mass,
        luminosity,
        gridsize = 200,
        bins = 'log',
        mincnt=1,
        extent = (2.1, 0.1, 0.1, 10)
    )
    pyplot.xlabel('Mass')
    pyplot.ylabel('luminosity')
    pyplot.pause(2)
    pyplot.draw()



def make_hr_diagram2(stars,color='blue'):
    #pyplot.figure(figsize = (8,8))
    #pyplot.title('Binary population', fontsize=12)
    mass = stars.mass.value_in(units.MSun)
    age = stars.age.value_in(units.Gyr)
    luminosity = stars.luminosity.value_in(units.LSun)

    filt=np.random.randint(0,len(stars),size=3000)


    pyplot.scatter(
        np.log10(mass[filt]),
        np.log10(luminosity)[filt],s=0.3,color=np.random.rand(3,1)
    )
    pyplot.xlabel('mass')
    pyplot.xlim([-4.0,2.05])
    pyplot.ylabel('luminosity')
    pyplot.ylim([-13,6])
    print "mean age", age.mean()
    pyplot.pause(3)
    pyplot.draw()
