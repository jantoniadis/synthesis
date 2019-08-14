import os
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, Float, Boolean, BigInteger, ForeignKey
from sqlalchemy.orm import sessionmaker

from sqlalchemy import or_, and_

IN_MEMORY = False
ECHO      = False
Base = declarative_base()

f = 'sqlite:///' + os.getcwd() + '/SimulationResults.db'






class StellarType(Base):
    __tablename__ = 'stellar_type'

    id = Column(Integer, primary_key =True,nullable=False)
    description = Column(String)


    def __repr__(self):
        return "<StellarType(id ='%s',description='%s')>" % (self.id, self.description)


class Simulation(Base):
    __tablename__ = 'simulation_log'

    id = Column(Integer, primary_key = True, nullable=False, autoincrement=True)
    date              = Column(String,   nullable=False)
    code              = Column(String,  nullable=False)
    stars_seed        = Column(Integer, nullable=True)
    positions_seed    = Column(Integer, nullable=True)
    eccentricity_seed = Column(Integer, nullable=True)
    q_seed            = Column(Integer, nullable=True)
    MaxTime           = Column(Float,   nullable=False)
    TimeStep          = Column(Float,   nullable=False)
    GalacticFraction  = Column(Float,   nullable=False)
    CEefficiency      = Column(Float,   nullable=False)
    BinaryFraction    = Column(Float,   nullable=False)
    MinimumMass       = Column(Float,   nullable=False)
    MaximumMass       = Column(Float,   nullable=False)
    



class Star(Base):
    __tablename__ = 'stars'
    id             = Column(
        Integer, 
        primary_key = True, 
        nullable=False, 
        autoincrement=True)
    key            = Column(
        Float,
        nullable=False,
        autoincrement=False)

    is_binary      = Column(Boolean, nullable=True)
    is_interacting = Column(Boolean, nullable=True)
    simulation_id  = Column(Integer, ForeignKey('simulation_log.id'))
    stellar_type_1 = Column(Integer, ForeignKey('stellar_type.id'))
    stellar_type_2 = Column(Integer, ForeignKey('stellar_type.id'))
    eccentricity   = Column(Float)
    semimajoraxis  = Column(Float)
    orbitalperiod  = Column(Float)
    m1             = Column(Float)
    m2             = Column(Float)
    temperature_1  = Column(Float)
    temperature_2  = Column(Float) 
    luminosity_1   = Column(Float)  
    luminosity_2   = Column(Float)  
    R_1            = Column(Float)
    R_2            = Column(Float)
    core_mass_1    = Column(Float)
    core_mass_2    = Column(Float)
    age            = Column(Float)






engine = create_engine(f, echo=ECHO)
Base.metadata.create_all(engine)
Session = sessionmaker(bind=engine)
session = Session()

try:
    add_stellar_type()
except:
    pass



def add_stellar_type(db = StellarType,session=session):

    session.add_all([
        db(id = 0,  description=' deeply convective low mass MS star'),
        db(id = 1,  description= 'main sequence star'),
        db(id = 2,  description= 'hertzsprung gap'),
        db(id = 3,  description= 'core helium burning'),
        db(id = 4,  description= 'core helium burning'),
        db(id = 5,  description= 'first asymptotic giant branch'),
        db(id = 6,  description= 'second asymptotic giant branch'),
        db(id = 7,  description= 'main sequence naked helium star'),
        db(id = 8,  description= 'hertzsprung gap naked helium star'),
        db(id = 9,  description= 'giant branch naked helium star'),
        db(id = 10, description= 'helium white dwarf'),
        db(id = 11, description= 'carbon/oxygen white dwarf'),
        db(id = 12, description= 'oxygen/neon white dwarf'),
        db(id = 13, description= 'neutron star'),
        db(id = 14, description= 'black hole'),
        db(id = 15, description= 'massless supernova'),
        db(id = 16, description= 'unknown stellar type'),
        db(id = 17, description= 'pre-main-sequence star'),
        db(id = -1, description= 'undefined/no star')])

    session.commit()


