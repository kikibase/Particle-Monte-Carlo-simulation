
"""
If not specified all units used are S.I units
all lenght based units are in metres(m)
all time based unnits are in seconds(s)
energy is in electron-volts(eV)
mass is in kilograms(kg)
"""

import matplotlib.pyplot as plt
import numpy as np


#constants

c = 299792458 #speed of light in m/s
h = 6.62607015 * pow(10, -34) # planck's constant in j/hz or js(joules second)
e = 2.71828 # Euler's constant
π = 3.141 #PI
ε = 8.854*pow(10,-12)# permitivity of free space in C/(Vm) coulombs per joules-metres
#MASS_OF_ELECTRON = 9.109*pow(10,-31) #mass of electron in kg
#CHARGE_OF_ELECTRON = 1.602 * pow(10,19)
Na = 6.022 * pow(10,23)# Avogdros number in /mol
"""Avogadros number """

class particle:

    def __init__(self,
        particle_name:str,
        charged:bool,
        mass:float,
        radius_of_particle,
        charge_of_particle:float=0,
    ):
        """
        parent class for Particle type

        :param float mass: mass in kilograms
        :type priority: integer or None
        :return: the message id
        :rtype: int
        :raises ValueError: if the message_body exceeds 160 characters
        :raises TypeError: if the message_body is not a basestring
        """
        self.particle_name:str = particle_name
        self.charged:bool = charged #if particle is charged
        self.mass:float = mass
        self.charge_value:float = charge_of_particle
        self.radius_of_particle:float = radius_of_particle


    def get_rest_mass_energy(self) -> float:
        e= self.mass*(c**2)
        return e

    def is_charged(self) -> bool:
        return self.charged

    def __str__(self) -> str:
        return f"particle name: {self.particle_name}, mass: {self.mass}, current velocity: {self.get_velocity()}"



class electron(particle):

    def __init__(self) -> None:
        self.particle_name = "electron"
        self.charged = True
        self.mass = 9.1093837 * pow(10, -31) #mass in kilograms
        self.charge_value = 1.602 * pow(10,19)
        self.radius_of_particle = 2.818 * pow(10,-15) # in metres


class material:


    def __init__(#default material is water
        self,
        name_of_material:str = "water",
        size_of_material:list = [1,1,0], #lenght height width of material
        density_of_material:float = 1,
        mass_of_material:float = 100,
        atomic_mass:float = 10, #this is an approximate value(Usually the atomic number exccept its a moelcue)
        molar_mass:float = 0.01801528, #in kg/mol
        mean_excitation_energy = 75 # of water in eV gotten from NIST database
    ):
        """
        Default material is water with mass of 1kg and
        """
        self.name:str = name_of_material
        #properties
        self.density:float = density_of_material
        self.size:list = size_of_material #the size in metres lenght, height, width
        self.mass:float = mass_of_material #the mass of the material in kg
        self.atomic_mass:float = atomic_mass
        self.molar_mass:float = molar_mass


class effects:
    """
    Effects from incident photons:
        Nuclear pair production
        Electronic pair production (Triplet production)
        Photo-electric effect
        Compton effect
        Rayleigh Scattering
        Photonuclear effect
    Effects from charged particles
        coulomb interaction
            elastic collsions
            Inelastic collsions
        Bremmstraung Radiation
    """

    def priority_table(self, particle:particle, material:material) -> dict:
        """
        Example priority table:
            {
            always: [bremstraung],
            random: [
            {hard-collision: 0.9},
            {soft_collision: 0},
            {process: priority:int}
            ]
            }
        """

        if particle.is_charged():
            datum:dict = {
            "always": [self.collision],
            "random": []
            }
        else:
            datum:dict = {
                "always": [],
                "random": []
            }
        return datum
    def collision(self,
                   # z:float #charge of incoming particle
                   ) -> tuple[float, float]:
        # CONSTANT = pow(e,2)/(4*e*π*ε)
        # CONST2 = pow(z,2)/(MASS_OF_ELECTRON*)# mulitply velocity and impact parameter
        # energy_loss:float = 2*pow((CONSTANT),2)* CONST2
        return 0.0, 0 # returns energyloss, angle deviated



    def bremsstrahlung(self, rest_energy:float, kinetic_energy:float, particle_radius:float, atomic_number:float, molar_mass:float) ->float:
        # Constants
        import math

        # Physical constants
        alpha = 1 / 137  # Fine-structure constant
        re = 2.818e-13  # Classical electron radius in cm
        electron_rest_mass_MeV = 0.511  # Electron rest mass energy in MeV
        water_molar_mass = 18.0153  # Molar mass of water in g/mol

        # Inputs for the problem
        electron_energy_MeV = 10  # Incident electron energy in MeV
        Z_eff = (2 * 1 + 8) / 3  # Effective Z for water (H2O)
        Na = 6.022e23  # Avogadro's number

        # Logarithmic factors
        log_factor1 = math.log(183 / (Z_eff ** (1 / 3)))  # ln(183 / Z^(1/3))
        log_factor2 = math.log(2 * electron_energy_MeV / electron_rest_mass_MeV)  # ln(2E / mc^2)

        # Calculate the total cross section formula
        pre_factor = 4 * alpha * re**2 * Z_eff**2 / 137  # Pre-factor (cm^2)
        total_cross_section_per_molecule = pre_factor * log_factor1 * log_factor2  # Cross-section in cm^2 per molecule

        # Convert to per atom by dividing by the number of atoms in a water molecule (3: 2 H + 1 O)
        total_cross_section_per_atom = total_cross_section_per_molecule / 3

        # Output results


        radiation_crossection:float = total_cross_section_per_molecule / 3
        total_energy:float = rest_energy+kinetic_energy
        stopping_power:float = (Na/molar_mass)*radiation_crossection*total_energy

        return stopping_power


    def pair_production(self):
        return None
    def triplet_production(self):
        return None
    def Photo_electric_effect(self):
        return None
    def Compton_effect(self):
        return None
    def Rayleigh_Scattering(self):
        return None
    def Photonuclear_effect(self):
        return None


class beam:

    def __init__(self,
        particle_type:particle,
        number_of_particles:int,
        is_monoenergetic:bool,
        polyenergetic_array:list = [],
        energy:float = 0):
        """
        This creates instances of the given particle to be used as a beam
        Energy should be in eV
        """
        self.is_beam_set:bool = False

        if is_monoenergetic and (energy == 0):
            print("You should give the energy level of the monoenergetic beam")
        elif not is_monoenergetic and (polyenergetic_array.__len__() == 0):
            print("You should give an array of energies for the polyenergetic beam")

        self.particle = particle_type
        self.number_of_particles:int = number_of_particles
        self.x_coords_of_particles = np.ones(number_of_particles)
        self.y_coords_of_particles =np.ones(number_of_particles)

        self.angle_of_particles = np.ones(number_of_particles)
        self.array_of_energies = np.multiply(np.ones(number_of_particles), energy)

    def beam_setup(self, point_of_entry:list, angle_of_entry:float):
        #sets the default position of all particles by the data given
        #current setup only works for monoenergetic beam
        self.is_beam_set:bool = True
        self.x_coords_of_particles = np.multiply(self.x_coords_of_particles, point_of_entry[0])
        self.y_coords_of_particles = np.multiply(self.y_coords_of_particles, point_of_entry[1])
        self.angle_of_particles = np.multiply(self.angle_of_particles, angle_of_entry)


    def get_velocity(self, index:int) -> float:
        energy_fraction:float = pow(self.particle.get_rest_mass_energy()/self.array_of_energies[index], 2)
        v:float = c * pow(1-energy_fraction , 0.5)
        return v

    def x_velocity(self,index:int) -> float:
        return self.get_velocity(index)* (np.cos(self.angle_of_particles[index]))
    def y_velocity(self,index:int) -> float:
        return self.get_velocity(index)* (np.sin(self.angle_of_particles[index]))


    def get_number_of_particles(self) -> int:
        return self.number_of_particles



class simulate:
    def __init__(self,
                beam_of_particle:beam,
                material_to_interact:material,
                point_of_entry:list,
                angle_of_entry:float,
                dicretization:int,
                energy_threshold:float = 0) -> None:
        if energy_threshold == 0:
            # if there is a zero energy threshold caclulate what it should be with the stopping power
            pass

        self.beam:beam = beam_of_particle
        self.material:material = material_to_interact
        self.energy_threshold:float = energy_threshold

        self.sample_x:float = material_to_interact.size[0]/dicretization
        self.sample_y:float = material_to_interact.size[1]/dicretization


        #sets the position of all particles by the data given
        self.beam.beam_setup(point_of_entry, angle_of_entry)
        self.simulation()


    def simulation(self) -> None:
        is_simulating = True
        rest_mass_energy:float = self.beam.particle.get_rest_mass_energy()
        particle_radius:float = self.beam.particle.radius_of_particle
        atomic_number:float = self.material.atomic_mass
        molar_mass:float = self.material.molar_mass
        self.history_of_x = np.array([[],[]])
        self.history_of_y = np.array([[],[]])
        history_of_stopping_power_rad = np.array([[],[]])
        history_of_stopping_power_col = np.array([[],[]])

        while is_simulating:
            i:int = 0#iterator for all the particles
            untracked:list = []

            while i < self.beam.number_of_particles:
                if (
                    (i not in untracked) and
                    self.beam.x_coords_of_particles[i] <=self.material.size[0] and
                    self.beam.y_coords_of_particles[i] <=self.material.size[1] and
                    self.beam.array_of_energies[i] >= self.energy_threshold
                ):
                    
                    mod_of_x = self.beam.x_coords_of_particles[i]% self.sample_x
                    old_mod_of_x = abs(self.beam.x_coords_of_particles[i]-self.beam.x_velocity(i)% self.sample_x)
                    mod_of_y = self.beam.y_coords_of_particles[i]% self.sample_x
                    old_mod_of_y = abs(self.beam.y_coords_of_particles[i]-self.beam.y_velocity(i)% self.sample_y)
                    if (
                        mod_of_x < old_mod_of_x or
                        mod_of_y < old_mod_of_y
                    ):
                        ke:float = self.beam.array_of_energies[i]
                        energy_loss:float
                        angle_deviated:float
                        energy_loss, angle_deviated = effects().collision()
                        energy_loss += effects().bremsstrahlung(rest_mass_energy, ke,particle_radius,atomic_number,molar_mass)

                        self.beam.array_of_energies[i] -= energy_loss
                        self.beam.angle_of_particles[i] += angle_deviated

                    self.beam.x_coords_of_particles[i] += self.beam.x_velocity(i)
                    self.beam.y_coords_of_particles[i] += self.beam.y_velocity(i)
                    self.history_of_x = np.append(self.history_of_x[i],self.beam.x_coords_of_particles[i])
                    self.history_of_y = np.append(self.history_of_y[i],self.beam.y_coords_of_particles[i])
                elif(i not in untracked):
                    untracked.append(i)
                elif len(untracked) >=self.beam.number_of_particles:
                    is_simulating = False


                i+=1




    def graph(self):#make graphs
        # Make a graph for the position of the electrons y against x

        fig, ax = plt.subplots()
        ax.plot(self.history_of_x, self.history_of_y)
        plt.show()
