
"""
If not specified all units used are S.I units
all lenght based units are in metres(m)
all time based unnits are in seconds(s)
energy is in electron-volts(eV)
mass is in kilograms(kg)
"""

#constants

c = 299792458 #speed of light in m/s
h = 6.62607015 * pow(10, -34) # planck's constant in j/hz or js(joules second)


class particle:

    def __init__(self, particle_name, charged, mass, coord:list):
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

        #positon
        self.coords:list = [0,0,0] #x,y,x
        self.velocity:list = [0,0,0] #x,y,x

    def setvelocity(self, x:float, y:float, z:float) -> None:
        self.velocity[0]=x
        self.velocity[1]=y
        self.velocity[2]=z

    def get_velocity(self) -> float:
        #calculate the vleocity of the particle
        return 0.0

    def get_rest_mass_energy(self) -> float:
        e= self.mass*(c**2)
        return e

    def is_charged(self) -> bool:
        return self.charged

    def __str__(self) -> str:
        return f"particle name: {self.particle_name}, mass: {self.mass}, current velocity: {self.get_velocity()}"

    def tick(self) -> None: #set of actions todo at each tick/interval
        self.velocity[0]+=self.coords[0]
        self.velocity[1]+=self.coords[1]
        self.velocity[2]+=self.coords[2]

class electron(particle):

    def __init__(self, coords:list) -> None:
        self.particle_name = "electon"
        self.charged = True
        self.mass = 9.1093837 * pow(10, -31) #mass in kilograms
        self.coords = [0,0,0]

        i:int = 0
        for coord in coords:
            self.coords[i] = coord

class material:


    def __init__(#default material is water
        self,
        name_of_material:str = "water",
        size_of_material:list = [1,1,0], #lenght height width of material
        density_of_material:float = 1,
        mass_of_material:float = 100
    ):
        """
        Default material is water with mass of 1kg and
        """
        self.name:str = name_of_material
        #properties
        self.density:float = density_of_material
        self.size:list = size_of_material #the size in metres lenght, height, width
        self.mass:float = mass_of_material #the mass of the material in kg


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
            "always": [self.bremmstraung],
            "random": []
            }
        else:
            datum:dict = {
                "always": [],
                "random": []
            }
        return datum

    def bremmstraung(self):

        return None

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

    def __init__(self, particle_type:particle, number_of_electrons:int, is_monoenergetic:bool, polyenergetic_array:list = [], energy:float = 0):
        """
        This creates instances of the given electron to be used as a beam
        Energy should be in eV
        """
        if not is_monoenergetic and (polyenergetic_array.__len__() == 0):
            print("You should give an array of energies for the polyenergetic beam")
        elif is_monoenergetic and (energy == 0):
            print("You should give the energy level of the monoenergetic beam")


class simulate:

    def __init__(self, beam_of_particle:beam, material_to_interact:material, point_of_entry:list, angle_of_entry:float) -> None:
        pass
