from modules import material, particle, electron, beam, simulate

#simulation setup

# create a 10MeV of monoenergetic electrons
#

particulate:particle = electron([0,0])# creates an electron
electron_beam:beam = beam(particulate, 500, True, energy=(10*pow(10,6)))# create an electron beam of 10MeV
water:material = material()

#simulation
simulate(electron_beam, water, [0.5,0.5,0], 90)

print(particulate)
