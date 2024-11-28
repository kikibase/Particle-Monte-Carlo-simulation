from modules import material, particle, electron, beam, simulate

#simulation setup


#CONSTANTS | Simulation parameters
electron_energy:int = (10*pow(10,6))
is_monoenergetic = True
number_of_electrons = 1 #500
enegy_threshold = pow(10,-4) #energy to stop tracking particles
steps = 1000 #discretization of values




# create a 10MeV of monoenergetic electrons
particulate:particle = electron([0,0])# creates an electron
electron_beam:beam = beam( # create a monoenergetic electron beam of 10MeV with 500 particles
    particulate, 
    number_of_electrons, 
    is_monoenergetic, 
    energy= electron_energy
    )

#Creates water as the material to 
water:material = material()


#simulation
simulation:simulate = simulate(
    electron_beam,
    water,
    [0.5,0.5,0],
    0,
    steps,
    enegy_threshold)




simulation.graph() #returns all graphs of the simulation
