import sparpy
import random


def test_simulation():
    N = 100
    D = 1.0
    lower_bound = [0,0]
    upper_bound = [1,1]
    periodic = [True,True]
    sim_time = 1.0
    number_of_observations = 100
    integrate_time = sim_time/number_of_observations
    dt = 0.001
    epsilon = 0.01
    cutoff = 0.1

    particles = sparpy.Particles2(N)
    for p in particles:
        p.position = [random.uniform(lower_bound[0],upper_bound[0]),random.uniform(lower_bound[1],upper_bound[1])]

    simulation = sparpy.Simulation2()
    simulation.set_domain(lower_bound,upper_bound,periodic)
    simulation.add_particles(particles,D)
    simulation.add_force(particles,particles,sparpy.exponential_force2(cutoff,epsilon))

    for i in range(number_of_observations):
        simulation.integrate(integrate_time,dt)
        vtk = particles.get_grid(True)

if __name__ == "__main__":
    test_simulation()
