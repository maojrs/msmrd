import MSMRD2 as m
import numpy as np

assert m.__version__ == '0.0.1'
assert m.add(1, 2) == 3
assert m.subtract(1, 2) == -1

particles = []
for i in range(10):
    particles.append(m.particle(i, 1, 2., 3., np.random.rand(4), np.random.rand(3)))

sim = m.simulation(particles)
sim.run(0.1, 100)