import rebound
import unittest
import math
import rebound.data
import warnings
   
def assert_back_after_one_orbit(self,sim):
    delta = 1e-8
    self.assertLess(math.fabs((sim.particles[1].x-1.)),delta)
    self.assertLess(math.fabs((sim.particles[1].y)),delta)
    self.assertLess(math.fabs((sim.particles[1].vx)),delta)
    self.assertLess(math.fabs((sim.particles[1].vy-1.)),delta)
    self.assertLess(math.fabs((sim.particles[0].x)),delta)
    self.assertLess(math.fabs((sim.particles[0].y)),delta)
    self.assertLess(math.fabs((sim.particles[0].vx)),delta)
    self.assertLess(math.fabs((sim.particles[0].vy)),delta)

    
class TestTestParticleTypeZeroOneOrbit(unittest.TestCase):
    def test_one_orbit_ias15(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(x=1.,vy=1.,m=1.)
        sim.N_active = 1
        sim.testparticle_type = 0
        sim.integrator = "ias15"
        sim.integrate(math.pi*2.,exact_finish_time=1)
        assert_back_after_one_orbit(self,sim)
    
    def test_one_orbit_leapfrog(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(x=1.,vy=1.,m=1.)
        sim.N_active = 1
        sim.testparticle_type = 0
        sim.integrator = "leapfrog"
        sim.dt=1e-5
        sim.integrate(math.pi*2.,exact_finish_time=1)
        assert_back_after_one_orbit(self,sim)
        
    def test_one_orbit_whfast_democraticheliocentric(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(x=1.,vy=1.,m=1.)
        sim.N_active = 1
        sim.testparticle_type = 0
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.integrate(math.pi*2.,exact_finish_time=1)
        assert_back_after_one_orbit(self,sim)
        
    def test_one_orbit_whfast_jacobi(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(x=1.,vy=1.,m=1.)
        sim.N_active = 1
        sim.testparticle_type = 0
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "jacobi" # default
        sim.integrate(math.pi*2.,exact_finish_time=1)
        assert_back_after_one_orbit(self,sim)
        
    def test_one_orbit_whfast_whds(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(x=1.,vy=1.,m=1.)
        sim.N_active = 1
        sim.testparticle_type = 0
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "whds"
        sim.integrate(math.pi*2.,exact_finish_time=1)
        assert_back_after_one_orbit(self,sim)
        



if __name__ == "__main__":
    unittest.main()
