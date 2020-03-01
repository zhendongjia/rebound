import rebound
import unittest
import warnings

class TestMercurana(unittest.TestCase):
    def test_encounter(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=1.14,f=-0.2)
        sim.dt = 0.1
        sim.move_to_com()
        E0 = sim.calculate_energy()
        sim2 = sim.copy()
        sim3 = sim.copy()

        sim.integrator = "whfast"
        sim.integrate(50.)
        E1 = sim.calculate_energy()
        self.assertGreater((abs((E1-E0)/E0)),.1)
    
        sim2.integrator = "mercurana"
        sim2.ri_mercurana.Nmaxshells = 30
        sim2.ri_mercurana.kappa0 = 0.001
        sim2.integrate(50.)
        E1 = sim2.calculate_energy()
        self.assertGreater(sim2.ri_mercurana.Nmaxshellused,10)
        self.assertLess((abs((E1-E0)/E0)),1e-5)
   
        sim3.integrator = "mercurana"
        sim3.ri_mercurana.Nmaxshells = 30
        sim3.ri_mercurana.N_dominant = 1
        sim3.ri_mercurana.n0=30
        sim3.ri_mercurana.kappa0 = 0.001
        sim3.integrate(50.)
        E1 = sim3.calculate_energy()
        self.assertGreater(sim3.ri_mercurana.Nmaxshellused,4)
        self.assertLess((abs((E1-E0)/E0)),2e-8)
   

if __name__ == "__main__":
    unittest.main()
