"""Contrôles des nouvelles interprétations, indépendants des grandes archives."""
import unittest
import numpy as np
from cu3au.cross_analysis import energy_channels, covariance_budget, phase_statistics
from cu3au.lattice import build_lattice
from cu3au.model import energy_int, SCALE, cowley, ordre, KB


class CrossAnalysisTests(unittest.TestCase):
    def test_energy_channels_against_direct_pair_energy(self):
        lat=build_lattice(3)
        n=len(lat.sublattice)
        rng=np.random.default_rng(1337)
        data=[]
        for _ in range(20):
            state=np.zeros(n,dtype=np.int8)
            state[rng.choice(n,n//4,replace=False)]=1
            a1,a2=cowley(state,lat.neighbors)
            data.append([energy_int(state,lat.neighbors)/SCALE,0,0,a1,a2])
        a=np.asarray(data)
        one,two=energy_channels(a,n)
        np.testing.assert_allclose(-3.813*n+one+two,a[:,0],atol=1e-10)

    def test_covariance_preserves_within_chain_variance(self):
        # Une chaîne piégée avec moyenne différente ne doit pas ajouter
        # Var(moyennes) à la chaleur thermique calculée dans chaque chaîne.
        rng=np.random.default_rng(912)
        a=np.zeros((2,200,2,9))
        a[...,3]=rng.normal(-.2,.02,(2,200,2))
        a[...,4]=.7+5*a[...,3]+rng.normal(0,.01,(2,200,2))
        a[1,...,3]+=.05
        n=108;t=np.array([900,1000])
        one,two=energy_channels(a,n)
        a[...,0]=-3.813*n+one+two
        parts=covariance_budget(a,t,n)
        expected=np.var(a[...,0],axis=1).mean(axis=0)/(n*(KB*t)**2)
        np.testing.assert_allclose(parts.sum(axis=0),expected,rtol=1e-11)
        self.assertTrue(np.all(parts[2]<0))
        pooled=np.var(a[...,0].reshape(-1,2),axis=0)/(n*(KB*t)**2)
        self.assertTrue(np.all(pooled>expected))

    def test_population_count_respects_chain_boundaries(self):
        a=np.zeros((2,4,1,9))
        a[0,:,0,2]=[.1,.8,.9,.2]
        a[1,:,0,2]=[.8,.8,.9,.9]
        fractions,changes=phase_statistics(a)
        np.testing.assert_array_equal(fractions[:,0],[.5,1])
        np.testing.assert_array_equal(changes[:,0],[2,0])


if __name__=='__main__':
    unittest.main()
