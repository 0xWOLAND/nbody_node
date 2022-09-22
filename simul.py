from init_config import init_config
import numpy as np

# LCDM = Cosmology(68.0, 0.31, 0.69)
# EdS = Cosmology(70.0, 1.0, 0.0)

class Cosmology:
    def __init__(self, H0, OmegaM, OmegaL):
        self._H0 = H0
        self._OmegaM = OmegaM
        self._OmegaL = OmegaL

    @property
    def OmegaK(self):
        return 1 - self._OmegaM - self._OmegaL

    @property
    def G(self):
        return 3./2 * self._OmegaM * self._H0**2

    def da(self, a):
        return self._H0 * a * np.sqrt(self._OmegaL + self._OmegaM * a**-3 + self._OmegaK * a ** -2) 
    
    def growing_mode(self, a):
        if isinstance(a, np.ndarray):
            return np.array([self.growing_mode(b) for b in a])
        elif a <= 0.001:
            return a
        else:
            return self.factor * self.adot(a)/a \
                * quad(lambda b: self.adot(b)**(-3), 0.00001, a)[0] + 0.00001