import random
import numpy as np

class particle:
    def __init__(self,x,y,z,u,v,w):
        self.x = x
        self.y = y
        self.z = z
        self.u = u
        self.v = v
        self.w = w


class init_config:
    def __init__(self, count):
        self._count = count
        self._arr = []
        
    
    def init_rand(self):
        for i in range(self._count):
            x = random.random() * 5
            y = random.random() * 5
            z = random.random() * 5
            u = v = w = 0
            self._arr.append(particle(x,y,z,u,v,w)) 
        
    def get_arr(self):
        return self._arr
    
    def xy(self):
        a = np.array([(u.x, u.y) for u in self._arr])
        return [a[:, 0], a[:, 1]]
    def yz(self):
        a = np.array([(u.y, u.z) for u in self._arr])
        return [a[:, 0], a[:, 1]]
    def xz(self):
        a = np.array([(u.x, u.z) for u in self._arr])
        return [a[:, 0], a[:, 1]]
