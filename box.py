import numpy as np
import matplotlib.pyplot as plt

class particle:
    def __init__(self, v=np.zeros(3), mass=1):
        self.v = v
        self.mass = 1

class point:
    def __init__(self, v=np.zeros(3)):
        self.v = v


class node:
    def __init__(self):
        self.com = np.zeros(3) 
        self.total_mass = 0
        self.point = point()
        self.dim = np.zeros(3)
        self.cnt = 0
        self.adj = []
        self.body = None
    
    def insertLevel(self):
        dim = self.dim / 2
        for i in [self.point.v[0], self.point.v[0] + self.dim[0] // 2]:
            for j in [self.point.v[1], self.point.v[1] + self.dim[1] // 2]:
                for k in [self.point.v[2], self.point.v[2] + self.dim[2] // 2]:
                    n = node()
                    n.point = np.array([i,j,k])
                    n.dim = dim
                    self.adj.append(n)

    def getQuadrant(self, p):
        if(len(self.adj) == 0):
            self.insertLevel()
        mx = self.point.v[0] + self.dim[0] // 2
        my = self.point.v[1] + self.dim[1] // 2
        mz = self.point.v[2] + self.dim[2] // 2
        idx = 0
        if(p.v[0] >= mx):
            idx += 4 
        if(p.v[1] >= my):
            idx += 2 
        if(p.v[2] >= mz):
            idx += 1 
        return self.adj[idx]
    
    def __str__(self, level=0):
        ret = "\t"*level+repr(self.total_mass)+"\n"
        for child in self.adj:
            ret += child.__str__(level+1)
        return ret

class box:
    def __init__(self, dim):
        self.dim = dim
        self.tree = None 

    def insertNode(self, ptc: particle, p_tree: node):
        s = []
        s.append((p_tree, ptc))
        while(len(s) > 0):
            tree,
            if(tree == None or tree.total_mass == 0):
                tree = node()
                tree.body = p
                tree.point.v = p.v 
                tree.cnt += 1
                tree.com = (tree.com * tree.total_mass + p.v * p.mass ) / (tree.total_mass + p.mass)
                tree.total_mass += p.mass
            elif(len(tree.adj) != 0):
                tree.cnt += 1
                tree.com = (tree.com * tree.total_mass + p.v * p.mass ) / (tree.total_mass + p.mass)
                tree.total_mass += p.mass
                s.append((p, tree.getQuadrant(p)))
                # self.insertNode(p, tree.getQuadrant(p))
            else:
                tree.insertLevel()
                s.append((tree.body, tree.getQuadrant(tree.body)))
                s.append((p, tree.getQuadrant(p)))
                tree.com = (tree.com * tree.total_mass + p.v * p.mass ) / (tree.total_mass + p.mass)
                tree.total_mass += p.mass
        print(tree)
    
    def createPlot(self, tree=None):
        ax = plt.gca()
        if tree is None:
            tree = self.tree
        if tree is None: return 
        if(len(tree.adj) == 0):
            print(tree.point[0], tree.point[1])
            ax.scatter(tree.point[0], tree.point[1])
        else:
            for x in tree.adj:
                self.createPlot(x)


        
        
    