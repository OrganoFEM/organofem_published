import numpy as np
import copy

class Organoid(object): 
 
 """Defines an organoid comprised of cells which can move, divide or die"""
    def __init__(self,cell_ids,lineage,topology,properties=None):
        """ Parameters:
        cell_ids: (N,) array ints
             unique id for each cell (N is number of cells)
        lineage: (N,) array ints
            id of mother for each cell (N is number of cells)
        topology: object
            defines cell locations and neighbour connections
        properties: dict or None
            dictionary available for any other cell properties
        """
        self.topology = topology
        N = len(topology)
        self.cell_ids = cell_ids
        self.lineage = lineage
        self.properties = properties or {}

    def __len__(self):
        return len(self.mesh)
    
    def copy(self):
        """create a copy of Organoid"""
        return Organoid(self.cell_ids.copy(),self.lineage.copy(),self.topology.copy(),self.properties.copy())
            
