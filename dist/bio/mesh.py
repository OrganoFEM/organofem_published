import numpy as np
import pandas as pd
from copy import deepcopy

class Mesh:
    """Base class defining a mesh representing the connective tissue."""

    def __init__(self, identifier, datasets, specs=None, coords=None, maxbackup=5):
        """Creates an epithelium
        Parameters
        ----------
        identifier : string
        datasets : dictionary of dataframes
            The keys correspond to the different geometrical elements
            constituting the epithelium:
            * `vert` contains a dataframe of vertices,
            * `edge` contains a dataframe of *oriented* half-edges between vertices,
            * `face` contains a dataframe of polygonal faces enclosed by half-edges,
            * `cell` contains a dataframe of polyhedral cells delimited by faces,
        specs : nested dictionnary of specifications
            The first key designs the element name: (`vert`, `edge`, `face`, `cell`),
            corresponding to the respective dataframes attribute in the dataset.
            The second level keys design column names of the above dataframes.
            For exemple:
            .. code::
                specs = {
                    "face": {
                        ## Face Geometry
                        "perimeter": 0.,
                        "area": 0.,
                        ## Coordinates
                        "x": 0.,
                        "y": 0.,
                        "z": 0.,
                        ## Topology
                        "num_sides": 6,
                        ## Masks
                        "is_alive": 1},
                    "vert": {
                        ## Coordinates
                        "x": 0.,
                        "y": 0.,
                        "z": 0.,
                        ## Masks
                        "is_active": 1},
                    "edge": {
                        ## Connections
                        "srce": 0,
                        "trgt": 1,
                        "face": 0,
                        "cell": 0,
                        ## Coordinates
                        "dx": 0.,
                        "dy": 0.,
                        "dz": 0.,
                        "length": 0.,
                        ## Normals
                        "nx": 0.,
                        "ny": 0.,
                        "nz": 1.}
                    "settings":
                        ## Custom values
                        "geometry": "flat"
                    }
