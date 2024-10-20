import sys,os
sys.path.append('../')
from bio.nodes import *

class Cells(dict):
    def resetElements(self):
        self.clear()

    def numberOfElements(self):
        return len(self)



class Cell:
    def __init__(self, dataset):
        self.datasets = datasets
        self.data_names = list(datasets.keys())
        self.element_names = ["srce", "trgt", "face", "cell"][: len(self.data_names)]

    def __repr__(self):
        return "Cell {} : {}".format(self.number(),  self.nodesConnectivity())

    def number(self):
        return self._numberCell

    def allNodes(self):
        return self._nodes

    def numberOfNodes(self):
        return 2

    def setActive(self):
        self._active = True

    def elStatus(self):
        """ Check if the element is active or not 
            In defaut all elements is inactive
        """
        return self._active


