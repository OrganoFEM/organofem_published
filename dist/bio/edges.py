import sys,os
sys.path.append('../')
from bio.nodes import *



class ElementsE(dict):
    def resetElements(self):
        self.clear()

    def numberOfElements(self):
        return len(self)



class ElementE:
    def __init__(self, numberElement, nodes, connectivity):

        if isinstance(nodes, Nodes):
            self._nodes = nodes
            self._connectivity = connectivity
            self._numberElement = numberElement
            self._active = False
        else:
            print ("Error : nodes must be an instance of class Nodes")
            # sys.exit("Error : nodes must be an instance of class Nodes")

    def __repr__(self):
        return "el {} : {}".format(self.number(),  self.nodesConnectivity())

    def number(self):
        return self._numberElement

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

    def firstNode(self):
        return self._nodes[self._connectivity[0]]

    def secondNode(self):
        return self._nodes[self._connectivity[1]]

    def nodesConnectivity(self):
        return self._connectivity[0], self._connectivity[1]

