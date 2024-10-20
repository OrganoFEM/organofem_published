
import sys,os
#sys.path.append(os.getcwd())
sys.path.append('../')
from bio.nodes import *


class ElementsTRI(dict):
    def resetElements(self):
        self.clear()

    def numberOfElements(self):
        return len(self)


class ElementsTET(dict):
    def resetElements(self):
        self.clear()

    def numberOfElements(self):
        return len(self)


class ElementTRI:
    def __init__(self, elementID, elementCoord, elementConnectivity):

        if isinstance(elementCoord, Nodes):
            self._nodes = elementCoord
            self._connectivity = elementConnectivity
            self._numberElement = elementID
            self._active = False
        else:
            #print ("Error : nodes must be an instance of class Nodes")
            sys.exit("Error : nodes must be an instance of class Nodes")

    def __repr__(self):
        return "el {} : {}".format(self.number(),  self.nodesConnectivity())

    def number(self):
        return self._numberElement

    def allNodes(self):
        return self._nodes

    def numberOfNodes(self):
        return 3

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

    def thirdNode(self):
        return self._nodes[self._connectivity[2]]

    def nodesConnectivity(self):
        return self._connectivity[0], self._connectivity[1], self._connectivity[2]


class ElementTET:
    def __init__(self, elementID, elementCoord, elementConnectivity):

        if isinstance(nodes, Nodes):
            self._nodes = elementCoord
            self._connectivity = elementConnectivity
            self._numberElement = elementID
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
        return 4

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

    def thirdNode(self):
        return self._nodes[self._connectivity[2]]

    def fourthNode(self):
        return self._nodes[self._connectivity[3]]

    def nodesConnectivity(self):
        return self._connectivity[0], self._connectivity[1], self._connectivity[2], self._connectivity[3]


if __name__ == '__main__':
    listNodes = [Node(0, 0, 0), Node(2, 0, 0),
                 Node(2, 4, 0), Node(0, 4, 0),
                 Node(3, 4, 0), Node(0, 3, 0)]
    j = 1
    nodes = Nodes()
    for i in listNodes:
        nodes[j] = i
        j += 1

    connectivity1 = [1, 2, 3, 6]
    connectivity2 = [1, 3, 4]
    elementri = ElementTRI(1, nodes, connectivity2)
    elementet = ElementTET(1, nodes, connectivity1)
    print(elementri.nodesConnectivity())
    print(elementet.nodesConnectivity())
    print(nodes)