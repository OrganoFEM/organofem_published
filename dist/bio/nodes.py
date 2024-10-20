import sys
import numpy as np
#sys.path.append('../..')
from geometry.ropoGeometryTools import  *
from geometry.ropoPoints import ropoPoint


class Node:
    def __init__(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z

    def __repr__(self):
        return ("{}, {}, {}".format(self.x(), self.y(), self.z()))

    def x(self):
        return self._x

    def y(self):
        return self._y

    def z(self):
        return self._z

    def point(self):
        return ropoPoint(self.x(), self.y(), self.z())


class Nodes(dict):
    def resetNodes(self):
        self.clear()

    def numberOfNodes(self):
        return len(self)


if __name__ == '__main__':
    nod = Node(2, 3, 4)
    nod1 = Node(2, 3, 4)
    print(ropoVector(nod.point(), nod1.point()) + ropoVector(nod.point(), nod1.point()))
