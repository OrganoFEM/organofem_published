import sys
sys.path.append('../..')
from geometry.ropoPoints import *
from geometry.ropoGeometryTools import *

class ropoVector:
    def __init__(self, Point1, Point2):
        pt1 = ropoPoint(0, 0, 0)
        self._x = Point2.x() - Point1.x()
        self._y = Point2.y() - Point1.y()
        self._z = Point2.z() - Point1.z()
        self._origine = pt1

    def __str__(self):
        return "({0},{1},{2})".format(self._x, self._y, self._z)

    def __add__(self, otherVector):
        x = self.x() + otherVector.x()
        y = self.y() + otherVector.y()
        z = self.z() + otherVector.z()
        pt = Point(x, y, z)
        return Vector(self.origine(), pt)

    def __sub__(self, otherVector):
        x = self.x() - otherVector.x()
        y = self.y() - otherVector.y()
        z = self.z() - otherVector.z()
        pt = ropoPoint()
        pt.setCoords(x, y, z)
        return ropoVector(self.origine(), pt)

    def __mul__(self, otherVector):
        """ Produit scalaire """
        prod = self.x() * otherVector.x() + self.y() * \
            otherVector.y() + self.z() * otherVector.z()
        return prod

    def scalarMultiplication(self, scalar):
        x = self.x() * scalar
        y = self.y() * scalar
        z = self.z() * scalar
        pt = ropoPoint(x, y, z)
        return ropoVector(self.origine(), pt)

    def origine(self):
        return self._origine

    def x(self):
        return self._x

    def y(self):
        return self._y

    def z(self):
        return self._z

    def setZ(self, z):
        self._z = z

    def norme(self):
        return sqrt(self.x()**2 + self.y()**2 + self.z()**2)
