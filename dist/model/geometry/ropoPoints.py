import sys
sys.path.append('../..')
from geometry.ropoGeometryTools import *
from geometry.ropoVector import *

from math import sqrt
import numpy as np
import copy
# import ropoGeometryTools


class ropoPoint:
    def __init__(self, x, y, z):

        self._x = x
        self._y = y
        self._z = z
        self._id = 0
        self._localLandMark = None

    def __repr__(self):
        """ Representation d'un point """
        return ("({},{},{})").format(self._x, self._y, self._z)

    def x(self):
        return self._x

    def y(self):
        return self._y

    def z(self):
        return self._z

    def setX(self, x):
        self._x = x

    def setY(self, y):
        self._y = y

    def setZ(self, z):
        self._z = z

    def getCoords(self):
        return self.x(), self.y(), self.z()

    def setId(self, id):
        self._id = id

    def id(self):
        return self._id

    def setLocalLandMark(self, vect1, vect2, vect3):
        self._localLandMark = np.array([vect1, vect2, vect3])

    def localLandMark(self):
        return self._localLandMark


class ropoPoints:
    def __init__(self):
        self._listOfPoints = []
        self._listOfPointsXCoordinate = []
        self._listOfPointsYCoordinate = []
        self._listOfPointsZCoordinate = []

    def __str__(self):
        return "{0}".format(self._listOfPoints[:])

    def __getitem__(self, key):
        return self._listOfPoints[key]

    def __setitem__(self, key, point):
        self._listOfPoints[key] = point
        self._listOfPointsXCoordinate[key] = point.x()
        self._listOfPointsYCoordinate[key] = point.y()
        self._listOfPointsZCoordinate[key] = point.z()

    def __len__(self):
        return len(self._listOfPoints)

    def addPoint(self, Point):
        self._listOfPoints.append(Point)
        self._listOfPointsXCoordinate.append(Point.x())
        self._listOfPointsYCoordinate.append(Point.y())
        self._listOfPointsZCoordinate.append(Point.z())

    def x(self):
        return self._listOfPointsXCoordinate

    def y(self):
        return self._listOfPointsYCoordinate

    def z(self):
        return self._listOfPointsZCoordinate

    def id(self):
        return [x.id() for x in self._listOfPoints]

    def getPointWithId(self, Id):
        return self._listOfPoints[Id - 1]

    def getPointXCoordWihtId(self, Id):
        return self.getPointWithId(Id).x()

    def getPointYCoordWihtId(self, Id):
        return self.getPointWithId(Id).y()

    def getPointZCoordWihtId(self, Id):
        return self.getPointWithId(Id).z()

    def distanceBetweenPoints(self):
        dist = []
        for i in range(len(self._listOfPoints) - 1):
            dist.append(ropo.modelling.surface.geometry.ropoGeometryTools.distanceBetweenPoints(
                self._listOfPoints[i], self._listOfPoints[i + 1]))
        return dist


class ropoTableOfPoint:
    def __init__(self):
        self._tableOfPoint = None

    def __str__(self):
        return "{}".format(self._tableOfPoint)

    def __getitem__(self, key):
        """key: tuple"""
        i, j = key
        return self._tableOfPoint[i, j]

    def __setitem__(self, key, Point):
        i, j = key
        self._tableOfPoint[i, j] = Point

    def shape(self):
        return np.shape(self._tableOfPoint)

    def x(self):
        nbline, nbcolum = self.shape()
        tablex = np.zeros([nbline, nbcolum])
        for i in range(nbline):
            for j in range(nbcolum):
                tablex[i, j] = self._tableOfPoint[i, j].x()
        return tablex

    def y(self):
        nbline, nbcolum = self.shape()
        tabley = np.zeros([nbline, nbcolum])
        for i in range(nbline):
            for j in range(nbcolum):
                tabley[i, j] = self._tableOfPoint[i, j].y()
        return tabley

    def z(self):
        nbline, nbcolum = self.shape()
        tablez = np.zeros([nbline, nbcolum])
        for i in range(nbline):
            for j in range(nbcolum):
                tablez[i, j] = self._tableOfPoint[i, j].z()
        return tablez

    def setTableOfPointWithArray(self, arr):
        self._tableOfPoint = arr

    def setTableOfPointWithPoints(self, Points, numberOfPointOnLayers):
        arr = np.array([])
        arr = np.append(arr, Points)
        arr = arr.reshape(
            int(len(Points) / numberOfPointOnLayers),  numberOfPointOnLayers)
        self._tableOfPoint = arr

    def reShape(self, a, b):
        return self._tableOfPoint.reshape(a, b)

    def distanceBetweenPoints(self):
        nbline, nbcolum = self.shape()
        tableDist = np.zeros([nbline - 1, nbcolum])
        for i in range(nbline - 1):
            for j in range(nbcolum):
                tableDist[i, j] = ropoGeometryTools.distanceBetweenPoints(
                    self._tableOfPoint[i, j], self._tableOfPoint[i + 1, j])
        return tableDist
