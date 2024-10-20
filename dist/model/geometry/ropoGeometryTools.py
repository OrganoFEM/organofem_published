import sys
sys.path.append('../..')
from geometry.ropoPoints import *
from geometry.ropoVector import *


def produitVectoriel(vect1, vect2):
    point = ropoPoint(vect1.y() * vect2.z() - vect1.z() * vect2.y(), vect2.x() *
                      vect1.z() - vect1.x() * vect2.z(), vect1.x() * vect2.y() - vect1.y() * vect2.x())
    vect = ropoVector(ropoPoint(0, 0, 0), point)
    return vect


def distanceBetweenPoints(point1, point2):
    dist = np.sqrt((point1.x() - point2.x())**2 + (point1.y() -
                                                   point2.y())**2 + (point1.z() - point2.z())**2)
    return dist


def squareDistanceBetweenPoints(point1, point2):
    dist = (point1.x() - point2.x())**2 + (point1.y() -
                                           point2.y())**2 + (point1.z() - point2.z())**2
    return dist


def distanceBetweenZCoordinatePoint(point1, point2):
    dist = np.abs(point2.z() - point1.z())
    return dist


def distanceBetweenXCoordinatePoint(point1, point2):
    dist = np.abs(point2.x() - point1.x())
    return dist


def distanceBetweenYCoordinatePoint(point1, point2):
    dist = np.abs(point2.y() - point1.y())
    return dist


def convertCoordToPoint(x, y, z):
    Points = ropoPoints()
    for ind, val in enumerate(x):
        pt = ropoPoint(val, y[ind], z[ind])
        Points.addPoint(pt)
    return Points


def convertToPointWithId(x, y, z):
    points = ropoPoints()
    id = 1
    for ind, val in enumerate(x):
        pt = ropoPoint(val, y[ind], z[ind])
        pt.setId(id)
        points.addPoint(pt)
        id += 1
    return points


def middleLayers(numLayBetCylScell, firtsPoints, lastPoints, partClose=False, numberOfLayersForJumpVanish=6):
    xInterieur = []
    yInterieur = []
    zInterieur = []
    sizeX = len(firtsPoints) / numberOfLayersForJumpVanish
    if partClose:
        for val1 in range(1, numLayBetCylScell, 1):
            for ind, val in enumerate(firtsPoints):
                if ind >= (numberOfLayersForJumpVanish - 1) * sizeX:
                    dist = distanceBetweenPoints(
                        firtsPoints[ind], lastPoints[0])
                    distToAdd = val1 * (dist / numLayBetCylScell)
                    pert = (distToAdd - dist) / dist
                    xInterieur.append(firtsPoints[ind].x(
                    ) + (lastPoints[0].x() - firtsPoints[ind].x()) * (1 + pert))
                    yInterieur.append(firtsPoints[ind].y(
                    ) + (lastPoints[0].y() - firtsPoints[ind].y()) * (1 + pert))
                    zInterieur.append(firtsPoints[ind].z(
                    ) + (lastPoints[0].z() - firtsPoints[ind].z()) * (1 + pert))

    else:
        for val1 in range(1, numLayBetCylScell, 1):
            for ind, val in enumerate(firtsPoints):
                if ind >= (numberOfLayersForJumpVanish - 1) * sizeX:
                    dist = distanceBetweenPoints(
                        firtsPoints[ind], lastPoints[ind - (numberOfLayersForJumpVanish - 1) * sizeX])
                    distToAdd = val1 * (dist / numLayBetCylScell)
                    pert = (distToAdd - dist) / dist
                    xInterieur.append(firtsPoints[ind].x(
                    ) + (lastPoints[ind - (numberOfLayersForJumpVanish - 1) * sizeX].x() - firtsPoints[ind].x()) * (1 + pert))
                    yInterieur.append(firtsPoints[ind].y(
                    ) + (lastPoints[ind - (numberOfLayersForJumpVanish - 1) * sizeX].y() - firtsPoints[ind].y()) * (1 + pert))
                    zInterieur.append(firtsPoints[ind].z(
                    ) + (lastPoints[ind - (numberOfLayersForJumpVanish - 1) * sizeX].z() - firtsPoints[ind].z()) * (1 + pert))
    return xInterieur, yInterieur, zInterieur


def distanceGeodesic(dictPointGeodesic, sommeDistance=False):
    "entree: prend un dictionnaire dont les cle vont de 0 a len(dictPointGeodesic)"
    dictDistGeodesic = {}
    sD = []
    for i in range(len(dictPointGeodesic)):
        dictDistGeodesic[i] = dictPointGeodesic[i].distanceBetweenPoints()
        sD.append(sum(dictPointGeodesic[i].distanceBetweenPoints()))
    if sommeDistance:
        return dictDistGeodesic, sD
    else:
        return dictDistGeodesic
