
#general functions
import os
import os.path
import sys
sys.path.append('../../../') #TO ADD /Functions/

import organofem
from organofem.dist.io.inp import save_inp_individual, save_inp
from organofem.dist.io.csv import eptm_to_csv,read_positions_csv
from organofem.dist.io import hdf5, obj
from organofem.dist.utils import *
from organofem.dist.utils.draw import vp_view as vip
from organofem.dist.utils.generation import *
from organofem.dist.bio.elements import *
from organofem.dist.bio.edges import *
from organofem.dist.bio.cells import *
from organofem.dist.bio.nodes import *

from tyssue import Sheet, ClosedMonolayerGeometry as geom
from tyssue.dynamics import model_factory, effectors
from tyssue.core.monolayer import Monolayer
from tyssue.core import Epithelium
from tyssue import SheetGeometry, ClosedSheetGeometry
from tyssue import config


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import numpy as np
import pandas as pd
import subprocess
import math
import random
from scipy.spatial import Voronoi
from itertools import combinations
from itertools import islice
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import distance as dist
from vedo import *
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from multiprocessing import Pool
os.chdir('../../test/tmp')

####################################
# Begin Global parameters
from_imaging = False
Nc = 200 # Define the total number of cells present in our oragnoid
apical = "in" # Define apical surface "in" for inside the organoid and "out" for outside.
radius = 1500
thick = 0.25  # 259 standard==0.374 define the thickness as a ratio of the radius (thickness = radius * thick  # pixels)
vsize=0.1
# select the method of distribution of the cells 
method = "uniform" 
SolverStep = True # solve the Rehologic vertex model defined in tyssue (exported from farhadifar 207)
C3D10Convert = False # add quradratic resolution to cell elements
IncreaseResolution = True # only compatible with C3D10Convert= True / to decimate the mesh
OutsidePart = False
IndividualCells = False # used to dissociate cells and create contacts between cells
AreaStep = False # option to create area positive and negatives starting from a unique cell and with growing cell radius distance

"""
Select the method of physical input for deformation used in abaqus:
- "PushArea" allow to create a gaussian field of displacement at the surface.
Only "pushArea" is defined at this time. 
"""
LoadType = "pushArea" 
# defining files names
GlobalPartName="ORGANOID"
timestr = time.strftime("%Y%m%d-%H%M%S") # defining a unique time-id for files saving used to save singular hdf5 file
#hdf5_file = str(timestr) + ".hdf5" # to create a unique file
hdf5_file = "organo.hdf5"
obj_file1 = str(Nc) + "cells_Organoid.obj"
obj_file2 = obj_file1.replace(".obj", "_2.obj")
new_inp = obj_file1.replace(".obj", str(method) + ".inp")
new_inp2 = obj_file1.replace(".obj", str(method) + "_2.inp")
new_inp3 = obj_file1.replace(".obj", "_divisions.inp")
###
# End Global parameters
#######################################

### BEGIN OF DEFINITIONS ####################################################

def distance3D(p0, p1):
    dx, dy, dz = p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]
    return math.hypot(math.hypot(dx, dy), dz)

def spheric3D(nb=200, dim=3, method="perfect"):
    
    if method == "perfect":
        print("perfect")

        indices = np.arange(0, nb, dtype=float) + 0.5

        phi = np.arccos(1 - 2 * indices / nb)
        theta = math.pi * (1 + 5 ** 0.5) * indices
        points = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
        points = np.array(points)

    if method == "normal":
        print("normal")
        nba= int(nb*0.8)
        nbb=nb-nba
        ####
        nba=nb
        nbb=nb
        ####
        norm = np.random.normal
        normal_deviates = norm(size=(dim, nba))
        radius = np.sqrt((normal_deviates ** 2).sum(axis=0))
        points = normal_deviates / radius

        r = 1 / np.sqrt(nbb)
        points2 = []
        normal_deviates2 = []
        i = 0
        last_success = 0
        while i < nbb:
            x = random.random()
            y = random.random()
            z = random.random()

            accept = True
            for p in normal_deviates2:
                if distance3D(p, (x, y, z)) < r:
                    accept = False
                    break
            if accept is True:
                normal_deviates2.append([x, y, z])
                i += 1
                if i - last_success > nbb:
                    break
                last_success = i
        normal_deviates2 = np.array(normal_deviates2).T
        radius2 = np.sqrt((normal_deviates2).sum(axis=0))
        points2 = normal_deviates2/radius2        
        points2=np.concatenate((points, points2), axis=1)


    if method == "uniform":
        print("uniform")
        r = 1 / np.sqrt(nb)
        points = []
        normal_deviates = []
        i = 0
        last_success = 0
        while i < nb:
            x = random.random()
            y = random.random()
            z = random.random()

            accept = True
            for p in normal_deviates:
                if distance3D(p, (x, y, z)) < r:
                    accept = False
                    break
            if accept is True:
                normal_deviates.append([x, y, z])
                i += 1
                if i - last_success > nb:
                    break
                last_success = i
        normal_deviates = np.array(normal_deviates).T
        radius = np.sqrt((normal_deviates).sum(axis=0))
        points = normal_deviates / radius
    return points

## Modelisation of a standard Organoid defined by
def spherical_sheet_julien(radius, centers, **kwargs):
    """Returns a spherical sheet with the given radius and (approximately)
    the given number of cells
    """
    # center the origin
    centers[:, 0] -= np.mean(centers[:, 0])
    centers[:, 1] -= np.mean(centers[:, 1])
    centers[:, 2] -= np.mean(centers[:, 2])
    # normalize organoid size
    centers[:, 0] /= np.max(centers[:, 0])
    centers[:, 1] /= np.max(centers[:, 1])
    centers[:, 2] /= np.max(centers[:, 2])
    
    eptm = sheet_from_cell_centers_julien(centers, **kwargs)    
    rhos = (eptm.vert_df[eptm.coords] ** 2).sum(axis=1).mean()
    ClosedSheetGeometry.scale(eptm, radius / rhos, eptm.coords)
    ClosedSheetGeometry.update_all(eptm)

    return eptm


def extrude(apical_datasets, method="homotecy", scale=0.3, vector=[0, 0, -1]):
    """Extrude a sheet to form a monlayer epithelium
    Parameters
    ----------
    * apical_datasets: dictionnary of three DataFrames,
    'vert', 'edge', 'face'
    * method: str, optional {'homotecy'|'translation'|'normals'}
    default 'homotecy'
    * scale: float, optional
    the scale factor for homotetic scaling, default 0.3.
    * vector: sequence of three floats, optional,
    used for the translation
    default [0, 0, -1]
    if `method == 'homotecy'`, the basal layer is scaled down from the
    apical one homoteticaly w/r to the center of the coordinate
    system, by a factor given by `scale`
    if `method == 'translation'`, the basal vertices are translated from
    the apical ones by the vector `vect`
    if `method == 'normals'`, basal vertices are translated from
    the apical ones along the normal of the surface at each vertex,
    by a vector whose size is given by `scale`
    """
    apical_vert = apical_datasets["vert"]
    apical_face = apical_datasets["face"]
    apical_edge = apical_datasets["edge"]

    apical_vert["segment"] = "apical"
    apical_face["segment"] = "apical"
    apical_edge["segment"] = "apical"

    coords = list("xyz")
    datasets = {}

    Nv = apical_vert.index.max() + 1
    Ne = apical_edge.index.max() + 1
    Nf = apical_face.index.max() + 1

    basal_vert = apical_vert.copy()

    basal_vert.index = basal_vert.index + Nv
    basal_vert["segment"] = "basal"

    cell_df = apical_face[coords].copy()
    cell_df.index.name = "cell"
    cell_df["is_alive"] = 1

    basal_face = apical_face.copy()
    basal_face.index = basal_face.index + Nf
    basal_face[coords] = basal_face[coords] * 1 / 3.0
    basal_face["segment"] = "basal"
    basal_face["is_alive"] = 1

    apical_edge["cell"] = apical_edge["face"]
    basal_edge = apical_edge.copy()
    # ## Flip edge so that normals are outward
    basal_edge[["srce", "trgt"]] = basal_edge[["trgt", "srce"]] + Nv
    basal_edge["face"] = basal_edge["face"] + Nf
    basal_edge.index = basal_edge.index + Ne
    basal_edge["segment"] = "basal"

    lateral_face = pd.DataFrame(
        index=apical_edge.index + 2 * Nf, columns=apical_face.columns
    )
    lateral_face["segment"] = "lateral"
    lateral_face["is_alive"] = 1

    lateral_edge = pd.DataFrame(
        index=np.arange(2 * Ne, 6 * Ne), columns=apical_edge.columns
    )

    lateral_edge["cell"] = np.repeat(apical_edge["cell"].values, 4)
    lateral_edge["face"] = np.repeat(lateral_face.index.values, 4)
    lateral_edge["segment"] = "lateral"

    lateral_edge.loc[np.arange(2 * Ne, 6 * Ne, 4), ["srce", "trgt"]] = apical_edge[
        ["trgt", "srce"]
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 1, 6 * Ne, 4), "srce"] = apical_edge[
        "srce"
    ].values
    lateral_edge.loc[np.arange(2 * Ne + 1, 6 * Ne, 4), "trgt"] = basal_edge[
        "trgt"
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 2, 6 * Ne, 4), ["srce", "trgt"]] = basal_edge[
        ["trgt", "srce"]
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 3, 6 * Ne, 4), "srce"] = basal_edge[
        "srce"
    ].values
    lateral_edge.loc[np.arange(2 * Ne + 3, 6 * Ne, 4), "trgt"] = apical_edge[
        "trgt"
    ].values

    if method == "homotecy":
        basal_vert[coords] = basal_vert[coords] * scale
    elif method == "translation":
        for c, u in zip(coords, vector):
            basal_vert[c] = basal_vert[c] + u
    elif method == "normals":
        field = apical_edge.groupby("srce")[["nx", "ny", "nz"]].mean()
        field = -field.values * scale / np.linalg.norm(field, axis=1)[:, None]
        basal_vert[coords] = basal_vert[coords] + field
    else:
        raise ValueError(
            """
`method` argument not understood, supported values are
'homotecy', 'translation' or 'normals'
        """
        )

    datasets["cell"] = cell_df
    datasets["vert"] = pd.concat([apical_vert, basal_vert])
    datasets["vert"]["is_active"] = 1
    datasets["edge"] = pd.concat([apical_edge, basal_edge, lateral_edge])
    datasets["face"] = pd.concat([apical_face, basal_face, lateral_face])
    datasets["edge"]["is_active"] = 1
    specs = bulk_spec.copy()

    for elem in ["vert", "edge", "face", "cell"]:
        datasets[elem].index.name = elem
        for col, value in specs[elem].items():
            if not col in datasets[elem]:
                datasets[elem][col] = value

    if (method == "normals") and (scale < 0):
        datasets["edge"][["srce", "trgt"]] = datasets["edge"][["trgt", "srce"]]
    return datasets

def appendRestrictedLoadDeformation(nodesA, nodesB):
    import operator

    distAB = {}
    distBA = {}
    for i, j in nodesA.items():
        c = list(map(operator.sub, nodesB[i], j))
        distAB[i] = c
    for i, j in nodesA.items():
        c = list(map(operator.sub, j, nodesB[i]))
        distBA[i] = c
    return distAB, distBA

def appendLoadDeformation(nodesA, nodesB, inputFileName):
    import operator

    distAB = {}
    for i, j in nodesA.items():
        c = list(map(operator.sub, nodesB[i], j))
        distAB[i] = c
    return distAB

def create_anchors(sheet):
    """Adds an edge linked to every vertices at the boundary
    and create anchor vertices
    """
    anchor_specs = {
        "face": {"at_border": 0},
        "vert": {"at_border": 0, "is_anchor": 0},
        "edge": {"at_border": 0, "is_anchor": 0},
    }

    sheet.update_specs(anchor_specs)
    # ## Edges with no opposites denote the boundary

    free_edge = sheet.edge_df[sheet.edge_df["opposite"] == -1]
    free_vert = sheet.vert_df.loc[free_edge["srce"]]
    free_face = sheet.face_df.loc[free_edge["face"]]

    sheet.edge_df.loc[free_edge.index, "at_border"] = 1
    sheet.vert_df.loc[free_vert.index, "at_border"] = 1
    sheet.face_df.loc[free_face.index, "at_border"] = 1

    # ## Make a copy of the boundary vertices
    anchor_vert_df = free_vert.reset_index(drop=True)
    anchor_vert_df[sheet.coords] = anchor_vert_df[sheet.coords] * 1.01
    anchor_vert_df.index = anchor_vert_df.index + sheet.Nv
    anchor_vert_df["is_anchor"] = 1
    anchor_vert_df["at_border"] = 0
    anchor_vert_df["is_active"] = 0

    sheet.vert_df = pd.concat([sheet.vert_df, anchor_vert_df])
    sheet.vert_df.index.name = "vert"
    anchor_edge_df = pd.DataFrame(
        index=np.arange(sheet.Ne, sheet.Ne + free_vert.shape[0]),
        columns=sheet.edge_df.columns,
    )

    anchor_edge_df["srce"] = free_vert.index
    anchor_edge_df["trgt"] = anchor_vert_df.index
    anchor_edge_df["line_tension"] = 0
    anchor_edge_df["is_anchor"] = 1
    anchor_edge_df["face"] = -1
    anchor_edge_df["at_border"] = 0
    sheet.edge_df = pd.concat([sheet.edge_df, anchor_edge_df], sort=True)
    sheet.edge_df.index.name = "edge"
    sheet.reset_topo()


def subdivide_faces(eptm, faces):
    """Adds a vertex at the center of each face, and returns a
    new dataset
    Parameters
    ----------
    eptm: a :class:`Epithelium` instance
    faces: list,
     indices of the faces to be subdivided
    Returns
    -------
    new_dset: dict
      a dataset with the new faces devided
    """

    face_df = eptm.face_df.loc[faces]

    remaining = eptm.face_df.index.delete(faces)
    untouched_faces = eptm.face_df.loc[remaining]
    edge_df = pd.concat([eptm.edge_df[eptm.edge_df["face"] == face] for face in faces])
    verts = set(edge_df["srce"])
    vert_df = eptm.vert_df.loc[verts]

    Nsf = face_df.shape[0]
    Nse = edge_df.shape[0]

    eptm.vert_df["subdiv"] = 0
    untouched_faces["subdiv"] = 0
    eptm.edge_df["subdiv"] = 0

    new_vs_idx = pd.Series(np.arange(eptm.Nv, eptm.Nv + Nsf), index=face_df.index)
    upcast_new_vs = new_vs_idx.loc[edge_df["face"]].values

    new_vs = pd.DataFrame(
        index=pd.Index(np.arange(eptm.Nv, eptm.Nv + Nsf), name="vert"),
        columns=vert_df.columns,
    )
    new_es = pd.DataFrame(
        index=pd.Index(np.arange(eptm.Ne, eptm.Ne + 2 * Nse), name="edge"),
        columns=edge_df.columns,
    )

    new_vs["subdiv"] = 1
    # new_fs['subdiv'] = 1
    new_es["subdiv"] = 1
    if "cell" in edge_df.columns:
        new_es["cell"] = np.concatenate([edge_df["cell"], edge_df["cell"]])
    new_vs[eptm.coords] = face_df[eptm.coords].values
    # eptm.edge_df.loc[edge_df.index, 'face'] = new_fs.index
    # new_es['face'] = np.concatenate([new_fs.index,
    #                                  new_fs.index])
    new_es["face"] = np.concatenate([edge_df["face"], edge_df["face"]])
    new_es["srce"] = np.concatenate([edge_df["trgt"].values, upcast_new_vs])
    new_es["trgt"] = np.concatenate([upcast_new_vs, edge_df["srce"].values])
    new_dset = {
        "edge": pd.concat([eptm.edge_df, new_es]),
        "face": eptm.face_df,  # pd.concat([untouched_faces, new_fs]),
        "vert": pd.concat([eptm.vert_df, new_vs]),
    }

    if "cell" in edge_df.columns:
        new_dset["cell"] = eptm.cell_df

    return new_dset

bulk_spec={
    "vert": {
        "z": 0.0,
        "x": 0.0,
        "is_active": 1,
        "y": 0.0
        },
    "face": {
        "z": 0.0,
        "x": 0.0,
        "num_sides": 6,
        "area": 0.0,
        "perimeter": 0.0,
        "is_alive": 1,
        "y": 0.0
        },
    "cell": {
        "z": 0.0,
        "x": 0.0,
        "area": 0.0,
        "vol": 0.0,
        "num_faces": 6,
        "is_alive": 1,
        "y": 0.0
        },
    "edge": {
        "dz": 0.0,
        "ny": 0.0,
        "dx": 0.0,
        "nx": 0.0,
        "length": 0.0,
        "sub_vol": 0.0,
        "sub_area": 0.0,
        "srce": 0,
        "face": 0,
        "cell": 0,
        "dy": 0.0,
        "trgt": 0,
        "nz": 0.0
        }
    }

def get_outer_sheet(eptm):
    """Return a Sheet object formed by all the faces w/o an opposite
    face.
    """
    eptm.get_opposite_faces() # tag -1 opposite face
    
    is_free_face = eptm.face_df["opposite"] == -1
    is_free_edge = eptm.upcast_face(is_free_face)
    edge_df = eptm.edge_df[is_free_edge].copy()
    face_df = eptm.face_df[is_free_face].copy()
    vert_df = eptm.vert_df.loc[edge_df["srce"].unique()].copy()

    datasets = {"edge": edge_df, "face": face_df, "vert": vert_df}
    specs = {k: eptm.specs.get(k, {}) for k in ["face", "edge", "vert", "settings"]}

    return Sheet(eptm.identifier + "outer", datasets, specs)



def get_sub_eptm(eptm, edges, copy=False):
    """
    Define sub-epithelium corresponding to the edges.
    Parameters
    ----------
    eptm: a :class:`Epithelium` instance
    edges: list of edges includes in the sub-epithelium
    Returns
    -------
    sub_eptm: a :class:`Epithelium` instance
    """

    datasets = {}
    edge_df = eptm.edge_df.loc[edges]
    if edge_df.empty:
        warnings.warn("Sub epithelium appears to be empty")
        return None
    datasets["edge"] = edge_df
    datasets["vert"] = eptm.vert_df.loc[np.unique(edge_df["srce"])]
    datasets["face"] = eptm.face_df.loc[np.unique(edge_df["face"])]
    if "cell" in eptm.datasets:
        datasets["cell"] = eptm.cell_df.loc[np.unique(edge_df["cell"])]

    if copy:
        for elem, df in datasets.items():
            datasets[elem] = df.copy()

    sub_eptm = Epithelium("sub", datasets, eptm.specs)
    sub_eptm.datasets["edge"]["edge_o"] = edges
    sub_eptm.datasets["edge"]["srce_o"] = edge_df["srce"]
    sub_eptm.datasets["edge"]["trgt_o"] = edge_df["trgt"]
    sub_eptm.datasets["edge"]["face_o"] = edge_df["face"]
    if "cell" in eptm.datasets:
        sub_eptm.datasets["edge"]["cell_o"] = edge_df["cell"]

    sub_eptm.datasets["vert"]["srce_o"] = np.unique(edge_df["srce"])
    sub_eptm.datasets["face"]["face_o"] = np.unique(edge_df["face"])
    if "cell" in eptm.datasets:
        sub_eptm.datasets["cell"]["cell_o"] = np.unique(edge_df["cell"])

    sub_eptm.reset_index()
    sub_eptm.reset_topo()
    return sub_eptm

def single_cell(eptm, cell, copy=False):
    """
    Define epithelium instance for all element to a define cell.
    Parameters
    -------
    eptm : a :class:`Epithelium` instance
    cell : identifier of a cell
    copy : bool, default `False`
    Returns
    -------
    sub_etpm: class:'Epithelium' instance corresponding to the cell
    """
    edges = eptm.edge_df[eptm.edge_df["cell"] == cell].index
    return get_sub_eptm(eptm, edges, copy)


def from_3d_voronoi(voro):

    specs3d = bulk_spec
    el_idx = []

    for f_idx, (rv, rp) in enumerate(zip(voro.ridge_vertices, voro.ridge_points)):

        if -1 in rv:
            continue
        face_verts = voro.vertices[rv]
        f_center = face_verts.mean(axis=0)
        c0 = voro.points[rp[0]]
        ctof = f_center - c0

        for rv0, rv1 in zip(rv, np.roll(rv, 1, axis=0)):
            fv0 = voro.vertices[rv0]
            fv1 = voro.vertices[rv1]
            edge_v = fv1 - fv0
            fto0 = fv0 - f_center
            normal = np.cross(fto0, edge_v)
            dotp = np.dot(ctof, normal)
            if np.sign(dotp) > 0:
                el_idx.append([rv0, rv1, f_idx, rp[0]])
                el_idx.append([rv1, rv0, f_idx, rp[1]])
            else:
                el_idx.append([rv1, rv0, f_idx, rp[0]])
                el_idx.append([rv0, rv1, f_idx, rp[1]])

    el_idx = np.array(el_idx)
    
    #creation of vertices
    nodes = Nodes()
    j=1
    for i in voro.vertices:
        nodes[j] = list(i)
        j += 1
    
    #creation of elements
    elements = ElementsE()
    j=1
    for i in el_idx:
        elements[j] = list(i)
        j += 1
    
    outedges={}
    for i, el in elements.items():
        if el[-1]==0:
            outedges[i-1]=el
    """
    print("outedges")
    print(outedges)
    print(len(outedges))
    """
    cells = Cells()
    j=1
    for i in voro.points:
        cells[j] = list(i)
        j += 1


    coords = ["x", "y", "z"]
    edge_idx = pd.Index(range(el_idx.shape[0]), name="edge")
    
    
    edge_df = make_df(edge_idx, specs3d["edge"])

    for i, elem in enumerate(["srce", "trgt", "face", "cell"]):
        edge_df[elem] = el_idx[:, i]

    vert_idx = pd.Index(range(voro.vertices.shape[0]), name="vert")

    vert_df = make_df(vert_idx, specs3d["vert"])
    vert_df[coords] = voro.vertices
    included_verts = edge_df["srce"].unique()
    vert_df = vert_df.loc[included_verts].copy()
    #print(vert_df) #### faire une sortie ici des vertices
    cell_idx = pd.Index(range(voro.points.shape[0]), name="cell")
    cell_df = make_df(cell_idx, specs3d["cell"])
    cell_df[coords] = voro.points
    included_cells = edge_df["cell"].unique()
    cell_df = cell_df.loc[included_cells].copy()
    #print(cell_df) #### faire une sortie ici des cells

    nfaces = len(voro.ridge_vertices)
    face_idx = pd.Index(np.arange(nfaces), name="face")
    face_df = make_df(face_idx, specs3d["face"])
    included_faces = edge_df["face"].unique()
    face_df = face_df.loc[included_faces].copy()
    #print(face_df) #### faire une sortie ici des faces

    edge_df.sort_values(by="cell", inplace=True)
    #print(edge_df) #### faire une sortie ici des edges

    datasets = {"vert": vert_df, "edge": edge_df, "face": face_df, "cell": cell_df}
    #return cells, elements, nodes
    return datasets

def ReadInp(fi, lineCount):
    s = fi.readline()
    if s != "":
        lineCount = lineCount + 1
    # endif
    return s, lineCount

def read_elsets(fin):
    elsetCellIds = []
    elsetOther = []
    elsets = {}
    elsetsO = {}
    iLine = 0
    with open(fin, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")
            if "*elset" in myLine:
                data = myLine.split(",")
                data = data[1].split("=")
                nsetName = data[1]
                nsetNdsIds = []
                try:
                    cellIds = int("".join(filter(str.isdigit, nsetName)))
                    elsetCellIds.append(int(cellIds))
                    elsets.setdefault(cellIds, [])
                    while True:
                        iLine += 1
                        try:
                            myLine = lines[iLine]
                        except:
                            break
                        myLine = myLine.strip()
                        myLine = myLine.replace(" ", "")
                        myLine = myLine.replace("\n", "")
                        data = myLine.split(",")
                        try:
                            for ndId in data:
                                elsets[cellIds].append(int(ndId))
                        except:
                            iLine -= 1
                            break
                except:
                    elsetOther.append(nsetName)
                    elsetsO.setdefault(nsetName, [])
                    while True:
                        iLine += 1
                        try:
                            myLine = lines[iLine]
                        except:
                            break
                        myLine = myLine.strip()
                        myLine = myLine.replace(" ", "")
                        myLine = myLine.replace("\n", "")
                        data = myLine.split(",")
                        try:
                            for ndId in data:
                                elsetsO[nsetName].append(int(ndId))
                        except:
                            iLine -= 1
                            break
            else:
                iLine += 1
    return elsets, elsetsO, elsetOther

def read_elems_types(fin):
    elemsSTRI65 = {}
    elemsS3RT = {}
    iLine = 0
    with open(fin, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")
            if "STRI65" in myLine:
                while True:
                    iLine += 1
                    myLine = lines[iLine]
                    myLine = myLine.strip()
                    myLine = myLine.replace(" ", "")
                    myLine = myLine.replace("\n", "")
                    data = myLine.split(",")
                    try:
                        elemsSTRI65[int(data[0])]=[int(data[1]),int(data[2]),int(data[3]),int(data[4]),int(data[5]),int(data[6])]
                    except:
                        break
            elif "S3RT" in myLine:
                while True:
                    iLine += 1
                    myLine = lines[iLine]
                    myLine = myLine.strip()
                    myLine = myLine.replace(" ", "")
                    myLine = myLine.replace("\n", "")
                    data = myLine.split(",")
                    try:
                        elemsS3RT[int(data[0])]=[int(data[1]),int(data[2]),int(data[3])]
                    except:
                        break
            else:
                iLine += 1
    return elemsSTRI65, elemsS3RT

def read_nsets(fin):

    nsets= {}
    iLine = 0
    with open(fin, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")
            
            if "*nset" in myLine:
                data = myLine.split(",")
                data = data[1].split("=")
                nsetName = data[1]
                nsetNdsIds = []
                while True:
                    iLine += 1
                    try:
                        myLine = lines[iLine]
                    except:
                        break
                    myLine = myLine.strip()
                    myLine = myLine.replace(" ", "")
                    myLine = myLine.replace("\n", "")
                    data = myLine.split(",")
                    try:
                        for ndId in data:
                            nsetNdsIds.append(int(ndId))
                    except:
                        iLine -= 1
                        break
                nsets[nsetName] = nsetNdsIds
            else:
                iLine += 1
    return nsets

def WriteIntermediatePoint(a, b, elements, nodes, median, x, fo, MaxNode, s2):
    ffout = open(fo,"a")
    medianK1 = str(s2[a]+"-"+s2[b])
    medianK2 = str(s2[b]+"-"+s2[a])
    if medianK1 in median:
        fid = int(median[medianK1][0])
        elements[x].append(int(fid))
    elif medianK2 in median:
        fid = int(median[medianK2][0])
        elements[x].append(int(fid))   
    else:
        MaxNode += 1
        elements[x].append(int(MaxNode))
        fx = float((nodes[int(s2[a])][0]) + (nodes[int(s2[b])][0])) / 2
        fy = float((nodes[int(s2[a])][1]) + (nodes[int(s2[b])][1])) / 2
        fz = float((nodes[int(s2[a])][2]) + (nodes[int(s2[b])][2])) / 2
        median.setdefault(medianK1, [])
        median[medianK1].append(int(MaxNode))
        median[medianK1].append(float(fx))
        median[medianK1].append(float(fy))
        median[medianK1].append(float(fz))
        median.setdefault(medianK2, [])
        median[medianK2].append(int(MaxNode))
        median[medianK2].append(float(fx))
        median[medianK2].append(float(fy))
        median[medianK2].append(float(fz))
        ffout.write("%d, %f, %f, %f\n" % (int(MaxNode), float(fx), float(fy), float(fz)))
        # endif
    ffout.close()
    return elements, median, MaxNode

def readNodesINP(inputFileName):
    print("read nodes from file %s" % inputFileName)
    nodes = {}
    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ","")
            myLine = myLine.replace("\n","")
            if myLine.lower() == "*node":
                break
            iLine += 1
        while True:
            iLine += 1
            try:
                myLine = lines[iLine]
                if len(myLine) == 1:
                    continue                     
            except:
                break
            myLine = myLine.strip()
            myLine = myLine.replace(" ","")
            myLine = myLine.replace("\n","")
            data = myLine.split(",")
            try:
                ndId = int(data[0])
                nodes.setdefault(ndId, [])
                nodes[ndId].append(float(data[1]))
                nodes[ndId].append(float(data[2]))
                nodes[ndId].append(float(data[3]))
            except:
                iLine -= 1
                break
    return nodes

def readElementsINP(inputFileName):
    print("read elements from file %s" % inputFileName)
    elements = {}
    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ","")
            myLine = myLine.replace("\n","")
            if myLine[:8].lower() == "*element":
                break
            iLine += 1
        while True:
            iLine += 1
            try:
                myLine = lines[iLine]
                if len(myLine) == 1:
                    continue                     
            except:
                break
            myLine = myLine.strip()
            myLine = myLine.replace(" ","")
            myLine = myLine.replace("\n","")
            data = myLine.split(",")
            try:
                elId = int(data[0])
                elements.setdefault(elId, [])
                elements[elId].append(int(data[1]))
                elements[elId].append(int(data[2]))
                elements[elId].append(int(data[3]))
                elements[elId].append(int(data[4]))
            except:
                iLine -= 1
                break
    return elements

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def DecimateINP(fin, fout, elemDim):
    # global lineCount
    lineCount = 0
    # ------------File----------------------------------------
    fi = open(fin, "r")
    fo = open(fout, "w+")
    # ----read  until '*Node---------------------------------
    s1 = ""
    while True:
        s, lineCount = ReadInp(fi, lineCount)
        if s == "":
            break
        s1 = s.strip()  # chomp
        if s1 == "*node":
            break
        # endif
        fo.write(s)
    # ----read  until '*Element---------------------------------
    # and collect nodes in dictionary "nodes"
    fo.write("*node\n")
    s1 = ""
    nodes = {}
    while True:
        s, lineCount = ReadInp(fi, lineCount)
        if s == "":
            break
        if len(s)==1:
            break
        s1 = s.strip()  # chomp
        s2 = s1.split(",")
        i = s1.find("*element")
        if i != -1:
            break
        # endif
        fo.write(
            "%d, %f, %f, %f\n" % (int(s2[0]), float(s2[1]), float(s2[2]), float(s2[3]))
        )
        x = int(s2[0])
        nodes.setdefault(x, [])
        nodes[x].append(float(s2[1]))
        nodes[x].append(float(s2[2]))
        nodes[x].append(float(s2[3]))
        MaxNode = int(s2[0])
    # end while
    #i = s2[1].find("C3D10M")
    fo.close()
    elset_cells, elsetsOs, elsetOthers=read_elsets(fin)
    nsets=read_nsets(fin)
    if elemDim ==4:
        elemsStri, elemsS3RT =read_elems_types(fin)



    apicalL = []
    basalL = []
    lateralL = []

    # don't touch to CF (center faces) but increment CN (cluster nodes)

    for i,j in nsets.items():
        if "apical" in i:
            apicalL=j
        if "basal" in i:
            basalL=j
        if "lateral" in i:
            lateralL=j



    oldMaxNode = MaxNode

    if elemDim == 3:
        while True:
            s2 = s1.split(",")
            if s == "":
                break
            for i in range(len(s2)):
                iT = s2[i].find("type=")
                if iT != -1:
                    iType = i
                    nType = iT
            # end For
            sType = s2[iType][nType + 5 :]
            iNodes = sType.find("D")
            # -----------------------------------------------------C3D4-------------------
            elements = {}
            median = {}
            if sType.find("C3D4") != -1:
                print("###3D-C3 found")
                while True:
                    s, lineCount = ReadInp(fi, lineCount)
                    if s == "":
                        print("empty")
                        break
                    #if len(s)==1:
                    #    continue
                    s1 = s.strip()  # read line
                    i = s1.find("elset")
                    if i != -1:
                        print("elset")
                        break
                    i = s1.find("nset")
                    if i != -1:
                        print("nset")
                        break
                    
                    s2 = s1.split(",")
                    x = int(s2[0])
                    elements.setdefault(x, [])
                    for i in range(4): # add the first 4 nodes to the element
                        elements[x].append(int(s2[i+1])) 
                    elements, median, MaxNode = WriteIntermediatePoint(
                        1, 2, elements, nodes, median, x, fout, MaxNode, s2
                    ) # create node 5
                    elements, median, MaxNode = WriteIntermediatePoint(
                        2, 3, elements, nodes, median, x, fout, MaxNode, s2
                    ) # create node 6
                    elements, median, MaxNode = WriteIntermediatePoint(
                        1, 3, elements, nodes, median, x, fout, MaxNode, s2
                    ) # create node 7
                    elements, median, MaxNode = WriteIntermediatePoint(
                        1, 4, elements, nodes, median, x, fout, MaxNode, s2
                    ) # create node 8
                    elements, median, MaxNode = WriteIntermediatePoint(
                        2, 4, elements, nodes, median, x, fout, MaxNode, s2
                    ) # create node 9
                    elements, median, MaxNode = WriteIntermediatePoint(
                        3, 4, elements, nodes, median, x, fout, MaxNode, s2
                    ) # create node 10
                sType = elemPreWord + sType[iNodes:] + elemPostWord
                foutt = open(fout,"a+")
                ss = "*element, type=C3D10M"
                foutt.write(ss)
                foutt.write("\n")
                for i in elements:
                    foutt.write(
                        "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n"
                        % (
                            int(i),
                            int(elements[i][0]),
                            int(elements[i][1]),
                            int(elements[i][2]),
                            int(elements[i][3]),
                            int(elements[i][4]),
                            int(elements[i][5]),
                            int(elements[i][6]),
                            int(elements[i][7]),
                            int(elements[i][8]),
                            int(elements[i][9]),
                        )
                    )
                foutt.close()
                break
    fi.close()
    fout2 = open(fout,"a")

    # write cell elsets
    for i, j in elset_cells.items():
        print(int(i))
        fout2.write("*elset,  elset=Cell%s\n" % int(i))
        for x in range(len(j)):
            if x + 1  == len(j):
                if x % 16 == 15:
                    fout2.write("\n")
                #print(j[x])
                fout2.write(str(j[x])+"\n")
            elif x % 16 == 15:
                #print(j[x])
                fout2.write(str(j[x])+"\n")
            else:
                #print(j[x])
                fout2.write(str(j[x])+",")
    apicalS6 = []
    basalS6 = []
    lateralS6 = []


    #CASE ElEM_DIM=3
    for i in elements:
        ca = intersection(elements[i], apicalL)
        cb = intersection(elements[i], basalL)
        cl = intersection(elements[i], lateralL)
        if len(ca):
            apicalS6.append(
                [
                    int(elements[i][1]),
                    int(elements[i][2]),
                    int(elements[i][3]),
                    int(elements[i][5]),
                    int(elements[i][9]),
                    int(elements[i][8]),
                ]
            )
        elif len(cb):
            basalS6.append(
                [
                    int(elements[i][1]),
                    int(elements[i][2]),
                    int(elements[i][3]),
                    int(elements[i][5]),
                    int(elements[i][9]),
                    int(elements[i][8]),
                ]
            )
        elif len(cl):
            lateralS6.append(
                [
                    int(elements[i][1]),
                    int(elements[i][2]),
                    int(elements[i][3]),
                    int(elements[i][5]),
                    int(elements[i][9]),
                    int(elements[i][8]),
                ]
            )


    apicalEl = []
    apicalS = []
    print("----------------------------------")
    print("# of apical Centroids -> %s " % len(apicalL))
    for i, j in elements.items():
        c = intersection(j, apicalL)
        if len(c):
            apicalEl.append(i)
            surf = [elements[i][1], elements[i][2], elements[i][3]]
            if surf in apicalS:
                pass
            else:
                apicalS.append(surf)


    basalEl = []
    basalS = []
    print("----------------------------------")
    print("# of basal Centroids -> %s " % len(basalL))
    for i, j in elements.items():
        c = intersection(j, basalL)
        if len(c):
            basalEl.append(i)
            surf = [elements[i][1], elements[i][2], elements[i][3]]
            if surf in basalS:
                pass
            else:
                basalS.append(surf)

    lateralEl = []
    lateralS = []
    print("----------------------------------")
    print("# of lateral Centroids -> %s " % len(lateralL))
    for i, j in elements.items():
        c = intersection(j, lateralL)
        if len(c):
            lateralEl.append(i)
            surf = [elements[i][1], elements[i][2], elements[i][3]]
            if surf in lateralS:
                pass
            else:
                lateralS.append(surf)
    


    MaxEl = len(elements)

    ss = "*element, type=STRI65\n"  # 6-node triangular thin shell, using five degrees of freedom per node equivalent of "triangle6" in meshio
    fout2.write(ss)


    MaxElAp6 = MaxEl + len(apicalS6)
    MaxElBa6 = MaxElAp6 + len(basalS6)
    MaxElS6 = MaxElBa6 + len(lateralS6)

    for i, j in enumerate(apicalS6):
        fout2.write(
            "%d, %d, %d, %d, %d, %d, %d\n"
            % (

                int(MaxEl + i + 1),
                int(j[0]),
                int(j[1]),
                int(j[2]),
                int(j[3]),
                int(j[4]),
                int(j[5]),
            )
        )


    for i, j in enumerate(basalS6):
        fout2.write(
            "%d, %d, %d, %d, %d, %d, %d\n"
            % (
                int(MaxElAp6 + i + 1),
                int(j[0]),
                int(j[1]),
                int(j[2]),
                int(j[3]),
                int(j[4]),
                int(j[5]),
            )
        )


    for i, j in enumerate(lateralS6):
        fout2.write(
            "%d, %d, %d, %d, %d, %d, %d\n"
            % (
                int(MaxElBa6 + i + 1),
                int(j[0]),
                int(j[1]),
                int(j[2]),
                int(j[3]),
                int(j[4]),
                int(j[5]),
            )
        )


    cnodeApical=[x for l in apicalS6 for x in l]
    cnodeBasal=[x for l in basalS6 for x in l]
    cnodeLateral=[x for l in lateralS6 for x in l]


    apicalT2 = []
    basalT2 = []
    lateralT2 = []
    for f in apicalS6:
        apicalT2.append([f[0],f[3],f[5]])
        apicalT2.append([f[3],f[1],f[4]])
        apicalT2.append([f[5],f[4],f[2]])
        apicalT2.append([f[3],f[4],f[5]])    
    for f in basalS6:
        basalT2.append([f[0],f[3],f[5]])
        basalT2.append([f[3],f[1],f[4]])
        basalT2.append([f[5],f[4],f[2]])
        basalT2.append([f[3],f[4],f[5]])   
    for f in lateralS6:
        lateralT2.append([f[0],f[3],f[5]])
        lateralT2.append([f[3],f[1],f[4]])
        lateralT2.append([f[5],f[4],f[2]])
        lateralT2.append([f[3],f[4],f[5]])   
    ss = "*element, type=S3RT\n"  # 6-node triangular thin shell, using five degrees of freedom per node equivalent of "triangle6" in meshio
    fout2.write(ss)
    MaxelT2A = int(MaxElS6)
    for i, j in enumerate(apicalT2):
        MaxelT2A += 1
        fout2.write("%d, %d, %d, %d\n" % (int(MaxelT2A), int(j[0]), int(j[1]), int(j[2])))
    MaxelT2B = int(MaxelT2A)
    for i, j in enumerate(basalT2):
        MaxelT2B += 1
        fout2.write("%d, %d, %d, %d\n" % (int(MaxelT2B), int(j[0]), int(j[1]), int(j[2])))
    MaxelT2L = int(MaxelT2B)
    for i, j in enumerate(lateralT2):
        MaxelT2L += 1
        fout2.write("%d, %d, %d, %d\n" % (int(MaxelT2L), int(j[0]), int(j[1]), int(j[2])))
    MaxelT2 = int(MaxelT2L)
    fout2.write("*nset, nset= CNapical\n")
    for x, a in enumerate(cnodeApical):
        if x + 1 == len(cnodeApical):
            if x % 16 == 15:
                fout2.write("\n")
            fout2.write(str(a) + "\n")
        elif x % 16 == 15:
            fout2.write(str(a) + "\n")
        else:
            fout2.write(str(a) + ",")
    fout2.write("*nset, nset= CNbasal\n")
    for x, a in enumerate(cnodeBasal):
        if x + 1 == len(cnodeBasal):
            if x % 16 == 15:
                fout2.write("\n")
            fout2.write(str(a) + "\n")
        elif x % 16 == 15:
            fout2.write(str(a) + "\n")
        else:
            fout2.write(str(a) + ",")
    fout2.write("*nset, nset= CNlateral\n")
    for x, a in enumerate(cnodeLateral):
        if x + 1 == len(cnodeLateral):
            if x % 16 == 15:
                fout2.write("\n")
            fout2.write(str(a) + "\n")
        elif x % 16 == 15:
            fout2.write(str(a) + "\n")
        else:
            fout2.write(str(a) + ",")
    ss = "*elset,  elset=ElsetApical\n"
    fout2.write(ss)
    for x in range(len(apicalEl)):
        if x + 1 == len(apicalEl):
            fout2.write(str(apicalEl[x]) + "\n")
        elif x % 16 == 15:
            fout2.write(str(apicalEl[x]) + "\n")
        else:
            fout2.write(str(apicalEl[x]) + ",")
    ss = "*elset,  elset=ElsetBasal\n"
    fout2.write(ss)
    for x in range(len(basalEl)):
        if x + 1 == len(basalEl):
            fout2.write(str(basalEl[x]) + "\n")
        elif x % 16 == 15:
            fout2.write(str(basalEl[x]) + "\n")
        else:
            fout2.write(str(basalEl[x]) + ",")
    ss = "*elset,  elset=ElsetLateral\n"
    fout2.write(ss)
    for x in range(len(lateralEl)):
        if x + 1 == len(lateralEl):
            fout2.write(str(lateralEl[x]) + "\n")
        elif x % 16 == 15:
            fout2.write(str(lateralEl[x]) + "\n")
        else:
            fout2.write(str(lateralEl[x]) + ",")
    alllEl = []
    alllEl = apicalEl + basalEl + lateralEl
    ss = "*elset,  elset=ElsetAll, generate\n"
    fout2.write(ss)
    ss = str(list(elements.keys())[0]) + ", " + str(list(elements.keys())[-1]) + ", 1\n"
    fout2.write(ss)

    """
    fo.write(ss)
    for x in range(len(alllEl)):
        if x + 1 == len(alllEl):
            fo.write(str(alllEl[x]) + "\n")
        elif x % 16 == 15:
            fo.write(str(alllEl[x]) + "\n")
        else:
            fo.write(str(alllEl[x]) + ",")
    """
    ss = "*elset,  elset=AllEl, generate\n"
    fout2.write(ss)
    ss = str(list(elements.keys())[0]) + ", " + str(list(elements.keys())[-1]) + ", 1\n"
    fout2.write(ss)
    ss = "*elset,  elset=surfApical, generate\n"
    fout2.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fout2.write(ss)
    ss = "*elset,  elset=surfBasal, generate\n"
    fout2.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fout2.write(ss)
    ss = "*elset,  elset=surfLateral, generate\n"
    fout2.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fout2.write(ss)

    ss = "*elset,  elset=CORTEXAPICAL, generate\n"
    fout2.write(ss)
    ss = str(MaxElS6 + 1) + ", " + str(MaxelT2A) + ", 1\n"
    fout2.write(ss)
    ss = "*elset,  elset=CORTEXBASAL, generate\n"
    fout2.write(ss)
    ss = str(MaxelT2A + 1) + ", " + str(MaxelT2B) + ", 1\n"
    fout2.write(ss)
    ss = "*elset,  elset=CORTEXLATERAL, generate\n"
    fout2.write(ss)
    ss = str(MaxelT2B + 1) + ", " + str(MaxelT2L) + ", 1\n"
    fout2.write(ss)

    ss = "*Elset, elset=APICALSPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fout2.write(ss)
    ss = "*Elset, elset=BASALSPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fout2.write(ss)
    ss = "*elset,  elset=LATERALSPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fout2.write(ss)

    ss = "*Elset, elset=_CPAPICAL_SPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElS6 + 1) + ", " + str(MaxelT2A) + ", 1\n"
    fout2.write(ss)
    ss = "*Elset, elset=_CPBASAL_SPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxelT2A + 1) + ", " + str(MaxelT2B) + ", 1\n"
    fout2.write(ss)
    ss = "*Elset, elset=_CPLATERAL_SPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxelT2B + 1) + ", " + str(MaxelT2L) + ", 1\n"
    fout2.write(ss)

    ss = "*Elset, elset=_SPAPICAL_SPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fout2.write(ss)
    ss = "*Elset, elset=_SNAPICAL_SNEG, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fout2.write(ss)


    ss = "*Elset, elset=_SPBASAL_SPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fout2.write(ss)
    ss = "*Elset, elset=_SNBASAL_SNEG, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fout2.write(ss)

    ss = "*Elset, elset=_SPLATERAL_SPOS, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fout2.write(ss)
    ss = "*Elset, elset=_SNLATERAL_SNEG, internal, generate\n"
    fout2.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fout2.write(ss)



    ss = "*Surface, type=ELEMENT, name=CPAPICAL\n"
    fout2.write(ss)
    ss = "_CPAPICAL_SPOS, SPOS\n"
    fout2.write(ss)
    ss = "*Surface, type=ELEMENT, name=CPBASAL\n"
    fout2.write(ss)
    ss = "_CPBASAL_SPOS, SPOS\n"
    fout2.write(ss)
    ss = "*Surface, type=ELEMENT, name=CPLATERAL\n"
    fout2.write(ss)
    ss = "_CPLATERAL_SPOS, SPOS\n"
    fout2.write(ss)

    ss = "*Surface, type=ELEMENT, name=SPAPICAL\n"
    fout2.write(ss)
    ss = "_SPAPICAL_SPOS, SPOS\n"
    fout2.write(ss)
    ss = "*Surface, type=ELEMENT, name=SNAPICAL\n"
    fout2.write(ss)
    ss = "_SNAPICAL_SNEG, SNEG\n"
    fout2.write(ss)    

    ss = "*Surface, type=ELEMENT, name=SPBASAL\n"
    fout2.write(ss)
    ss = "_SPBASAL_SPOS, SPOS\n"
    fout2.write(ss)
    ss = "*Surface, type=ELEMENT, name=SNBASAL\n"
    fout2.write(ss)
    ss = "_SNBASAL_SNEG, SNEG\n"
    fout2.write(ss)

    ss = "*Surface, type=ELEMENT, name=SPLATERAL\n"
    fout2.write(ss)
    ss = "_SPLATERAL_SPOS, SPOS\n"
    fout2.write(ss)
    ss = "*Surface, type=ELEMENT, name=SNLATERAL\n"
    fout2.write(ss)
    ss = "_SNLATERAL_SNEG, SNEG\n"
    fout2.write(ss)


    ss = "*Shell Section, elset=CORTEXAPICAL, material=cortexa\n"
    fout2.write(ss)
    ss = "1., 5\n"
    fout2.write(ss)
    ss = "*Shell Section, elset=CORTEXBASAL, material=cortexb\n"
    fout2.write(ss)
    ss = "1., 5\n"
    fout2.write(ss)
    ss = "*Shell Section, elset=CORTEXLATERAL, material=cortexl\n"
    fout2.write(ss)
    ss = "1., 5\n"
    fout2.write(ss)
    ss = "*Solid Section, elset=ELSETALL, material=cytoplasm1\n"
    fout2.write(ss)
    ss = ",\n"
    fout2.write(ss)
    ss = "*Shell Section, elset=SURFAPICAL, material=membranea\n"
    fout2.write(ss)
    ss = "1., 5\n"
    fout2.write(ss)
    ss = "*Shell Section, elset=SURFBASAL, material=membraneb\n"
    fout2.write(ss)
    ss = "1., 5\n"
    fout2.write(ss)
    ss = "*Shell Section, elset=SURFLATERAL, material=membranel\n"
    fout2.write(ss)
    ss = "1., 5\n"



    fout2.write(ss)
    fout2.write("*End Part\n")
    fout2.write("**\n")
    fout2.close()
    return oldMaxNode, MaxNode, cnodeApical, cnodeBasal, cnodeLateral 


def Convert_element(fin, fout, elemDim, elemPreWord=None, elemPostWord=None):
    print("processing %s -> %s  for conversion with dimension %s" % (fin, fout, elemDim))
    C3D4nodes, MaxNode, cnodeApical, cnodeBasal, cnodeLateral = DecimateINP(fin, fout, elemDim)
    return C3D4nodes, cnodeApical, cnodeBasal, cnodeLateral


def subdivide_edges(subdivisions, points, facets, facet_markers=None):
    def intermediate_points(pa, pb, n):
        for i in range(1, n):
            tau = i / n
            yield [pai * (1 - tau) + tau * pbi for pai, pbi in zip(pa, pb)]

    if isinstance(subdivisions, int):
        from itertools import repeat

        subdiv_it = repeat(subdivisions, len(facets))
    else:
        assert len(facets) == len(subdivisions)
        subdiv_it = subdivisions.__iter__()

    new_points = points[:]
    new_facets = []

    if facet_markers is not None:
        assert len(facets) == len(facet_markers)
        new_facet_markers = []

    for facet_idx, ((pidx_a, pidx_b), subdiv) in enumerate(zip(facets, subdiv_it)):
        facet_points = [pidx_a]
        for p in intermediate_points(points[pidx_a], points[pidx_b], subdiv):
            facet_points.append(len(new_points))
            new_points.append(p)
        facet_points.append(pidx_b)

        for i, p1 in enumerate(facet_points[:-1]):
            p2 = facet_points[i + 1]
            new_facets.append((p1, p2))

            if facet_markers is not None:
                new_facet_markers.append(facet_markers[facet_idx])

    if facet_markers is not None:
        return new_points, new_facets, new_facet_markers
    else:
        return new_points, new_facets

def segment_edge(a, b):
    # list coords point a and list coords point b
    cx = float((a[0] + b[0]) / 2)
    cy = float((a[1] + b[1]) / 2)
    cz = float((a[2] + b[2]) / 2)
    c = [cx, cy, cz]
    return c  # list coords point c

def C3D4_c3d10(elements, coord):
    coords = coord.copy()
    ips = {}
    intermediate_points_coords = {}
    new_elements = []
    for e in elements:
        ei5 = (e[0], e[1])
        if ei5 not in ips.keys():
            e5 = segment_edge(coords[e[0] - 1], coords[e[1] - 1])
            intermediate_points_coords[ei5] = e5
            coords.append(e5)
            node5 = coords.index(e5) + 1
            ips[ei5] = node5
            ips[(e[1], e[0])] = node5
        else:
            node5 = ips[ei5]

        ei6 = (e[1], e[2])
        if ei6 not in ips:
            e6 = segment_edge(coords[e[1] - 1], coords[e[2] - 1])
            intermediate_points_coords[ei6] = e6
            coords.append(e6)
            node6 = coords.index(e6) + 1
            ips[ei6] = node6
            ips[(e[2], e[1])] = node6
        else:
            node6 = ips[ei6]

        ei7 = (e[0], e[2])
        if ei7 not in ips:
            e7 = segment_edge(coords[e[0] - 1], coords[e[2] - 1])
            intermediate_points_coords[(e[0], e[2])] = e7
            coords.append(e7)
            node7 = coords.index(e7) + 1
            ips[ei7] = node7
            ips[(e[2], e[0])] = node7
        else:
            node7 = ips[ei7]

        ei8 = (e[0], e[3])
        if ei8 not in ips:
            e8 = segment_edge(coords[e[0] - 1], coords[e[3] - 1])
            intermediate_points_coords[(e[0], e[3])] = e8
            coords.append(e8)
            node8 = coords.index(e8) + 1
            ips[ei8] = node8
            ips[(e[3], e[0])] = node8
        else:
            node8 = ips[ei8]

        ei9 = (e[1], e[3])
        if ei9 not in ips:
            e9 = segment_edge(coords[e[1] - 1], coords[e[3] - 1])
            intermediate_points_coords[(e[1], e[3])] = e9
            coords.append(e9)
            node9 = coords.index(e9) + 1
            ips[ei9] = node9
            ips[(e[3], e[1])] = node9
        else:
            node9 = ips[ei9]

        ei10 = (e[2], e[3])
        if (ei10) not in ips:
            e10 = segment_edge(coords[e[2] - 1], coords[e[3] - 1])
            intermediate_points_coords[(e[2], e[3])] = e10
            coords.append(e10)
            node10 = coords.index(e10) + 1
            ips[ei10] = node10
            ips[(e[3], e[2])] = node10
        else:
            node10 = ips[ei10]

        elm = [e[0], e[1], e[2], e[3], node5, node6, node7, node8, node9, node10]
        new_elements.append(elm)

    return new_elements, coords

def c3d10_C3D4s(nodes, tetra):
    new_nodes = nodes.copy()
    new_tetra = []

    for tet in tetra:
        if len(tet)==11:
            tet.pop(0)
        new_tetra.append([tet[0], tet[4], tet[6], tet[7]])
        new_tetra.append([tet[4], tet[5], tet[6], tet[7]])
        new_tetra.append([tet[1], tet[5], tet[4], tet[8]])
        new_tetra.append([tet[5], tet[7], tet[4], tet[8]])
        new_tetra.append([tet[3], tet[7], tet[9], tet[8]])
        new_tetra.append([tet[7], tet[5], tet[9], tet[8]])
        new_tetra.append([tet[2], tet[6], tet[5], tet[9]])
        new_tetra.append([tet[6], tet[7], tet[5], tet[9]])
    return new_nodes, new_tetra

def subdivide_tetra(ref_nodes, nodes, tetra):
    new_nodes = nodes.copy()
    new_tetra = tetra.copy()
    new_tetra, new_nodes = C3D4_c3d10(tetra, nodes)
    new_nodes, new_tetra = c3d10_C3D4s(new_nodes, new_tetra)
    return new_nodes, new_tetra

def subdivide_tetra_elem(ref_nodes, nodes, tetra):
    new_nodes = nodes.copy()
    new_tetra = tetra.copy()
    new_tetra, new_nodes = C3D4_c3d10(tetra, nodes)
    new_nodes, new_tetra = c3d10_C3D4s(new_nodes, new_tetra)
    return new_nodes, new_tetra

def sheet_from_cell_centers_julien(points, noise=0):
    """Returns a Sheet object from the Voronoï tessalation
    of the cell centers.
    """
    points2 = points.copy()

    if noise:
        print("NOISE request")
        points += np.random.normal(0, scale=noise, size=points.shape)
    points -= points.mean(axis=0)
    bbox = np.ptp(points, axis=0)
    points /= bbox
    
    rhos = np.linalg.norm(points, axis=1)
    thetas = np.arcsin(points[:, 2] / rhos)
    phis = np.arctan2(points[:, 0], points[:, 1])
    
    sphere_rad = rhos.max() * 1.1
    
    #spherisation
    points_sphere = np.vstack(
        (
            sphere_rad * np.cos(thetas) * np.cos(phis),
            sphere_rad * np.cos(thetas) * np.sin(phis),
            sphere_rad * np.sin(thetas),
        )
    ).T
    
    points_sphere = np.concatenate(([[0, 0, 0]], points_sphere))

    vor3D = Voronoi(points_sphere)
    dsets = from_3d_voronoi(vor3D)
    eptm_ = Epithelium("v", dsets)
    eptm_ = single_cell(eptm_, 0)
    eptm = get_outer_sheet(eptm_)
    eptm.reset_index()
    eptm.reset_topo()
    eptm.vert_df["rho"] = np.linalg.norm(eptm.vert_df[eptm.coords], axis=1)
    mean_rho = eptm.vert_df["rho"].mean()
    eptm.vert_df["basal_shift"] = 0
    SheetGeometry.scale(eptm, sphere_rad / mean_rho, ["x", "y", "z"])
    SheetGeometry.update_all(eptm)
    eptm.vert_df["rho"] = np.linalg.norm(eptm.vert_df[eptm.coords], axis=1)
    
    eptm.face_df["phi"] = np.arctan2(eptm.face_df.y, eptm.face_df.x)
    eptm.face_df["rho"] = np.linalg.norm(eptm.face_df[["x", "y", "z"]], axis=1)
    eptm.face_df["theta"] = np.arcsin(eptm.face_df.z / eptm.face_df["rho"])
    eptm.face_df["x"] = eptm.face_df.eval("rho * cos(theta) * cos(phi)")
    eptm.face_df["y"] = eptm.face_df.eval("rho * cos(theta) * sin(phi)")
    eptm.face_df["z"] = eptm.face_df.eval("rho * sin(theta)")

    eptm.edge_df[["fx", "fy", "fz"]] = eptm.upcast_face(eptm.face_df[["x", "y", "z"]])
    eptm.vert_df[["x", "y", "z"]] = eptm.edge_df.groupby("srce")[
        ["fx", "fy", "fz"]].mean()
    #vip(eptm)
    for i, c in enumerate("xyz"):
        eptm.vert_df[c] *= bbox[i]

    SheetGeometry.update_all(eptm)
    #vip(eptm)
    eptm.sanitize(trim_borders=True)

    eptm.reset_index()
    eptm.reset_topo()
    SheetGeometry.update_all(eptm)
    null_length = eptm.edge_df.query("length == 0")

    while null_length.shape[0]:
        type1_transition(eptm, null_length.index[0])
        SheetGeometry.update_all(eptm)
        null_length = eptm.edge_df.query("length == 0")
        
    
    return eptm


def extrude(apical_datasets, method="homotecy", scale=0.3, vector=[0, 0, -1]):
    """Extrude a sheet to form a monlayer epithelium
    Parameters
    ----------
    * apical_datasets: dictionnary of three DataFrames,
    'vert', 'edge', 'face'
    * method: str, optional {'homotecy'|'translation'|'normals'}
    default 'homotecy'
    * scale: float, optional
    the scale factor for homotetic scaling, default 0.3.
    * vector: sequence of three floats, optional,
    used for the translation
    default [0, 0, -1]
    if `method == 'homotecy'`, the basal layer is scaled down from the
    apical one homoteticaly w/r to the center of the coordinate
    system, by a factor given by `scale`
    if `method == 'translation'`, the basal vertices are translated from
    the apical ones by the vector `vect`
    if `method == 'normals'`, basal vertices are translated from
    the apical ones along the normal of the surface at each vertex,
    by a vector whose size is given by `scale`
    """
    apical_vert = apical_datasets["vert"]
    apical_face = apical_datasets["face"]
    apical_edge = apical_datasets["edge"]

    apical_vert["segment"] = "apical"
    apical_face["segment"] = "apical"
    apical_edge["segment"] = "apical"

    coords = list("xyz")
    datasets = {}

    Nv = apical_vert.index.max() + 1
    Ne = apical_edge.index.max() + 1
    Nf = apical_face.index.max() + 1

    basal_vert = apical_vert.copy()

    basal_vert.index = basal_vert.index + Nv
    basal_vert["segment"] = "basal"

    cell_df = apical_face[coords].copy()
    cell_df.index.name = "cell"
    cell_df["is_alive"] = 1

    basal_face = apical_face.copy()
    basal_face.index = basal_face.index + Nf
    basal_face[coords] = basal_face[coords] * 1 / 3.0
    basal_face["segment"] = "basal"
    basal_face["is_alive"] = 1

    apical_edge["cell"] = apical_edge["face"]
    basal_edge = apical_edge.copy()
    # ## Flip edge so that normals are outward
    basal_edge[["srce", "trgt"]] = basal_edge[["trgt", "srce"]] + Nv
    basal_edge["face"] = basal_edge["face"] + Nf
    basal_edge.index = basal_edge.index + Ne
    basal_edge["segment"] = "basal"

    lateral_face = pd.DataFrame(
        index=apical_edge.index + 2 * Nf, columns=apical_face.columns
    )
    lateral_face["segment"] = "lateral"
    lateral_face["is_alive"] = 1

    lateral_edge = pd.DataFrame(
        index=np.arange(2 * Ne, 6 * Ne), columns=apical_edge.columns
    )

    lateral_edge["cell"] = np.repeat(apical_edge["cell"].values, 4)
    lateral_edge["face"] = np.repeat(lateral_face.index.values, 4)
    lateral_edge["segment"] = "lateral"

    lateral_edge.loc[np.arange(2 * Ne, 6 * Ne, 4), ["srce", "trgt"]] = apical_edge[
        ["trgt", "srce"]
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 1, 6 * Ne, 4), "srce"] = apical_edge[
        "srce"
    ].values
    lateral_edge.loc[np.arange(2 * Ne + 1, 6 * Ne, 4), "trgt"] = basal_edge[
        "trgt"
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 2, 6 * Ne, 4), ["srce", "trgt"]] = basal_edge[
        ["trgt", "srce"]
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 3, 6 * Ne, 4), "srce"] = basal_edge[
        "srce"
    ].values
    lateral_edge.loc[np.arange(2 * Ne + 3, 6 * Ne, 4), "trgt"] = apical_edge[
        "trgt"
    ].values

    if method == "homotecy":
        basal_vert[coords] = basal_vert[coords] * scale
    elif method == "translation":
        for c, u in zip(coords, vector):
            basal_vert[c] = basal_vert[c] + u
    elif method == "normals_JL":
        field = apical_edge.groupby("srce")[["nx", "ny", "nz"]].mean()
        field2 = apical_edge.groupby("srce")[["cell"]].mean()
        field = -field.values * field2.values / np.linalg.norm(field, axis=1)[:, None]
        apical_vert[coords] = apical_vert[coords] + field *0.3
        basal_vert[coords] = basal_vert[coords] + field *0.7
    elif method == "normals":
        field = apical_edge.groupby("srce")[["nx", "ny", "nz"]].mean()
        field = -field.values * scale / np.linalg.norm(field, axis=1)[:, None]
        basal_vert[coords] = basal_vert[coords] + field

    else:
        raise ValueError(
            """
`method` argument not understood, supported values are
'homotecy', 'translation' or 'normals'
        """
        )

    datasets["cell"] = cell_df
    datasets["vert"] = pd.concat([apical_vert, basal_vert])
    datasets["vert"]["is_active"] = 1
    datasets["edge"] = pd.concat([apical_edge, basal_edge, lateral_edge])
    datasets["face"] = pd.concat([apical_face, basal_face, lateral_face])
    datasets["edge"]["is_active"] = 1
    specs = bulk_spec

    for elem in ["vert", "edge", "face", "cell"]:
        datasets[elem].index.name = elem
        for col, value in specs[elem].items():
            if not col in datasets[elem]:
                datasets[elem][col] = value

    if (method == "normals") and (scale < 0):
        datasets["edge"][["srce", "trgt"]] = datasets["edge"][["trgt", "srce"]]
    return datasets


def create_anchors(sheet):
    """Adds an edge linked to every vertices at the boundary
    and create anchor vertices
    """
    anchor_specs = {
        "face": {"at_border": 0},
        "vert": {"at_border": 0, "is_anchor": 0},
        "edge": {"at_border": 0, "is_anchor": 0},
    }

    sheet.update_specs(anchor_specs)
    # ## Edges with no opposites denote the boundary

    free_edge = sheet.edge_df[sheet.edge_df["opposite"] == -1]
    free_vert = sheet.vert_df.loc[free_edge["srce"]]
    free_face = sheet.face_df.loc[free_edge["face"]]

    sheet.edge_df.loc[free_edge.index, "at_border"] = 1
    sheet.vert_df.loc[free_vert.index, "at_border"] = 1
    sheet.face_df.loc[free_face.index, "at_border"] = 1

    # ## Make a copy of the boundary vertices
    anchor_vert_df = free_vert.reset_index(drop=True)
    anchor_vert_df[sheet.coords] = anchor_vert_df[sheet.coords] * 1.01
    anchor_vert_df.index = anchor_vert_df.index + sheet.Nv
    anchor_vert_df["is_anchor"] = 1
    anchor_vert_df["at_border"] = 0
    anchor_vert_df["is_active"] = 0

    sheet.vert_df = pd.concat([sheet.vert_df, anchor_vert_df])
    sheet.vert_df.index.name = "vert"
    anchor_edge_df = pd.DataFrame(
        index=np.arange(sheet.Ne, sheet.Ne + free_vert.shape[0]),
        columns=sheet.edge_df.columns,
    )

    anchor_edge_df["srce"] = free_vert.index
    anchor_edge_df["trgt"] = anchor_vert_df.index
    anchor_edge_df["line_tension"] = 0
    anchor_edge_df["is_anchor"] = 1
    anchor_edge_df["face"] = -1
    anchor_edge_df["at_border"] = 0
    sheet.edge_df = pd.concat([sheet.edge_df, anchor_edge_df], sort=True)
    sheet.edge_df.index.name = "edge"
    sheet.reset_topo()

def write_properties_inp(inp_file,load,center_node,GlobalPartName,allEl,OutsidePart=False):
    print("\n------- Creation of properties and steps for the model -------")
    print("Creating properties in file %s ..." % inp_file)
    file1 = open(inp_file, "a+")
    file1.write("**\n")
    file1.write("** ASSEMBLY\n")
    file1.write("**\n")
    file1.write("*Assembly, name=Assembly\n")
    file1.write("**\n")
    file1.write("*Instance, name=%s-1, part=%s\n" % (GlobalPartName,GlobalPartName))
    file1.write("*End Instance\n")
    file1.write("**\n")
    file1.write("*Instance, name=InsidePart-1, part=InsidePart\n")
    file1.write("*End Instance\n")
    file1.write("**\n")
    if OutsidePart == True:
        file1.write("*Instance, name=OutsidePart-1, part=OutsidePart\n")
        file1.write("*End Instance\n")
        file1.write("**\n")
    file1.write("*Nset, nset=InsidePart_center, instance=InsidePart-1\n")
    file1.write("%d, \n" % int(center_node))
    """
    if method == "perfect":
        file1.write("*Nset, nset=maxload, instance=%s-1\n" % GlobalPartName)
        file1.write("%d, \n" % int(zmax_node))
    """
    ss = "*Surface, type=NODE, name=%s-CNAPICAL, internal\n" % GlobalPartName
    file1.write(ss)
    ss = "%s-1.CNAPICAL, 1.\n" % GlobalPartName
    file1.write(ss)
    ss = "*Surface, type=NODE, name=%s-CNBASAL, internal\n" % GlobalPartName
    file1.write(ss)
    ss = "%s-1.CNBASAL, 1.\n" % GlobalPartName
    file1.write(ss)
    file1.write("*Elset, elset=PRESSURE, internal, generate, instance=%s-1\n" % GlobalPartName)
    file1.write(str(allEl[0])+", "+str(allEl[1])+", "+str(allEl[2])+"\n" )
    file1.write("*Surface, type=ELEMENT, name=AREAPRESSURE\n")
    file1.write("PRESSURE, S1\n")
    file1.write("PRESSURE, S2\n")
    file1.write("PRESSURE, S4\n")
    file1.write("PRESSURE, S3\n")
    file1.write("** Constraint: CP-1-%s-1-InsidePart-1\n" % GlobalPartName)
    file1.write("*Tie, name=CP-1-%s-1-InsidePart-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName)
    file1.write("InsidePart-1.InsidePartContact, %s-1.SPAPICAL\n" % GlobalPartName)
    if OutsidePart == True:
        file1.write("** Constraint: CP-1-%s-1-OutsidePart-1\n" % GlobalPartName)
        file1.write(
            "*Tie, name=CP-1-%s-1-OutsidePart-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
        )
        file1.write("OutsidePart-1.OutsidePartContact, %s-1.SPBASAL\n" % GlobalPartName)
    file1.write("** Constraint: CN-1-%s-CENTERS-APICAL-1\n" % GlobalPartName)
    file1.write(
        "*Tie, name=CN-1-%s-APICAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
    )
    file1.write("%s-1.SPAPICAL, %s-1.CPAPICAL\n" % (GlobalPartName, GlobalPartName))

    file1.write("** Constraint: CN-1-%s-CENTERS-BASAL-1\n" % GlobalPartName)
    file1.write(
        "*Tie, name=CN-1-%s-BASAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
    )
    file1.write(" %s-1.SPBASAL, %s-1.CPBASAL\n" % (GlobalPartName, GlobalPartName))
    
    file1.write("** Constraint: CN-1-%s-CENTERS-LATERAL-1\n" % GlobalPartName)
    file1.write("*Tie, name=CN-1-%s-LATERAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName)
    file1.write(" %s-1.SPLATERAL, %s-1.CPLATERAL\n" % (GlobalPartName, GlobalPartName))
    file1.write("*End Assembly\n")
    file1.write("**\n")
    file1.write("** MATERIALS\n")
    file1.write("**\n")
    
    file1.write("*Material, name=cytoplasm1\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("0.0004, 0.4\n")
    
    file1.write("*Material, name=cytoplasm2\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("0.002, 0.44\n")
    
    file1.write("*Material, name=membranea\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("0.02, 0.35\n")
    
    file1.write("*Material, name=membraneb\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("0.08, 0.38\n")
    
    file1.write("*Material, name=membranel\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("5e-08, 0.35\n")
    
    file1.write("*Material, name=cortexa\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("0.004, 0.4\n")
    file1.write("*Expansion\n")
    file1.write("0.03\n") 
    file1.write("*Specific Heat\n")
    file1.write("200\n") 
    file1.write("*Conductivity\n")
    file1.write("0\n")
    
    file1.write("*Material, name=cortexb\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("0.09, 0.38\n")
    file1.write("*Expansion\n")
    file1.write("1.5e-05\n") 
    file1.write("*Specific Heat\n")
    file1.write("400\n") 
    file1.write("*Conductivity\n")
    file1.write("0\n")
    
    file1.write("*Material, name=cortexl\n")
    file1.write("*Density\n")
    file1.write("7.8e-09,\n")
    file1.write("*Elastic\n")
    file1.write("6e-08, 0.35\n")
    file1.write("*Expansion\n")
    file1.write("6e-08\n") 
    file1.write("*Specific Heat\n")
    file1.write("400\n") 
    file1.write("*Conductivity\n")
    file1.write("0\n")
    
    file1.write("*Material, name=inside\n")
    file1.write("*Density\n")
    file1.write("1e-09,\n")
    file1.write("*Elastic\n")
    file1.write("1e-10, 0.3\n")
    
    file1.write("**\n")    
    file1.write("** BOUNDARY CONDITIONS\n")
    file1.write("** \n")
    file1.write("** Name: Disp-BC-1 Type: Symmetry/Antisymmetry/Encastre\n")
    file1.write("*BOUNDARY\n")
    file1.write("INSIDEPART_CENTER, ENCASTRE\n")
    file1.write("**--------------------------------------------------------------\n")
    file1.write("** \n")
    file1.write("** STEP: Step-1\n")
    file1.write("** \n")
    file1.write("*Step, name=Step-1, nlgeom=NO, inc=1000000\n")
    file1.write("*Coupled Temperature-displacement, creep=none, steady state\n")
    file1.write("0.01, 1.,\n")
    file1.write("**\n")    
    file1.write("** BOUNDARY CONDITIONS\n")
    file1.write("**\n")
    file1.write("** Name: Temp-BC-1 Type: Temperature\n")
    file1.write("*BOUNDARY\n")
    file1.write("AREA1, 11, 11, -8.\n")
    file1.write("** Name: Temp-BC-2 Type: Temperature\n")  
    file1.write("*BOUNDARY\n")
    file1.write("AREA2, 11, 11, -3.\n")
    file1.write("** Name: Temp-BC-3 Type: Temperature\n") 
    file1.write("*BOUNDARY\n")
    file1.write("AREA3, 11, 11, -12.\n")
    file1.write("**\n")
    file1.write("** LOADS\n")
    file1.write("**\n")
    file1.write("** Name: SURFFORCE-1   Type: Pressure\n")
    file1.write("*Dsload\n")
    file1.write("ORGANOID-1.CPAPICAL, P, -0.0005\n")
    file1.write("** \n")
    file1.write("** Name: SURFFORCE-2   Type: Pressure\n")
    file1.write("*Dsload\n")
    file1.write("ORGANOID-1.CPBASAL, P, -0.0002\n")
    file1.write("** \n")
    
    file1.write("** OUTPUT REQUESTS\n")
    file1.write("** \n")
    file1.write("*Restart, write, number interval=1, time marks=NO\n")
    file1.write("** \n")
    file1.write("** FIELD OUTPUT: F-Output-1\n")
    file1.write("** \n")
    file1.write("*Output, field\n")
    file1.write("*Node Output\n")
    file1.write("CF, COORD, MOT, RF, TF, U, UR, UT\n")
    file1.write("*element Output, directions=YES\n") 
    file1.write("ELASE, ELSE, EVOL, IVOL, PRESS, TEMP\n")
    file1.write("** \n")
    file1.write("** HISTORY OUTPUT: H-Output-1\n")
    file1.write("** \n")
    file1.write("***Output, history, frequency=0\n")
    file1.write("***FILE FORMAT, ASCII\n")
    file1.write("***NODE FILE\n")        
    file1.write("**COORD\n")
    file1.write("**U\n")
    file1.write("*End Step\n")
    file1.close()
    print("properties added to inp file %s " % inp_file)

def surface_obj(surface, connectivity, objfile):
    # export surfaces mesh
    fileobj = open(objfile, "w")
    fileobj.write("# OBJ file\n")
    for i in surface:
        fileobj.write("v %f %f %f\n" % (float(i[0][0]), float(i[0][1]), float(i[0][2])))
    fileobj.write("s off\n")
    for c in connectivity:
        connect = []
        for r, n in enumerate(surface):
            if c[0] == n[1]:
                connect.append(r + 1)
        for r, n in enumerate(surface):
            if c[1] == n[1]:
                connect.append(r + 1)
        for r, n in enumerate(surface):
            if c[2] == n[1]:
                connect.append(r + 1)
        fileobj.write(
            "f %d %d %d\n" % (int(connect[0]), int(connect[1]), int(connect[2]))
        )
    fileobj.close()


def readInputFile(inputFileName):
    from copy import deepcopy
    nodes = {}
    elements = {}
    nsets = {}
    elsetCellIds = []
    elsetOther = []
    elsets = {}
    elsetsO = {}
    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")
            if myLine.lower() == "*node":
                while True:
                    iLine += 1
                    try:
                        myLine = lines[iLine]
                    except:
                        break
                    myLine = myLine.strip()
                    myLine = myLine.replace(" ", "")
                    myLine = myLine.replace("\n", "")
                    data = myLine.split(",")
                    try:
                        ndId = int(data[0])
                        nodes.setdefault(ndId, [])
                        nodes[ndId].append(float(data[1]))
                        nodes[ndId].append(float(data[2]))
                        nodes[ndId].append(float(data[3]))
                    except:
                        iLine -= 1
                        lastNdId = deepcopy(ndId)
                        break
            elif myLine[:8].lower() == "*element":
                while True:
                    iLine += 1
                    try:
                        myLine = lines[iLine]
                    except:
                        break
                    myLine = myLine.strip()
                    myLine = myLine.replace(" ", "")
                    myLine = myLine.replace("\n", "")
                    data = myLine.split(",")
                    try:
                        elId = int(data[0])
                        ndsIds = []
                        for ndId in data[1:]:
                            ndsIds.append(int(ndId))
                        elements[elId] = ndsIds
                    except:
                        iLine -= 1
                        break
            elif myLine[:6].lower() == "*elset":
                data = myLine.split(",")
                data = data[1].split("=")
                nsetName = data[1]
                nsetNdsIds = []
                try:
                    cellIds = int("".join(filter(str.isdigit, nsetName)))
                    elsetCellIds.append(int(cellIds))
                    elsets.setdefault(cellIds, [])
                    while True:
                        iLine += 1
                        try:
                            myLine = lines[iLine]
                        except:
                            break
                        myLine = myLine.strip()
                        myLine = myLine.replace(" ", "")
                        myLine = myLine.replace("\n", "")
                        data = myLine.split(",")
                        try:
                            for ndId in data:
                                elsets[cellIds].append(int(ndId))
                        except:
                            iLine -= 1
                            break
                except:
                    elsetOther.append(nsetName)
                    elsetsO.setdefault(nsetName, [])
                    while True:
                        iLine += 1
                        try:
                            myLine = lines[iLine]
                        except:
                            break
                        myLine = myLine.strip()
                        myLine = myLine.replace(" ", "")
                        myLine = myLine.replace("\n", "")
                        data = myLine.split(",")
                        try:
                            for ndId in data:
                                elsetsO[nsetName].append(int(ndId))
                        except:
                            iLine -= 1
                            break

            elif myLine[:5].lower() == "*nset":
                data = myLine.split(",")
                data = data[1].split("=")
                nsetName = data[1]
                nsetNdsIds = []
                while True:
                    iLine += 1
                    try:
                        myLine = lines[iLine]
                    except:
                        break
                    myLine = myLine.strip()
                    myLine = myLine.replace(" ", "")
                    myLine = myLine.replace("\n", "")
                    data = myLine.split(",")
                    try:
                        for ndId in data:
                            nsetNdsIds.append(int(ndId))
                    except:
                        iLine -= 1
                        break

                nsets[nsetName] = nsetNdsIds

            elif myLine[:5].lower() == "*part":
                data = myLine.split(",")
                data = data[1].split("=")
                partName = data[1]

            iLine += 1
    ndim = len(nodes[lastNdId])
    return (
        ndim,
        lastNdId,
        partName,
        nsetNdsIds,
        elements,
        nodes,
        elsetCellIds,
        elsets,
        elsetsO,
        nsets,
    )


def export_obj_surfaces(infile, outfile):
    # load Organoid informations from inp
    apicalOut = outfile.replace(".obj", "_apical.obj")  # generated during export
    basalOut = outfile.replace(".obj", "_basal.obj")  # generated during export
    lateralOut = outfile.replace(".obj", "_lateral.obj")  # generated during export

    surfaceA = []
    surfaceB = []
    surfaceL = []

    surfaceAconnect = []
    surfaceBconnect = []
    surfaceLconnect = []
    (
        ndim,
        lastNdId,
        partName,
        nsetNdsIds,
        elements2,
        nodes,
        elsetCellIds,
        elsets,
        elsetsO,
        nsets,
    ) = readInputFile(infile)

    elements = {}
    for key, value in elements2.items():
        if len(value) > 4:
            elements[key] = value

    for i in range(elsetsO["surfApical"][0], elsetsO["surfApical"][1] + 1, 1):
        if any(elements[i][0] in sublist for sublist in surfaceA) == False:
            surfaceA.append([nodes[elements[i][0]], elements[i][0]])
        if any(elements[i][1] in sublist for sublist in surfaceA) == False:
            surfaceA.append([nodes[elements[i][1]], elements[i][1]])
        if any(elements[i][2] in sublist for sublist in surfaceA) == False:
            surfaceA.append([nodes[elements[i][2]], elements[i][2]])
        surfaceAconnect.append([elements[i][0], elements[i][1], elements[i][2]])
    for i in range(elsetsO["surfBasal"][0], elsetsO["surfBasal"][1] + 1, 1):
        if any(elements[i][0] in sublist for sublist in surfaceB) == False:
            surfaceB.append([nodes[elements[i][0]], elements[i][0]])
        if any(elements[i][1] in sublist for sublist in surfaceB) == False:
            surfaceB.append([nodes[elements[i][1]], elements[i][1]])
        if any(elements[i][2] in sublist for sublist in surfaceB) == False:
            surfaceB.append([nodes[elements[i][2]], elements[i][2]])
        surfaceBconnect.append([elements[i][0], elements[i][1], elements[i][2]])
    for i in range(elsetsO["surfLateral"][0], elsetsO["surfLateral"][1] + 1, 1):
        if any(elements[i][0] in sublist for sublist in surfaceL) == False:
            surfaceL.append([nodes[elements[i][0]], elements[i][0]])
        if any(elements[i][1] in sublist for sublist in surfaceL) == False:
            surfaceL.append([nodes[elements[i][1]], elements[i][1]])
        if any(elements[i][2] in sublist for sublist in surfaceL) == False:
            surfaceL.append([nodes[elements[i][2]], elements[i][2]])
        surfaceLconnect.append([elements[i][0], elements[i][1], elements[i][2]])

    surface_obj(surfaceA, surfaceAconnect, apicalOut)
    surface_obj(surfaceB, surfaceBconnect, basalOut)
    surface_obj(surfaceL, surfaceLconnect, lateralOut)

    return (
        surfaceA,
        surfaceAconnect,
        surfaceB,
        surfaceBconnect,
        surfaceL,
        surfaceLconnect,
        apicalOut,
        basalOut,
        lateralOut,
    )

def write_inside_part_inp(inp_file, nodes, tetra, external_nodes):
    print("\n------- Creation of Lumen part into the model -------")
    print("Creating internal part in file %s ..." % inp_file)
    
    file2 = open(inp_file, "a+")

    file2.write("*PART,  name=InsidePart\n")
    file2.write("*NODE, nset=InsidePart\n")
    internal_volume = []
    for i, n in enumerate(nodes):
        file2.write(
            "%d, %f, %f, %f\n" % (int(i + 1), float(n[0]), float(n[1]), float(n[2]))
        )
    if len(tetra[1]) == 4:
        file2.write("*element, TYPE=C3D4\n")
        for k, x in enumerate(tetra):
            file2.write(
                "%d, %d, %d, %d, %d\n"
                % (int(k + 1), int(x[0]), int(x[1]), int(x[2]), int(x[3]))
            )
            try:
                extt = [nodes[x[0] + 1], nodes[x[1] + 1], nodes[x[2] + 1]]
                for r, j in enumerate(extt):
                    if j in external_nodes:
                        internal_volume.append(int(k + 1))
            except:
                pass
    elif len(tetra[1]) == 10:
        file2.write("*element, TYPE=C3D10\n")
        for k, x in enumerate(tetra):
            file2.write(
                "%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n"
                % (
                    int(k + 1),
                    int(x[0]),
                    int(x[1]),
                    int(x[2]),
                    int(x[3]),
                    int(x[4]),
                    int(x[5]),
                    int(x[6]),
                    int(x[7]),
                    int(x[8]),
                    int(x[9]),
                )
            )
    
    file2.write("*Elset, elset=all_InsidePart, generate\n")
    file2.write("%d, %d, %d\n" % (int(1), int(k + 1), int(1)))
    
    # ajouter les surfaces -->
    file2.write("*Elset, elset=_InsidePartContact_S1, internal\n")
    for x in range(len(internal_volume)):
        if x + 1 == len(internal_volume):
            file2.write(str(internal_volume[x]) + "\n")
        elif x % 16 == 15:
            file2.write(str(internal_volume[x]) + "\n")
        else:
            file2.write(str(internal_volume[x]) + ",")
    file2.write("*Surface, type=ELEMENT, name=InsidePartContact\n")
    file2.write("_InsidePartContact_S1, S1\n")
    file2.write("** Section: internal\n")
    file2.write("*Solid Section, elset=all_InsidePart, material=inside\n")
    file2.write(",\n")
    file2.write("*End Part\n")
    file2.write("**\n")
    file2.close()
    print("Checking informations in file -> %s " % inp_file)
    nsets=read_nsets(inp_file)
    apicalL = []
    basalL = []
    lateralL = []
    for i,j in nsets.items():
        if "apical" in i:
            apicalL=j
        if "basal" in i:
            basalL=j
        if "lateral" in i:
            lateralL=j
    print("----------------------------------")
    print("# of apical Centroids -> %s " % len(apicalL))
    print("----------------------------------")
    print("# of basal Centroids -> %s " % len(basalL))
    print("----------------------------------")
    print("# of lateral Centroids -> %s " % len(lateralL))  
    print("internal part added to inp file %s " % inp_file)
    
def punch(input, center, radius):
    try:
        nodes = readNodesINP(inpfile)
        #nodes = readNodesINP(input)
    except:
        nodes = input
    punchnodes = {}
    punchmode = {}
    p = 0.1
    h = 0.6
    before = nodes.copy()
    x = []
    y = []
    z = []
    xx = []
    zz = []
    zzz = []
    zmin=200
    idzmin=0
    for i, j in nodes.items():
        if j[2] > center[2]:
            x.append(j[0])
            xx.append(j[0] + radius * 3)
            y.append(j[1])
            z.append(j[2])
            ppush = (1 + p * (j[0] ** 2 + j[1] ** 2)) / (radius + p + h ** 2)
            punchme = 1 / ppush
            zcoord = j[2] - punchme
            punchnodes[i] = [j[0], j[1], zcoord]
            punchmode[i] = punchme * 2
            zz.append(zcoord)
            zzz.append(punchme)
        else:
            if j[2]< zmin:
                zmin=j[2]
                idzmin= int(i)
            punchnodes[i] = [j[0], j[1], j[2]]

    opposite = idzmin
    variances = zzz
    vmin, vmax = np.min(variances), np.max(variances)
    return punchnodes, punchmode, opposite

def write_header_inp(fin, fout, GlobalPartName):
    print("\n------- Adding new header into the model -------")
    print("Creating header in file %s ..." % fout)
    lineCount = 0
    modelname = fin.replace(".inp", "")

    # ------------File----------------------------------------

    fii = open(fin, "r")
    foo = open(fout, "w")
    foo.write("*Heading\n")
    foo.write("** Job name: %s\n" % str(modelname))
    foo.write("** Generated by: Abaqus/CAE 2018\n")
    foo.write("** Copyright by: Julien Laussu\n")
    foo.write("**\n")
    foo.write("** PARTS\n")
    foo.write("**\n")

    while True:
        s, lineCount = ReadInp(fii, lineCount)
        s1 = s.strip()  # chomp
        if s1 == "*part,  name=mono":
            print("replace name by:  %s in %s file" % (GlobalPartName, fout))
            foo.write("*part,  name=%s\n" % GlobalPartName)
        #elif s1 == "*node":
        #    fo.write("*node, nset=nall\n")
        else:
            if s1 == "*End Step":
                foo.write(s)
                break
            foo.write(s)
    
    fii.close()
    foo.close()

def create_InsidePart(insurfacePTS, insurfaceCON, outfile):
    # organoid barycenter
    inner_points = [[row[i] for row in insurfacePTS] for i in range(2)]
    xcentroid = np.mean([row[0] for row in inner_points[0]])
    ycentroid = np.mean([row[1] for row in inner_points[0]])
    zcentroid = np.mean([row[2] for row in inner_points[0]])
    InternalPart_center = [xcentroid, ycentroid, zcentroid]
    coord_in = inner_points[0].copy()
    labels_in = inner_points[1].copy()
    connect_in = insurfaceCON.copy()
    center_label = len(labels_in) + 1
    # create first elements tetra
    tetra_in = []
    for c in connect_in:
        triang = []
        for cc in c:
            triang.append(labels_in.index(cc) + 1)
        triang.append(center_label)
        tetra_in.append(triang)
    ref_nodes = [*range(1, len(inner_points[1]) + 1, 1)]

    nodes_in = coord_in.copy()
    nodes_in.append(InternalPart_center)
    new_nodes, new_tetra = subdivide_tetra(ref_nodes, nodes_in, tetra_in)
    write_inside_part_inp(outfile, new_nodes, new_tetra, nodes_in)
    center_node = new_nodes.index(InternalPart_center) + 1
    return center_node, InternalPart_center, new_nodes, new_tetra 


def angle_triangle(t1,t2,t3):
    [x1, y1, z1], [x2, y2, z2], [x3, y3, z3] = t1,t2,t3
    num = (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+(z2-z1)*(z3-z1)
    den = math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)*math.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
    angle = math.degrees(math.acos(num / den))
    return round(angle, 3)

def triangle_area(face):
    a = distance3D(face[0],face[1])
    b = distance3D(face[1],face[2])
    c = distance3D(face[2],face[0])
    s = (a + b + c) / 2
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    angle_A = angle_triangle(face[0],face[1],face[2])
    angle_B = angle_triangle(face[1],face[2],face[0])
    angle_C = angle_triangle(face[2],face[0],face[1])
    minangle=min([angle_A, angle_B, angle_C])
    minedge=min([a,b,c])
    return area, minedge, minangle

def tet_integrity(tet):
    areas=[]
    min_edges=[]
    min_angle=[]

    for face in (  (tet[0], tet[1], tet[2]), 
                   (tet[0], tet[2], tet[3]), 
                   (tet[0], tet[3], tet[2]),
                   (tet[1], tet[3], tet[2]) ):
        area, minedge, minangle = triangle_area(face)
        areas.append(area)
        min_edges.append(minedge)
        min_angle.append(minangle)

    ratio=min(areas)/max(areas)
    if min(min_edges)<0.1:
        ratio=0.001
    if min(min_angle)<5:
        ratio=0.001
    return ratio

def test_plane(tet):
    [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3], [x, y, z]]=tet
    Coplanar=False
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * x1 - b * y1 - c * z1)
    if(a * x + b * y + c * z + d < 100):
        #print("coplanar")
        Coplanar=True
    return Coplanar

def min_edge_face(elem, tet):
    #print(tet)
    edge_lenght=[]
    edges=[[elem[0],elem[1]], [elem[1],elem[2]], [elem[2],elem[0]], 
                   [elem[0],elem[2]],[elem[2],elem[3]], [elem[3],elem[0]], 
                   [elem[0],elem[3]], [elem[3],elem[1]], [elem[1],elem[0]],
                   [elem[1],elem[3]], [elem[3],elem[2]], [elem[2],elem[1]]]
    for face in (  (tet[0], tet[1], tet[2]), 
                   (tet[0], tet[2], tet[3]), 
                   (tet[0], tet[3], tet[1]),
                   (tet[1], tet[3], tet[2]) ):    
        a = distance3D(face[0],face[1])
        b = distance3D(face[1],face[2])
        c = distance3D(face[2],face[0])
        edge_lenght.append(a)
        edge_lenght.append(b)
        edge_lenght.append(c)
    min_min = min(edge_lenght)
    edge_idx = edge_lenght.index(min_min)
    return edges[edge_idx]
    
def identify_poorElements(inp_file, organo):
    geom.update_all(organo)

    scores={}
    planar=[]
    nodes = readNodesINP(inp_file)
    elements = readElementsINP(inp_file)
    for elId, elems in elements.items():
        tet=[nodes[elems[0]],nodes[elems[1]],nodes[elems[2]],nodes[elems[3]]]
        score=tet_integrity(tet)
        scores[elId]=score
        planarity=test_plane(tet)
        if planarity==True:
            planar.append(elId)

    elongated=[]
    for el, sc in scores.items():
        if sc < 0.05:
            elongated.append(el)
            
    elongated = [x for x in elongated if x not in planar]
    toremove= elongated + planar
    #toremove= []
    

    #save_inp(inp_file, organo,elemtosup=toremove, add_spaces=True)
    
    return elongated, planar

from tyssue.topology import *


### END OF DEFINITIONS ####################################################



####################
###### START #######
####################

# Creation of cells centroids in normalized space
centers = spheric3D(Nc, 3, method)
#print(centers)
if from_imaging == True:
    centers = read_positions_csv("positions.csv")
    
# organoid barycenter
xmean = np.mean(centers[:, 0])
ymean = np.mean(centers[:, 1])
zmean = np.mean(centers[:, 2])
organo_center = np.array([[xmean, ymean, zmean]])
# transpose centroid array
centers = centers.T
# measure mean distance of centroids to the organoid barycenter
dist_radius = dist.cdist(centers, organo_center)
#radius = np.mean(dist_radius[:, 0])  # pixels
#radius = Nc / 10 * radius
#radius = 1350

thickness = radius * thick  # pixels


## ---- STEP 1 ---- ##
## General properties of the Organoid
R_in = int(radius) - thickness / 2  # internal diameter
R_out = int(radius) + thickness / 2  # external diameter

## Modelisation of a standard Organoid defined by
sheet = spherical_sheet_julien(R_in, centers)

delta_R = R_out - R_in
#print(delta_R)


mono = Monolayer("ORGANOID", extrude(sheet.datasets, method="normals", scale=-delta_R))


#vip(mono)
#sys.exit()
if apical == "out":
    swap_apico_basal(mono)
else:
    mono.settings["lumen_side"] = "apical"
geom.update_all(mono)
organo = mono
print("\n##############################")
print("Orgonoid in numbers: ")
print(str(organo.Nc) + " cells")
print(str(organo.Nf) + " faces")
print(str(organo.Nv) + " vertices")
print(str(organo.Ne) + " edges")
print("##############################\n")
## Saving informations of the Organoid before deformation
# save geometrical informations for abaqus

if IndividualCells == True:
    save_inp_individual(new_inp, organo, add_spaces=True)
else:
    save_inp(new_inp, organo, add_spaces=True)

# creation of area sets for abaqus
inputeditFileName = new_inp  # input inp filename
outputeditFileName = inputeditFileName.replace(
    ".inp", "_area.inp")  # generated during export
# optional, to visualized organo panda dataframes in csv
eptm_to_csv(organo, "csv_export/")
# Creation of an archive file for panda organo informations in .hdf5 / to extend for a complete archive of data and metadata
geom.update_all(organo)


### Sheet for visualisation
specs = config.geometry.spherical_sheet()
sheet = Sheet("emin", organo.datasets, specs)
# saving in Object format ".obj"
obj.save_triangulated(obj_file1, sheet)
## extracting vertices and faces from Mesh model of the organoid
vertices, faces = sheet.triangular_mesh(sheet.coords, False)

# VISUALISATION of Organoid before deformations using Vispy
#plt.show()

geom.update_all(organo)

import pandas as pd
vols=organo.datasets["cell"]["vol"].values
Stats=pd.DataFrame(vols*vsize*vsize*vsize,columns=["volume"])

#calcul des surfaces et barycentres from the mesh

from sklearn.neighbors import KDTree
bary=organo.datasets["cell"].values[:,[0,1,2]]
tree= KDTree(bary)
nearst_dist, nearest_ind = tree.query(bary, k=5)
nearst_dist = nearst_dist

apical_surf={}
basal_surf={}
apical_lateral_surf={}
basal_lateral_surf={}
neibhoors_number={}
for cell in organo.datasets["cell"].index:
    cell_edges = organo.edge_df[organo.edge_df["cell"] == cell]
    apical_edges = cell_edges[cell_edges["segment"] == "apical"]
    basal_edges = cell_edges[cell_edges["segment"] == "basal"]
    lateral_edges = cell_edges[cell_edges["segment"] == "lateral"]
    
    srce_segment = organo.upcast_srce(organo.vert_df["segment"]).loc[
        lateral_edges.index
    ]
    trgt_segment = organo.upcast_trgt(organo.vert_df["segment"]).loc[
        lateral_edges.index
    ]
    srce_segment1 = organo.upcast_srce(organo.vert_df["segment"]).loc[
        apical_edges.index
    ]
    trgt_segment1 = organo.upcast_trgt(organo.vert_df["segment"]).loc[
        apical_edges.index
    ]
    srce_segment2 = organo.upcast_srce(organo.vert_df["segment"]).loc[
        basal_edges.index
    ]
    trgt_segment2 = organo.upcast_trgt(organo.vert_df["segment"]).loc[
        basal_edges.index
    ]
    
    ab_edges = lateral_edges[
        (srce_segment == "apical") & (trgt_segment == "basal")
    ]["sub_area"]
    ba_edges = lateral_edges[
        (trgt_segment == "apical") & (srce_segment == "basal")
    ]["sub_area"]
    cd_edges = apical_edges["sub_area"]
    ef_edges = basal_edges["sub_area"]
    #cell_edges_dic.setdefault(cell, [])

    apical_surface=np.sum(cd_edges.values)
    basall_surface=np.sum(ef_edges.values)
    lateral_surface=np.sum(ab_edges.values)

    apical_surf[cell]=apical_surface
    basal_surf[cell]=basall_surface
    apical_lateral_surf[cell]=apical_surface/lateral_surface
    basal_lateral_surf[cell]=basall_surface/lateral_surface
    neibhoors_number[cell]=len(ab_edges)

Stats["spacing"]=list(nearst_dist[:,2]*vsize)
Stats["apical"]=apical_surf.values()
Stats["apical"]=Stats["apical"]*vsize*vsize
Stats["basal"]=basal_surf.values()
Stats["basal"]=Stats["basal"]*vsize*vsize
Stats["apical-lateral"]=apical_lateral_surf.values()
Stats["basal-lateral"]=basal_lateral_surf.values()
Stats["neibhoors_number"]=neibhoors_number.values()


new_csv=obj_file1.replace(".obj","_stats_pre.csv")
new_csv2=obj_file1.replace(".obj","_stats_pre_2.csv")
Stats.to_csv(new_csv)
Stats = Stats[Stats.apical > 1]
Stats = Stats[Stats.basal > 1]
Stats["basal_apical"]=Stats.basal/Stats.apical
Stats.to_csv(new_csv2)
print("Before deformation stats saved in: "+str(new_csv2))
save_inp(new_inp, organo, add_spaces=True)



## ---- STEP 2 ---- ####
# optional: to compute organoid morphochanges in a model of equilibrium
## Initial model settings
dyn_specs = {
    "settings": {
        "lumen_prefered_vol": organo.settings["lumen_vol"],#*(1/Nc),
        "lumen_vol_elasticity": 10e-1, # organo.Nc,  # 1e-1
        "threshold_length": 1e-3,
    },
    "cell": {
        "prefered_vol": organo.cell_df.vol.mean(),
        "vol_elasticity": 10000, #10000
        "area_elasticity": 1000.0,
        "prefered_area": organo.cell_df.area.mean(),
    },
    "face": {
        "surface_tension": 10,
        "contractility": 0.1,
    },
    "edge": {"line_tension": 50},
}

if SolverStep == True:
    model = model_factory(
        [
            effectors.LineTension,
            effectors.LumenVolumeElasticity,
            effectors.CellAreaElasticity,
            effectors.CellVolumeElasticity,
            effectors.SurfaceTension,
        ]
    )
    ## Equilibrium model to deform Organoid
    solver = QSSolver()
    organo.update_specs(dyn_specs, reset=True)

    res = solver.find_energy_min(organo, geom, model)
    if res["success"] == True:
        print("Min Energy Solver step successfully completed")

    organo.validate()
    organo.reset_index()
    organo.reset_topo()

    #vip(organo) # vispy view after deformation

    ## Saves the informations of the Organoid after deformation
    inputeditFileName = new_inp2  # input inp filename
   
    ## update properties of the Organoid to be sure that all informations are complete
    geom.update_all(organo)
   
    import pandas as pd
    vols=organo.datasets["cell"]["vol"].values
    Stats=pd.DataFrame(vols*vsize*vsize*vsize,columns=["volume"])
    
    #calcul des surfaces et barycentres from the mesh
    
    from sklearn.neighbors import KDTree
    bary=organo.datasets["cell"].values[:,[0,1,2]]
    tree= KDTree(bary)
    nearst_dist, nearest_ind = tree.query(bary, k=5)
    nearst_dist = nearst_dist
    
    apical_surf={}
    basal_surf={}
    apical_lateral_surf={}
    basal_lateral_surf={}
    neibhoors_number={}
    for cell in organo.datasets["cell"].index:
        cell_edges = organo.edge_df[organo.edge_df["cell"] == cell]
        apical_edges = cell_edges[cell_edges["segment"] == "apical"]
        basal_edges = cell_edges[cell_edges["segment"] == "basal"]
        lateral_edges = cell_edges[cell_edges["segment"] == "lateral"]
        
        srce_segment = organo.upcast_srce(organo.vert_df["segment"]).loc[
            lateral_edges.index
        ]
        trgt_segment = organo.upcast_trgt(organo.vert_df["segment"]).loc[
            lateral_edges.index
        ]
        srce_segment1 = organo.upcast_srce(organo.vert_df["segment"]).loc[
            apical_edges.index
        ]
        trgt_segment1 = organo.upcast_trgt(organo.vert_df["segment"]).loc[
            apical_edges.index
        ]
        srce_segment2 = organo.upcast_srce(organo.vert_df["segment"]).loc[
            basal_edges.index
        ]
        trgt_segment2 = organo.upcast_trgt(organo.vert_df["segment"]).loc[
            basal_edges.index
        ]
        
        ab_edges = lateral_edges[
            (srce_segment == "apical") & (trgt_segment == "basal")
        ]["sub_area"]
        ba_edges = lateral_edges[
            (trgt_segment == "apical") & (srce_segment == "basal")
        ]["sub_area"]
        cd_edges = apical_edges["sub_area"]
        ef_edges = basal_edges["sub_area"]
        #cell_edges_dic.setdefault(cell, [])
    
        apical_surface=np.sum(cd_edges.values)
        basall_surface=np.sum(ef_edges.values)
        lateral_surface=np.sum(ab_edges.values)
    
        apical_surf[cell]=apical_surface
        basal_surf[cell]=basall_surface
        apical_lateral_surf[cell]=apical_surface/lateral_surface
        basal_lateral_surf[cell]=basall_surface/lateral_surface
        neibhoors_number[cell]=len(ab_edges)

    Stats["spacing"]=list(nearst_dist[:,2]*vsize)
    Stats["apical"]=apical_surf.values()
    Stats["apical"]=Stats["apical"]*vsize*vsize
    Stats["basal"]=basal_surf.values()
    Stats["basal"]=Stats["basal"]*vsize*vsize
    Stats["apical-lateral"]=apical_lateral_surf.values()
    Stats["basal-lateral"]=basal_lateral_surf.values()
    Stats["neibhoors_number"]=neibhoors_number.values()
    
    new_csv=obj_file1.replace(".obj","_stats_post.csv")
    new_csv2=obj_file1.replace(".obj","_stats_post_2.csv")
    
    Stats.to_csv(new_csv)
    Stats = Stats[Stats.apical > 1]
    Stats = Stats[Stats.basal > 1]
    Stats["basal_apical"]=Stats.basal/Stats.apical
    Stats.to_csv(new_csv2)
    print("After deformation stats saved in: "+str(new_csv2)+"\n")
    save_inp(new_inp2, organo, add_spaces=True)
    dir_path = os.path.dirname(os.path.realpath(new_inp2))
    print("All inp files are saved in: "+str(dir_path)+"\n")

