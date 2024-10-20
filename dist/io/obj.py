import numpy as np
import pandas as pd

try:
    from vispy.io import write_mesh
except ImportError:
    print("You need vispy to use the .OBJ export")

import logging
from .inp import readInputFile

logger = logging.getLogger(name=__name__)


def save_triangulated(filename, eptm):

    vertices, faces = eptm.triangular_mesh(eptm.coords, False)
    write_mesh(
        filename,
        vertices=vertices,
        faces=faces,
        normals=None,
        texcoords=None,
        overwrite=True,
    )
    logger.info("Saved %s as a trianglulated .OBJ file", eptm.identifier)


def save_junction_mesh(filename, eptm):

    vertices, faces, normals = eptm.vertex_mesh(eptm.coords, vertex_normals=True)

    write_mesh(
        filename,
        vertices=vertices,
        faces=faces,
        normals=normals,
        texcoords=None,
        overwrite=True,
        reshape_faces=False,
    )  # GH 1155
    logger.info("Saved %s as a junction mesh .OBJ file", eptm.identifier)


def write_splitted_cells(*args, **kwargs):
    logger.warning("Deprecated, use `save_splitted_cells` instead")
    save_splitted_cells(*args, **kwargs)


def save_splitted_cells(fname, sheet, epsilon=0.1):

    coords = sheet.coords
    up_srce = sheet.upcast_srce(sheet.vert_df[coords])
    up_trgt = sheet.upcast_trgt(sheet.vert_df[coords])
    up_face = sheet.upcast_face(sheet.face_df[coords])
    up_srce = (up_srce - up_face) * (1 - epsilon) + up_face
    up_trgt = (up_trgt - up_face) * (1 - epsilon) + up_face

    cell_faces = pd.concat([sheet.face_df[coords], up_srce, up_trgt], ignore_index=True)
    Ne, Nf = sheet.Ne, sheet.Nf

    triangles = np.vstack(
        [sheet.edge_df["face"], np.arange(Ne) + Nf, np.arange(Ne) + Ne + Nf]
    ).T
    write_mesh(
        fname,
        cell_faces.values,
        triangles,
        normals=None,
        texcoords=None,
        overwrite=True,
    )


def surface_obj(surface, connectivity, objfile):
    with open(objfile, "w") as f:
        f.write("# OBJ file\n")
        for i in surface:
            f.write("v %f %f %f\n" % (float(i[0][0]), float(i[0][1]), float(i[0][2])))
        f.write("s off\n")
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
            f.write(
                "f %d %d %d\n" % (int(connect[0]), int(connect[1]), int(connect[2]))
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
