
# define the files requiered:
nodefilfile = "node_fil_organofem.fil"
elfilfile = "elem_fil_organofem.fil"
original_file="node_inp_organofem.inp"
refined_file="node_inp_organofem.inp"
beforedef_file="node_inp_organofem.inp"
hdf5save= original_file.replace(".inp",".h5")


def readNodesINP(inputFileName):
    precision = 6
    from copy import deepcopy
    nodes = {}
    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")

            myLine_1 = lines[iLine-1]
            myLine_1 = myLine_1.strip()
            myLine_1 = myLine_1.replace(" ", "")
            myLine_1 = myLine_1.replace("\n", "")

            if myLine[:5].lower() == "*node" and "ORGANOID" in myLine_1:
                myLine_1
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
                        nodes[ndId].append("{:.{}f}".format( float(data[1]), precision ))
                        nodes[ndId].append("{:.{}f}".format( float(data[2]), precision ))
                        nodes[ndId].append("{:.{}f}".format( float(data[3]), precision ))
                    except:
                        iLine -= 1
                        lastNdId = deepcopy(ndId)
                        break
            iLine += 1
    return nodes


def readElementsINP(inputFileName):
    precision = 6
    from copy import deepcopy
    elements = {}
    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")

            if myLine[:8].lower() == "*element" and "C3D4" not in myLine:
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
            iLine += 1
    return elements


def readNsetsINP(inputFileName):
    precision = 6
    from copy import deepcopy
    nsets = {}
    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")

            if myLine[:5].lower() == "*nset":
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

            iLine += 1
    return nsets


def readElsetsINP(inputFileName):
    precision = 6
    from copy import deepcopy
    elsets = {}
    elsetsO = {}
    elsetOther=[]

    with open(inputFileName, "r") as myFile:
        lines = myFile.readlines()
        nLines = len(lines)
        iLine = 0
        while iLine < nLines:
            myLine = lines[iLine]
            myLine = myLine.strip()
            myLine = myLine.replace(" ", "")
            myLine = myLine.replace("\n", "")

            if myLine[:6].lower() == "*elset":
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
                    if "_S" not in nsetName:
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

            iLine += 1
    return elsets,elsetsO




######## TRIANGLE AREA ############

def heron(a,b,c):
    s = (a + b + c) / 2
    area = (s*(s-a) * (s-b)*(s-c)) ** 0.5
    return area

def distance3d(x1,y1,z1,x2,y2,z2):
    a=(float(x1)-float(x2))**2+(float(y1)-float(y2))**2 + (float(z1)-float(z2))**2
    d= a ** 0.5
    return d

def areatriangle3d(pointlist):
    point1, point2, point3=pointlist[0],pointlist[1],pointlist[2]
    x1, y1, z1=point1[0],point1[1],point1[2]
    x2, y2, z2=point2[0],point2[1],point2[2]
    x3, y3, z3=point3[0],point3[1],point3[2]
    a=distance3d(x1,y1,z1,x2,y2,z2)
    b=distance3d(x2,y2,z2,x3,y3,z3)
    c=distance3d(x3,y3,z3,x1,y1,z1)
    A = heron(a,b,c)
    #print("area of triangle is %r " %A)
    return A

def areatriangleelement(element,elements,nodes):
    node1, node2, node3 = elements[element][0], elements[element][1],elements[element][2]
    B=areatriangle3d([nodes[node1],nodes[node2],nodes[node3]])
    return B

######## TETRAHEDRON VOLUME ############

def determinant_3x3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
            m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1]) +
            m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1]))


def subtract(a, b):
    return (float(a[0]) - float(b[0]),
            float(a[1]) - float(b[1]),
            float(a[2]) - float(b[2]))

def tetrahedron_calc_volume(a, b, c, d):
    A= abs(determinant_3x3((subtract(a, b),
                                 subtract(b, c),
                                 subtract(c, d),
                                 ))) / 6.0
    #print("volume of tetrahedron is %r " %A)
    return(A)

def volumetetrahedronelement(element,elements,nodes):
    node1, node2, node3, node4 = elements[element][0], elements[element][1],elements[element][2],elements[element][3]
    B=tetrahedron_calc_volume(nodes[node1],nodes[node2],nodes[node3],nodes[node4])
    return B


def volumetetrahedronelementD10(element,elements,nodes):
    node1, node2, node3, node4, node5, node6, node7, node8, node9, node10 = elements[element][0], elements[element][1],elements[element][2],elements[element][3], elements[element][4],elements[element][5],elements[element][6], elements[element][7],elements[element][8],elements[element][9]
    B=tetrahedron_calc_volume(nodes[node1],nodes[node2],nodes[node3],nodes[node4])
    return B

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
                if c[1] == n[1]:
                    connect.append(r + 1)
                if c[2] == n[1]:
                    connect.append(r + 1)
            if len(connect)==3:
                f.write(
                    "f %d %d %d\n" % (int(connect[0]), int(connect[1]), int(connect[2]))
            )


def export_obj_surf(infile, outfile, newposition, deleted_positions):
    # load Organoid informations from inp
    apicalOut = outfile.replace(".obj", "_apical.obj")  # generated during export
    basalOut = outfile.replace(".obj", "_basal.obj")  # generated during export
    lateralOut = outfile.replace(".obj", "_lateral.obj")  # generated during export

    SurfAp=[]
    SurfBas=[]
    SurfLat=[]
    surfaceA = []
    surfaceB = []
    surfaceL = []
    surfaceAconnect = []
    surfaceBconnect = []
    surfaceBconnectOBJ = []
    surfaceLconnect = []
     
    surfaceAD = []
    surfaceBD = []
    surfaceLD = []
    surfaceAconnectD = []
    surfaceBconnectD = []
    surfaceBconnectDOBJ = []
    surfaceLconnectD = []
    elements = readElementsINP(infile)
    nodes = readNodesINP(infile)
    nsets = readNsetsINP(infile)
    elsets,elsetsO = readElsetsINP(infile)
    
    
    for i in range(elsetsO["MEMAPICAL"][0],elsetsO["MEMAPICAL"][1]):
        if [nodes[elements[i][0]], elements[i][0]] not in surfaceA:
            surfaceA.append([nodes[elements[i][0]], elements[i][0]])
            surfaceAD.append([newposition[elements[i][0]], elements[i][0]])
        if [nodes[elements[i][1]], elements[i][1]] not in surfaceA:
            surfaceA.append([nodes[elements[i][1]], elements[i][1]])
            surfaceAD.append([newposition[elements[i][1]], elements[i][1]])
        if [nodes[elements[i][2]], elements[i][2]] not in surfaceA:
            surfaceA.append([nodes[elements[i][2]], elements[i][2]])
            surfaceAD.append([newposition[elements[i][2]], elements[i][2]])
        surfaceAconnect.append([elements[i][0], elements[i][1], elements[i][2]])
        surfaceAconnectD.append([elements[i][0], elements[i][1], elements[i][2]])
    for i in range(elsetsO["MEMBASAL"][0],elsetsO["MEMBASAL"][1]):
        if [nodes[elements[i][0]], elements[i][0]] not in surfaceB:
            surfaceB.append([nodes[elements[i][0]], elements[i][0]])
            surfaceBD.append([newposition[elements[i][0]], elements[i][0]])
        if [nodes[elements[i][1]], elements[i][1]] not in surfaceB:
            surfaceB.append([nodes[elements[i][1]], elements[i][1]])
            surfaceBD.append([newposition[elements[i][1]], elements[i][1]])
        if [nodes[elements[i][2]], elements[i][2]] not in surfaceB:
            surfaceB.append([nodes[elements[i][2]], elements[i][2]])
            surfaceBD.append([newposition[elements[i][2]], elements[i][2]])
        surfaceBconnect.append([elements[i][0], elements[i][1], elements[i][2]])
        surfaceBconnectD.append([elements[i][0], elements[i][1], elements[i][2]])    
    for i in range(elsetsO["MEMLATERAL"][0],elsetsO["MEMLATERAL"][1]):
        if [nodes[elements[i][0]], elements[i][0]] not in surfaceL:
            surfaceL.append([nodes[elements[i][0]], elements[i][0]])
            surfaceLD.append([newposition[elements[i][0]], elements[i][0]])
        if [nodes[elements[i][1]], elements[i][1]] not in surfaceL:
            surfaceL.append([nodes[elements[i][1]], elements[i][1]])
            surfaceLD.append([newposition[elements[i][1]], elements[i][1]])
        if [nodes[elements[i][2]], elements[i][2]] not in surfaceL:
            surfaceL.append([nodes[elements[i][2]], elements[i][2]])
            surfaceLD.append([newposition[elements[i][2]], elements[i][2]])
        surfaceLconnect.append([elements[i][0], elements[i][1], elements[i][2]])
        surfaceLconnectD.append([elements[i][0], elements[i][1], elements[i][2]])

    apicalOutDeformed = outfile.replace(".obj", "_apical_Def.obj")  # generated during export
    basalOutDeformed = outfile.replace(".obj", "_basal_Def.obj")  # generated during export
    lateralOutDeformed = outfile.replace(".obj", "_lateral_Def.obj")  # generated during export
    
    from os.path import exists as file_exists
    if file_exists(apicalOut)==False:
        print("export surface OBJ Step 0")
        surface_obj(surfaceA, surfaceAconnect, apicalOut)
        surface_obj(surfaceB, surfaceBconnect, basalOut)
        surface_obj(surfaceL, surfaceLconnect, lateralOut)
    print("export surface OBJ Step 0")
    surface_obj(surfaceA, surfaceAconnect, apicalOut)
    surface_obj(surfaceB, surfaceBconnect, basalOut)
    surface_obj(surfaceL, surfaceLconnect, lateralOut)
    print("export surface OBJ Step 1")
    surface_obj(surfaceAD, surfaceAconnect, apicalOutDeformed)
    surface_obj(surfaceBD, surfaceBconnect, basalOutDeformed)
    surface_obj(surfaceLD, surfaceLconnect, lateralOutDeformed)

    return surfaceAconnect,surfaceAconnectD,surfaceBconnect,surfaceBconnectD,surfaceLconnect,surfaceLconnectD


import numpy as np
def misesfrom(stress):
    S1=float(stress[0])
    S2=float(stress[1])
    S3=float(stress[2])
    mises= np.sqrt(((S1-S2)**2+(S2-S3)**2+(S3-S1)**2)*1/2)
    return mises


import sys
import hickle as hkl
import os

if os.path.isfile('./FIL_read.hkl'):
    inpNodes,inpElements,inpElsets, inpElsets0, inpNsets, inpNodes1, inpNodes2, inpElements2, inpElsets2,inpElsets02,inpNsets2 =  hkl.load('analysis_read.hkl')
    
else:
    ### Import informations from inp and fil files
    precision = 6 #precision to uniform floats from inp and fil files 
    
    nodefilfil = open(nodefilfile, "r")
    import time
    start = time.time()
    results=""
    myfile = open('concat.txt', 'w')
    with open(nodefilfile) as file:
        for l in file:
            l = l.replace("\n", "")
            l = l.replace(" ", "")
            myfile.write("%s" % l)
    myfile.close()
    myfile = open("concat.txt", "r")
    results = myfile.readline()
    myfile.close()
    end =  time.time()
    
    
    print("Read .fil abaqus output file and it concatenate in a simple line of ",len(results), "characters")
    print("Execution time in seconds: ",(end-start), "\n")
    steps=results.split("I223I42000D") # search I 223I 42000D
    print("Number of steps increment identified: ",len(steps)-1)
    print("Start analyzing for the last step id",len(steps)-1,"\n")
    
    NodesID={}
    NodesID1={}
    if "16I3107I" in steps[1]:
        #print("3D coordinate properties found as 16I3107I ID")
        string=steps[1].replace("I19I41901I", "")
        node_string=string.split("I16I3107I")
        #print(node_string[2])
        for i in node_string:
            try:
                elem=i.split("*")
                elem=elem[0].split("D")
                NodesID1[int(elem[0][1:])]=[
                    "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                    "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                    "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                    ]
                
            except:
                pass
    ### FIRST STEP
    NodesFS={}
    NodesFS1={}
    if "16I3107I" in steps[11]:
        #print("3D coordinate properties found as 16I3107I ID for step 1")
        string=steps[11].replace("I12I42001", "")
        node_string=string.split("I16I3107I")
        #print(node_string[2])
        for i in node_string:
            try:
                elem=i.split("*")
                elem=elem[0].split("D")
                NodesFS1[int(elem[0][1:])]=[
                    "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                    "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                    "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                    ]
            except:
                pass
    #print(len(NodesFS))
    ### LAST STEP
    NodesLS={}
    NodesLS1={}
    if "16I3107I" in steps[33]:
        print("3D coordinate properties found as 16I3107I ID for last step")
        string=steps[33].replace("I12I42001", "")
        node_string=string.split("I16I3107I")
        #print(node_string[2])
        for i in node_string:
            try:
                elem=i.split("*")
                elem=elem[0].split("D")
                #print("element of the last step")
                NodesLS1[int(elem[0][1:])]=[
                    "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                    "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                    "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                    ]        
            except:
                pass

    UID={}
    if "19I3101I" in steps[1]:
        #print("3D coordinate properties found as 16I3107I ID")
        string=steps[1].replace("I19I41901I", "")
        node_string=string.split("19I3101I")
        for i in node_string:
            try:
                elem=i.replace("*","")
                elem=elem.split("D")
                UID[int(elem[0][1:])]=[
                    "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                    "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                    "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                    ]
            except:
                pass
    UFS={}    # FIRST STEP
    if "19I3101I" in steps[11]:
        #print("3D coordinate properties found as 16I3107I ID for step 1")
        string=steps[11].replace("I12I42001", "")
        node_string=string.split("19I3101I")
        for i in node_string:
            try:
                elem=i.split("*")
                elem=elem[0].split("D")
                UFS[int(elem[0][1:])]=[
                    "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                    "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                    "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                    ]
            except:
                pass
    ULS={}  # LAST STEP
    if "19I3101I" in steps[33]:
        #print("3D coordinate properties found as 16I3107I ID for step [-1]")
        string=steps[33].replace("I12I42001", "")
        node_string=string.split("19I3101I")
        
        for i in node_string:
            try:
                node=i.split("*")
                node=node[0].split("D")
                #print("element of the last step")
                ULS[int(node[0][1:])]=[
                    "{:.{}f}".format( float(node[1])*10**int(node[2]), precision ),
                    "{:.{}f}".format( float(node[3])*10**int(node[4]), precision ),
                    "{:.{}f}".format( float(node[5])*10**int(node[6]), precision )
                    ]        
            except:
                pass
    #print(ULS.keys())

    if len(UID) >= 1:
        for n, node in NodesID1.items():
                try:
                    NodesID[n]=[
                        float(node[0])+float(UID[n][0]),
                        float(node[1])+float(UID[n][1]),
                        float(node[2])+float(UID[n][2]),                
                        ]
                except:
                    NodesID[n]=[
                        float(node[0]),
                        float(node[1]),
                        float(node[2]),                
                        ]
    if len(UFS) >= 1:
        for n, node in NodesFS1.items():
                try:
                    NodesFS[n]=[
                        float(node[0])+float(UFS[n][0]),
                        float(node[1])+float(UFS[n][1]),
                        float(node[2])+float(UFS[n][2]),                
                        ]
                    
                except:
                    NodesFS[n]=[
                        float(node[0]),
                        float(node[1]),
                        float(node[2]),                
                        ]
    if len(ULS) >= 1:
        for n, node in NodesLS1.items():
                try:
                    NodesLS[n]=[
                        float(node[0])+float(ULS[n][0]),
                        float(node[1])+float(ULS[n][1]),
                        float(node[2])+float(ULS[n][2]),                
                        ]
                except:
                    NodesLS[n]=[
                        float(node[0]),
                        float(node[1]),
                        float(node[2]),                
                        ]
    
    print("first point to test:")
    print(list(NodesID.keys())[0])
    print("test position at increment 11")
    #print(NodesID[list(ULS.keys())[0]])
    print("test position at increment 22")
    #print(NodesFS[list(ULS.keys())[0]])
    print("test position at increment 33")
    #print(NodesLS[list(ULS.keys())[0]])

    ################ ELEMENTS ######## FIL OUTPUT
    from collections import defaultdict
    print("Research element porperties")
    precision=9

    import time
    results=""
    start = time.time()
    count = 0
    myfile = open('concat.txt', 'w')
    with open(elfilfile) as file:
        for l in file:
            l = l.replace("\n", "")
            l = l.replace(" ", "")
            count = count + 1
            #results += l
            myfile.write("%s" % l)
    myfile.close()
    myfile = open("concat.txt", "r")
    results = myfile.readline()
    myfile.close()
    end =  time.time()
    print("Execution time in seconds: ",(end-start))

    print("Read .fil abaqus output file and it concatenate in a simple line of ",len(results), "characters\n")
    steps=results.split("I223I42000D") # search I 223I 42000D
    ### FIRST STEP
    ElementsStrainFS= defaultdict(list)
    ElementsStressFS= defaultdict(list)
    ElementsElseFS={}
    step2=""
    try:
        steps2=steps[1].split("*I15I41911I10A") # split element type
    except:
        print("No element informations")
    
    for e in steps2:
        if "*I15I291D" in e: #*I 15I 291D for Strain
            node_string=e.split("*I211I11I")
            for i in node_string:
                try:
                    elem1=i.split("*I")
                    elem=elem1[1].split("D")
                    element=elem1[0].replace("I11I11I10AI12I11I10I16", "")
                    element=element.replace("I11I15I10AI12I11I10I16", "")
                    Els=[
                        "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                        "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                        "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                        ]
                    for el in Els:
                        ElementsStrainFS[int(element[1:])].append(el)
                except:
                    pass
        elif "*I15I211D" in e: #*I 15I 211D for Stress
            node_string=e.split("*I211I11I")
            for i in node_string:
                try:
                    elem1=i.split("*I")
                    elem=elem1[1].split("D")
                    element=elem1[0].replace("I11I11I10AI12I11I10I16", "")
                    element=element.replace("I11I15I10AI12I11I10I16", "")
                    Els=[
                        "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                        "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                        "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                        ]
                    for el in Els:
                        ElementsStressFS[int(element[1:])].append(el)
                except:
                    pass
        elif "*I212I219D" in e: #*I 15I 211D for Stress
            node_string=e.split("*I211I11I")
            for i in node_string:
                try:
                    elem1=i.split("*I")
                    elem=elem1[1].split("D")
                    element=elem1[0].replace("I10I10I15AI12I11I10I16", "")
                    ElementsElseFS[int(element[1:])]=float("{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ))
                except:
                    pass
    ### LAST STEP
    ElementsStrainLS= defaultdict(list)
    ElementsStressLS= defaultdict(list)
    ElementsElseLS={}
    step2=""
    try:
        steps2=steps[-1].split("*I15I41911I10A") # split element type
    except:
        print("No element informations")
    for e in steps2:
        if "*I15I291D" in e: #*I 15I 291D for Strain
            node_string=e.split("*I211I11I")
            for i in node_string:
                try:
                    elem1=i.split("*I")
                    elem=elem1[1].split("D")
                    element=elem1[0].replace("I11I11I10AI12I11I10I16", "")
                    element=element.replace("I11I15I10AI12I11I10I16", "")
                    Els=[
                        "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                        "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                        "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                        ]
                    for el in Els:
                        ElementsStrainLS[int(element[1:])].append(el)
                except:
                    pass
        elif "*I18I211D" in e: #*I 15I 211D for Stress
            if "C3D10M" in e:
                node_string=e.split("*I211I11I")
                for i in node_string:
                    try:
                        elem1=i.split("*I")
                        elem=elem1[1].split("D")
                        element=elem1[0]
                        element=elem1[0].replace("I11I10I10AI13I13I10I10", "")
                        element=element.replace("I12I10I10AI13I13I10I10", "")
                        element=element.replace("I13I10I10AI13I13I10I10", "")
                        element=element.replace("I14I10I10AI13I13I10I10", "")
                        Els=[
                            "{:.{}f}".format( float(elem[1])*10**int(elem[2]), precision ),
                            "{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ),
                            "{:.{}f}".format( float(elem[5])*10**int(elem[6]), precision )
                            ]
                        for el in Els:
                            ElementsStressLS[int(element[1:])].append(el)
                    except:
                        pass
        elif "*I212I219D" in e: #*I 15I 211D for Stress
            if "C3D10M" in e:
                node_string=e.split("*I211I11I")
                for i in node_string:
                    try:
                        elem1=i.split("*I")
                        elem=elem1[1].split("D")
                        element=elem1[0].replace("I10I10I15AI12I11I10I16", "")
                        element=element.replace("I10I10I15AI13I13I10I10", "")
                        ElementsElseLS[int(element[1:])]=float("{:.{}f}".format( float(elem[3])*10**int(elem[4]), precision ))
                    except:
                        pass


#print(ElementsStressLS)

inpNodes=readNodesINP(original_file)
inpElements=readElementsINP(original_file)
inpElsets,inpElsets0=readElsetsINP(original_file) # Elsets0 correspond à elsets cells
inpNsets=readNsetsINP(original_file)
inpNodes1=readNodesINP(beforedef_file)
inpNodes2=readNodesINP(refined_file)
inpElements2=readElementsINP(refined_file)
inpElsets2,inpElsets02=readElsetsINP(refined_file) # Elsets0 correspond à elsets cells
inpNsets2=readNsetsINP(refined_file)
    
print("Number of nodes in inp file",len(inpNodes))
print("Number of nodes in fil file",len(NodesID))

CellELSE={}
try:
    
    delta=inpElsets02['ALL_INSIDEPART'][1]
except:
    delta=inpElsets02['ALL_WATER'][-1]
"""
print("delta elements")
print(delta)
"""
for cells, elements in inpElsets02.items():
    if "Cell" in cells:
        elese=0
        nbe=0
        try:
            c = int(cells.replace("Cell", ""))
        except:
            pass
        
        for e in elements:
            try:
                elese+=ElementsElseLS[e+delta]
                nbe+=1
            except:
                pass
        elese=elese/nbe
        CellELSE[c]=elese
    elif "CELL" in cells:
        elese=0
        nbe=0
        try:
            c = int(cells.replace("CELL", ""))
        except:
            pass
        
        for e in elements:
            try:
                elese+=ElementsElseLS[e+delta]
                nbe+=1
            except:
                pass
        elese=elese/nbe
        CellELSE[c]=elese



CellMISES={}
for cells, elements in inpElsets02.items():
    if "Cell" in cells:
        smises=0
        nbe=0
        try:
            c = int(cells.replace("Cell", ""))
        except:
            pass
        for e in elements:
            try:
                smises+=misesfrom(ElementsStressLS[e+delta])
                nbe+=1
            except:
                pass
        smises=smises/nbe
        CellMISES[c]=smises
    elif "CELL" in cells:
        smises=0
        nbe=0
        try:
            c = int(cells.replace("CELL", ""))
        except:
            pass
        for e in elements:
            try:
                smises+=misesfrom(ElementsStressLS[e+delta])
                nbe+=1
            except:
                pass
        #print(cells)
        if nbe != 0:
            smises=smises/nbe
            CellMISES[c]=smises

try:
    Delta=inpNsets2['INSIDEPART'][1] #Delta=5990 pour les analyses 401cells du 8 juin

except:
    Delta=inpNsets2['ALL_WATER'][1]

"""    
print("Delta nodes")
print(Delta)

newposition={}
oldposition = readNodesINP(original_file)
for i, j in oldposition.items():
    try:
        newposition[i]=NodesLS[i+Delta] #+Delta
    except:
        print("node is missing ", i+Delta)

"""
oldposition = readNodesINP(original_file)

newposition={}
for i, n in inpNodes.items():
    #newpos=[float( NodesLS[i+Delta][0]),float( NodesLS[i+Delta][1]),float( NodesLS[i+Delta][2])]
    newpos=[float(NodesLS[i+Delta][0]),float(NodesLS[i+Delta][1]),float(NodesLS[i+Delta][2])]
    newposition[i]= newpos 

NodesID
deleted_positions=[]
#print(newposition)
#print(i)


outputobjFileName = original_file.replace(".inp", ".obj")  # generated during export
triApi,triApiDef,triBas,triBasDef,triLat,triLatDef=export_obj_surf(original_file, outputobjFileName, newposition, deleted_positions)

apicalOBJ = outputobjFileName.replace(".obj", "_apical_Def.obj")  # generated during export
basalOBJ = outputobjFileName.replace(".obj", "_basal_Def.obj")  # generated during export

########################  Apico-basal and Basal-lateral AREA RATIO  #########################

Api_Bas_ratio={}
Bas_Lat_ratio={}
xapi=[]
ybas=[]
zlat=[]

# and after deformation
Api_Bas_ratiod={}
Bas_Lat_ratiod={}
xapid=[]
ybasd=[]
zlatd=[]

apical_surf={}
basal_surf={}
apical_lateral_surf={}
basal_lateral_surf={}
neibhoors_number={}


for cell, elems in inpElsets0.items():
    if "CELL" in cell:
        apiArea=0
        basArea=0
        latArea=0
    
        apiAread=0
        basAread=0  
        latAread=0
        neighbor=0
        for el in elems:
            if len(inpElements[el]) == 10:
                trinodes=[inpElements[el][1],inpElements[el][5],inpElements[el][8]]
                trinodes2=[inpElements[el][8],inpElements[el][5],inpElements[el][1]]
                if trinodes in triApi or trinodes2 in triApi:
                    neighbor+=1
                    try:
                       apiArea+=float(areatriangleelement(el,inpElements,oldposition))
                       apiAread+=float(areatriangleelement(el,inpElements,newposition))
                    except:
                        print("warning! error with area element apical ",el," on cell ", cell.replace("CELL",""))
                elif trinodes in triBas or trinodes2 in triBas:
                    try:
                        basArea+=float(areatriangleelement(el,inpElements,oldposition))
                        basAread+=float(areatriangleelement(el,inpElements,newposition))
                    except:
                        print("warning! error with area element basal ",el," on cell ", cell.replace("CELL",""))
                elif trinodes in triLat or trinodes2 in triLat:
                    try:
                        latArea+=float(areatriangleelement(el,inpElements,oldposition))
                        latAread+=float(areatriangleelement(el,inpElements,newposition))
                    except:
                        print("warning! error with area element lateral ",el," on cell ", cell.replace("CELL",""))
                
        if apiArea != 0 and basArea!= 0 and latAread != 0:
            Api_Bas_ratio[cell.replace("CELL","")]=float(apiArea/basArea)
            Api_Bas_ratiod[cell.replace("CELL","")]=float(apiAread/basAread)
            Bas_Lat_ratio[cell.replace("CELL","")]=float(basArea/latArea)
            Bas_Lat_ratiod[cell.replace("CELL","")]=float(basAread/latAread)
            apical_surf[cell.replace("CELL","")]=float(apiAread)
            basal_surf[cell.replace("CELL","")]=float(basAread)
            apical_lateral_surf[cell.replace("CELL","")]=float(apiAread/latAread)
            basal_lateral_surf[cell.replace("CELL","")]=float(basAread/latAread)
            neibhoors_number[cell.replace("CELL","")]=neighbor
            
            
            xapi.append(apiArea)
            ybas.append(basArea)
            zlat.append(latArea)
        
            xapid.append(apiAread)
            ybasd.append(basAread)
            zlatd.append(latAread)
        else:
            print("no area for cell ",cell.replace("CELL",""))
        
#### Volumes of the cells ####
print("recording cell volumes")
"""
apical_surf["367"]=0
basal_surf["367"]=0
apical_lateral_surf["367"]=0
basal_lateral_surf["367"]=0
neibhoors_number["367"]=0

apical_surf["360"]=0
basal_surf["360"]=0
apical_lateral_surf["360"]=0
basal_lateral_surf["360"]=0
neibhoors_number["360"]=0

apical_surf["223"]=0
basal_surf["223"]=0
apical_lateral_surf["223"]=0
basal_lateral_surf["223"]=0
neibhoors_number["223"]=0
"""
Volume_cell={}
for cell, elems in inpElsets0.items():
    if "CELL" in cell:
        Volume=0
        for el in elems:
            try:
                Volume+=float(volumetetrahedronelementD10(el,inpElements,newposition))
            except:
                print("warning! error with volume element ",el ," on cell ", cell.replace("CELL",""))
        
        Volume_cell[cell.replace("CELL","")]=float(Volume)


DifVolume_cell={}
for cell, elems in inpElsets0.items():
    if "CELL" in cell:
        Volume1=0
        Volume2=0
        for el in elems:
            try:
                Volume1+=float(volumetetrahedronelementD10(el,inpElements,newposition))
                Volume2+=float(volumetetrahedronelementD10(el,inpElements,oldposition))
            except:
                print("warning! error with volume element ",el ," on cell ", cell.replace("CELL",""))

        Volume3=Volume1-Volume2
        DifVolume_cell[cell.replace("CELL","")]=float(Volume3)


# 2D scatter plot for cells inforamtions
from matplotlib import pyplot as plt
from vedo import *
import numpy as np
import pandas as pd

barycenters={}
x=1

print("recording cell barycenters")
for cell, elems in inpElsets02.items():
    if "CELL" in cell:
        #print(x)
        barycenters[x]=newposition[x+Delta]
        x+=1
        """
        try:
            barycenters[x]=newposition[x]
            x+=1
        except:
            barycenters[x]=oldposition[x]
            x+=1
        """
######################## AREA and Volume inforamtions plot and export #########################

print("preparing individual cell plot ")


# load mesh from .obj to plot with vedo
mesh = load(apicalOBJ).compute_normals().linewidth(1)
#mesh = load(apicalOBJ).computeNormals().lineWidth(1)#old vedo function
mesh2 = load(basalOBJ).compute_normals().linewidth(1)
mesh3 = load(apicalOBJ).compute_normals().linewidth(1)
mesh4 = load(basalOBJ).compute_normals().linewidth(1)
mesh5 = load(apicalOBJ).compute_normals().linewidth(1)
mesh6 = load(apicalOBJ).compute_normals().linewidth(1)
mesh7 = load(apicalOBJ).compute_normals().linewidth(1)

# generate a numpy array for cell ratio map on vtk format (vedo)

print("generating list for cell aspect ratio Api/Bas")

ratiolist=[]
for tri in triApiDef:
    for cell, elems in inpElsets02.items():
        if "CELL" in cell:
            for e in elems:
                if tri == [inpElements[e][0],inpElements[e][1],inpElements[e][2]]:
                    ratiolist.append(Api_Bas_ratio[cell.replace("CELL","")])
                    break

print("generating list for cell aspect ratio Bas/Lat")

ratiolist2=[]
for tri in triApiDef:
    for cell, elems in inpElsets02.items():
        if "CELL" in cell:
            for e in elems:
                if tri == [inpElements[e][0],inpElements[e][1],inpElements[e][2]]:
                    ratiolist2.append(Bas_Lat_ratio[cell.replace("CELL","")])
                    break

print("generating list for cell strain ELSE")

eleselist=[]
for tri in triApiDef:
    for cell, elems in inpElsets02.items():
        if "Cell" in cell:
            for e in elems:
                if tri == [inpElements[e][0],inpElements[e][1],inpElements[e][2]]:
                    eleselist.append(CellELSE[int(cell.replace("CELL",""))])
                    break

print("generating list for cell volume")

volumelist=[]
for tri in triApiDef:
    for cell, elems in inpElsets02.items():
        if "CELL" in cell:
            for e in elems:
                if tri == [inpElements[e][1],inpElements[e][2],inpElements[e][3]]:
                    volumelist.append(Volume_cell[cell.replace("CELL","")])
                    break

print("generating list for cell volume variation")

difvolumelist=[]
for tri in triApiDef:
    for cell, elems in inpElsets02.items():
        if "CELL" in cell:
            for e in elems:
                if tri == [inpElements[e][1],inpElements[e][2],inpElements[e][3]]:
                    difvolumelist.append(DifVolume_cell[cell.replace("CELL","")])
                    break

print("generating list for cell MISES stress")

miseslist=[]
for tri in triApiDef:
    for cell, elems in inpElsets02.items():
        if "CELL" in cell:
            for e in elems:
                if tri == [inpElements[e][1],inpElements[e][2],inpElements[e][3]]:
                    miseslist.append(CellMISES[int(cell.replace("CELL",""))])
                    break

scals = np.array(ratiolist)
scals2 = np.array(eleselist)
scals3 = np.array(volumelist)
scals4 = np.array(miseslist)
scals5 = np.array(ratiolist2)
scals6 = np.array(difvolumelist)

# define color map using cells as determinant
#mesh.cmap('jet', scals5, on='cells').addScalarBar3D() # old vedo version
mesh.cmap('jet', scals5, on='cells').add_scalarbar3d() #mesh.cmap('RdYlBu', scals, on='cells')
mesh2.cmap('jet', scals5, on='cells').add_scalarbar3d()
mesh3.cmap('jet', scals6, on='cells').add_scalarbar3d() #mesh.cmap('RdYlBu', scals, on='cells')
mesh4.cmap('jet', scals6, on='cells').add_scalarbar3d() #mesh.cmap('RdYlBu', scals, on='cells')

# define the camera
cam = dict(pos=(59.8, -191, 78.9),
           focalPoint=(27.9, -2.94, 3.33),
           viewup=(-0.0170, 0.370, 0.929),
           distance=205,
           clippingRange=(87.8, 355))
plt = Plotter(N=6, axes=0)
# VEDO plot based on vtk
#show(mesh, __doc__, bg='w', camera=cam, axes=11) # apical
#cmesh = mesh.cutWithPlane(normal=(50,-1,-1)) # old vedo version
cmesh = mesh.cut_with_plane(normal=(50,-1,-1)) # halfcut
cmesh2 = mesh2.cut_with_plane(normal=(50,-1,-1))# halfcut
cmesh3 = mesh3.cut_with_plane(normal=(50,-1,-1)) # halfcut
cmesh4 = mesh4.cut_with_plane(normal=(50,-1,-1))# halfcut

plt.show(cmesh3,cmesh4, "Volume changes", at=3) # both apical and basal
#show(cmesh,cmesh2, __doc__, bg='w', camera=cam, axes=11) # both apical and basal
plt.show(cmesh,cmesh2, "Basal/Lateral ratio", at=0, interactive=1).close()


print("recording cell dispersion from barycenter distance")



from sklearn.neighbors import KDTree
bary = np.array(list(barycenters.values()))
tree= KDTree(bary)
nearst_dist, nearest_ind = tree.query(bary, k=5)
nearst_dist = nearst_dist

import pandas as pd
vols= np.fromiter(Volume_cell.values(), dtype=float)/100
Stats=pd.DataFrame(vols,columns=['volume'])

Stats["spacing"]=np.array(list(nearst_dist[:,2]))/2
Stats["apical"]=np.fromiter(apical_surf.values(), dtype=float)/40
Stats["basal"]=np.fromiter(basal_surf.values(), dtype=float)/40
Stats["basal-lateral"]=basal_lateral_surf.values()
Stats["apical-lateral"]=apical_lateral_surf.values()
Stats["neibhoors_number"]=neibhoors_number.values()

new_csv=refined_file.replace(".inp","_statistiques.csv")
new_csv2=refined_file.replace(".inp","_statistiques_post.csv")

#Stats.to_csv(new_csv)

Stats = Stats[Stats.apical > 1]
Stats = Stats[Stats.basal > 1]
Stats["basal_apical"]=Stats.basal/Stats.apical
Stats.to_csv(new_csv2)
print("statistiques saved in: "+str(new_csv2))