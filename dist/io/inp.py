import numpy as np
import pandas as pd
import sys, os, inspect
import subprocess
from multiprocessing import Pool


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))
    
def save_inp_individual(inp_file, mono,elemtosup=None, add_spaces=None, coords=None):
    """
    Write tetrahedra in inp format for Finite Element Modelisation in Abaqus
    """    
    if coords is None:
        coords = list("xyz")
    all_nodes = pd.concat(
        [mono.cell_df[coords], mono.face_df[coords], mono.vert_df[coords]],
        axis=0,
        ignore_index=True,
    )
    tetrahedra = mono.edge_df[["cell", "face", "srce", "trgt"]].copy()
    tetrahedra["face"] += mono.Nc
    tetrahedra[["srce", "trgt"]] += mono.Nc + mono.Nf
    tetrahedra.index += 1
    tetrahedra += 1
    all_nodes.index += 1
    array = range(1,len(mono.cell_df)+1)
    side = mono.face_df['segment']
    nsetside= side.to_frame().groupby('segment').groups
    apical = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[0]])].copy()
    basal = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[1]])].copy()
    lateral = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[2]])].copy()
    elcells = tetrahedra['cell'].loc[tetrahedra['cell'].isin(array)].to_frame().copy()
    elcells = elcells.groupby('cell').groups

    # "individual segmentation"
    nodes_tocreate=[]
    Cell_nodes={}
    nodeMax=len(all_nodes)
    tetrahedra2=tetrahedra.copy()
    for cells_id, elems_ids in elcells.items():
        if cells_id not in Cell_nodes.keys():
            Cell_nodes[cells_id]=[]
        for el in elems_ids:
            for node in tetrahedra2.values[el-1]:
                if any(node in i for i in Cell_nodes.values()) == False or node in Cell_nodes[cells_id]:
                    Cell_nodes[cells_id].append(node)
                else:
                    nodeMax+=1
                    nodes_tocreate.append(list(all_nodes.values[node-1]))
                    column=list(tetrahedra.values[el-1]).index(node)
                    tetrahedra.iloc[el-1, column]=nodeMax

    df2=pd.DataFrame(nodes_tocreate, columns=['x','y','z'])
    all_nodes=pd.concat([all_nodes, df2], axis=0, ignore_index=True)
    all_nodes.index += 1

    # deplacer les nodes en direction du centroid de la cellule
    with open(inp_file, "w+") as inph:
        inph.write(f"*part, name={mono.identifier}\n")
        inph.write("*node\n")
    with open(inp_file, "a+") as inph:
        all_nodes.to_csv(inph, header=False, index=True, sep=",", line_terminator='\n')
    with open(inp_file, "a+") as inph:
        inph.write("*element, type=C3D4\n")
    with open(inp_file, "a+") as inph:
        tetrahedra.to_csv(inph, header=False, index=True, sep=",", line_terminator='\n')
    for i in range(len(elcells)):
        with open(inp_file, "a+") as inph:
            inph.write("*elset, elset=Cell"+str(i+1)+"\n")
            for x, a in enumerate(elcells[i+1]):
                if a == elcells[i+1][-1]:
                    if x % 16 == 15:
                        inph.write("\n")
                    inph.write(str(a)+"\n")
                elif x % 16 == 15:
                    inph.write(str(a)+"\n")
                else:
                    inph.write(str(a)+",")
            inph.write("*Surface, type=ELEMENT, name=Cell"+str(i+1)+"_S3\n")
            inph.write("Cell"+str(i+1)+", S3\n")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFapicalN\n")
        for x, a in enumerate(apical.index):
            if x+1 == len(apical.index):
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFbasalN\n")
        for x, a in enumerate(basal.index):
            if x+1 == len(basal.index):
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFlateralN\n")
        for x, a in enumerate(lateral.index):
            if x+1 == len(lateral.index):
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

def ReadInp(fi, lineCount):
    s = fi.readline()
    if s != "":
        lineCount = lineCount + 1
    # endif
    return s, lineCount



def readNodesINP(inputFileName):
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
                        break
            iLine += 1
    return nodes

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
        fo.write("%d, %f, %f, %f\n" % (int(MaxNode), float(fx), float(fy), float(fz)))
        # endif
    return elements, median, MaxNode

def decimate_inp(fin, fout, elemDim):
    # global lineCount
    lineCount = 0
    # ------------File----------------------------------------
    fi = open(fin, "r")
    fo = open(fout, "w")
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
            continue
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
    i = s2[1].find("C3D10M")


    elset_cells, elsetsOs, elsetOthers=read_elsets(fin)
    nsets=read_nsets(fin)
    if elemDim ==4:
        elemsStri, elemsS3RT =read_elems_types(fin)

    apicalL = []
    basalSL = []
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
                if elemDim == 3:
                    print("###3D-C3 found")
                    while True:
                        s, lineCount = ReadInp(fi, lineCount)
                        if s == "":
                            break
                        if len(s)==1:
                            continue
                        s1 = s.strip()  # read line
                        i = s1.find("elset")
                        if i != -1:
                            break
                        i = s1.find("nset")
                        if i != -1:
                            break
                        s2 = s1.split(",")
                        x = int(s2[0])
                        elements.setdefault(x, [])
                        for i in range(4): # add the first 4 nodes to the element
                            elements[x].append(int(s2[i+1])) 
                        elements, median, MaxNode = WriteIntermediatePoint(
                            1, 2, elements, nodes, median, x, fo, MaxNode, s2
                        ) # create node 5
                        elements, median, MaxNode = WriteIntermediatePoint(
                            2, 3, elements, nodes, median, x, fo, MaxNode, s2
                        ) # create node 6
                        elements, median, MaxNode = WriteIntermediatePoint(
                            1, 3, elements, nodes, median, x, fo, MaxNode, s2
                        ) # create node 7
                        elements, median, MaxNode = WriteIntermediatePoint(
                            1, 4, elements, nodes, median, x, fo, MaxNode, s2
                        ) # create node 8
                        elements, median, MaxNode = WriteIntermediatePoint(
                            2, 4, elements, nodes, median, x, fo, MaxNode, s2
                        ) # create node 9
                        elements, median, MaxNode = WriteIntermediatePoint(
                            3, 4, elements, nodes, median, x, fo, MaxNode, s2
                        ) # create node 10
                    sType = sType[iNodes:]
                    ss = "*element, type=C3D10M"
                    fo.write(ss)
                    fo.write("\n")
                    for i in elements:
                        fo.write(
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
                    break
                fo.write("%s\n" % s1)




    elif elemDim == 4:
        while True:
            elements3 = {}
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
            if sType.find("C3D10M") != -1:
                fint=fi
                lineCount2=lineCount
                while True:
                    ss, lineCount2 = ReadInp(fint, lineCount2)
                    if ss == "":
                        break
                    ss1 = ss.strip()  # read line
                    i = ss1.find("elset")
                    if i != -1:
                        break
                    i = ss1.find("nset")
                    if i != -1:
                        break
                    i = ss1.find("element")
                    if i != -1:
                        break
                    ss2 = ss1.split(",")
                    x = int(ss2[0])
                    elements3.setdefault(x, [int(i) for i in ss2])
                break
        Nnodes, Nelements = c3d10_C3D4s(nodes, elements3.values())
        elements = {}
        median = {}   

        for x,y in enumerate(Nelements):
            x+=1
            elements.setdefault(x, [])  
            strr=[str(x+1),str(y[0]),str(y[1]),str(y[2]),str(y[3])]
            for i in range(4): # add the first 4 nodes to the element
                elements[x].append(int(y[i])) 

            elements, median, MaxNode = WriteIntermediatePoint(
                1, 2, elements, nodes, median, x, fo, MaxNode, strr)
            elements, median, MaxNode = WriteIntermediatePoint(
                2, 3, elements, nodes, median, x, fo, MaxNode, strr)
            elements, median, MaxNode = WriteIntermediatePoint(
                1, 3, elements, nodes, median, x, fo, MaxNode, strr)
            elements, median, MaxNode = WriteIntermediatePoint(
                1, 4, elements, nodes, median, x, fo, MaxNode, strr)
            elements, median, MaxNode = WriteIntermediatePoint(
                2, 4, elements, nodes, median, x, fo, MaxNode, strr)
            elements, median, MaxNode = WriteIntermediatePoint(
                3, 4, elements, nodes, median, x, fo, MaxNode, strr)
 
        ss = "*element, type=C3D10M"
        fo.write(ss)
        fo.write("\n")
        for i in elements:
            fo.write(
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


    # when SUBTETRA creation: C3D10 --> C3D4s increase the number of final elements
    # add 8x elements to cell elsets 
    if elemDim == 4:
        Cortex={}
        ElsetSS={}
        for i, j in elset_cells.items():
            new_list=[]
            for x in j:
                new_list.append(x*8-7)
                new_list.append(x*8-6)
                new_list.append(x*8-5)
                new_list.append(x*8-4)
                new_list.append(x*8-3)
                new_list.append(x*8-2)
                new_list.append(x*8-1)
                new_list.append(x*8)
            elset_cells[i]=new_list
        for i,j in elsetsOs.items():
            if "CORTEX" in i: 
                Cortex[i.replace("CORTEX","")]=list(range(j[0],j[1]+1,1))



    # write cell elsets
    for i, j in elset_cells.items():
        fo.write("*elset,  elset=Cell%s\n" % int(i))
        for x, a in enumerate(j):
            if a == j[-1]:
                if x % 16 == 15:
                    fo.write("\n")
                fo.write(str(a)+"\n")
            elif x % 16 == 15:
                fo.write(str(a)+"\n")
            else:
                fo.write(str(a)+",")
        fo.write("*Surface, type=ELEMENT, name=Cell"+str(i)+"_S3\n")
        fo.write("Cell"+str(i)+", S3\n")

    apicalS6 = []
    basalS6 = []
    lateralS6 = []

    if elemDim == 4:
        start_time = time.time()
        for tri in Cortex["APICAL"]:
            for el in elsetsOs["ElsetApical"]:
                for i in list(range(el*8-7,el*8+1,1)):
                    ca = intersection([elements[i][1],elements[i][2],elements[i][3],elements[i][5],elements[i][9],elements[i][8]], elemsS3RT[tri])
                    if len(ca)==3:
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
                        break
                    cb = intersection([elements[i][0],elements[i][1],elements[i][3],elements[i][4],elements[i][8],elements[i][7]], elemsS3RT[tri])
                    if len(cb)==3:
                        apicalS6.append(
                            [
                                int(elements[i][0]),
                                int(elements[i][1]),
                                int(elements[i][3]),
                                int(elements[i][4]),
                                int(elements[i][8]),
                                int(elements[i][7]),
                            ]
                        )
                        break
                    cc = intersection([elements[i][0],elements[i][1],elements[i][2],elements[i][4],elements[i][5],elements[i][6]], elemsS3RT[tri])
                    if len(cc)==3:
                        apicalS6.append(
                            [
                                int(elements[i][0]),
                                int(elements[i][1]),
                                int(elements[i][2]),
                                int(elements[i][4]),
                                int(elements[i][5]),
                                int(elements[i][6]),
                            ]
                        )
                        break
                    cd = intersection([elements[i][2],elements[i][0],elements[i][3],elements[i][6],elements[i][7],elements[i][9]], elemsS3RT[tri])
                    if len(cd)==3:
                        apicalS6.append(
                            [
                                int(elements[i][2]),
                                int(elements[i][0]),
                                int(elements[i][3]),
                                int(elements[i][6]),
                                int(elements[i][7]),
                                int(elements[i][9]),
                            ]
                        )
                        break
                else:
                    continue
                break
        print("--- %s seconds ---" % (time.time() - start_time))
        #exit()
        for tri in Cortex["BASAL"]:
            for el in elsetsOs["ElsetBasal"]:
                for i in list(range(el*8-7,el*8+1,1)):
                    ca = intersection([elements[i][1],elements[i][2],elements[i][3],elements[i][5],elements[i][9],elements[i][8]], elemsS3RT[tri])
                    cb = intersection([elements[i][0],elements[i][1],elements[i][3],elements[i][4],elements[i][8],elements[i][7]], elemsS3RT[tri])
                    cc = intersection([elements[i][0],elements[i][1],elements[i][2],elements[i][4],elements[i][5],elements[i][6]], elemsS3RT[tri])
                    cd = intersection([elements[i][2],elements[i][0],elements[i][3],elements[i][6],elements[i][7],elements[i][9]], elemsS3RT[tri])
                    if len(ca)==3:
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
                        break
                    elif len(cb)==3:
                        basalS6.append(
                            [
                                int(elements[i][0]),
                                int(elements[i][1]),
                                int(elements[i][3]),
                                int(elements[i][4]),
                                int(elements[i][8]),
                                int(elements[i][7]),
                            ]
                        )
                        break
                    elif len(cc)==3:
                        basalS6.append(
                            [
                                int(elements[i][0]),
                                int(elements[i][1]),
                                int(elements[i][2]),
                                int(elements[i][4]),
                                int(elements[i][5]),
                                int(elements[i][6]),
                            ]
                        )
                        break
                    elif len(cd)==3:
                        basalS6.append(
                            [
                                int(elements[i][2]),
                                int(elements[i][0]),
                                int(elements[i][3]),
                                int(elements[i][6]),
                                int(elements[i][7]),
                                int(elements[i][9]),
                            ]
                        )       
                        break
                else:
                    continue
                break
        print("--- %s seconds ---" % (time.time() - start_time))
        for tri in Cortex["LATERAL"]:
            for el in elsetsOs["ElsetLateral"]:
                for i in list(range(el*8-7,el*8+1,1)):
                    ca = intersection([elements[i][1],elements[i][2],elements[i][3],elements[i][5],elements[i][9],elements[i][8]], elemsS3RT[tri])
                    if len(ca)==3:
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
                        break
                    cb = intersection([elements[i][0],elements[i][1],elements[i][3],elements[i][4],elements[i][8],elements[i][7]], elemsS3RT[tri])
                    if len(cb)==3:
                        lateralS6.append(
                            [
                                int(elements[i][0]),
                                int(elements[i][1]),
                                int(elements[i][3]),
                                int(elements[i][4]),
                                int(elements[i][8]),
                                int(elements[i][7]),
                            ]
                        )
                        break
                    cc = intersection([elements[i][0],elements[i][1],elements[i][2],elements[i][4],elements[i][5],elements[i][6]], elemsS3RT[tri])
                    if len(cc)==3:
                        lateralS6.append(
                            [
                                int(elements[i][0]),
                                int(elements[i][1]),
                                int(elements[i][2]),
                                int(elements[i][4]),
                                int(elements[i][5]),
                                int(elements[i][6]),
                            ]
                        )
                        break
                    cd = intersection([elements[i][2],elements[i][0],elements[i][3],elements[i][6],elements[i][7],elements[i][9]], elemsS3RT[tri])
                    if len(cd)==3:
                        lateralS6.append(
                            [
                                int(elements[i][2]),
                                int(elements[i][0]),
                                int(elements[i][3]),
                                int(elements[i][6]),
                                int(elements[i][7]),
                                int(elements[i][9]),
                            ]
                        )  
                        break
                else:
                    continue
                break
        print("--- %s seconds ---" % (time.time() - start_time))
    else:
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
    fo.write(ss)

    MaxElAp6 = MaxEl + len(apicalS6)
    MaxElBa6 = MaxElAp6 + len(basalS6)
    MaxElS6 = MaxElBa6 + len(lateralS6)

    for i, j in enumerate(apicalS6):
        fo.write(
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
        fo.write(
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
        fo.write(
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
    fo.write(ss)
    MaxelT2A = int(MaxElS6)
    for i, j in enumerate(apicalT2):
        MaxelT2A += 1
        fo.write("%d, %d, %d, %d\n" % (int(MaxelT2A), int(j[0]), int(j[1]), int(j[2])))
    MaxelT2B = int(MaxelT2A)
    for i, j in enumerate(basalT2):
        MaxelT2B += 1
        fo.write("%d, %d, %d, %d\n" % (int(MaxelT2B), int(j[0]), int(j[1]), int(j[2])))
    MaxelT2L = int(MaxelT2B)
    for i, j in enumerate(lateralT2):
        MaxelT2L += 1
        fo.write("%d, %d, %d, %d\n" % (int(MaxelT2L), int(j[0]), int(j[1]), int(j[2])))
    MaxelT2 = int(MaxelT2L)

    fo.write("*nset, nset= CNapical\n")
    for x, a in enumerate(cnodeApical):
        if x + 1 == len(cnodeApical):
            if x % 16 == 15:
                fo.write("\n")
            fo.write(str(a) + "\n")
        elif x % 16 == 15:
            fo.write(str(a) + "\n")
        else:
            fo.write(str(a) + ",")
    fo.write("*nset, nset= CNbasal\n")
    for x, a in enumerate(cnodeBasal):
        if x + 1 == len(cnodeBasal):
            if x % 16 == 15:
                fo.write("\n")
            fo.write(str(a) + "\n")
        elif x % 16 == 15:
            fo.write(str(a) + "\n")
        else:
            fo.write(str(a) + ",")
    fo.write("*nset, nset= CNlateral\n")
    for x, a in enumerate(cnodeLateral):
        if x + 1 == len(cnodeLateral):
            if x % 16 == 15:
                fo.write("\n")
            fo.write(str(a) + "\n")
        elif x % 16 == 15:
            fo.write(str(a) + "\n")
        else:
            fo.write(str(a) + ",")

    ss = "*elset,  elset=ElsetApical\n"
    fo.write(ss)
    for x in range(len(apicalEl)):
        if x + 1 == len(apicalEl):
            fo.write(str(apicalEl[x]) + "\n")
        elif x % 16 == 15:
            fo.write(str(apicalEl[x]) + "\n")
        else:
            fo.write(str(apicalEl[x]) + ",")
    ss = "*elset,  elset=ElsetBasal\n"
    fo.write(ss)
    for x in range(len(basalEl)):
        if x + 1 == len(basalEl):
            fo.write(str(basalEl[x]) + "\n")
        elif x % 16 == 15:
            fo.write(str(basalEl[x]) + "\n")
        else:
            fo.write(str(basalEl[x]) + ",")
    ss = "*elset,  elset=ElsetLateral\n"
    fo.write(ss)
    for x in range(len(lateralEl)):
        if x + 1 == len(lateralEl):
            fo.write(str(lateralEl[x]) + "\n")
        elif x % 16 == 15:
            fo.write(str(lateralEl[x]) + "\n")
        else:
            fo.write(str(lateralEl[x]) + ",")
    alllEl = []
    alllEl = apicalEl + basalEl + lateralEl
    ss = "*elset,  elset=ElsetAll, generate\n"
    fo.write(ss)
    ss = str(list(elements.keys())[0]) + ", " + str(list(elements.keys())[-1]) + ", 1\n"
    fo.write(ss)

    ss = "*elset,  elset=AllEl, generate\n"
    fo.write(ss)
    ss = str(list(elements.keys())[0]) + ", " + str(list(elements.keys())[-1]) + ", 1\n"
    fo.write(ss)
    ss = "*elset,  elset=surfApical, generate\n"
    fo.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fo.write(ss)
    ss = "*elset,  elset=surfBasal, generate\n"
    fo.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fo.write(ss)
    ss = "*elset,  elset=surfLateral, generate\n"
    fo.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fo.write(ss)

    ss = "*elset,  elset=CORTEXAPICAL, generate\n"
    fo.write(ss)
    ss = str(MaxElS6 + 1) + ", " + str(MaxelT2A) + ", 1\n"
    fo.write(ss)
    ss = "*elset,  elset=CORTEXBASAL, generate\n"
    fo.write(ss)
    ss = str(MaxelT2A + 1) + ", " + str(MaxelT2B) + ", 1\n"
    fo.write(ss)
    ss = "*elset,  elset=CORTEXLATERAL, generate\n"
    fo.write(ss)
    ss = str(MaxelT2B + 1) + ", " + str(MaxelT2L) + ", 1\n"
    fo.write(ss)

    ss = "*Elset, elset=APICALSPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fo.write(ss)
    ss = "*Elset, elset=BASALSPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fo.write(ss)
    ss = "*elset,  elset=LATERALSPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fo.write(ss)

    ss = "*Elset, elset=_CPAPICAL_SPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElS6 + 1) + ", " + str(MaxelT2A) + ", 1\n"
    fo.write(ss)
    ss = "*Elset, elset=_CPBASAL_SPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxelT2A + 1) + ", " + str(MaxelT2B) + ", 1\n"
    fo.write(ss)
    ss = "*Elset, elset=_CPLATERAL_SPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxelT2B + 1) + ", " + str(MaxelT2L) + ", 1\n"
    fo.write(ss)

    ss = "*Elset, elset=_SPAPICAL_SPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fo.write(ss)
    ss = "*Elset, elset=_SNAPICAL_SNEG, internal, generate\n"
    fo.write(ss)
    ss = str(MaxEl + 1) + ", " + str(MaxElAp6) + ", 1\n"
    fo.write(ss)

    ss = "*Elset, elset=_SPBASAL_SPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fo.write(ss)
    ss = "*Elset, elset=_SNBASAL_SNEG, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElAp6 + 1) + ", " + str(MaxElBa6) + ", 1\n"
    fo.write(ss)

    ss = "*Elset, elset=_SPLATERAL_SPOS, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fo.write(ss)
    ss = "*Elset, elset=_SNLATERAL_SNEG, internal, generate\n"
    fo.write(ss)
    ss = str(MaxElBa6 + 1) + ", " + str(MaxElS6) + ", 1\n"
    fo.write(ss)

    ss = "*Surface, type=ELEMENT, name=CPAPICAL\n"
    fo.write(ss)
    ss = "_CPAPICAL_SPOS, SPOS\n"
    fo.write(ss)
    ss = "*Surface, type=ELEMENT, name=CPBASAL\n"
    fo.write(ss)
    ss = "_CPBASAL_SPOS, SPOS\n"
    fo.write(ss)
    ss = "*Surface, type=ELEMENT, name=CPLATERAL\n"
    fo.write(ss)
    ss = "_CPLATERAL_SPOS, SPOS\n"
    fo.write(ss)

    ss = "*Surface, type=ELEMENT, name=SPAPICAL\n"
    fo.write(ss)
    ss = "_SPAPICAL_SPOS, SPOS\n"
    fo.write(ss)
    ss = "*Surface, type=ELEMENT, name=SNAPICAL\n"
    fo.write(ss)
    ss = "_SNAPICAL_SNEG, SNEG\n"
    fo.write(ss)    

    ss = "*Surface, type=ELEMENT, name=SPBASAL\n"
    fo.write(ss)
    ss = "_SPBASAL_SPOS, SPOS\n"
    fo.write(ss)
    ss = "*Surface, type=ELEMENT, name=SNBASAL\n"
    fo.write(ss)
    ss = "_SNBASAL_SNEG, SNEG\n"
    fo.write(ss)

    ss = "*Surface, type=ELEMENT, name=SPLATERAL\n"
    fo.write(ss)
    ss = "_SPLATERAL_SPOS, SPOS\n"
    fo.write(ss)
    ss = "*Surface, type=ELEMENT, name=SNLATERAL\n"
    fo.write(ss)
    ss = "_SNLATERAL_SNEG, SNEG\n"
    fo.write(ss)

    ss = "*Shell Section, elset=CORTEXAPICAL, material=cortexa\n"
    fo.write(ss)
    ss = "1., 5\n"
    fo.write(ss)
    ss = "*Shell Section, elset=CORTEXBASAL, material=cortexb\n"
    fo.write(ss)
    ss = "1., 5\n"
    fo.write(ss)
    ss = "*Shell Section, elset=CORTEXLATERAL, material=cortexl\n"
    fo.write(ss)
    ss = "1., 5\n"
    fo.write(ss)
    ss = "*Solid Section, elset=ELSETALL, material=cytoplasm1\n"
    fo.write(ss)
    ss = ",\n"
    fo.write(ss)
    ss = "*Shell Section, elset=SURFAPICAL, material=membranea\n"
    fo.write(ss)
    ss = "1., 5\n"
    fo.write(ss)
    ss = "*Shell Section, elset=SURFBASAL, material=membraneb\n"
    fo.write(ss)
    ss = "1., 5\n"
    fo.write(ss)
    ss = "*Shell Section, elset=SURFLATERAL, material=membranel\n"
    fo.write(ss)
    ss = "1., 5\n"



    fo.write(ss)
    fo.write("*End Part\n")
    fo.write("**\n")
    fi.close()
    fo.close()
    return oldMaxNode, MaxNode, cnodeApical, cnodeBasal, cnodeLateral 



def save_inp(inp_file, mono,elemtosup=None, add_spaces=None, coords=None):
    """
    Write tetrahedra in inp format for Finite Element Modelisation in Abaqus
    """
    import pprint
    #pprint.pprint(mono.cell_df)
    
    if coords is None:
        coords = list("xyz")
    all_nodes = pd.concat(
        [mono.cell_df[coords], mono.face_df[coords], mono.vert_df[coords]],
        axis=0,
        ignore_index=True,
    )
    tetrahedra = mono.edge_df[["cell", "face", "srce", "trgt"]].copy()
    tetrahedra["face"] += mono.Nc
    tetrahedra[["srce", "trgt"]] += mono.Nc + mono.Nf
    #pprint.pprint(tetrahedra)
   
    tetrahedra.index += 1
    tetrahedra += 1
    if elemtosup!=None:
        tetrahedra.drop(elemtosup, axis=0, inplace=True)


    all_nodes.index += 1

    array = range(1,len(mono.cell_df)+1)
   
    side = mono.face_df['segment']
    nsetside= side.to_frame().groupby('segment').groups

    apical = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[0]])].copy()
    basal = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[1]])].copy()
    lateral = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[2]])].copy()
    #pprint.pprint(tetrahedra)
    elcells = tetrahedra['cell'].loc[tetrahedra['cell'].isin(array)].to_frame().copy()

    elcells = elcells.groupby('cell').groups
    #pprint.pprint(elcells)

    with open(inp_file, "w+") as inph:
        inph.write(f"*part, name={mono.identifier}\n")
        inph.write("*node\n")
    with open(inp_file, "a+") as inph:
        all_nodes.to_csv(inph, header=False, index=True, sep=",", line_terminator='\n')
    with open(inp_file, "a+") as inph:
        inph.write("*element, type=C3D4\n")
    with open(inp_file, "a+") as inph:
        tetrahedra.to_csv(inph, header=False, index=True, sep=",", line_terminator='\n')
    """
    for i in range(len(elcells)):
        #print("*elset, elset=Cell"+str(i+1))
        with open(inp_file, "a+") as inph:
            inph.write("*elset, elset=Cell"+str(i+1)+"\n")
            for x, a in enumerate(elcells[i+1]):
                if a == elcells[i+1][-1]:
                    if x % 16 == 15:
                        inph.write("\n")
                    inph.write(str(a)+"\n")
                elif x % 16 == 15:
                    inph.write(str(a)+"\n")
                else:
                    inph.write(str(a)+",")
    """
    for i in range(len(elcells)):
        with open(inp_file, "a+") as inph:
            inph.write("*elset, elset=Cell"+str(i+1)+"\n")
            for x, a in enumerate(elcells[i+1]):
                if a == elcells[i+1][-1]:
                    if x % 16 == 15:
                        inph.write("\n")
                    inph.write(str(a)+"\n")
                elif x % 16 == 15:
                    inph.write(str(a)+"\n")
                else:
                    inph.write(str(a)+",")
            inph.write("*Surface, type=ELEMENT, name=Cell"+str(i+1)+"_S3\n")
            inph.write("Cell"+str(i+1)+", S3\n")
    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFapicalN\n")
        #print (elcells[i+1])
        for x, a in enumerate(apical.index):
            #print(x)
            #print (len(apical))
            if x+1 == len(apical.index):
            #print (i+1)
                #if x % 16 == 15:
                    #print(x)
                    #inph.write("\n")
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFbasalN\n")
        #print (elcells[i+1])
        for x, a in enumerate(basal.index):
            #print(x)
            #print (len(basal))
            if x+1 == len(basal.index):
                #print (i+1)
                #if x % 16 == 15:
                    #print(x)
                    #inph.write("\n")
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFlateralN\n")
        #print (lateral.index)
        for x, a in enumerate(lateral.index):
            
            if x+1 == len(lateral.index):
                #print (i+1)
                #if x % 16 == 15:
                    #print(x)
                    #inph.write("\n")
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")
    #if add_spaces:
      # If you need the space after the coma in the file
    #    subprocess.run(["sed", "-i", "-e", "s/,/, /g", inp_file])   


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

def save_inp_job(inp_file, mono, add_spaces=None, coords=None):
    """
    Write tetrahedra in inp format for Finite Element Modelisation in Abaqus
    """
    if coords is None:
        coords = list("xyz")
    all_nodes = pd.concat(
        [mono.cell_df[coords], mono.face_df[coords], mono.vert_df[coords]],
        axis=0,
        ignore_index=True,
    )
    tetrahedra = mono.edge_df[["cell", "face", "srce", "trgt"]].copy()
    tetrahedra["face"] += mono.Nc
    tetrahedra[["srce", "trgt"]] += mono.Nc + mono.Nf

    tetrahedra.index += 1
    tetrahedra += 1
    all_nodes.index += 1

    array = range(len(mono.cell_df))
    side = mono.face_df['segment']
    nsetside= side.to_frame().groupby('segment').groups

    apical = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[0]])].copy()
    basal = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[1]])].copy()
    lateral = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[2]])].copy()
 
    elcells = tetrahedra['cell'].loc[tetrahedra['cell'].isin(array)].to_frame().copy()
    elcells = elcells.groupby('cell').groups

    with open(inp_file, "w+") as inph:
        inph.write(f"*Heading\n")
        inph.write(f"** Job name: Job-5 Model name: test\n")
        inph.write(f"** Generated by: Abaqus/CAE 2018\n")
        inph.write(f"*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        inph.write(f"**\n")
        inph.write(f"** PARTS\n")
        inph.write(f"*part, name={mono.identifier}\n")
        inph.write(f"**\n")
    with open(inp_file, "a+") as inph:
        all_nodes.to_csv(inph, header=False, index=True, sep=",", line_terminator='\n')
    with open(inp_file, "a+") as inph:
        inph.write("*element, type=C3D4\n")
    with open(inp_file, "a+") as inph:
        tetrahedra.to_csv(inph, header=False, index=True, sep=",", line_terminator='\n')
    for i in range(len(elcells)):
        #print (i)
        with open(inp_file, "a+") as inph:
            inph.write("*elset, elset=Cell"+str(i+1)+"\n")
            #print (elcells[i+1])
            for x, a in enumerate(elcells[i+1]):
                if a == elcells[i+1][-1]:
                    #print (i+1)
                    if x % 16 == 15:
                        inph.write("\n")
                    inph.write(str(a)+"\n")
                elif x % 16 == 15:
                    inph.write(str(a)+"\n")
                else:
                    inph.write(str(a)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFapical\n")
        #print (elcells[i+1])
        for x, a in enumerate(apical.index):
            #print(x)
            #print (len(apical))
            if x+1 == len(apical.index):
            #print (i+1)
                if x % 16 == 15:
                    inph.write("\n")
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFbasal\n")
        #print (elcells[i+1])
        for x, a in enumerate(basal.index):
            #print(x)
            #print (len(basal))
            if x+1 == len(basal.index):
                #print (i+1)
                if x % 16 == 15:
                    inph.write("\n")
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    with open(inp_file, "a+") as inph:
        inph.write("*nset, nset= CFlateral\n")
        #print (elcells[i+1])
        for x, a in enumerate(lateral.index):
            #print(x)
            #print (len(lateral))
            if x+1 == len(lateral.index):
                #print (i+1)
                if x % 16 == 15:
                    inph.write("\n")
                inph.write(str(a+mono.Nc+1)+"\n")
            elif x % 16 == 15:
                inph.write(str(a+mono.Nc+1)+"\n")
            else:
                inph.write(str(a+mono.Nc+1)+",")

    if add_spaces:
        # If you need the space after the coma in the file
        subprocess.run(["sed", "-i", "-e", "s/,/, /g", inp_file])   
        print(len(tetrahedra))

def write_inside_part_inp(inp_file, nodes, tetra, external_nodes):

    with open(inp_file, "a+") as f:

        f.write("*PART,  name=InsidePart\n")
        f.write("*NODE, nset=InsidePart\n")
        internal_volume = []

        for i, n in enumerate(nodes):
            f.write(
                "%d, %f, %f, %f\n" % (int(i + 1), float(n[0]), float(n[1]), float(n[2]))
            )
        if len(tetra[1]) == 4:
            f.write(f"*element, TYPE=C3D4\n")
            for k, x in enumerate(tetra):
                f.write(
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
            f.write(f"*element, TYPE=C3D10\n")
            for k, x in enumerate(tetra):
                f.write(
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

        f.write(f"*Elset, elset=all_InsidePart, generate\n")
        f.write("%d, %d, %d\n" % (int(1), int(k + 1), int(1)))

        # ajouter les surfaces -->
        f.write(f"*Elset, elset=_InsidePartContact_S1, internal\n")
        for x in range(len(internal_volume)):
            if x + 1 == len(internal_volume):
                f.write(str(internal_volume[x]) + "\n")
            elif x % 16 == 15:
                f.write(str(internal_volume[x]) + "\n")
            else:
                f.write(str(internal_volume[x]) + ",")
        f.write(f"*Surface, type=ELEMENT, name=InsidePartContact\n")
        f.write(f"_InsidePartContact_S1, S1\n")
        f.write(f"** Section: internal\n")
        f.write(f"*Solid Section, elset=all_InsidePart, material=inside\n")
        f.write(f",\n")
        f.write("*End Part\n")
        f.write("**\n")
        print("internal part added to inp file %s " % inp_file)

def write_properties_inp(inp_file, load, center_node, opposite_node, GlobalPartName,allEl, OutsidePart=False):
    with open(inp_file, "a+") as f:
        if method == "perfect":
            a = 0
            for i, j in load.items():
                if j > a:
                    a = j
                    zmax_node = i
            # calculating standard deviation using np.std
            std = np.std(list(load.values()))
            mean = np.percentile(list(load.values()), 85)
        f.write("**\n")
        f.write("** ASSEMBLY\n")
        f.write("**\n")
        f.write("*Assembly, name=Assembly\n")
        f.write("**\n")
        f.write("*Instance, name=%s-1, part=%s\n" % (GlobalPartName,GlobalPartName))
        f.write("*End Instance\n")
        f.write("**\n")
        f.write("*Instance, name=InsidePart-1, part=InsidePart\n")
        f.write("*End Instance\n")
        f.write("**\n")
        if OutsidePart == True:
            f.write("*Instance, name=OutsidePart-1, part=OutsidePart\n")
            f.write("*End Instance\n")
            f.write("**\n")
        f.write("*Nset, nset=InsidePart_center, instance=InsidePart-1\n")
        f.write("%d, \n" % int(center_node))
        f.write("*Nset, nset=opposite, instance=%s-1\n" % GlobalPartName)
        f.write("%d, \n" % int(opposite_node))        
        if method == "perfect":
            f.write("*Nset, nset=maxload, instance=%s-1\n" % GlobalPartName)
            f.write("%d, \n" % int(zmax_node))
        ss = "*Surface, type=NODE, name=%s-CNAPICAL, internal\n" % GlobalPartName
        f.write(ss)
        ss = "%s-1.CNAPICAL, 1.\n" % GlobalPartName
        f.write(ss)
        ss = "*Surface, type=NODE, name=%s-CNBASAL, internal\n" % GlobalPartName
        f.write(ss)
        ss = "%s-1.CNBASAL, 1.\n" % GlobalPartName
        f.write(ss)

        f.write("*Elset, elset=PRESSURE, internal, generate, instance=%s-1\n" % GlobalPartName)
        f.write(str(allEl[0])+", "+str(allEl[1])+", "+str(allEl[2])+"\n" )
        f.write("*Surface, type=ELEMENT, name=AREAPRESSURE\n")
        f.write("PRESSURE, S1\n")
        f.write("PRESSURE, S2\n")
        f.write("PRESSURE, S4\n")
        f.write("PRESSURE, S3\n")
        
        f.write("** Constraint: CP-1-%s-1-InsidePart-1\n" % GlobalPartName)
        f.write("*Tie, name=CP-1-%s-1-InsidePart-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName)
        f.write("InsidePart-1.InsidePartContact, %s-1.SPAPICAL\n" % GlobalPartName)
        if OutsidePart == True:
            f.write("** Constraint: CP-1-%s-1-OutsidePart-1\n" % GlobalPartName)
            f.write(
                "*Tie, name=CP-1-%s-1-OutsidePart-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
            )
            f.write("OutsidePart-1.OutsidePartContact, %s-1.SPBASAL\n" % GlobalPartName)
        f.write("** Constraint: CN-1-%s-CENTERS-APICAL-1\n" % GlobalPartName)
        f.write(
            "*Tie, name=CN-1-%s-APICAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
        )
        f.write("%s-1.SPAPICAL, %s-1.CPAPICAL\n" % (GlobalPartName, GlobalPartName))

        f.write("** Constraint: CN-1-%s-CENTERS-BASAL-1\n" % GlobalPartName)
        f.write(
            "*Tie, name=CN-1-%s-BASAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
        )
        f.write(" %s-1.SPBASAL, %s-1.CPBASAL\n" % (GlobalPartName, GlobalPartName))
        
        f.write("** Constraint: CN-1-%s-CENTERS-LATERAL-1\n" % GlobalPartName)
        f.write("*Tie, name=CN-1-%s-LATERAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName)
        f.write(" %s-1.SPLATERAL, %s-1.CPLATERAL\n" % (GlobalPartName, GlobalPartName))


        f.write("*End Assembly\n")
        f.write("**\n")
        f.write("** MATERIALS\n")
        f.write("**\n")
        f.write("*Material, name=cytoplasm1\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.002, 0.44\n")
        f.write("*Material, name=cytoplasm2\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.002, 0.44\n")
        f.write("*Material, name=membranea\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.009, 0.33\n")
        f.write("*Material, name=membraneb\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.0009, 0.33\n")
        f.write("*Material, name=membranel\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("4.5e-09, 0.42\n")
        f.write("*Material, name=cortexa\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.003, 0.15\n")
        f.write("*Expansion\n")
        f.write("0.005\n") 
        f.write("*Specific Heat\n")
        f.write("400\n") 
        f.write("*Conductivity\n")
        f.write("0\n")
        f.write("*Material, name=cortexb\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.01, 0.3\n")
        f.write("*Expansion\n")
        f.write("0.0006\n") 
        f.write("*Specific Heat\n")
        f.write("400\n") 
        f.write("*Conductivity\n")
        f.write("0\n")
        f.write("*Material, name=cortexl\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.0001, 0.3\n")
        f.write("*Expansion\n")
        f.write("0.000006\n") 
        f.write("*Specific Heat\n")
        f.write("400\n") 
        f.write("*Conductivity\n")
        f.write("0\n")

        f.write("*Material, name=inside\n")
        f.write("*Density\n")
        f.write("1e-09,\n")
        f.write("*Elastic\n")
        f.write("1e-10, 0.32\n")
        f.write("** BOUNDARY CONDITIONS\n")
        f.write("** \n")
        f.write("** Name: BC-5 Type: Symmetry/Antisymmetry/Encastre\n")
        f.write("*Boundary\n")
        f.write("opposite, ENCASTRE\n")
        f.write("**--------------------------------------------------------------\n")
        f.write("** \n")
        f.write("** STEP: Step-1\n")
        f.write("** \n")
        f.write("*Step, name=Step-1, nlgeom=NO, inc=1000000\n")
        f.write("*Coupled Temperature-displacement, creep=none, steady state\n")
        f.write("0.01, 1.,\n")
        f.write("*BOUNDARY\n")
        f.write("maxload, 1,1\n")
        f.write("maxload, 2,2\n")
        f.write("maxload, 6,6\n")




        f.write("**\n")
        f.write("** LOADS\n")
        f.write("**\n")
        f.write("** Name: SURFFORCE-1   Type: Pressure\n")
        f.write("*Dsload\n")
        f.write("AREAPRESSURE, P, -0.01\n")

        """
        if LoadType == "pushArea": 
            f.write("**\n")
            f.write("** LOADS\n")
            f.write("**\n")
            f.write("** Name: Load-1   Type: Concentrated force\n")
            f.write("*Cload, follower\n")
            mean = np.percentile(list(load.values()), 85)
            for i, j in load.items():
                if j > mean:
                    f.write("%s."% GlobalPartName + str(i) + ", 3,-" + str(j) + "e-03\n" )
        else:
            for i, j in load.items():
                f.write("*BOUNDARY, type= DISP\n")
                f.write("%s-1."% GlobalPartName + str(i) + ", 1,1," + str(j[0]) + "\n" )
                f.write("*BOUNDARY, type= DISP\n")
                f.write("%s-1."% GlobalPartName + str(i) + ", 2,2," + str(j[1]) + "\n" )
                f.write("*BOUNDARY, type= DISP\n")
                f.write("%s-1."% GlobalPartName + str(i) + ", 3,3," + str(j[2]) + "\n" )
        """
        f.write("** \n")
        f.write("** OUTPUT REQUESTS\n")
        f.write("** \n")
        f.write("*Restart, write, number interval=1, time marks=NO\n")
        f.write("** \n")
        f.write("** FIELD OUTPUT: F-Output-1\n")
        f.write("** \n")
        f.write("*Output, field\n")
        f.write("*Node Output\n")
        f.write("CF, COORD, MOT, RF, TF, U, UR, UT\n")
        f.write("*element Output, directions=YES\n") 
        f.write("ELASE, ELSE, EVOL, IVOL, PRESS, TEMP\n")
        f.write("** \n")
        f.write("** HISTORY OUTPUT: H-Output-1\n")
        f.write("** \n")
        f.write("***Output, history, variable=PRESELECT\n")
        f.write("***FILE FORMAT, ASCII\n")
        f.write("***NODE FILE\n")        
        f.write("**COORD\n")
        f.write("**U\n")
        f.write("*End Step\n")
        print("properties added to inp file %s " % inp_file)

def write_individual_properties_inp(inp_file, mono, load, center_node, GlobalPartName, OutsidePart=False):
    a,b,c=read_elsets(inp_file)
    allEl=b["ElsetAll"]
    coords = list("xyz")
    all_nodes = pd.concat(
        [mono.cell_df[coords], mono.face_df[coords], mono.vert_df[coords]],
        axis=0,
        ignore_index=True,
    )
    tetrahedra = mono.edge_df[["cell", "face", "srce", "trgt"]].copy()
    tetrahedra["face"] += mono.Nc
    tetrahedra[["srce", "trgt"]] += mono.Nc + mono.Nf
    tetrahedra.index += 1
    tetrahedra += 1
    all_nodes.index += 1
    array = range(len(mono.cell_df))
    side = mono.face_df['segment']
    nsetside= side.to_frame().groupby('segment').groups
    apical = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[0]])].copy()
    basal = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[1]])].copy()
    lateral = mono.face_df.loc[mono.face_df.index.isin(nsetside[list(nsetside)[2]])].copy()
    elcells = tetrahedra['cell'].loc[tetrahedra['cell'].isin(array)].to_frame().copy()
    elcells = elcells.groupby('cell').groups


    with open(inp_file, "a+") as f:

        f.write("**\n")
        f.write("** ASSEMBLY\n")
        f.write("**\n")
        f.write("*Assembly, name=Assembly\n")
        f.write("**\n")
        f.write("*Instance, name=%s-1, part=%s\n" % (GlobalPartName,GlobalPartName))
        f.write("*End Instance\n")
        f.write("**\n")
        f.write("*Instance, name=InsidePart-1, part=InsidePart\n")
        f.write("*End Instance\n")
        f.write("**\n")
        """
        for i in range(len(elcells)):
            with open(inp_file, "a+") as inph:
                inph.write("*elset, elset=Cell"+str(i+1)+"\n")

                f.write("*Instance, name=Cell-%s, part=%s\n"% (i, GlobalPartName))
                f.write("*End Instance\n")
                f.write("**\n")
        """

        f.write("*Nset, nset=InsidePart_center, instance=InsidePart-1\n")
        f.write("%d, \n" % int(center_node))
        ss = "*Surface, type=NODE, name=%s-CNAPICAL, internal\n" % GlobalPartName
        f.write(ss)
        ss = "%s-1.CNAPICAL, 1.\n" % GlobalPartName
        f.write(ss)
        ss = "*Surface, type=NODE, name=%s-CNBASAL, internal\n" % GlobalPartName
        f.write(ss)
        ss = "%s-1.CNBASAL, 1.\n" % GlobalPartName
        f.write(ss)
        for i in range(mono.Nc):
            f.write("*Surface, type=ELEMENT, name=Cell"+str(i+1)+"_S3\n")
            f.write("ORGANOID-1.Cell"+str(i+1)+", S3\n")

        f.write("** Constraint: CP-1-%s-1-InsidePart-1\n" % GlobalPartName)
        f.write("*Tie, name=CP-1-%s-1-InsidePart-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName)
        f.write("InsidePart-1.InsidePartContact, %s-1.SPAPICAL\n" % GlobalPartName)
        f.write("** Constraint: CN-1-%s-CENTERS-APICAL-1\n" % GlobalPartName)
        f.write(
            "*Tie, name=CN-1-%s-APICAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
        )
        f.write("%s-1.SPAPICAL, %s-1.CPAPICAL\n" % (GlobalPartName, GlobalPartName))

        f.write("** Constraint: CN-1-%s-CENTERS-BASAL-1\n" % GlobalPartName)
        f.write(
            "*Tie, name=CN-1-%s-BASAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName
        )
        f.write(" %s-1.SPBASAL, %s-1.CPBASAL\n" % (GlobalPartName, GlobalPartName))
        f.write("** Constraint: CN-1-%s-CENTERS-LATERAL-1\n" % GlobalPartName)
        f.write("*Tie, name=CN-1-%s-LATERAL-1, adjust=yes, type= SURFACE TO SURFACE\n" % GlobalPartName)
        f.write(" %s-1.SPLATERAL, %s-1.CPLATERAL\n" % (GlobalPartName, GlobalPartName))
        f.write("*End Assembly\n")
        f.write("**\n")
        f.write("** MATERIALS\n")
        f.write("**\n")
        f.write("*Material, name=cytoplasm1\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.002, 0.44\n")
        f.write("*Material, name=cytoplasm2\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.002, 0.44\n")
        f.write("*Material, name=membranea\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.009, 0.33\n")
        f.write("*Material, name=membraneb\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.0009, 0.33\n")
        f.write("*Material, name=membranel\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("4.5e-09, 0.42\n")
        f.write("*Material, name=cortexa\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.003, 0.15\n")
        f.write("*Expansion\n")
        f.write("0.005\n") 
        f.write("*Specific Heat\n")
        f.write("400\n") 
        f.write("*Conductivity\n")
        f.write("0\n")
        f.write("*Material, name=cortexb\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.01, 0.3\n")
        f.write("*Expansion\n")
        f.write("0.0006\n") 
        f.write("*Specific Heat\n")
        f.write("400\n") 
        f.write("*Conductivity\n")
        f.write("0\n")
        f.write("*Material, name=cortexl\n")
        f.write("*Density\n")
        f.write("7.8e-09,\n")
        f.write("*Elastic\n")
        f.write("0.0001, 0.3\n")
        f.write("*Expansion\n")
        f.write("0.000006\n") 
        f.write("*Specific Heat\n")
        f.write("400\n") 
        f.write("*Conductivity\n")
        f.write("0\n")
        f.write("*Material, name=inside\n")
        f.write("*Density\n")
        f.write("1e-09,\n")
        f.write("*Elastic\n")
        f.write("1e-10, 0.32\n")



        f.write("**--------------------------------------------------------------\n")
        f.write("** STEP: Step-1\n")
        f.write("** \n")
        f.write("*Step, name=Step-1, nlgeom=NO, inc=1000000\n")
        f.write("*Coupled Temperature-displacement, creep=none, steady state\n")
        f.write("0.01, 1.,\n")
        f.write("** BOUNDARY CONDITIONS\n")
        f.write("** \n")
        f.write("** Name: BC-1 Type: Symmetry/Antisymmetry/Encastre\n")
        f.write("* BOUNDARY\n")
        f.write("INSIDEPART_CENTER, ENCASTRE\n")
        f.write("**\n")
        f.write("** LOADS\n")
        f.write("**\n")
        f.write("** Name: Load-apical   Type: Pressure\n")
        f.write("*Dsload\n")
        f.write("ORGANOID-1.CPAPICAL, P, -0.0005\n")
        f.write("** Name: Load-basal   Type: Pressure\n")
        f.write("*Dsload\n")
        f.write("ORGANOID-1.CPBASAL, P, 0.0005\n")


        for i in range(mono.Nc):
            f.write("** Name: Cell_"+str(i+1)+"  Type: Pressure\n")
            f.write("*Dsload\n")
            f.write("Cell"+str(i+1)+"_S3 , P, -1E-10\n")
        f.write("** \n")
        """
        f.write("** INTERACTION PROPERTIES\n")
        f.write("**\n")
        f.write("*Surface Interaction, name=IntProp-1\n")
        f.write("1.,\n")
        f.write("*Surface Behavior, pressure-overclosure=HARD\n")
        f.write("**\n")
        f.write("** INTERACTIONS\n")
        f.write("** \n")
        f.write("** Interaction: Int-2\n")
        f.write("*Contact\n")
        f.write("*Contact Inclusions\n")
        for i in range(mono.Nc):
            f.write(", ORGANOID-1.CELL"+str(i+1)+"_S3\n")
        f.write("*Contact Property Assignment\n")
        f.write(",  , IntProp-1\n")
        """
        f.write("** OUTPUT REQUESTS\n")
        f.write("** \n")
        f.write("*Restart, write, number interval=1, time marks=NO\n")
        f.write("** \n")
        f.write("** FIELD OUTPUT: F-Output-1\n")
        f.write("** \n")
        f.write("*Output, field\n")
        f.write("*Node Output\n")
        f.write("CF, COORD, MOT, RF, TF, U, UR, UT\n")
        f.write("*element Output, directions=YES\n") 
        f.write("ELASE, ELSE, EVOL, IVOL, PRESS, TEMP\n")
        f.write("** \n")
        f.write("** HISTORY OUTPUT: H-Output-1\n")
        f.write("** \n")
        f.write("***Output, history, variable=PRESELECT\n")
        f.write("***FILE FORMAT, ASCII\n")
        f.write("***NODE FILE\n")        
        f.write("**COORD\n")
        f.write("**U\n")
        f.write("*End Step\n")


        f.write("**--------------------------------------------------------------\n")
        f.write("**--------------------------------------------------------------\n")
        f.write("** STEP: Step-2\n")
        f.write("** \n")
        f.write("*Step, name=Step-2, nlgeom=NO, inc=1000000\n")
        f.write("*Coupled Temperature-displacement, creep=none, steady state\n")
        f.write("0.01, 1.,\n")
        f.write("** BOUNDARY CONDITIONS\n")
        f.write("** \n")
        f.write("** Name: BC-1 Type: Symmetry/Antisymmetry/Encastre\n")
        f.write("* BOUNDARY\n")
        f.write("INSIDEPART_CENTER, ENCASTRE\n")
        f.write("**\n")
        f.write("** LOADS\n")
        f.write("**\n")
        f.write("** Name: Load-apical   Type: Pressure\n")
        f.write("*Dsload\n")
        f.write("ORGANOID-1.CPAPICAL, P, -5E-09\n")
        f.write("** Name: Load-basal   Type: Pressure\n")
        f.write("*Dsload\n")
        f.write("ORGANOID-1.CPBASAL, P, 5E-09\n")
        f.write("** \n")
        f.write("** OUTPUT REQUESTS\n")
        f.write("** \n")
        f.write("*Restart, write, number interval=1, time marks=NO\n")
        f.write("** \n")
        f.write("** FIELD OUTPUT: F-Output-1\n")
        f.write("** \n")
        f.write("*Output, field\n")
        f.write("*Node Output\n")
        f.write("CF, COORD, MOT, RF, TF, U, UR, UT\n")
        f.write("*element Output, directions=YES\n") 
        f.write("ELASE, ELSE, EVOL, IVOL, PRESS, TEMP\n")
        f.write("** \n")
        f.write("** HISTORY OUTPUT: H-Output-1\n")
        f.write("** \n")
        f.write("***Output, history, variable=PRESELECT\n")
        f.write("***FILE FORMAT, ASCII\n")
        f.write("***NODE FILE\n")        
        f.write("**COORD\n")
        f.write("**U\n")
        f.write("*End Step\n")

        print("properties added to inp file %s " % inp_file)