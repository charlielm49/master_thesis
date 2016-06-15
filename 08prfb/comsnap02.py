#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""""

#XXXXXXXXX-BEGIN: REQUIREMENTS & SPECS

#from __future__ import print_function # to run in python2 & 3
import os
import sys
import struct # Packing and unpacking of heterogeneous binary data
import array
import math
#import networkx as nx  # graph library
#import snap # graph library # actually it is not needed; only AGM
import numpy as np
import matplotlib.pyplot as plt
import community  # for community detection to form hypergraphs
import operator  # for helping sort a dictionary
from scipy.sparse import *

#opt
import pprint

#XXXXXXXXX-END: REQUIREMENTS & SPECS


# community Graph / Community creation
def detectComm(ingraph, layout):
    parts = community.best_partition(ingraph)
    values = [parts.get(node) for node in ingraph.nodes()]

    plt.axis("off")
    nx.draw_networkx(ingraph, pos=layout, cmap=plt.get_cmap("jet"),
                     node_color=values,
                     node_size=35,
                     with_labels=False
                     )
    plt.show()


# construct incidence matrix (and check that there are no isolated 
# vertices in a community) get list of unique communities
# IN: list of values (str, int) - not doubles
# OUT: list of unique (de-duplicated) values 
def getUniqueComm(vertexCommunityL):
    commUniqList = []
    for comm in vertexCommunityL:
        if comm not in commUniqList:
            commUniqList.append(comm)

    lenUL = len(commUniqList)
    #print commUniqList
    return commUniqList, lenUL
# This code obtains list of all unique nodes in ALL the file of edges
# It scans both lists, source and destination for node ids
# RESTRICT: doesn't work well if origin-source are the same
# IN: list of 2-values (edges - [1, 3]). str, int but not doubles
# OUT: list of unique (de-duplicated) values 
def getUniqueCommLL(vertexCommunityLL):
    commUniqList = []
    for comm in vertexCommunityLL:
        if comm[0] not in commUniqList:
            commUniqList.append(comm[0])
        if comm[1] not in commUniqList:
            commUniqList.append(comm[1])
            
    lenUL = len(commUniqList)
    #print commUniqList
    return commUniqList, lenUL

# Create format requisito for csr from in the Klein format
# THIS RETURNS INCIDENCE MATRIX OF HYPERGRAPH (H)
# Iterate over the list and substitute
# IN01: (indexList: list of 2-tuples node-edge)
# [(0, 0), (1, 0), (2, 0), ... (9, 0), ...]
# IN02: number of edges (hyperedges)
# OUT: (matrix in CSR format)
# for this format, we need:
#  rows = []
#  cols = []
#  data = []
#  of pairs {(a, b), (c, d), ...}
# <- first coord of each pair (a, ...) (c, )goes in rows of csr
# <- second coord (..., b) ( , d)goes in cols of csr
def createPrevCSR(indexList, numEdges):
    #global csrMat # OUTPUT
    print "\nEnter createPrevCSR for incidence matr"
    print "numEdges", numEdges

    matSize = len(indexList)
    print "matSize", matSize

    rows = []
    cols = []
    data = []
    for tupleItem in indexList:
        # look for item in indexTitles
        rows.append(tupleItem[0])
        cols.append(tupleItem[1])
        data.append(1)

    csrMat = csr_matrix( (data, (rows, cols)), shape = (matSize, numEdges) ).toarray()
    #ok print("csrMat created")
    #ok print (csrMat)
    #ok-butchangetoCSR return rows, cols, data
    return csrMat
# second version of CSR creation that takes (column index, value)
# because this comes from inverse vertex degree
# returns square matrix
def createPrevCSR2(invDegreePairs):
    #global csrMat # OUTPUT
    print "\nEnter createPrevCSR2 for inv Degree"
    matSize = len(invDegreePairs)
    print "matSize", matSize

    rows = []
    cols = []
    data = []
    for tupleItem in invDegreePairs:
        # look for item in indexTitles
        rows.append(tupleItem[0])
        cols.append(tupleItem[0])
        data.append(tupleItem[1])

    csrMat = csr_matrix( (data, (rows, cols)), shape = (matSize, matSize) ).toarray()
    #ok print("csrMat created")
    #ok print (csrMat)
    #ok-butchangetoCSR return rows, cols, data
    return csrMat

# IN01: list of 2-tuples node-edge
# OUT01: dictionary with node degree (how many edges a vertex belongs to)
# just counts number of edges each vertex belongs to
def getInvDvMatr(nodeEdgeL):
    # type: (list) -> dictionary
    nodeEdgeCount = {}

    for nodeEdge in nodeEdgeL:
        if nodeEdge[0] in nodeEdgeCount:
            #ok: print "found", nodeEdgeCount[nodeEdge[0]]
            nodeEdgeCount.update({nodeEdge[0]: nodeEdgeCount[nodeEdge[0]] + 1})
        else:
            nodeEdgeCount.update({nodeEdge[0] : 1})

    print "node-edge count:", nodeEdgeCount.items()[:10]
    # check sum for one-class vertices

    # get inverse of each degree
    invNodeEdgeCount = []
    for node in nodeEdgeCount.keys():
        #print nodeCount
        invNodeEdgeCount.append((node, 1./nodeEdgeCount[node]))

    print "Inverse of vertex degree:"
    print "first", invNodeEdgeCount[:5]
    print "last", invNodeEdgeCount[-5:]

    #make sparse CSR matrix
    return createPrevCSR2(invNodeEdgeCount)


# Matrix = A
# initVector = v
def iteration(A, v):
    MAX_ITER = 100 # works ok for simple cases of PageRank
    MAX_ITER = 100000
    EPS = 1.0e-6 # works ok for simple cases of PageRank
    EPS = 1.0e-24
    for i in range(MAX_ITER):
        # guardamos el vector anterior
        vOld = v.copy()
        # calculamos el vector nuevo
        z = np.dot(A, v)
        # print("z", z)
        # z2 = mulMatVecCUDA(A, v)
        # print("z2", z2)
        # calculamos su magnitud
        zMag = math.sqrt(np.dot(z.T, z))
        # print("zMag", zMag)
        # dividimos vector entre la norma para normalizar
        v = z / zMag
        if np.dot(vOld.T, v) < 0.0:
            sign = -1.0
            v = -v
        else:
            sign = 1.0

        # Condición de terminación
        vecDif = vOld - v
        if math.sqrt(np.dot(vecDif.T, vecDif)) < EPS:
            break

    print("num iters = ", i)

    eigVal = sign * zMag
    # Another way to get eigVal
    # eigValAnother = np.dot(v.T, np.dot(A, v))
    # print("eigValAnother = ", eigValAnother)

    # eigVec = v
    return eigVal, v


def makeInitVec(listSize):
    # pdb.set_trace()
    matRows = listSize  # size of the matrix (assuming square)
    # 1. inicializamos la matriz con ceros
    newList = [1.0 for i in range(matRows)]
    # print newMat

    return newList



# XXXXXXXXXXXXXXX BEG: extra processing for SNAP XXXXXXXXXXXXXXXXXXXX



# entra: un archivo de texto con renglones con \n al final de cada renglón
# sale: una lista con todo el contenido del archivo en una sola lista
# Nota1: no se hace ningún procesamiento whatsoever, sólo se quitan los
#  \n y se hace split por espacios (con el fp.read() no se quitan los \n)
def leer_texto(fp):
    tmpFull = []
    # procesamiento de archivo de entrada completo
    for line in fp:
        tmp = line.split()
        #ok print tmp
        tmpFull.append([int(tmp[0]), int(tmp[1])])

    # regresa una lista con todas las palabras del texto en esa lista
    return tmpFull
# entra: un archivo de texto con renglones con \n al final de cada renglón
# IN: (tab-separated file - each line are the members of a group)
# 1912	2347	2543 ...
# 107	1888	1800 ...
# OUT: [(1912, 0), (2347, 0), (2543, 0), (2266, 0),...]
def leer_texto2(fp):
    tmpFull = []
    initialGroup = 0
    # procesamiento de archivo de entrada completo
    for line in fp:
        tmp = line.split()
        #ok print tmp
        tmp2 = [int(item) for item in tmp]
        for item2 in tmp2:
            tmpFull.append((item2, initialGroup))
        initialGroup += 1

    # regresa una lista con todas las palabras del texto en esa lista
    return tmpFull


# toma lista de "cosas" (lista ([key, val]), etc.)
# en este caso, lista de listas: [[1, 2], [3, 4], ...]
# guarda una "cosa" por linea
# Version para cosas de mas de 1 elemento (p.ej lista ([key, val]), etc.)
# outputs to tab-separated file,(input for BigCLAM):
# 1\t2
# 3\t4
def escribir_texto(listaL, fp):
    # type: (object, object) -> object
    separator = "\t"

    for item in listaL:
        item = separator.join(str(x) for x in item) + "\n"
        fp.write(item)


# obtain nodes and_edges_for graph
# 2016MAY26: I don't need to build the graph because:
# a) I need to run the agm algorthm in C++ which takes the edge list
# b) agm algorithm returns groups and grp membership of nodes
# I used it because I though the alg was in the interface (it isn't)
# DISCONTINUED
def obtainGraph():
    # listEdges is list of lists [[1, 2], [3, 4], ...]
    #ok foirst 10: listEdges = leer_texto(abrir_archivos())[:10]
    listEdges = leer_texto(abrir_archivos())
    gfb = snap.TNGraph.New() # directed
    #gfb = snap.TUNGraph.New() # undirected

    # get unique nodes to add them to graph
    uniqCommL, lenUL = getUniqueCommLL(listEdges)
    #ok: print uniqCommL, lenUL
    # add nodes to graph
    for node in uniqCommL:
        gfb.AddNode(node)

    # add edges to graph
    for edge in listEdges:
        gfb.AddEdge(edge[0], edge[1])

    '''
    # ok: test to see if graph was well created
    for node in gfb.Nodes():
        print "node: %d, outdegree %d, in-degree %d" % (node.GetId(), \
                                                        node.GetOutDeg(), \
                                                        node.GetInDeg(),)
    '''

# Process to run BigCLAM:
# 1) get list of edges (from main())
# 2) out list of edges for BC (tab-sep) (to:
#  outFilePointer01/facebook_combined_tab)
# 3) call command with parameters for BC:
#  a) numCommunities
#  b) inFileNameBC (data/facebook_combined_tab.txt)
#  c) outFileNameBC (defined in this file: data/communities.txt)
# ---IN01: (list of edges - space)
# listEdges is list of lists [[1, 2], [3, 4], ...]
# ---IN02: outFilePointer01: ptr to outfile to write tab-separated file 
# of edges
# This format (tab-separated is necessary for BIGCLAM)
# defined in abrir_archivos(): outFilePointer01 / 
# data/facebook_combined_tab.txt
# ---IN03 inFileNameBC is the name of the processed file with tab-separated 
# edges it is called inFile because it is the input for BigCLAM)
# defined in abrir_archivos(): data/facebook_combined_tab.txt
def commPartition(listEdges, outFilePointer01, inFileNameBC):
    pathToBC = '/ul/bin/agm-package/bigclam/' # for CIC server
    #pathToBC = './' # for MTA PC
    numCommunities = 10
    outFileNameBC = "data/communities.txt"  # out File name for BigCLAM
    outFileNameCy = "data/commCyto.txt"     # out File name for Cytoscape
    # save edge list in file
    escribir_texto(listEdges, outFilePointer01)
    # Run snap algorithm in C++ with previous file
    cmd = pathToBC \
            + "bigclam" \
            + " -c:" + str(numCommunities) \
            + " -i:" + inFileNameBC \
            + " -o:" + outFileNameBC
    print "Running:", cmd
    os.system(cmd)

    # outfile format data/communities.txt (each line is a community)
    # 1912\t  2347\t    2543 ...
    # 1684\t  2839\t    3363 ...
    # ...

    # read grps and grp membership of nodes from file and put in a dict
    # This is necessary because BC outputs to file and we need that info
    inFilePointer02 = open(outFileNameBC, "r")
    communitiesL = leer_texto2(inFilePointer02)
    print communitiesL[:10]
    print len(communitiesL)
    #ok full print: pprint.pprint communitiesL

    # communitiesL is a list of 2-tuples, having pairs: (node, community)
    # [(0, 0), (1, 0), (2, 0), (3, 0), ..., (4037, 6), (4038, 6)]
    # export commList to file for visualization with cytoscape
    outFilePointer02 = open(outFileNameCy, "w")
    escribir_texto(communitiesL, outFilePointer02)
    
    # close file opened here:
    inFilePointer02.close()
    
    return communitiesL


# Abre archivos entrada/salida / Regresa: apuntadores a archivos
# 2015FEB21: ya no se ponen parámetros en línea de comandos porque
# ya hay muchos archivos que leer - sólo el del archivo a procesar
# file_in01: archivo en raw text a procesar
# file_out01: archivo de salida ya procesado
def abrir_archivos():
    # print sys.argv # 0: nombre de programa 1: nombre de archivo
    try:
        filename_in01 = "data/facebook_combined.txt"  # este sí cambia
        fpIN01 = open(filename_in01, 'r')
        # outfile for SNAP C++ INPUT (output of read process)
        filename_OUT01 = "data/facebook_combined_tab.txt"  
        fpOUT01 = open(filename_OUT01, 'w')
    except IOError as e:
        print "Error de E/S ({0}): {1}".format(e.errno, e.strerror)
        sys.exit()

    # regresa apuntadores a archivos
    return fpIN01, fpOUT01, filename_OUT01

# XXXXXXXXXXXXXXX END: extra processing for SNAP XXXXXXXXXXXXXXXXXXXX



def main():
    #obtainGraph() # DROPPED because not needed, in favor of commPart()
    # inFilePointer01: ptr to original file with edges (sep: " ")
    # outFilePointer01: ptr to outfle to write tab-separated file of edges
    # inFileName: name of the processed file with tab-separated edges
    # (it is called inFile because it is the input for BigCLAM)
    inFilePointer01, outFilePointer01, inFileNameBC = abrir_archivos()

    # listEdges is list of lists [[1, 2], [3, 4], ...]
    #ok first 10: listEdges = leer_texto(inFilePointer01)[:10]
    listEdges = leer_texto(inFilePointer01)

    '''
    gfb = nx.read_edgelist("data/facebook_combined.txt",
                           create_using=nx.Graph(),
                           nodetype=int
                           )
    '''
    
    '''
    gfb = nx.read_edgelist("data/facebook_combined_small.txt",
                           create_using=nx.Graph(),
                           nodetype=int
                           )

    print "\n--Reading matrix from file"
    print "first gfb.nodes", gfb.nodes()[:5]
    print "first gfb.edges", gfb.edges()[:5]
    print nx.info(gfb)
    '''
    print "\n--Community creation"
    # drawing of communities disabled for PR code
    #okButDisabled: detectComm(gfb, None)


    #ok but disabled becuase it is networkx
    # ok original networkx: parts is a dict having pairs of node-community
    # ok original networkx: parts = community.best_partition(gfb)
    parts = commPartition(listEdges, outFilePointer01, inFileNameBC)
    #print type(parts)
    print "communityInfo"
    # parts is a list of 2-tuples, having pairs of (node, community)
    # [(0, 0), (1, 0), (2, 0), (3, 0), ..., (4037, 6), (4038, 6)]
    #ok_first10:
    print "first parts", parts[:5]
    #ok_all_items: pprint.pprint(parts)
    #ok_last10:
    print "last parts ", parts[-5:]

    # Check communities
    # send only the values of communities to get the list of unique communit.
    # ok original networkx: uniqCommL, lenUL = getUniqueComm(parts.values())
    commValues = [y for (x, y) in parts]
    uniqCommL, lenUL = getUniqueComm(commValues)
    print "uniqComm", uniqCommL
    pprint.pprint(uniqCommL)
    print lenUL

    # obtain incidence matrix vertex-hyperedge (H)
    # en networkx basta con meter la lista que salio de las comunidades
    # so, it is not necessary to obtain unique communities
    # 1) generate the graph (NOT NEEDED anymore)
    """"
    print "communityGraphInfo"
    GComm=nx.Graph()
    GComm.add_edges_from(parts)
    print GComm.nodes()[:10]
    #ok_4039_nodes: print GComm.nodes()
    print GComm.edges()[:10]
    print GComm.edges()
    print nx.info(GComm)
    """


    print "\n--Creating Incidence Matrix from node-community pairs"
    # 2) send 1) node-community 2-tuple list & number of unique communities
    # get the incidence matrix
    # sparse format
    H = createPrevCSR(parts, lenUL)
    print H[:10]

    print "\n--Creating Inverse Matrix from node-community pairs"
    # 3) get invDvMatr

    ################# TEST ##########
    copy = parts
    #ok print copy
    print len(copy)
    copy.extend([(0,1), (1,1)])
    print len(copy)
    invDvMatr = getInvDvMatr(copy)
    ################# TEST ##########
    #ok_GOOD: invDvMatr = getInvDvMatr(parts)
    print invDvMatr

    dv_1H = invDvMatr.dot(H)
    print "dv-1H", dv_1H


    # 4) get invDeMatr
    # es lo mismo que las Dv, pero volteamos los n'umeros y aplicamos lo 
    # mismo que Dv. volteamos indices:
    newList = [(item[1], item[0]) for item in copy]
    # ok newlist: print "newList for De", newList
    invDeMatr = getInvDvMatr(newList)
    print invDeMatr

    de_1Ht = invDeMatr.dot(H.T)
    print "de-1Ht", de_1Ht


    # --- Ahora la lmultiplicacion completa:

    MH = dv_1H.dot(de_1Ht)
    print "MH", MH

    # --- Matriz de teletransportacion
    matSize = len(MH)
    teleMat = [[1. / matSize for j in range(matSize)] for i in range(matSize)]
    #check: print "teleMat", teleMat[:5][:5]

    teleMat = np.matrix(teleMat)
    print "teleMat", teleMat[:5, :5]

    MH2 = 0.85 * MH + 0.15 * teleMat
    print "MH2", MH2[:5,:5]

    print type(MH2)


    # eigenvectors
    '''
    evalM, evecM = np.linalg.eig(MH2)
    print "eigenvals", evalM[:5], "\n"
    print "eigenvec", evecM[:, 1]
    '''
    # iteraciones para sacar lso demas eigenvectores
    # necesitamos el primer eigenvector (eigenvec)
    # Define initial vector v
    v = np.matrix(makeInitVec(matSize))
    print("v", v.T)

    # Obtenemos los eigenvec y eigenval con A y v
    # NOTE: we need to send the INVERSE of MH2 to make it work with algorithm
    eigVal, eigVec = iteration(MH2, v.T)  # v needs to be a matrix
    print("eigVal = ", eigVal, type(eigVal))
    print("eigVec = ", eigVec, type(eigVec))

    nodeEigVec = [(i, j) for i, j in enumerate(np.matrix.tolist(eigVec))]
    print nodeEigVec[:5]
    pprint.pprint(nodeEigVec)
    print "sorted", sorted(nodeEigVec, key=operator.itemgetter(1), reverse=True)[:5]


    '''
    print "\n--Obtaining PageRank for Graph"
    # PR of graph - dictionary {node: pr value}
    pr = nx.pagerank(gfb)
    #print(type(pr))
    print pr.items()[-3:] # print last 10 elements
    pprint.pprint(pr.items() )
    #print top 10 items w higher pr score
    print sorted(pr.items(), key=operator.itemgetter(1), reverse = True)[:5]
    '''
        
    # close files    
    inFilePointer01.close()
    outFilePointer01.close()

if __name__ == '__main__':
    main()
