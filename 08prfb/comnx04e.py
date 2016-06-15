#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""""

#XXXXXXXXX-BEGIN: REQUIREMENTS & SPECS

#from __future__ import print_function # to run in python2 & 3
import os.path
import sys
import struct # Packing and unpacking of heterogeneous binary data
import array
import math
import networkx as nx  # graph library
import numpy as np
import matplotlib.pyplot as plt
import community  # for community detection to form hypergraphs
import operator  # for helping sort a dictionary
from scipy.sparse import *

#opt
import pprint

#XXXXXXXXX-END: REQUIREMENTS & SPECS



# construct incidence matrix (and check that there are no isolated vertices in a community)
# get list of unique communities
def getUniqueComm(vertexCommunityL):
    commUniqList = []
    for comm in vertexCommunityL:
        if comm not in commUniqList:
            commUniqList.append(comm)

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
#  of pairs [(a, b), (c, d), ...]
# <- first coord of each pair (a, ...) (c, )goes in rows of csr
# <- second coord (..., b) ( , d)goes in cols of csr
def createPrevCSR(indexList, numVertices, numEdges):
    #global csrMat # OUTPUT
    print "\nEnter createPrevCSR for incidence matr"
    print "numEdges", numEdges

    # get the
    matSize = numVertices # assignment for downstream code consistency
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
    MAX_ITER = 10000
    EPS = 1.0e-6 # works ok for simple cases of PageRank
    EPS = 1.0e-12
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

# NOTE: ONLY works for non-overlapping communities
# if node-comm pair not in original node-pair list, add it
# edgeList: [[0, 0], [1, 2], [2, 1], ..., [5, 1], [5, 4]]
# oldComm (partsList): [[0, 0], [1, 0], ..., [4037, 6], [4038, 6]]
# NOTE: partsList can have many groups for the same node, so we can't 
#       have a dict for this list
def newEdgeComm(listEdges, oldComm): #ok
    lista = {} #node-comms
    #{1:[1, 2, 3], 2:[3, 4], ...}
    
    # change oldComm to dir for efficiently finding of node
    oldCommD = dict(oldComm)
    #ok: print "oldCommD", oldCommD
    
    for edgePair in listEdges: #[0, 0], [0, 1]
        #ok: print edgePair
        Ed1 = edgePair[0] # first edge
        Ed2 = edgePair[1] # second edge
        newCommVal = oldCommD[Ed2] # group of the second edge
        #ok: print "newCommVal", newCommVal
        # if first edge in dict, add more 
        if edgePair[0] in lista:
            # prev list of grp values
            oldCommVal = lista[Ed1]
            #ok: print "oldCommVal", oldCommVal
            #paste new value to old list
            oldCommVal.extend([newCommVal])
            # add group of second edge to previous grp value/values
            lista.update({Ed1: oldCommVal})
        else:
            currentVal = [oldCommD[Ed1]]
            #ok: print "currentVal", currentVal
            currentVal.extend([newCommVal])
            lista.update({Ed1: currentVal})
            #ok: print "lista", lista.items()
    
    return lista
#version 2 -> returns a list of tuples, instead of a dictionary
# this format is necessary for building sparse matrices
# edgeList = [[0, 0], [1, 2], [2, 1], [3, 0], [3, 1], [4, 1], \
# [4, 3], [4, 5], [5, 1], [5, 4]]
# partsList = [[0, 0], [1, 1], [2, 2], [3, 0], [4, 1], [5, 1]]
# OUT01: [[0, 0], [1, 0], ..., [4037, 6], [4038, 6]] <- w/new grps
def newEdgeComm2(listEdges, oldComm):  # ok
    lista = []  # node-comms
    # OUT: [[0, 0], [0, 1], [1, 1], [2, 3], [3, 4], ...]

    # change oldComm to dir for efficiently finding of node
    oldCommD = dict(oldComm)
    print "oldCommD", oldCommD

    for edgePair in listEdges:  # [0, 0], [0, 1] loop over EDGES
        #ok: print "edgePair", edgePair
        Ed1 = edgePair[0]  # first edge
        Ed2 = edgePair[1]  # second edge
        oldCommVal = oldCommD[Ed1]  # group of the first edge
        #ok: print "oldCommVal", oldCommVal
        newCommVal = oldCommD[Ed2]  # group of the second edge
        #ok: print "newCommVal", newCommVal
        # build edge-comm pair
        oldPair = (Ed1, oldCommVal) # need tuples to deduplicate fast
        # build edge-comm pair
        newPair = (Ed1, newCommVal)
        # if edge-comm pair in list, add more
        lista.append(oldPair)
        if newPair not in lista:
            # add group of second edge to previous grp value/values
            lista.append(newPair)
            #ok: print "lista", lista

    #deduplicate list:
    lista = set(lista)

    return lista


def leer_texto(fp):
    tmpFull = []
    # procesamiento de archivo de entrada completo
    for line in fp:
        tmp = line.split()
        #ok print tmp
        tmpFull.append([int(tmp[0]), int(tmp[1])])

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
    separator = " "

    for item in listaL:
        item = separator.join(str(x) for x in item) + "\n"
        fp.write(item)

# Abre archivos entrada/salida / Regresa: apuntadores a archivos
# 2015FEB21: ya no se ponen parámetros en línea de comandos porque
# ya hay muchos archivos que leer - sólo el del archivo a procesar
# file_in01: archivo en raw text a procesar
# file_out01: archivo de salida ya procesado
def abrir_archivos():
    # print sys.argv # 0: nombre de programa 1: nombre de archivo
    try:
        #filename_IN01 = "data/facebook_combined.txt"  # este sí cambia
        # no blank lines at the end of file
        filename_IN01 = "data/smallnet.txt"
        fpIN01 = open(filename_IN01, 'r')
        # outfile for Cytoscape visualization
        filename_IN02 = "data/smallgrp.txt"
        fpIN02 = open(filename_IN02, 'r')
    except IOError as e:
        print "Error de E/S ({0}): {1}".format(e.errno, e.strerror)
        sys.exit()

    # regresa apuntadores a archivos
    return fpIN01, fpIN02




def main():

    #XXXXXXXXXXXXXXX BEG: FILE & GRAPH READING XXXXXXXXXXXXXXX

    # Read the network to make the hypergraphs
    # inFilePointer01 is ptr (read) to network file
    # inFilePointer02 is ptr (read) to groups file
    inFilePointer01, inFilePointer02 = abrir_archivos()

    # This is for automatic dtection of groups with NX
    # NOTE: for the Klein example, the graph needs to be directed
    '''
    gfb = nx.read_edgelist("data/facebook_combined.txt",
                           create_using=nx.Graph(),
                           nodetype=int
                           )


    gfb = nx.read_edgelist("data/smallnet.txt",
                           create_using=nx.Graph(),
                           nodetype=int
                           )
    '''
    # Use DiGraph for directed graph
    gfb = nx.read_edgelist("data/smallnet.txt",
                           create_using=nx.DiGraph(),
                           nodetype=int
                           )

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

    # Network reading from file
    edgeList = leer_texto(inFilePointer01)
    print "\nEdge/Graph Info"
    print "edgeList (data read, first 10):"
    pprint.pprint(edgeList[:10])

    #XXXXXXXXXXXXXXX END: FILE & GRAPH READING XXXXXXXXXXXXXXX



    #XXXXXXXXXXXXXXX BEG: COMMUNITY STUFF XXXXXXXXXXXXXXX
    '''
    print "\n--Community creation"
    # community Graph
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


    # drawing of communities disabled for PR code
    #okButDisabled: detectComm(gfb, None)
    '''
    #original automatic community detection: parts = community.best_partition(gfb)
    # parts is a dictionary having pairs of node-community
    # print type(parts)
    #partsList = parts.items() # create a new list var to manipulate

    # Manual reading of node-Community file (artificial example)
    partsList = leer_texto(inFilePointer02)

    print "\nCommunity Info"
    # partsList is a list of 2-tuples, having pairs of (node, community)
    # [(0, 0), (1, 0), (2, 0), (3, 0), ..., (4036, 6), (4037, 6), (4038, 6)]
    #ok_first10:
    print "first partsList "
    pprint.pprint(partsList[:10])
    #ok_all_items: pprint.pprint(partsList)
    #ok_last10:
    print "last partsList "
    pprint.pprint(partsList[-10:])
    # send node-comm values to file
    #escribir_texto(partsList, inFilePointer01)

    #XXXXXXXXXXXXXXX END: COMM STUFF XXXXXXXXXXXXXXX


    '''
    #XXXXXXXXXXXXXXX BEG: HYPERG STUFF XXXXXXXXXXXXXXX

    # construct new overlapping hyperedges from non-overlapping ones
    # get new edge-community
    # edgeList = [[0, 0], [1, 2], [2, 1], [3, 0], [3, 1], [4, 1], \
    # [4, 3], [4, 5], [5, 1], [5, 4]]
    # partsList = [[0, 0], [1, 1], [2, 2], [3, 0], [4, 1], [5, 1]]
    #ok for adjacency-list format:
    # x = newEdgeComm(edgeList, partsList) # send edge list & old comm list
    #ok: print "x", x
    # This returns one list of lists with all node-comm pairs: (overlapping)
    overlapComms = newEdgeComm2(edgeList, partsList) # send edge list & old comm list
    print "overlapComms", overlapComms

    # Var reassignment to avoid disruption in the downstream code
    partsList = overlapComms

    #XXXXXXXXXXXXXXX END: HYPERG STUFF XXXXXXXXXXXXXXX
    '''


    #XXXXXXXXXXXXXXX BEG: MATRIX STUFF FOR HG XXXXXXXXXXXXXXX

    # send only the values of commununities to get the list of unique communit.
    commVals = [operator.itemgetter(1)(item) for item in partsList]
    uniqCommL, lenUComL = getUniqueComm(commVals)
    print "uniqCommL, lenUComL", uniqCommL[:5], lenUComL

    # send only the values of edges to get the list of unique edges.
    edgeVals = [operator.itemgetter(0)(item) for item in partsList]
    uniqEdgeL, lenUEdgL = getUniqueComm(edgeVals)
    print "uniqEdgeL, lenUEdgL", uniqEdgeL[:5], lenUEdgL


    # obtain incidence matrix vertex-hyperedge (H)
    # en networkx basta con meter la lista que salio de las comunidades
    print "\n--Creating Incidence Matrix from node-community pairs\n"
    # 2) send 1) node-community 2-tuple list & number of unique communities
    # get the incidence matrix
    # sparse format
    H = createPrevCSR(partsList, lenUEdgL, lenUComL)
    print H[:10]

    print "\n--Creating Inverse Matrix from node-community pairs\n"
    # 3) get invDvMatr
    #ok_GOOD: 
    invDvMatr = getInvDvMatr(partsList)
    print "inv vertex Degree Matrix:\n", invDvMatr[:10]

    dv_1H = invDvMatr.dot(H)
    print "Mult of inv vertex Degree Matrix by Incidence:)\n", \
        dv_1H[:10]


    # 4) get invDeMatr
    # es lo mismo que las Dv, pero volteamos los numeros y aplicamos
    # lo mismo que Dv, solo intercambiamos indices i -> j:
    newList = [(item[1], item[0]) for item in partsList]
    # ok newlist: 
    print "newList for De (changing the indexes)", newList
    invDeMatr = getInvDvMatr(newList)
    print "inv h-edge Degree Matrix:\n", invDeMatr[:10]

    de_1Ht = invDeMatr.dot(H.T)
    print "Mult of inv h-edge Degr Matr by Incidence.T:)\n", \
        de_1Ht[:10]


    # --- Ahora la lmultiplicacion completa:

    MH = dv_1H.dot(de_1Ht)
    print "\nMH\n", MH[:10]

    # --- Matriz de teletransportacion
    matSize = len(MH)
    teleMat = [[1. / matSize for j in range(matSize)] for i in range(matSize)]
    #check: print "teleMat", teleMat[:3][:3]

    teleMat = np.matrix(teleMat)
    print "teleMat\n", teleMat[:3, :3]

    MH2 = 0.85 * MH + 0.15 * teleMat
    print "MH2\n", MH2[:10,:10]

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
    print "v", v.T[:3]

    # Obtenemos los eigenvec y eigenval con A y v
    # NOTE: the code works correctly for column stochastic matrices
    # so, if the matrix is row-stochastic, we need to transpose it
    eigVal, eigVec = iteration(MH2.T, v.T)  # need to be a matrix
    print "eigVal = ", eigVal, type(eigVal)
    print "eigVec = \n", eigVec[:10], type(eigVec)

    nodeEigVec = [(i, j) for i, j in enumerate(np.matrix.tolist(eigVec))]
    #ok: print nodeEigVec[:10]
    pprint.pprint(nodeEigVec[:10])
    print "sorted", sorted(nodeEigVec, key=operator.itemgetter(1), \
        reverse=True)[:10]

    #XXXXXXXXXXXXXXX END: MATRIX STUFF FOR HG XXXXXXXXXXXXXXX


    #XXXXXXXXXXXXXXX BEG: MATRIX STUFF FOR G XXXXXXXXXXXXXXX

    # We need to build the vertex-edge list [(), (), ...]
    # We already have the graph node-node information in edgeList
    
    # obtain incidence matrix vertex-vertex (H)
    # en networkx we only need to add the list of (edge, edge)
    print "\n--Creating Incidence Matrix from edge-edge pairs\n"
    # sparse format
    G = createPrevCSR(edgeList, lenUEdgL, lenUEdgL)
    print G[:10]
    
    print "\n--Creating Inverse Matrix from edge-edge pairs\n"
    # get invDvMatr
    invDvMatrG = getInvDvMatr(edgeList)
    print "inv vertex Degree Matrix G:\n", \
        invDvMatrG[:10]

    dv_1G = invDvMatrG.dot(G)
    print "Mult of inv vertex Degree Matrix by Incidence \n", \
        dv_1G[:10]
    
    # --- Matriz de teletransportacion
    matSize = len(dv_1G)
    teleMat = [[1. / matSize for j in range(matSize)] for i in range(matSize)]
    #check: print "teleMat", teleMat[:3][:3]

    teleMat = np.matrix(teleMat)
    print "teleMat\n", teleMat[:3, :3]

    # up to here, everything is fine with example from Klein
    # correct results ae obtained
    #ok simple graph of Klein: MG = 0.85 * dv_1G + 0.15 * teleMat
    #ok: print "MG\n", MG[:10,:10]
    
    # but now we add the HG matrix to this (telemat already included)
    # THINK about whre to add the teletransp matrix
    # CAREFUL: add two rw-stochastic or two column-stochastic matrices
    MG = 0.85 * dv_1G + 0.15 * MH2
    print "MG\n", MG[:10,:10]

    print type(MG)


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
    print "v", v.T[:3]

    # Obtenemos los eigenvec y eigenval con A y v
    # NOTE: the code works correctly for column stochastic matrices
    # so, if the matrix is row-stochastic, we need to transpose it
    eigVal, eigVec = iteration(MG.T, v.T)  # need to be a matrix
    print "eigVal = ", eigVal, type(eigVal)
    print "eigVec = \n", eigVec[:10], type(eigVec)

    # printing node-eigenvector value
    nodeEigVec = [(i, j) for i, j in enumerate(np.matrix.tolist(eigVec))]
    #ok: print nodeEigVec[:10]
    pprint.pprint(nodeEigVec[:10])
    print "sorted", sorted(nodeEigVec, key=operator.itemgetter(1), \
        reverse=True)[:10]
    # (Sorted)Results for original Klein example of simple graph
    #  (1, [0.618189433764439]), (2, [0.5737952417102831]),
    # (0, [0.522037084149513]), (4, [0.07830556262242697]),
    # (3, [0.07052079908686405]), (5, [0.07052079908686405])

    #XXXXXXXXXXXXXXX END: MATRIX STUFF FOR G XXXXXXXXXXXXXXX


    print "\n--Obtaining PageRank for Graph"
    # PR of graph - dictionary {node: pr value}
    pr = nx.pagerank(gfb)
    # print(type(pr))
    print pr.items()[-3:]  # print last 10 elements
    # ok print:
    pprint.pprint(pr.items())
    # print top 10 items w higher pr score
    print "sorted PR", sorted(pr.items(), key=operator.itemgetter(1), \
                              reverse=True)[:10]

if __name__ == '__main__':
    main()
