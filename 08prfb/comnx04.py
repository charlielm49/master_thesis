#!/home/carlos/bin/anaconda2/bin/python
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
        #filename_in01 = "data/facebook_combined.txt"  # este sí cambia
        #fpIN01 = open(filename_in01, 'r')
        # outfile for Cytoscape visualization
        filename_OUT01 = "data/fb_nx_modular.txt"
        fpOUT01 = open(filename_OUT01, 'w')
    except IOError as e:
        print "Error de E/S ({0}): {1}".format(e.errno, e.strerror)
        sys.exit()

    # regresa apuntadores a archivos
    return fpOUT01


def main():
    outFilePointer01 = abrir_archivos()


    gfb = nx.read_edgelist("data/facebook_combined.txt",
                           create_using=nx.Graph(),
                           nodetype=int
                           )
    '''
    '''
    #test-begin
    #gfb = nx.Graph()
    #gfb.add_edges_from()
    #test-end
    '''
    gfb = nx.read_edgelist("data/facebook_combined_small.txt",
                           create_using=nx.Graph(),
                           nodetype=int
                           )
    '''
    print "\n--Reading matrix from file"
    print "first gfb.nodes", gfb.nodes()[:5]
    print "first gfb.edges", gfb.edges()[:5]
    print nx.info(gfb)

    #XXXXXXXXXXXXXXX BEG: COMMUNITY STUFF XXXXXXXXXXXXXXX
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

    # parts is a dictionary having pairs of node-community
    parts = community.best_partition(gfb)
    # print type(parts)
    partsList = parts.items() # create a new list var to manipulate

    print "communityInfo"
    # partsList is a list of 2-tuples, having pairs of (node, community)
    # [(0, 0), (1, 0), (2, 0), (3, 0), ..., (4036, 6), (4037, 6), (4038, 6)]
    #ok_first10:
    print "first partsList ", partsList[:5]
    #ok_all_items: pprint.pprint(partsList)
    #ok_last10:
    print "last partsList ", partsList[-5:]
    # send node-comm values to file
    escribir_texto(partsList, outFilePointer01)

    # send only the values of commununities to get the list of unique communit.
    uniqCommL, lenUL = getUniqueComm(parts.values())

    #XXXXXXXXXXXXXXX END: COMM STUFF XXXXXXXXXXXXXXX


    #XXXXXXXXXXXXXXX BEG: HYPERG STUFF XXXXXXXXXXXXXXX

    #


    #XXXXXXXXXXXXXXX END: HYPERG STUFF XXXXXXXXXXXXXXX


    #XXXXXXXXXXXXXXX BEG: MATRIX STUFF XXXXXXXXXXXXXXX

    # obtain incidence matrix vertex-hyperedge (H)
    # en networkx basta con meter la lista que salio de las comunidades
    print "\n--Creating Incidence Matrix from node-community pairs"
    # 2) send 1) node-community 2-tuple list & number of unique communities
    # get the incidence matrix
    # sparse format
    H = createPrevCSR(partsList, lenUL)
    print H[:10]

    print "\n--Creating Inverse Matrix from node-community pairs"
    # 3) get invDvMatr

    ################# TEST ##########
    copy = partsList
    #ok print copy
    print len(copy)
    copy.extend([(0,1), (1,1)])
    print len(copy)
    invDvMatr = getInvDvMatr(copy)
    ################# TEST ##########
    #ok_GOOD: invDvMatr = getInvDvMatr(partsList)
    print invDvMatr

    dv_1H = invDvMatr.dot(H)
    print "dv-1H", dv_1H


    # 4) get invDeMatr
    # es lo mismo que las Dv, pero volteamos los n'umeros y aplicamos lo mismo que Dv
    # volteamos indices:
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
    eigVal, eigVec = iteration(MH2, v.T)  # need to be a matrix
    print("eigVal = ", eigVal, type(eigVal))
    print("eigVec = ", eigVec, type(eigVec))

    nodeEigVec = [(i, j) for i, j in enumerate(np.matrix.tolist(eigVec))]
    print nodeEigVec[:5]
    pprint.pprint(nodeEigVec)
    print "sorted", sorted(nodeEigVec, key=operator.itemgetter(1), reverse=True)[:5]



    print "\n--Obtaining PageRank for Graph"
    # PR of graph - dictionary {node: pr value}
    pr = nx.pagerank(gfb)
    #print(type(pr))
    print pr.items()[-3:] # print last 10 elements
    pprint.pprint(pr.items() )
    #print top 10 items w higher pr score
    print sorted(pr.items(), key=operator.itemgetter(1), reverse = True)[:5]




if __name__ == '__main__':
    main()
