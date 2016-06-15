#!/usr/bin/env python
# -*- coding: utf-8 -*-
""""""
import sys
import pprint
import operator  # for helping sort a dictionary

def abrir_archivos():
    # print sys.argv # 0: nombre de programa 1: nombre de archivo
    try:
        file_in01 = "communities.txt"  # este sí cambia
        fpIN01 = open(file_in01, 'r')
    except IOError as e:
        print "Error de E/S ({0}): {1}".format(e.errno, e.strerror)
        sys.exit()

    # regresa apuntadores a archivos
    return fpIN01

# OUT: list
# ['1912', '2347', '2543', '2266', ...] (some are repeated)
def leer_texto(fp):
    tmpFull = []
    # procesamiento de archivo de entrada completo
    for line in fp:
        #ok prints each line of file: print line
        tmp = line.split()
        #ok print tmp
        tmpFull.extend(tmp) # adds all the elements in each line
    
    # regresa una lista con todas las palabras del texto en esa lista
    return tmpFull


# counts the number of groups each node belongs to
# OUT01: dict
#{'4023': 1, '1868': 1, '345': 1, ,... }
def contar(lista):
    conteo = {}
    
    for item in lista:
        if item in conteo:
            conteo[item] += 1
        else:
            conteo.update({item: 1})
    
    return conteo


def getUniqueComm(vertexCommunityL):
    commUniqList = []
    for comm in vertexCommunityL:
        if comm not in commUniqList:
            commUniqList.append(comm)

    lenUL = len(commUniqList)
    #print commUniqList
    return commUniqList, lenUL


def buscar_faltantes(lista_uniq_ori):
    # sacar lista de unicos de todos los originales (ya)
    pass
    # sacar lista de unicos de lista con comunidades
    
    # obtener diferencia
    
    # ver esa diferencia, a que otros nodos esta asociado
    # (buscar en la lista original de edges)
    
    # asignarle el grupo (o LOS grupos de aquellos gupos a los que 
    # pertenezca el nodo al cual esta conectado)


def main():
    # lista has the list of all nodes (some repeated)
    lista = leer_texto(abrir_archivos())
    #ok: print lista
    # ['1912', '2347', '2543', '2266', ...]
    conteo = contar(lista)
    pprint.pprint(sorted(conteo.items(), \
        key=operator.itemgetter(1) ) \
        )
    #    [('4023', 1),     ('1868', 1),     ('345', 1),     ('346', 1),     ('341', 1),
    #  ('1166', 2), ('414', 3), ('917', 3), ('107', 5)]
    print len(conteo.items())
    
    #buscar faltantes en lista de grupos
    # conteo.keys() has unique nodes de lista de nodos-comunidad
    buscar_faltantes(conteo.keys())
    

if __name__ == '__main__':
    main()
