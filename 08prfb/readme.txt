comnx01.py - ok
first version, only communities, works ok


comnx02.py - ok
adding with networkx:
PR w/graphs (nx algorithm)
communities
Incidence Matrix with CSR format from node-community
Inverse vertex degree Matrix


comnx03.py - ok
continuing 02 for security
FB Graph info from networkx:
first 5 gfb.nodes [0, 1, 2, 3, 4]
first 5 gfb.edges [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
Name: 
Type: Graph
Number of nodes: 4039
Number of edges: 88234
Average degree:  43.6910

PR w/hypergraphs: ok
 sparse Multiplication of Dv-1 * H
 sparse Multiplication of De-1 * Ht
 sparse Multiplication of (Dv-1 * H) * (De-1 * Ht)

non-overlapping comm detection ok (Modularity)

HERE obtained results for pagerank


comnx04.py - ok
continuing 03 for security
-adding variable partsList to avoid usage of dict.items() - ok
-adding export of node-comm file for vis in cyto - ok
-constructing hg by hand complementing modularity


comnx04a.py - small branch ONLY to MANUALLY test all hypergraph things- ok
continuing 04 for security
-changed code to external file testing - ok

comnx04b.py - ok
continuing 04a for security
-constructing hg by hand complementing modularity - ok algorithm
need to chamge the output

comnx04c.py - ok
continuing 04b for security
-constructing hg by hand complementing modularity - ok algorithm
working og but output is adjacency list, need to chamge the output - ok
- chamging size of matrix in new schema - ok
(now is unique nodes x unique comms size)

comnx04d.py - ok
continuing 04c for security
add the matrix for the network to the hypergraphs - ok


comnx05.py - 
continuing 04d for security
example of facebook, obtaining ranks, compared to PR-

obtain results
1) overlapping (ol) hg only vs graph PR
2) (graph + ol hg) vs graph PR




Note on the matrix addition of all matrices:
1) isolated hypergraph (non-overlapping grps)
2) overlapping hg (overlapping grps)
3) graph matrix
4) teletransp matrix

maybe we can't add the overlapping hg to the graph matrix because we already have some structure of the graph in the overlapping matrix so, some effect would be duplicated

Do we need to add the graph matrix??


--------------------------------------------------


comsnap01.py - not ok
All is the same as comnx03.py (starting point), only 1 change: overlapping 
community detection with SNAP.

numEdges 10
matSize 2808
Traceback (most recent call last):
  File "/media/KING16-3/maestria/TESIS/CODE/08prfb/comm03a.py", line 510, in <module>
    main()
  File "/media/KING16-3/maestria/TESIS/CODE/08prfb/comm03a.py", line 422, in main
    H = createPrevCSR(parts, lenUL)
  File "/media/KING16-3/maestria/TESIS/CODE/08prfb/comm03a.py", line 99, in createPrevCSR
    csrMat = csr_matrix( (data, (rows, cols)), shape = (matSize, numEdges) ).toarray()
  File "/home/carlos/bin/anaconda2/lib/python2.7/site-packages/scipy/sparse/compressed.py", line 48, in __init__
    other = self.__class__(coo_matrix(arg1, shape=shape))
  File "/home/carlos/bin/anaconda2/lib/python2.7/site-packages/scipy/sparse/coo.py", line 182, in __init__
    self._check()
  File "/home/carlos/bin/anaconda2/lib/python2.7/site-packages/scipy/sparse/coo.py", line 236, in _check
    raise ValueError('row index exceeds matrix dimensions')
ValueError: row index exceeds matrix dimensions


this maybe is happening because when obtaining the communities, not all nodes are assigned a community number. so this causes the problme in scipy:
there are less nodes asigned to communities and the error occurs

WILL TRY: paste the nodes to whichever goup they are inked to

NOTE: there was a bug in the code: I used the space-separated file as input for big clam. However, theprogram didn't complain and threw output which looked correct!!! (but wasn't)


comsnap02.py - DISCONTINUED
continuing from comsnap01 for safety
Info from BigClam process:s
Graph: 3984 Nodes 87995 Edges [clm: ????]
rearrage nodes
...
[(1912, 0), (2347, 0), (2543, 0), (2266, 0), (1985, 0), (2233, 0), (2142, 0), (2206, 0), (2410, 0), (2229, 0)]
2815 [clm: nodes in output]
output unique nodes to file so it can be read from readComm.py - ok
FIX: paste the nodes to whichever goup they are inked to - DISCONT.





