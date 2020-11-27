#!/usr/bin/env python3
"""
Here is a sample program to decompose a hydrogen bond network into cycles.
It is as close as possible to the one used in the paper, but it is slightly different from the code used in the actual analysis.
"""

import numpy as np
import networkx as nx
import pairlist as pl

def read_water(file):
    """
    A sample of the loader

    (Assume an MDVIEW-type file)
    """
    # Assume that the first line is the cell dimension
    cell = [float(x) for x in file.readline().split()[1:]]
    # cell shape matrix (unit is AA)
    cell = np.diag(cell)
    while True:
        line = file.readline()
        if line[0] != "-":
            break
    N = int(line)
    waters = []
    H = []
    O = []
    for i in range(N):
        cols = file.readline().split()  # name X Y Z
        if cols[0][0] == "O":
            O += [float(x) for x in cols[1:]]
        elif cols[0][0] == "H":
            H += [float(x) for x in cols[1:]]
        if len(O) == 3 and len(H) == 6:
            waters.append(O+H) # nine numbers
            H = []
            O = []
    # inverse of the cell matrix
    celli = np.linalg.inv(cell)
    # fractional coordinate
    waters = np.array(waters).reshape([-1, 3]) # OHH order
    rpos = np.array(waters) @ celli
    return rpos, cell




def all_cycles(g, size):
    """
    List the cycles of the size only. No shortcuts are allowed during the search.
    """

    def find_cycle(history):
        """
        Recursively find a homodromic cycle.

        The label of the first vertex in the history (head) must be the smallest.
        """
        head = history[0]
        last = history[-1]
        if len(history) == size:
            for next in g.successors(last):
                if next == head:
                    # test the dipole moment of a cycle.
                    d = np.zeros(3)
                    for i in range(len(history)):
                        a,b = history[i-1], history[i]
                        d += g[a][b]["vec"]
                    if np.allclose(d, np.zeros(3)):
                        yield history
        else:
            for next in g.successors(last):
                if next < head:
                    # members must be greater than the head
                    continue
                if next not in history:
                    yield from find_cycle(history+[next])

    for head in g.nodes():
        yield from find_cycle([head])


def AtoX(A):
    # 環の個数のほうが十分多い場合の近似
    # 特異値分解
    Q1, S, Q2T = np.linalg.svd(A)
    # print(Q1.shape, S.shape, Q2T.shape)
    # ほぼ0の要素を0にする
    S[np.abs(S)<1e-12] = 0
    rank = np.sum(np.abs(S) > 0)
    # print(S,r)
    # SS = Q1.T @ A @ Q2T.T
    # 対角行列S†の準備
    Sd = np.zeros_like(A).T
    Sd[:rank, :rank] = np.diag(1/S[:rank])
    # print(SS)
    # print(Sd.shape)
    # A†
    Ad = Q2T.T @ Sd @ Q1.T
    #print(Ad.shape)
    b = np.ones(A.shape[0])
    x = Ad@b
    # print(A@x)
    # print(x)
    return x, rank


def weight(g):
    # label the edges
    HBcode = {(d,a): x for x, (d,a) in enumerate(dg.edges())}

    A = []
    cycles = []
    for s in range(4, 16):
        lastA = len(A)
        for cycle in all_cycles(g,s):
            cycles.append(cycle)
            row = np.zeros(len(HBcode))
            for j in range(len(cycle)):
                edge = HBcode[cycle[j-1], cycle[j]]
                row[edge] = 1.0
            A.append(row)
        if lastA==len(A):
            continue
        lastA = len(A)
        AT = np.array(A).T
        # Quasi-inverse matrix
        x, rank = AtoX(AT)
        print(f"Largest cycle size {s}")
        print(f"Number of cycles {len(cycles)}")
        b = AT @ x
        print(f"Sum weight: {b}")
        if np.allclose(b, np.ones_like(b)):
            return cycles, x



# read the atomic positions
# e.g. genice 1c -r 4 4 4 -f mdview > 1c.mdv
with open("1c.mdv") as file:
    rpos, cell = read_water(file)

# define the directed graph
Opos = rpos[::3].copy()
H1pos = rpos[1::3].copy()
H2pos = rpos[2::3].copy()
dg = nx.DiGraph()
# add the edge if intermolecular O-H distance is less than 2.45 AA
for i,j,d in pl.pairs_iter(H1pos, 2.45, cell=cell, pos2=Opos, distance=True):
    if d > 1.5:
        vec = Opos[j] - Opos[i]
        vec -= np.floor(vec + 0.5) # PBC
        dg.add_edge(i,j, vec=vec)
for i,j,d in pl.pairs_iter(H2pos, 2.45, cell=cell, pos2=Opos, distance=True):
    if d > 1.5:
        vec = Opos[j] - Opos[i]
        vec -= np.floor(vec + 0.5) # PBC
        dg.add_edge(i,j, vec=vec)


# find all the cycles and determine the weights
cycles, weights = weight(dg)
for cycle, weight in zip(cycles, weights):
    print(cycle, weight)
