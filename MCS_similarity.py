import time
import functools
from collections import defaultdict
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        return run_time, value
    return wrapper_timer

def get_atom_envs(mol):
    info = {}
    fp = AllChem.GetMorganFingerprint(mol, 1, bitInfo=info)
    return info
    
def dfs(node, graph, subgraph, visited, nodes):
    visited.add(node)
    for neighbor in graph[node]:
        if neighbor in nodes and neighbor not in visited:
            subgraph.append(neighbor)
            dfs(neighbor, graph, subgraph, visited, nodes)

def search_subgraphs(graph, nodes):
    subgraphs = []
    visited = set()
    for node in nodes:
        if node not in visited:
            subgraph = [node]
            dfs(node, graph, subgraph, visited, nodes)
            subgraphs.append(subgraph)
    return subgraphs

def search_frags_with_atoms(mol, atoms_to_use):
    if len(atoms_to_use) == 0:
        return [[]]
    graph = {atom.GetIdx(): [n.GetIdx() for n in atom.GetNeighbors()] for atom in mol.GetAtoms()}
    expanded_atoms_to_use = set(atoms_to_use)
    [expanded_atoms_to_use.update(graph[idx]) for idx in atoms_to_use]
    subgraphs = search_subgraphs(graph, expanded_atoms_to_use)
    return subgraphs

@timer
def Tanimoto_Sim(mol1, mol2, radius = 2):
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

@timer
def rdkit_MCS_Sim(mol1, mol2):
    res = rdFMCS.FindMCS([mol1, mol2], ringMatchesRingOnly=True,completeRingsOnly=True)
    return (2*res.numAtoms)/(mol1.GetNumAtoms()+mol2.GetNumAtoms())

@timer
def fast_MCS_Sim(mol1, mol2, radius = 1, print_frag = False):
    info1, info2 = get_atom_envs(mol1), get_atom_envs(mol2)
    common_atoms1, common_atoms2 = [], []
    for fp, irs in info1.items():
        if fp in info2 and irs[0][1] == radius:
            common_atoms1 += [ir[0] for ir in irs]
            common_atoms2 += [ir[0] for ir in info2[fp]]
    frags1, frags2 = search_frags_with_atoms(mol1, common_atoms1), search_frags_with_atoms(mol2, common_atoms2)
    if print_frag: print (frags1, frags2)
    max_frag1, max_frag2 = max([len(frag) for frag in frags1]), max([len(frag) for frag in frags2])
    return 2*min([max_frag1, max_frag2])/(mol1.GetNumAtoms()+mol2.GetNumAtoms())

class sim_recorder():
    def __init__(self, smi, topk = 3):
        self.smi = smi
        self.topk = topk
        self.times = []
        self.similarities = []
        
    def record(self, time, similarity):
        self.times.append(time)
        self.similarities.append(similarity)
        
    def summerize(self, smiles_list):
        sorted_similarities = sorted(enumerate(self.similarities), key=lambda x:-x[1])[:self.topk]
        return np.mean(self.times), [(smiles_list[i], s) for i, s in sorted_similarities]
    
        
        
        