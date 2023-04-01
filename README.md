# MCS Similarity
Molecular similarity based on fast estimation of maximum common substructure (MCS).

## Definition

Let the MCS of two molecules mol<sub>1</sub> and mol<sub>2</sub> as mol<sub>MCS</sub>, here the similarity is defined as 

$$ MCS\_Sim = {2n(mol_{MCS})\over n(mol_1)+ n(mol_2)} $$

Where n(mol) denotes the number of atoms in the molecule.

## Fast MCS esmimation

Noticing the `rdkit.Chem.rdFMCS.FindMCS` is slow in some cases, such as large molecules or highly similar molecules, I use molecular fingerprint and deep first search (DFS) to quickly estimate the atoms in MCS in three steps.

1. Get the atom environment of each atom in mol1 and mol2 by MorganFingerprint and **get the atoms having same environments**.
2. Apply DFS to make the substructures by connecting the atoms having common environment, and **choose the substructure including the most atoms as MCS**.
3. When the size of MCS obtain from mol1 and mol2 are different, **pick the MCS having smaller size** due to the possibility of same environment on different atoms.

## Experiment
First, I use the products in USPTO_50K dataset and apply [CReM](https://github.com/DrrDom/crem) to generate 5 mutations for each molecule. These mutations are supposed to be "similar" with the parent molecules.

Next, I calculate the similarities between the parent molecules and their mutations and find the top 3 molst similar mutations using 
1. Tanimoto Similarity
2. rdkit MCS similarity
3. fast MCS similarity
4. Hybrid similarity ($\sqrt{(1)*(3)}$)

To reduce the experimental time, I random chose 200 molecules from 50K molecules from USPTO_50K.

## Results

### 1. Similarity

Both RDKit_MCS and fast_MCS similairty show high similarity scores, and Tanimoto shows lower similarity scores between parents mutations.

![](https://i.imgur.com/ZnTIWC3.png)

### 2. Computational time

The computational time using Tanimoto similarity as unit. While rdkit_MCS requires 60K loger time than Tanimoto, fast_MCS only needs around double time of Tanimoto.

![](https://i.imgur.com/yve1aHm.png)















