import sys
import allel
import numpy as np
from Bio import SeqIO
from itertools import combinations

#converts an aligned fasta file to a sci-kit genotype array
def fasta_to_genotype(fasta):
    records = list(SeqIO.parse(fasta, "fasta"))
    ids = [r.id for r in records]
    seqs = np.array([list(str(r.seq).upper()) for r in records])

    BASES = np.array(list("ACGT"))
    rows = []
    for col in seqs.T:
        valid = np.isin(col, BASES)
        alleles, counts = np.unique(col[valid], return_counts=True)
        if len(alleles) != 2:
            continue
        major, minor = alleles[np.argsort(-counts)]
        geno = np.full(len(col), -1, dtype="int8")
        geno[col == major] = 0
        geno[col == minor] = 1
        rows.append(geno)

    geno = np.array(rows)
    g = np.stack([geno, geno], axis=2)
    return ids, allel.GenotypeArray(g)

#Function to calculate fst value given a genotype array and a list of each indice of the populations
def calc_fst(g, group1, group2):
    a, b, c = allel.weir_cockerham_fst(g, [group1, group2], max_allele=None, blen=None)
    fst = np.nansum(a) / (np.nansum(a) + np.nansum(b) + np.nansum(c))
    return fst

#Function to generate a list of possible groupings using X and Y as group identifiers
def permutation_maker(sample_ids):
    n = len(sample_ids)
    indexes = set(range(n))
    for size in range(2, n // 2 + 1):
        for x in combinations(range(n), size):
            if size * 2 == n and 0 not in x:
                continue
            group1 = list(x)
            group2 = sorted(indexes - set(x))
            yield group1, group2

#Store the fasta file
fasta_file = sys.argv[1]

#generate a list of ids and the genotype array
ids, g = fasta_to_genotype(fasta_file)

#calculate the fst for the true groupings
east = [i for i, name in enumerate(ids) if name.endswith("E")]
west = [i for i, name in enumerate(ids) if name.endswith("W")]
true_fst = calc_fst(g, east, west)

num_extremes = 0
num_perms = 0

#for each possible permuation calculate the fst to see if greater than true value
for group1, group2 in permutation_maker(ids):
    num_perms += 1
    fst = calc_fst(g, group1, group2)
    
    if fst >= true_fst:
        num_extremes += 1

#Calculate a p-value from the fsts
p_value = (num_extremes / num_perms)

#Create output file
output = f"True FST: {true_fst}\n"
output += f"# permutations with FST >= true FST: {num_extremes}\n"
output += f"p-value: {p_value}\n"

with open(f"fst_output.txt", "w") as outfile:
    outfile.write(output)