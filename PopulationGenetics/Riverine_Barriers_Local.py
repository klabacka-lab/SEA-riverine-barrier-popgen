import itertools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import tempfile
import pathlib
from pathlib import Path

#grab input fasta... to be replaced with a loop through all of the files later
directory = Path("/Users/christianfurness/Desktop/Meta-Analysis/Trofimets et al., 2024/Edited_Fastas")
for input_file in directory.glob("*.fasta"):
    #make a dictionary of the IDs and sequences from fasta
    IDs = {}
    for record in SeqIO.parse(input_file.resolve().as_posix(), "fasta"):
        IDs[record.id]=record.seq

    #import R packages
    hierfstat = importr("hierfstat")
    ape = importr("ape")
    adegenet = importr("adegenet")
    rownames = r['rownames']
    sub = r['sub']
    table = r['table']

    #function to calculate FST values
    def calc_FST(file_input, group1_char, group2_char):
            dna = ape.read_dna(file_input, format="fasta")
            hdrs = rownames(dna)
            pop_code = sub(".*([A-Za-z])$", "\\1", hdrs)
            gen = adegenet.DNAbin2genind(dna)
            pop_factor = r['factor'](pop_code)
            gen = r['pop<-'](gen, pop_factor)
            hf = hierfstat.genind2hierfstat(gen)
            pwFst = hierfstat.pairwise_WCfst(hf)
            FST_value = list(pwFst.rx(group1_char, group2_char))[0]
            return FST_value


    #generate a true FST value for the original grouping of the data (E and W)
    true_FST = calc_FST(input_file.resolve().as_posix(), "E", "W")
    output = f"True FST: {true_FST}\n"

    #function that returns a list of all possible groupings split into two groups
    def permutationmaker(IDlist):
        groupings = []
        length = len(IDlist)
        for i in range(2, length//2 + 1):
            for group in itertools.combinations(IDlist, i):
                group1 = set(group)
                group2 = IDlist - group1
                groupings.append([group1, group2])
        return groupings

    #create a list of all possible groupings of the data into two groups
    Possible_combos = permutationmaker(set(IDs.keys()))

    #output variables
    greater_than_true = 0
    greater_than_list = []

    #create list of fasta data for each grouping
    for groupings in Possible_combos:
        fastaseq = []
        for ID in groupings[0]:
            newrecord = SeqRecord(
                IDs[ID],
                id=ID + "_X",
                description=""
            )
            fastaseq.append(newrecord)
        for ID in groupings[1]:
            newrecord = SeqRecord(
                IDs[ID],
                id=ID + "_Y",
                description=""
            )
            fastaseq.append(newrecord)

        #for each grouping, make a temp fasta and run the fst calculation on it
        with tempfile.NamedTemporaryFile("w", suffix = ".fasta") as tmp:
            tmp.close()

            SeqIO.write(fastaseq, tmp.name, "fasta")

            #calculate FST for possible permutation
            FST_value = calc_FST(tmp.name, "X", "Y")

            #compare with true FST
            if FST_value > true_FST:
                greater_than_true += 1
                greater_than_list.append((groupings[0], groupings[1], FST_value))


    #create output string
    output += f"# permutations with FST > true FST: {greater_than_true}\n"
    output += f"Total permutations: {len(Possible_combos)}\n"
    output += f"p-value: {greater_than_true/len(Possible_combos)}\n\n"
    output += "Group1\tGroup2\tFST\n"

    for group1, group2, FST_value in greater_than_list:
        group1_str = ",".join(sorted(group1))
        group2_str = ",".join(sorted(group2))
        output += f"{group1_str}\t|\t{group2_str}\t|\t{FST_value}\n"

    #write output to file
    with open(f"Output_{input_file.name[0:-6]}.txt", "w") as out_file:
        out_file.write(output)



