#Hello! If you're reading this, I want you to know that I really enjoyed working on this project
#and I find this subject matter realy interesting. I used this as an opportunity to learn some more about python and numpy, 
#neither of which are languages very familiar to me. I think this may make some of this code somewhat
#more bulky and less pretty/efficient than I would have made if it were in another language, 
#bug I still think it does what it's supposed to do. 

import numpy as np
import time as t
from global_alignment import *
from local_alignment import *
from sub_matrix import *
from sub import *


def generate_sub_matrix(filename):
    fn = np.loadtxt(open("sub/"+filename), dtype=np.str, delimiter=",")
    return sub_matrix(fn)

def generate_pam_n(n):
    smat = generate_sub_matrix("sub/pam1.txt")
    diagarr = []
    pw = smat.sub_scores.astype(float) / 10000
    pw = np.linalg.matrix_power(pw,n) * 100

    freq_array = np.loadtxt(open("sub/freq.txt"), dtype=np.str, delimiter="\n")
    
    
    print(freq_array[0].astype(float))

    for idx,row in enumerate(pw):
            row = (row / freq_array.astype(float)[idx])

    print(np.log10(row))
    
    return pw

def read_file(filename):
    f = open("seq/"+filename, "r") 
    seq = f.read()    
    seq = seq.replace("\n", "")  
    seq = seq.replace("\r", "") 
    return seq

def read_file_as_protein(filename):
    f = read_file("seq/"+filename)
    p = translate(f)
    return p


def print_formatted(tbl,indent):
    for row in tbl:
        print("".join('{:^{indent}}'.format(str(elem),indent=indent) for elem in row))


def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X', 
        'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    print_formatted(protein,3)
    return protein 
def print_seq(arr1,arr2):
    match = np.array([" "] * len(arr1))
    match[np.where(arr1!=arr2) and np.where(arr1!="-") and np.where(arr2!="-")] = "."
    match[np.where(arr1=="-") or np.where(arr2=="-")]=" "
    match[np.where(np.char.equal(arr1,arr2))] = "|"
    print_arr = [arr1,match,arr2]
    for x in range(0,((len(arr1)%50))):
        for y in range(0,3):
            print("".join('{:^1}'.format(str(elem)) for elem in print_arr[y][(50*x):50*(x)+50]))
    return print_arr


def main():
    
    print("Hello!")
    print("First, let's pick a substitution matrix.")
    sm_fn = input("please enter a filename: ")
    smat = generate_sub_matrix(sm_fn)
    print(smat.sub_scores)
    nuc_seq = input("Is this a nucleotide sequence? Y/N: ")
    filename1 = input("Please enter filename of first sequence: ")
    if nuc_seq == "Y" or nuc_seq == "y":
        cseq = read_file_as_protein(filename1)
    else:
        cseq = read_file(filename1)
        
    nuc_seq2 = input("Is this a nucleotide sequence? Y/N: ")
    filename2 = input("Please enter filename of second sequence: ")
    if nuc_seq2 == "Y" or nuc_seq2 == "y":
        rseq = read_file_as_protein(filename2)
    else:
        rseq = read_file(filename2)
    print("If this is a global alignment, enter 1")
    print("If this is a semiglobal alignment, enter 2")
    print("If this is a local alignment, enter 3")
    aln = input("Please enter 1, 2, or 3: ")
    gap_open = input("Please enter a gap score: ")
    if aln=="1":
        alignment = global_alignment(smat, cseq, rseq, False, int(gap_open))
    elif aln=="2":
        alignment = global_alignment(smat, cseq, rseq, True, int(gap_open))
    else:
        alignment = local_alignment(smat,cseq,rseq,int(gap_open))

    print_formatted(alignment.out,4)
    print_formatted(alignment.dirct,3)
    print_seq(alignment.arr1,alignment.arr2)
    
if __name__ == "__main__":
    main()

    