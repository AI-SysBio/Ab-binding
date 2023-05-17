'''
Created on 6 Feb 2013

@author: dunbar

@change: Added the new version of anarci. This overwrites the old anarci function which is available as anarci_old.
@change: anarci must be intalled separately now as it will be distributed as a separate package. 
'''

# python
import os, sys
import tempfile
import subprocess
import re

from Bio.PDB.Polypeptide import aa1, aa3 # to allow me to return "X" if not found. 
to_one_letter_code = dict(zip(aa3,aa1))

# ANARCI imported as a result of the setup options
from anarci.anarci import number as anarci_algorithm
    
def anarci( seq, scheme="c" ):
    """
    Use the ANARCI program to number the sequence. 
    
    The anarci submodule has more functionality. Here we restrict the program to number only those sequences that are from an antibody. 
    
    @param seq: An amino acid sequence that you wish to number.
    @type seq: C{str} 
    
    @param scheme: The scheme that should be applied. Choose either c, k, a, imgt for chothia, kabat, martin or imgt respectively. 
    @type scheme: C{str} 

    @return: numbering, chain type
    
    """
    numbering , chain_type = anarci_algorithm(seq, scheme=scheme, allow=set( ["H","K","L"]) )
    
    if numbering:
        # replace the Kappa annotation with a light annotation (will come back as L for a lambda chain already).
        chain_type.replace("K","L")
        return [ (_,aa) for _, aa in numbering if aa != "-"] , chain_type
    else:
        return False, False
        

def annotate(chain, method="user", force=False, scheme="c"):
    sequence_list, sequence_str = extract_sequence(chain) 
    numbering, chain_type = anarci(sequence_str, scheme=scheme)
    
    # align the original residue id's to the numbering    
    aligned_numbering = align_numbering(numbering, sequence_list)

    # aligned numbering is a dictionary of the original residue ids and the new numbering
    return aligned_numbering, chain_type       

def extract_sequence(chain,selection=False,return_warnings=False, ignore_hets=False,backbone=False):
    """
    Get the amino acid sequence of the chain.
    @change:    Residues containing HETATOMs are skipped -->  Residues containing HETATOMs are checked as an amino acid.
    
    Residues containing HETATOMs are checked  to be amino acids and the single letter returned.
    
    This works provided the residues in the chain are in the correct order.
    
    @param selection: a selection object to select certain residues
    @param return_warnings: Flag to return a list of warnings or not
    @param backbone: Flag whether to only show residues with a complete backbone (in the structure) or not.
    @return: The sequence in a resid:aa tuple list and the sequence as a string.
    
    """
    sequence_list = []
    warnings=[]
    for residue in chain.get_list():
        if residue.id[0] != " ": # skip HETATOMs - this is not necesserily a good idea, flag to the user that is has been done.
            if residue.get_resname() in to_one_letter_code:
                if ignore_hets:
                    if return_warnings:
                        warnings.append("Warning: HETATM residue %s at position %s (PDB numbering) found in chain %s. Not including it in structure's sequence."%(residue.get_resname(), str(residue.id[1])+residue.id[2].strip(), residue.parent.id ))
                    else:
                        print("Warning: HETATM residue %s position %s (PDB numbering) found in chain %s. Not including it in structure's sequence."%(residue.get_resname(), str(residue.id[1])+residue.id[2].strip(), residue.parent.id ), file=sys.stderr)
                    continue
            else:
                continue
        if selection:
            if not selection.accept(residue): continue

        atoms_of_residue = residue.child_dict.keys()
        backboneCondition = ('N' in atoms_of_residue and 'C' in atoms_of_residue and 'CA' in atoms_of_residue and 'O' in atoms_of_residue) # Boolean to hold if residue has a full backbone

		# CASE 1: backbone = True, and residue has a full backbone; convert a.a into single letter
        if backbone and backboneCondition:
            sequence_list.append( (residue.id, to_one_letter_code.get(residue.get_resname(), 'X') ) )
        # CASE 2: backbone = True, but residue does not have a full backbone; use a gap in sequence annotation
        elif backbone and not backboneCondition:
            sequence_list.append( (residue.id, '-' ) )
        # CASE 0 (default): don't care about backbone, just write it to sequence if it's found in structure.
        elif not backbone:
            sequence_list.append( (residue.id, to_one_letter_code.get(residue.get_resname(), 'X') ) ) # i am 

    sequence_str = "".join([r[1] for r in sequence_list])
    if not return_warnings:
        return sequence_list, sequence_str
    else:
        return sequence_list, sequence_str, warnings
    
def interpret(x):
    """
    Function to interpret an annotation in the form H100A into the form ( 100, 'A' )
    """
    assert x[0] == "H" or x[0] == "L", x
    try:
        return ( int(x[1:]), ' ')
    except ValueError:
        return ( int(x[1:-1]), x[-1] )


def align_numbering(numbering, sequence_list,alignment_dict={}):
    """
    Align the sequence that has been numbered to the sequence you input.
    The numbered sequence should be "in" the input sequence.
    If not, supply an alignment dictionary.(align sequences and use get_alignment_dict(ali1,ali2))
    """
    if numbering:
        numbered_sequence = "".join( [r[1] for r in numbering])
        input_sequence = "".join( [r[1] for r in sequence_list])
        if not alignment_dict:
            if numbered_sequence in input_sequence:
                numbered_sequence_ali,input_sequence_ali=easy_alignment(numbered_sequence, input_sequence)
                alignment_dict = get_alignment_dict(input_sequence_ali,numbered_sequence_ali)
            else:
                raise Exception("Could not align numbered sequence to aligned sequence"+"\n"+str(numbered_sequence)+"\n"+str(input_sequence))

        aligned_numbering = {}
        n=-1
        after_flag=False
        for i in range(len(input_sequence)):
            if i in alignment_dict:
                #during
                assert after_flag is False, "Extra residue in structure than expected from provided sequence"
                assert input_sequence[i] == numbered_sequence[alignment_dict[i]], "alignment dictionary failed"
                aligned_numbering[sequence_list[i][0]] = numbering[alignment_dict[i]][0]
                n = numbering[-1][0][0] +1
            elif n > -1:
                # after
                after_flag=True
                aligned_numbering[sequence_list[i][0]] = (n,' ')
                n+=1
            else:
                # before numbering
                aligned_numbering[sequence_list[i][0]] = '' 

        return aligned_numbering
    else:
        return False


def get_alignment_dict(ali1,ali2):
    """
    Get a dictionary which tells you the index in sequence 2 that should align with the index in sequence 1 (key)
    
    ali1:  ----bcde-f---        seq1: bcdef
    ali2:  ---abcd--f---        seq2: abcdf

    alignment_dict={
        0:1,
        1:2,
        2:3,
        4:4
        }
    
    If the index is aligned with a gap do not include in the dictionary.
    e.g  1 in alignment_dict  --> True
    e.g  3 in alignment_dict  --> False
    """
    assert len(ali1)==len(ali2), "aligned sequences must be same lengths (including gaps)"
    alignment_dict={}
    p1=-1
    p2=-1
    for ap in range( len(ali1) ):
        if ali1[ap] != "-" and ali2[ap] != "-":
            p1+=1
            p2+=1
            alignment_dict[p1] = p2
        elif ali1[ap] != "-": 
            p1+=1 
        elif ali2[ap] != "-": 
            p2+=1    
    return alignment_dict

def easy_alignment(seq1, seq2):
    """
    Function to align two sequences by checking if one is in the other.
    This function will conserve gaps.
    """
    assert type(seq1) is str and type(seq2) is str, "Sequences must be strings for easy_alignment" 
    if seq1 in seq2:
        start = seq2.index(seq1)
        seq1_ali = "-"*start + seq1 + "-"*(len(seq2) - start - len(seq1) )
        return seq1_ali, seq2
    elif seq2 in seq1:
        start = seq1.index(seq2)
        seq2_ali = "-"*start + seq2 + "-"*(len(seq1) - start - len(seq2) )
        return seq1, seq2_ali
    else:
        # Can't align them # I return just one value here. 
        return False
       

def write_pir(seq,filename=""):
    """
    Create a .pir file
    This is the format used as input for abnum
    """
    if not filename:
        pirfd, filename=tempfile.mkstemp('.pir' )
        f = os.fdopen( pirfd, 'w' )
    else:
        f = open(filename,'w')
    # Output to pir file.
    f.write( ">P1;abchain\n" )
    f.write("pir file generated by %s\n"%os.path.split(__file__)[-1])
    f.write( seq )
    f.write("*\n")
    f.close()
    return filename # return the path to the pir file

if __name__ =="__main__":
    pass
    # tests needed

