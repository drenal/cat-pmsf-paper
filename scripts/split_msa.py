#!/usr/bin/env python3
"""
Split phylip or fasta MSAs into specified number of chunks

Author:
- Lenard Szantho <lenard@drenal.eu>

Version:
- 0.1 (2020-08-30)
"""

import argparse
import re
import math

from itertools import tee, islice, zip_longest

#
# Helper functions
#

# from stackoverflow.com/questions/4197805
def get_next(iterable, window=1):
    """Function to look ahead of @param window in an @param iterable

    """
    items, nexts = tee(iterable, 2)
    nexts = islice(nexts, window, None)
    return zip_longest(items,nexts)

# from stackoverflow/question/18854620
def chunkstring(string, length, offset=0):
    return (string[offset+i: length+i] for i in range(0, len(string), length))


#
# Input handlers
#
def parse_fasta(inputfile):
    """
    Parse fasta formatted MSA
    
    Parameters:
    -----------
    - inputfile: path to file to read
    
    Returns:
    --------
    - dict of sequence ID and sequences
    """
    sequences={}

    with open(inputfile, 'r') as finput:
        for i,line in enumerate(finput):
            if line.strip():
                if line.strip()[0] == '>':
                    gene_name=line.strip('>').strip()
                    sequences[gene_name] = ''
                else:
                    sequences[gene_name] += line.strip()

    return sequences



def parse_phylip(inputfile):
    """
    Parse phylip formatted MSA
    Special care is taken for different 
    non-standard but widely used variations:
    - interleaved: the parts of a sequence is
                    split in multiple lines
                    preserving the original
                    order of sequence IDs
    - sequential: the parts of a sequence is 
                    following continously 
    
    Parameters:
    -----------
    - inputfile: path to file to read
    
    Returns:
    --------
    - (array of sequences, array of sequence IDs) 
    """
    reported_genes_num = 0
    reported_seq_length = 0

    sequences={}
    genes=[]

    interleaved_sequences={}
    interleaved_genes =[]

    with open(inputfile, 'r') as finput:
        first_line = True
        first_block = True
        reading_window = True
        interleaved_counter = 0
        last_line = False
        # we suppose the phylip is a sequential
        # and not an interleaved one 
        sequential = True
        line_parts_num = 0

        for line,next_line in get_next(finput):
            clean_line = line.strip()
            line_parts = re.split(r"\s+", clean_line)

            clean_next_line = ''
            next_line_splited = []
            next_line_parts = 0

            if not next_line is None:
                clean_next_line = next_line.strip()
                next_line_splited = re.split(r"\s+", clean_next_line)
                next_line_parts = sum([len(l) for l in next_line_splited])
            else:
                last_line = True
                
            curr_line_parts = sum([len(l) for l in line_parts])

            # Parsing phylip header
            if first_line:
                #print("Line %s is a header" % clean_line)
                first_line = False
                if len(line_parts) != 2:
                    print("Error: phylip header not found or wrong in file %s" % inputfile)
                    print("Line is: %s" % line_parts)
                    exit(5)
                reported_genes_num = int(line_parts[0])
                reported_seq_length = int(line_parts[1])
                continue

            # Empty lines are only in the interleaved format
            if curr_line_parts == 0:
                first_block = False
                sequential = False
                interleaved_counter = 0
                # pointer copy is enough as we don't want to 
                # access the old, false data
                genes = interleaved_genes
                sequences = interleaved_sequences
                if len(genes) != reported_genes_num:
                    print("Error: we left the first block but didn't yet find all the genes in file %s" % inputfile)
                    print("Genes: %s" % genes)
                    exit(4)
                continue
            
            # If format is interleaved and we're after the first block
            # the job is simple
            if not sequential:
                #print("Added line %s to %s" % (clean_line, genes[interleaved_counter]))
                sequences[genes[interleaved_counter]] += ''.join(line_parts)
                interleaved_counter += 1
                continue

            # Collect data as it would be an interleaved format
            # as long as the genes are not yet collected
            if len(genes) < reported_genes_num:
                interleaved_genes.append(line_parts[0])
                interleaved_sequences.update({line_parts[0]: ''.join(line_parts[1:])})

            # Special case: one-line-sequences
            if reported_seq_length == sum([len(part) for part in line_parts[1:]]):
                #print("Line %s is a one-liner" % line)
                gene_name = line_parts[0]
                sequences.update({gene_name: ''.join(line_parts[1:])})
                genes.append(gene_name)
                reading_window=False
                gene_name = ''
                continue

            # Normal case, the sequences are at least two lines long
            # Normal case A: the sequence IDs are in a separate line
            # Normal case B: the sequences begin right after the sequence ID


            # Bootstrap step, first line
            if line_parts_num == 0:
                #print("Line %s is the firs data line" % clean_line)
                # Bootstrap value, we are on the first data line
                # thus the first part has to be a seqID
                gene_name = line_parts[0]
                sequences.update({gene_name: ''})
                genes.append(gene_name)
                if curr_line_parts > 1:
                    sequences[gene_name]+=''.join(line_parts[1:])
                reading_window = True
                line_parts_num = curr_line_parts
                continue

            if last_line and reading_window:
                #print("Adding last line %s to %s " % (clean_line, gene_name))
                sequences[gene_name] += ''.join(line_parts)
                reading_window=False
                gene_name=''
                line_parts_num = curr_line_parts
                continue

            print("Curr line {}, next line {}, parts {}".format(curr_line_parts, next_line_parts, line_parts_num))

            if not reading_window and gene_name == '':
                #print("New sequence id found on line %s" % clean_line)
                gene_name = line_parts[0]
                sequences.update({gene_name: ''})
                genes.append(gene_name)
                if curr_line_parts > 1:
                    sequences[gene_name] += ''.join(line_parts[1:])
                reading_window=True
                line_parts_num = curr_line_parts
            elif curr_line_parts > line_parts_num:
                #print("Line %s is longer than previous" % clean_line)
                line_parts_num = curr_line_parts
                # here is something more than before
                # - maybe an ID with sequence blocks next to it
                # - maybe the beginning of a sequence after a oneliner ID
                if next_line_parts == curr_line_parts:
                    # this is a sequence
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=True
                elif next_line_parts > curr_line_parts:
                    # next line is a new ID, this is the last line of the sequence
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=False
                    gene_name=''
                elif len(next_line_splited) == 1:
                    # next line is a seqID line without sequences, this is the end of current seq
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=False
                    gene_name=''
                else:
                    # next line will be the end of this sequence
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=True
            elif curr_line_parts < line_parts_num:
                #print("Line %s is shorter than previous" % clean_line)
                line_parts_num = curr_line_parts
                # here is something less than before
                # - maybe the end of a sequence
                # - maybe the next line after an ID
                # - maxbe an ID on a single line
                if next_line_parts == curr_line_parts:
                    # sequence continues
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=True
                elif next_line_parts > curr_line_parts:
                    # next line is a new ID, this is the last line of the sequence
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=False
                    gene_name=''
                elif len(next_line_splited) == 1:
                    # next line is a seqID line without sequences, this is the end of current seq
                    sequences[gene_name]+=''.join(line_parts)
                    reading_window=False
                    gene_name=''
                else:
                    #print(next_line_parts)
                    #print(curr_line_parts)
                    # new ID here without sequence
                    # ID's length is smaller than the sequence length
                    print("Error: This line is shorter than previous, and next is even shorter in file %s, this is not possible." % inputfile)
                    exit(4)
                    #gene_name = line_parts[0]
                    #reading_window=True
                    #sequences.update({gene_name: ''})
                    #genes.append(gene_name)
            else:
                #print("Line %s is the same length as before" % clean_line)
                # - maybe a continous sequence
                sequences[gene_name] += ''.join(line_parts)
                if reported_seq_length == sum([len(part) for part in line_parts]):
                    reading_window = False
                    gene_name = ''
        
        # Some checks
        if len(genes) != reported_genes_num:
            print("Error: the collected genes number %s is not the same as the stated one %s in phylip header in file %s" % (len(genes), reported_genes_num, inputfile))
            print("Collected genes: %s" % genes)
            exit(6)

        for seq_key in sequences:
            if len(sequences[seq_key]) != reported_seq_length:
                print("Error: the sequence length %s doesn't match the stated one %s in phylip header in file %s for sequence id %s" % (len(sequences[seq_key]), reported_seq_length, inputfile, seq_key))
                print("Sequence: %s" % sequences)
                exit(7)

    return sequences

#
# Output handlers
#
def output_as_nexus(sequences, out, order_of_families=[], seq_length_of_families=[]):
    with open(out+".nexus", "w") as output:
        output.write("#nexus\n")

        output.write("begin data;\n")
        output.write("dimensions ntax=%s nchar=%s;\n" % (len(sequences), len(list(sequences.values())[0]) ) )
        output.write("format datatype=protein missing=? gap=-;\n")
        output.write("matrix\n")

        for seq in sequences:
            output.write("%s  %s" % (seq, sequences[seq]))
            output.write('\n')

        output.write("end;\n")
        
        
def output_as_fasta(sequences, out):
    with open(out+".fasta", "w") as output:
        for seq in sequences:
            chunks = chunkstring(sequences[seq], 71)
            output.write(">%s\n" % seq)
            for chunk in chunks:
                output.write(chunk+'\n')
            output.write('\n')

def output_as_phylip(sequences, out):
    max_length = 0
    chunked_sequences = {}

    for seq in list(sequences.keys()):
        length = len(sequences[seq])
        if length > max_length:
            if max_length != 0:
                print("Warning: the %s concatenated sequence is not the same length as the others. It's %s long instead of %s" % (seq, length, max_length))
            max_length = length
        #if len(seq) > 10:
            #new_seq_name = seq[:10]
            #sequences[new_seq_name] = sequences.pop(seq)
            #print("Notice: The sequence id %s is too long for phylip format, it was reduced to %s" % (seq, new_seq_name))
    # Optimize line length
    optimal_line_length = 0
    min_remainder = 50
    #for i in reversed(range(40,50)):
    #    if max_length % i < min_remainder:
    #        min_remainder = max_length % i
    #        optimal_line_length = i
    #        if min_remainder == 0:
    #            break
    # Do not optimize, simply 5 columns and the rest
    if max_length < 50:
        optimal_line_length = max_length
        min_remainder = 0
    else:
        optimal_line_length = 50
        min_remainder = max_length % 50
        

    with open(out+".phylip", "w") as output:
        output.write("%s   %s" % (len(sequences), max_length))

        for i in range(0, math.floor(max_length/optimal_line_length)):
            output.write('\n')

            for seq in sequences:
                if i == 0:
                    #output.write(seq + ' '*(12-len(seq)-1))
                    output.write(seq + ' ')
                else:
                    output.write(' '*11)

                counter = 0
                for j in range(i*optimal_line_length, i*optimal_line_length + optimal_line_length):
                    if counter != 0 and counter % 10 == 0:
                        output.write(' ')
                    output.write(sequences[seq][j])
                    counter += 1
                output.write('\n')

        if min_remainder != 0:
            output.write('\n')
            
            for seq in sequences:
                output.write(' '*11)
                
                for j in range(max_length - min_remainder, max_length):
                    if j != max_length-min_remainder and j % 10 == 0:
                        output.write(' ')
                    output.write(sequences[seq][j])
                output.write('\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-i', '--input', required=True, 
            help="path to input MSA")
    parser.add_argument('-t', '--input-type', required=True, choices=["phylip", "nexus", "fasta"],
            help="format of the input file")
    parser.add_argument('-o', '--output', required=True, help="output file name base (without extension)")
    parser.add_argument('-s', '--output-type', required=True, choices=["phylip", "nexus", "fasta"], help="output file format")
    parser.add_argument('-p', '--parts', required=True, help="Number of parts to split the file", type=int)
    
    args = parser.parse_args()
    
    # read in sequences
    sequences = {}

    if args.input_type == "fasta":
        sequences = parse_fasta(args.input)
    elif args.input_type == "phylip":
        sequences = parse_phylip(args.input)
    else:
        print("Nexus input not yet supported")
        exit(-1)

    # calculate how many splits will be made
    seq_len = len(sequences[list(sequences.keys())[0]])
    
    split_len = int(seq_len / args.parts)
    if seq_len % args.parts != 0:
        split_len += 1
        
    # get file path without extension
    #basename="".join(args.output.split(".")[:-1])
    #print(basename)
    basename = args.output
    
    for i in range(0,args.parts):
        seq_split = {}
        
        if (i+1)*split_len > seq_len:
            for seq_k in sequences:
                seq_split[seq_k] = sequences[seq_k][i*split_len:] 
        else:            
            for seq_k in sequences:
                seq_split[seq_k] = sequences[seq_k][i*split_len:(i+1)*split_len]
        
        if args.output_type == 'fasta':
            output_as_fasta(seq_split, basename+"_"+str(i))
        elif args.output_type == "nexus":
            output_as_nexus(seq_split, basename+"_"+str(i))
        else:
            output_as_phylip(seq_split, basename+"_"+str(i))    
    
    
if __name__ == "__main__":
    main()
