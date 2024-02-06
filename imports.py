# VAPiD is an extremely lightweight virus genome annotator that takes any number of viral genomes and annotates them
# producing files suitable for NCBI submission

# Vapid Version
import shutil
import time
from Bio import SeqIO
from Bio import Entrez
import sys
import platform
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
import os
import timeit
import argparse
import re
import subprocess

from arg_parse import arg_init, check_os

VERSION = 'v1.6.7'

Entrez.email = 'uwvirongs@gmail.com'

args = arg_init() 
SLASH = check_os()


# quick check to make sure slashes go the right way on both Windows and
# Mac/Linux
def check_os():
    if platform.system() == 'Linux' or platform.system() == 'Darwin':
        return '/'
    else:
        return '\\'

SLASH = check_os()

# Spell checking functionality provided by Entrez
# takes an input string and returns the Entrez corrected string as long as
# it exists
def spell_check(query_string):
    # new entrez rules limit requests to no more than 3 a second, this will
    # ensure we don't send more than two
    time.sleep(0.5)
    handle = Entrez.espell(term=query_string)
    record = Entrez.read(handle)
    corrected_query = record["CorrectedQuery"]
    # Entrez returns blank strings for numerals or things that are spelled correctly
    # Since this is based on NCBI's spell checking protein names are included and correct
    # However this won't correct SUPER messed up words or made up words
    if corrected_query != '':
        # print('Checking spelling on ' + query_string)
        # print(query_string + ' was corrected to: ' + corrected_query)
        return corrected_query
    else:
        return query_string

# takes a strain name and a genome and writes and saves a fasta to the
# correct directory
def write_fasta(strain, genome):
    w = open(strain + SLASH + strain + '.fasta', 'w')
    w.write('>' + strain + '\n')
    w.write(genome)
    w.close()

# Writes an fsa file based of the name, strain and genome, honestly we should allow for much more flexibility
# and automation here
def write_fsa(strain, virus_genome, full_name, nucleic_acid_type):
    if nucleic_acid_type == 'RNA':
        fsa = open(strain + SLASH + strain + '.fsa', 'w')
        fsa.write(
            '>' +
            full_name.strip() +
            ' [moltype=genomic] [molecule=' +
            nucleic_acid_type +
            ']' +
            '\n')
        fsa.write(virus_genome)
        fsa.write('\n')
        fsa.close()
    else:  # nucleic_acid_type =='DNA'
        # table2asn uses "genomic DNA" as the default mol_type
        fsa = open(strain + SLASH + strain + '.fsa', 'w')
        fsa.write('>' + full_name.strip() + '\n')
        fsa.write(virus_genome)
        fsa.write('\n')
        fsa.close()

# this takes in all of our information and makes a feature table that contains correct annotations for for ribosomal slippage and RNA editing
# - as well as creation of a .pep file for rna editing -- Now we also pass two possibly empty lists to write tbl so we can write gene annotations
def write_tbl(
        strain,
        gene_product_list,
        gene_locations,
        genome,
        gene_of_intrest,
        note,
        name_o_vir,
        all_loc_list,
        all_product_list,
        full_name,
        name_of_the_feature_list):
    # covers the nipah situation where there's RNA editing on more than 1 protein - if this happens for more viruses I'll need to code a more
    # robust sollution, but for now this works
    if 'nipah' in name_o_vir.lower():
        pep = open(strain + SLASH + strain + '.pep', 'w')

    tbl = open(strain + SLASH + strain + '.tbl', 'w')
    tbl.write('>Feature ' + full_name)

    # This block should write all gene annotations to tbl file as long as we got passed genes, and the only way that will ever happen is if the
    # User put the -all flag
    if len(all_product_list) > 0:
        for x in range(0, len(all_product_list)):
            print(all_product_list[x] + str(all_loc_list[x]))

            e_flag = ''
            s_flag = ''
            s_all = all_loc_list[x][0]
            e_all = all_loc_list[x][1]
            p_all = all_product_list[x]
            if int(e_all) >= len(genome):
                e_flag = ''
            if int(s_all) < 1:
                s_flag = '<'
                s_all = '1'

            if int(e_all) < 1:
                e_all = len(genome)
                e_flag = '>'
            if p_all == 'inverted terminal repeat':
                if int(s_all) < (len(genome) / 2):
                    s_all = 1
                else:
                    e_all = len(genome)
                tbl.write(
                    '\n' +
                    s_flag +
                    str(s_all) +
                    '\t' +
                    e_flag +
                    str(e_all) +
                    '\trepeat_region\n')
                tbl.write('\t\t\tnote\t' + p_all + '\n')
                tbl.write('\t\t\trpt_type\tinverted')
            else:
                tbl.write(
                    '\n' +
                    s_flag +
                    str(s_all) +
                    '\t' +
                    e_flag +
                    str(e_all) +
                    '\t' +
                    name_of_the_feature_list[x] +
                    '\n')
                if 'UTR' in name_of_the_feature_list[x] or 'gene' in name_of_the_feature_list[x]:
                    feat_des = 'gene'
                elif 'mat_peptide' in name_of_the_feature_list[x] or 'CDS' in name_of_the_feature_list:
                    feat_des = 'product'
                else:
                    feat_des = 'rpt_type'
                tbl.write('\t\t\t' + feat_des + '\t' + p_all)

    for x in range(0, len(gene_product_list)):
        print(gene_product_list[x] + ' ' + str(gene_locations[x]))
        flag = ''
        xtra = ''
        sflag = ''
        product = gene_product_list[x]

        if gene_of_intrest in product:
            xtra = note

        if 'nipah' in name_o_vir.lower():
            nts_of_gene = genome[int(
                gene_locations[x][0]) - 1:int(gene_locations[x][1]) - 1]

            if product.lower() == 'v protein':
                xtra = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds non templated ' \
                    'G\n\t\t\tprotein_id\tn_1' + strain
                start_of_poly_g = nts_of_gene.find('AAAAAGG')
                nts_of_gene = nts_of_gene[0:start_of_poly_g + \
                    1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
                new_translation = str(Seq(nts_of_gene).translate())
                pep.write('>n_1' + strain + '\n' + new_translation)
                pep.write('\n')
            elif product.lower() == 'w protein':
                xtra = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds non templated ' \
                    'G\n\t\t\tprotein_id\tn_2' + strain
                start_of_poly_g = nts_of_gene.find('AAAAAGG')
                nts_of_gene = nts_of_gene[0:start_of_poly_g + \
                    1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
                new_translation = str(Seq(nts_of_gene).translate())
                pep.write('>n_2' + strain + '\n' + new_translation)
                pep.write('\n')

        if 'HIV' in name_o_vir and (
                'Pol polyprotein' == product or 'Pol' == product):
            sflag = '<'

        location_info = gene_locations[x]
        if len(location_info) == 4:

            start_1 = str(location_info[0])
            end_1 = str(location_info[1])
            start_2 = str(location_info[2])
            end_2 = str(location_info[3])

            tbl.write('\n' + start_1 + '\t' + end_1 + '\tCDS\n')
            tbl.write(start_2 + '\t' + end_2 + '\n')
            tbl.write('\t\t\tproduct\t' + product + '\n')

            if 'HEPATITIS B' not in name_o_vir and 'BK polyamavirus' not in name_o_vir:
                tbl.write('\t\t\texception\tRibosomal Slippage\n')

        else:
            start = int(location_info[0])
            end = int(location_info[1])

            it_count = 0
            modifid_orf = False
            # won't execute this block of code for complemented genes
            # print(genome[end - 3:end].upper())
            # this makes sure that our end is in frame, and will adjust hopefully this doesn't break everything
            # added a check to only do this in the absence of RNA editing or
            # ribosomal sliippage
            if xtra == '' and end != len(genome):
                if ((end - start) + 1) % 3 != 0:
                    end += 1
                if ((end - start) + 1) % 3 != 0:
                    end += 1

            if end > start and 'IIIA' not in product.upper():
                if (genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG') and (end < len(
                        genome) - 3) and not re.search('[MRWSYKVHDBN]', genome[end - 3:end].upper()):
                    if re.search('[MRWSYKVHDBN]', genome[end - 3:end].upper()):
                        print(
                            'Ambiguous base detected in a putative stop codon, this can cause problems with VAPiD annotations')
                    print('Modifying ORF length for ' + str(product))
                    end = find_end_stop(genome, start, end)
            # This should now correctly annotate assemblies that come in with the very edges chopped off
            # print(genome[end - 3:end].upper())
            pie = ''
            die = ''
            if int(start) < 1:
                sflag = '<'
                pie = str((int(end) % 3) + 1)
                start = '1'

            if int(end) < 1:
                end = len(genome)
                flag = '>'

            if 'HPIV-1' in name_o_vir or 'human parainfluenza virus 1' in name_o_vir.lower():
                if 'C\'' in product or 'Y2' in product:
                    die = '\n\t\t\ttransl_except\t(pos:' + str(
                        start) + '..' + str(int(start) + 2) + ',aa:Met)'
            tbl.write(
                '\n' +
                sflag +
                str(start) +
                '\t' +
                flag +
                str(end) +
                '\tCDS\n')
            tbl.write('\t\t\tproduct\t' + product + xtra)
            if pie != '':
                tbl.write('\n\t\t\tcodon_start\t' + pie)
            if die != '':
                tbl.write(die)

    tbl.write('\n')
    tbl.close()
    if 'nipah' in name_o_vir.lower():
        pep.close()


# Take the name of a virus sample, and write the .cmt file for it using supplied coverage information
# NOTE: only writes coverage length - so now if we want to say our sequencing platform we have to edit this code
# Now also writes in the comment the reference that this subission was annotated off - this should provide some more
# accountability
def write_cmt(sample_name, coverage, ref_gb, did_we_rc):
    cmt = open(sample_name + SLASH + 'assembly.cmt', 'w')
    cmt.write('##Assembly-Data-START##\n')
    if coverage != '':
        cmt.write('Coverage\t' + coverage + '\n')
    if did_we_rc:
        cmt.write(
            'Original input sequence was reverse complemented by MAFFT during the alignment phase')
    cmt.write(
        'Created with VAPiD' +
        VERSION +
        ' Reference annotations were pulled from ' +
        ref_gb +
        '.\n')
    if args.revica:
        cmt.write('Sequence generated by Revica ' + args.revica +
                  ' (github.com/greninger-lab/revica).' + '\n')
    cmt.write('##Assembly-Data-END##\n')
    cmt.close()

# Build the metadata for every virus that's been submitted
def do_meta_data(strain, sheet_exists, full_name):
    first = True
    s = ''
    coverage = ''
    metadata_sheet_location = ''

    if sheet_exists:
        metadata_sheet_location = args.metadata_loc
        for line in open(metadata_sheet_location):
            if first:
                names = line.split(',')
                first = False
            elif line.split(',')[0] == strain:
                for dex in range(0, len(names)):
                    if names[dex].strip().lower == 'coverage':
                        coverage = line.split(',')[dex].strip()
                    elif names[dex].strip().lower() == 'full_name':
                        if line.split(',')[dex].strip() != '':
                            full_name = line.split(',')[dex].strip()
                    else:
                        s = s + ' [' + names[dex].strip() + '=' + \
                            line.split(',')[dex].strip() + ']'
                break

    if s == '':
        print(
            'metadata not found in provided .csv or .csv not created -  time for minimal manual entry for sequence - ' +
            strain)
        col = ' [collection-date=' + input(
            'Enter collection date in the format (23-Mar-2005, Mar-2005, or 2005): ').strip() + ']'
        con = ' [country=' + \
            input('Enter country sample was collected in (example - USA): ').strip() + ']'
        st = ' [strain=' + \
            input('Enter strain name - if unknown just put ' + strain + ': ').strip() + ']'
        cov = input(
            'Enter coverage as a number (example 42.3), if unknown just leave this blank and hit enter: ')
        meta_data = col + con + st
        coverage = cov
        # Here's one line of code to unilaterally standardize defualt naming
        # scheme
        if full_name == '':
            full_name = strain + \
                ' (' + con.split('=')[1][:-1] + '/' + col.split('=')[1][:-1] + ')'

    else:
        meta_data = s
        if full_name == '':
            full_name = strain
            print('Automatic strain naming failed but submission will proceed without metadata appended to the fasta header.')
    return meta_data, coverage, full_name

# create an sbt file to include DBlink info in the work directory and
# return the full path of that sbt file.
def create_new_sbt(sbt, dblink_meta_dict, strain):
    if dblink_meta_dict['bioproject'] is not None or dblink_meta_dict['biosample'] is not None or dblink_meta_dict['sra'] is not None:
        sbt_content = open(sbt).read()
        sbt_add_text = "Seqdesc ::= user {type str \"DBLink\",data {"
        if dblink_meta_dict['bioproject'] is not None:
            sbt_add_text += "{label str \"BioProject\",num 1,data strs {\"" + \
                dblink_meta_dict['bioproject'] + "\"}},"
        if dblink_meta_dict['biosample'] is not None:
            sbt_add_text += "{label str \"BioSample\",num 1,data strs {\"" + \
                dblink_meta_dict['biosample'] + "\"}},"
        if dblink_meta_dict['sra'] is not None:
            sbt_add_text += "{label str \"Sequence Read Archive\",num 1,data strs {\"" + \
                dblink_meta_dict['sra'] + "\"}},"
        sbt_add_text += "}}\n"
        sbt_content += sbt_add_text

        if not os.path.isdir("work"):
            os.mkdir("work")
        file_name = "work" + SLASH + strain + ".sbt"
        new_sbt_file = open(file_name, "w")
        new_sbt_file.write(sbt_content)
        new_sbt_file.close()

        return os.path.abspath(file_name)


# Takes in two sequences with gaps inserted inside of them and returns arrays that have a -1 in the gap locations and
# count up from 1 in the nucleotide areas - This data structure allows for extremely rapid conversion between relative
# locations in the two sequences although does assume that these genes are of uniform length
# NOTE: This means that when we have reads that like don't have the start codons of the first gene or something we'll
# get a -1 for the start location on our annotation
def build_num_arrays(our_seq, ref_seq):
    ref_count = 0
    our_count = 0
    ref_num_array = []
    our_num_array = []

    for x in range(0, len(ref_seq)):
        if ref_seq[x] != '-':
            ref_count += 1
            ref_num_array.append(ref_count)
        else:
            ref_num_array.append(-1)

        if our_seq[x] != '-':
            our_count += 1
            our_num_array.append(our_count)
        else:
            our_num_array.append(-1)

    return our_num_array, ref_num_array

# Takes a nucleotide sequence and a start and end position [1 indexed] and search for a stop codon from the start
# to the end + 60 so every codon in the provided gene and then 3 after it. Return the first stop codon found or if no
# stop codon is found return the original end value and print a warning
def find_end_stop(genome, start, end):
    # save the provided end
    start -= 1
    old_end = end
    end = start + 3
    # Search for stop codons in DNA and RNA space until 3 codons after the provided end.
    # Turns out 3 codons isn't enough
    while genome[end -
                 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG' and end <= (old_end +
                                                                             60):
        end += 3
    if end == old_end + 60:
        # print('WARNING no stop codon found, examine reference and original sequence')
        return old_end
    else:
        return end

