import argparse
import platform
import sys

# This file creates values needed for global variables in vapid4.py and imports.py. 
## Arguments from argparse
## SLASH from check_os()

VERSION = 'v1.6.7'

#### return args to be set globally in other scripts
def arg_init(): 
    parser = argparse.ArgumentParser(
        description='Version ' +
        VERSION +
        '\nPrepares FASTA file for NCBI Genbank submission '
        'through local or online blastn-based annotation of viral sequences. '
        'In default mode, VAPiD searches this folder for our viral databases.')
    parser.add_argument(
        'fasta_file',
        help='Input file in .fasta format containing complete or near complete '
        'genomes for all the viruses that you want to have annotated')
    # parser.add_argument(
    #         '--work_dir',
    #         help='select directory where input .fasta is located')
    parser.add_argument(
        'author_template_file_loc',
        help='File path for the NCBI-provided sequence author template file'
        ' (should have a .sbt extension)\n https://submit.ncbi.nlm.nih.gov/genbank/template/submission/')

    ### DEPRECATING THIS FEATURE


    # parser.add_argument(
    #     '--metadata_loc',
    #     help='If you\'ve input the metadata in the provided csv, specify the location '
    #     'with this optional argument. Otherwise all metadata will be manually prompted for.')

    ####
    
    parser.add_argument(
        '--src_file', help='A tab-delimited source modifiers table (suffix .src). '
        'https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html')
    parser.add_argument(
        '--dblink_metadata_loc',
        help='Metadata csv file containing bioproject, biosample, or sra accessions.')
    parser.add_argument(
        '--r', help='If you want to specify a specific NCBI reference, put the accession number here '
        '- must be the exact accession number - note: feature forces all sequences in FASTA to be this viral species.')
    parser.add_argument(
        '--f',
        help='specify a custom gbf file that you would like to annotate off of')
    parser.add_argument(
        '--db',
        help='specify the local blast database name.  You MUST have blast+ with blastn'
        'installed correctly on your system path for this to work.')
    parser.add_argument(
        '--online',
        action='store_true',
        help='Force VAPiD to blast against online database.  This is good for machines that don\'t '
        'have blast+ installed or if the virus is really strange.'
        'Warning: this can be EXTREMELY slow, up to ~5-25 minutes a virus')
    parser.add_argument(
        '--spell_check',
        action='store_true',
        help='Turn on spellchecking for protein annoations ')
    parser.add_argument(
        '--all',
        action='store_true',
        help='Use this flag to transfer ALL annotations from reference, this is largely untested')
    # parser.add_argument(
    #     '--slashes', action='store_true', help='Use this flag to allow any characters in the name of your virus - This allows '
    #     'you to submit with a fasta file formated like >Sample1 (Human/USA/2016/A) Complete CDS'
    #     ' make sure that your metadata file only contains the first part of your name \'Sample1\' in the example above. '
    #     'You can also submit names with slashes by specifying in the metadata sheet under the header full_name, if you do that '
    #     'you do not need to use this flag')
    parser.add_argument(
        '--dna',
        action='store_true',
        help='Make all files annotated by this run be marked as DNA instead of the default (RNA)')
    parser.add_argument(
        '--output_loc',
        help='Specifies an output location for all files.')
    parser.add_argument(
            '--align_dir',
            help='Specifies an output location for fasta header-specific working files')
    parser.add_argument(
        '--revica',
        help='Add Revica version used for consensus calling and github url.')

    try: 
        return parser.parse_args()
    except BaseException: 
        parser.print_help()
        sys.exit(0)

#return appropriate slash to use, to be set to SLASH in main scripts
def check_os():
    if platform.system() == 'Linux' or platform.system() == 'Darwin':
        return '/'
    else:
        return '\\'


