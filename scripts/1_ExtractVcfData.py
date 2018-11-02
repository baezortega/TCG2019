#!/usr/bin/env python

# Analysis scripts for Baez-Ortega et al., 2019
# Step 1. Extract variant information from VCF files

# Adrian Baez-Ortega, 2018


# 1_ExtractVcfData.py
# Extracts metadata, NR and NV values from a Platypus VCF (or gzipped VCF) into three text files.
#
# ARGUMENTS
#  Path to VCF or gzipped VCF file (.vcf or .vcf.gz)
#  Path to output directory
#
# USAGE
#  1_ExtractVcfData.py input_variants.vcf[.gz] path/to/output_dir

# This script was written for Python 2.7 on Unix; it may not work on Python 3 or Windows.


"""
This script is used to extract the metadata (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO), 
the total number of reads (NR) and the number of reads supporting the variant (NV), for 
every variant in a VCF or gzipped VCF file (with Platypus output format), into three text files. 
"""

import sys
import gzip
import os


# If not 2 arguments: print help
if len(sys.argv) != 3:
    print '\n1_ExtractVcfData.py: Extracts metadata (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO), total number'
    print '                     of reads (NR), and number of reads supporting the variant (NV), for every variant'
    print '                     in a VCF or gzipped VCF file (with Platypus output format), into three text files.'
    print '          Arguments: Path to a VCF or gzipped VCF file (.vcf[.gz] extension, with Platypus output format).'
    print '                     Path to output directory.'
    print '              Usage: 1_ExtractVcfData.py input_variants.vcf[.gz] path/to/output_dir\n'
    sys.exit(0)


script, vcf_file, out_dir = sys.argv


# Helper function: open a VCF or gzipped VCF file
def open_vcf(vcf_path, gzipped):
    if gzipped:
        return gzip.open(vcf_path, 'r')
    else:
        return open(vcf_path, 'r')


# Compose output file paths
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

if vcf_file[-7:] == '.vcf.gz':
    gzipped = True
    ext_len = 7
elif vcf_file[-4:] == '.vcf':
    gzipped = False
    ext_len = 4
else:
    print '\nERROR: Invalid file extension. Only extensions .vcf and .vcf.gz are accepted.\n'
    sys.exit(1)

out_file_NR = os.path.join(out_dir, os.path.basename(vcf_file)[:-ext_len] + '_NR.tsv')
out_file_NV = os.path.join(out_dir, os.path.basename(vcf_file)[:-ext_len] + '_NV.tsv')
out_file_MD = os.path.join(out_dir, os.path.basename(vcf_file)[:-ext_len] + '_Metadata.tsv')

print '\nInput file:          ', vcf_file
print 'Output NR file:      ', out_file_NR
print 'Output NV file:      ', out_file_NV
print 'Output metadata file:', out_file_MD


# Process VCF
with open_vcf(vcf_file, gzipped) as vcf, open(out_file_NR, 'w') as out_NR, open(out_file_NV, 'w') as out_NV, open(out_file_MD, 'w') as out_MD:
    
    for line in vcf:
        
        # Skip VCF header
        if not line.startswith('##'):
            
            # Write column names to output files
            if line.startswith('#CHROM'):
                MD_head = line[1:].strip().split('\t')[:8]
                NR_head = line.strip().split('\t')[9:]
                out_MD.write('\t'.join(MD_head) + '\n')
                out_NR.write('\t'.join(NR_head) + '\n')
                out_NV.write('\t'.join(NR_head) + '\n')
            
            else:
                # Extract metadata
                metadata = line.strip().split('\t')[:8]
                out_MD.write('\t'.join(metadata) + '\n')
            
                # Extract NR and NV from sample data
                NR = []
                NV = []
                data = line.strip().split('\t')[9:]
            
                for record in data:
                    # Extract total reads (nr) and supp. reads (nv)
                    nr = record.split(':')[4]
                    nv = record.split(':')[5]
                    # Add coverage and VAF values to list
                    NR.append(nr)
                    NV.append(nv)
                
                out_NR.write('\t'.join(NR) + '\n')
                out_NV.write('\t'.join(NV) + '\n')
               
print '\nDone\n'

