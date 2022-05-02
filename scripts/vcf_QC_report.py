#!/usr/bin/env python3
PROGRAM_VERSION = 0.1

# include standard modules
import sys, os, argparse, csv, re

# require cyvcf2 for VCF parsing
from cyvcf2 import VCF

Transitions = ('G>A', 'C>T','A>G','T>C')
Transversions = ('A>C','A>T','C>G','G>T',
                 'C>A','T>A','G>C','T>G')

info="This script will QC VCF file to help identify errors in run or incompleteness."
usage="vcf_QC_report.py --vcf in.vcf --output report.txt --plot report.pdf"

parser = argparse.ArgumentParser(description = info, usage = usage)

parser.add_argument("-V", "--version", help="show program version",
                    action="store_true")

parser.add_argument("--output", "-o", help="Write report output to this file",
                    type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument("--plot", "-p", help="Draw plots to this file")
parser.add_argument("--vcf", "-v", help="snpEff VCF file as bgzip VCF or BCF file")

args = parser.parse_args()

if args.version:
    print("Version %s"%(PROGRAM_VERSION))
    sys.exit(0)

csvout = csv.writer(args.output,delimiter="\t")
stats = {
        'Ts':     0,
        'Tv':     0,
        'SNP':    0,
        'FILTER': {},
        'INDEL':  0,
        'PASS':   0,
        'CpG':    0,
        'nCpG':   0,
        'Total':  0,
        'Length': 0,
        }

# deal with PASS/FILTER?
chromstats = {}
variants = VCF(args.vcf)
for chrom in variants.seqnames:
    chromstats[chrom] = stats.copy()

hdr = variants.raw_header
mtch = re.finditer(r'##contig=\<ID=([^,]+),length=(\d+)',hdr)
for m in mtch:
    chrname = m.group(1)
    chrlen  = int(m.group(2))
    chromstats[chrname]['Length'] = chrlen

for variant in variants:
    #print(variant.gt_types)
    #if not all(variant.gt_types == [VCF.HET, VCF.HOM_REF, VCF.HOM_REF]): continue
    if variant.CHROM not in chromstats:
        chromstats[variant.CHROM] = stats.copy()
    chromstats[variant.CHROM]['Total'] += 1
    if variant.FILTER not in chromstats[variant.CHROM]['FILTER']:
        chromstats[variant.CHROM]['FILTER'][variant.FILTER] = 0
    chromstats[variant.CHROM]['FILTER'][variant.FILTER] += 1
    chromstats[variant.CHROM]['Total'] += 1

    vartype = 'SNP'
    if len(variant.REF) > 1:
        vartype = 'INDEL'
    else:
        for altAllele in variant.ALT:
            if altAllele == '<NON_REF>':
                continue
            elif altAllele != '<NON_REF>' and len(altAllele) > 1:
                vartype = 'INDEL'
                break
            fromTo = "{}>{}".format(variant.REF,altAllele)
            fromTo = fromTo.upper()
            if fromTo in Transitions:
                chromstats[variant.CHROM]['Ts'] += 1
            elif fromTo in Transversions:
                chromstats[variant.CHROM]['Tv'] += 1
#            else:
#                print("unknown allele combo '{}'".format(fromTo))
    chromstats[variant.CHROM][vartype] += 1

allstats = stats.copy()

csvout.writerow(['Chromosome','Length','Total Variants','VarPerKb',
                 'SNP','INDEL','Ts','Tv','Ts/Tv'])
for chrom in sorted(chromstats):
    cstats = chromstats[chrom]
    for t in cstats:
        if type(cstats[t]) is dict:
            for filter,fval in cstats[t].items():
                if filter not in allstats[t]:
                    allstats[t][filter] = fval
                else:
                    allstats[t][filter] += fval
        else:
            allstats[t] += cstats[t]

    ratio = 0
    if cstats['Tv'] > 0:
        ratio = '%.2f'%(cstats['Ts']/cstats['Tv'])

    varPerKb = 0
    if cstats['Length'] > 0:
        varPerKb = "%.2f"%(1000 * cstats['Total'] / cstats['Length'])

    csvout.writerow([chrom,cstats['Length'], cstats['Total'],
                    varPerKb,
                    cstats['SNP'],cstats['INDEL'],
                    cstats['Ts'],cstats['Tv'],ratio])
ratio = 0
if allstats['Tv'] > 0 :
    ratio = "%.2f"%(allstats['Ts'] / allstats['Tv'])

varPerKb = 0
if allstats['Length'] > 0:
    varPerKb = "%.2f"%(1000 * allstats['Total'] / allstats['Length'])

csvout.writerow( ['Totals', allstats['Length'], allstats['Total'],varPerKb,
                  allstats['SNP'], allstats['INDEL'], allstats['Ts'],
                  allstats['Tv'], ratio])
