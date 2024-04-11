import argparse
import pysam
from typing import List
from .models import GenomePosition,GenomeRange
from collections import Counter
from tqdm import tqdm
import tabulate
import json
from . import __version__

def get_alleles(
    read:pysam.libcalignedsegment.AlignedSegment,
    positions=List[GenomePosition]
):

    alleles = {}
    for read_pos, ref_pos, read_nt in read.get_aligned_pairs(with_seq=True):
        if read_nt is None:
            continue
        if read_pos is None:
            continue
        p = GenomePosition(chrom=read.reference_name,pos=ref_pos+1)
        if p not in positions:
            continue
        if read.query_qualities[read_pos] < 12:
            continue
        if read_nt.islower():
            read_nt = read.query_sequence[read_pos].upper()
        
        alleles[p] = read_nt

    return [alleles.get(p) for p in positions]

def get_haplotype_counts(
    bam:pysam.AlignmentFile,
    positions:List[GenomePosition]
):
    allele_combinations = []
    chrom = positions[0].chrom
    start = min(positions).pos
    end = max(positions).pos
    for read in bam.fetch(contig=chrom,start=start,end=end):
        if read.mapping_quality < 10:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.is_unmapped:
            continue
        alleles = get_alleles(read,positions=positions)
        if all(alleles):
            allele_combinations.append(''.join(alleles))
    counts = Counter(allele_combinations)
    return counts

def get_num_haplotype(
    bam:pysam.AlignmentFile,
    positions:List[GenomePosition],
    min_count:float=10
):
    counts = get_haplotype_counts(bam,positions)
    if min_count < 1:
        temp_count = sum(counts.values()) * min_count
        print(f"Filtering haplotypes with less than {temp_count} reads")
        filtered_counts = {k:v for k,v in counts.items() if v > temp_count}
    else:
        filtered_counts = {k:v for k,v in counts.items() if v > min_count}
    return len(filtered_counts)

def get_triplets(
    vcf:pysam.VariantFile,
    maxdist:int=100
):
    snp_positions = []
    for var in vcf:
        snp_positions.append(GenomePosition(chrom=var.contig,pos=var.pos))  
    triplets = []
    for p in snp_positions:
        gr = GenomeRange(chrom=p.chrom,start=p.pos,end=p.pos+maxdist)
        positions = [p for p in snp_positions if p in gr][:3]
        if len(positions) == 3:
            triplets.append(positions) 
    return triplets

def main(
    bam:pysam.AlignmentFile,
    vcf:pysam.VariantFile,
    outfile: str,
    min_count:int=10,
    maxdist:int=500
):
    triplets = get_triplets(
        vcf,
        maxdist=maxdist
    )

    haplotype_counts = []       
    for positions in tqdm(triplets):
        num_haps = get_num_haplotype(
            bam=bam,
            positions=positions,
            min_count=min_count
        )
        if num_haps > 0:
            haplotype_counts.append(num_haps)
    
    haplotype_counts.sort()
    
    count = Counter(haplotype_counts)
    table = [[k,v] for k,v in sorted(count.items(),key=lambda x: x[0])]
    print(tabulate.tabulate(table,headers=["Num haplotypes","Sites"]))

    moi = haplotype_counts[int(len(haplotype_counts)*0.9)]
    print()
    print(f"Estimated MOI: {moi}")
    result = {
        "moi": moi,
        "haplotype_counts": count,
        "software": "pymoi",
        "version": __version__
    }
    json.dump(result,open(outfile,"w"))





def cli():
    parser = argparse.ArgumentParser(description='A simple command line tool for Moi', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam', type=pysam.AlignmentFile, help='BAM (or cram) file',required = True)
    parser.add_argument('--vcf', type=pysam.VariantFile, help='VCF file',required = True)
    parser.add_argument('--outfile', type=str, help='Name of output file',required = True)
    parser.add_argument('--maxdist', type=int, default=500, help='Maximum distance between the first and last SNP')
    parser.add_argument('--min_count', type=float, default=10, help='Minimum count of haplotype. If less than 1, it is interpreted as a fraction of the total number of reads at the site')
    args = parser.parse_args()
    
    main(
        bam=args.bam,
        vcf=args.vcf,
        outfile=args.outfile,
        min_count=args.min_count,
        maxdist=args.maxdist
    )