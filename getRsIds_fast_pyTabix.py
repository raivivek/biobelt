#! /usr/bin/env python3

import sys
import time
import argparse

import numpy
import tabix
import pandas


def getOpts():
    parser = argparse.ArgumentParser(
        description="For an input file with chrom and pos fields for SNPs, fetch rsids and output in a dataframe. "
        "Uses dbsnp 150 vcf by default "
        "!!!IMPORTANT!!!: The script will return the first rsID that overlaps the SNP position, this might not be what is needed if there are overlapping indels. "
        "If you only need SNPs, subset your vcf first using vcftools --remove-indels and then use this script."
        "This script uses chunks and is fast for, say n hundred thousand queries, but could take hours for example working with GWAS summary data with many millions of queries.",
        usage="python getRsIds_fast_pyTabix.py --inputfile <inputfile> --outputfile <outputfile> --vcfFile <filepath> --chrom chromosome --pos position",
    )
    parser.add_argument(
        "--inputfile",
        required=True,
        type=str,
        help="""Tab delimited input file. Should contain at least chrom and pos information. See other options to specify format""",
    )
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        default="/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz",
        help="""The vcf file. Default = `/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz.`""",
    )
    parser.add_argument(
        "--header-pos",
        type=int,
        help="""Row number (0 indexed) which is to be taken as header. If no header, specify the --header-names argement instead.""",
    )
    parser.add_argument(
        "--header-names",
        type=str,
        nargs="+",
        help="""If input file is headerless, provide column names. """,
    )
    parser.add_argument(
        "--subset-col",
        type=str,
        help="""If need to subset input data to snps that pass a certain p value threshold, provide the column name to subset""",
    )
    parser.add_argument(
        "--subset-val",
        type=float,
        help="""Retain snps with value of subset_col column less than or equal to this threshold""",
    )
    parser.add_argument(
        "-c",
        "--chrom",
        type=str,
        default="chrom",
        help="""Column name for chromosome. Default considered = chrom """,
    )
    parser.add_argument(
        "-p",
        "--pos",
        type=str,
        default="pos",
        help="""Column name for SNP position. Default considered = pos """,
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=300000,
        help="""Number of rows in a chunk as inputfile is read in chunks """,
    )
    parser.add_argument(
        "--keep_best_in",
        nargs="+",
        help="""If need to keep snps with min p value per gene, provide space separated gene_col and p_value_col""",
    )
    parser.add_argument(
        "--keep-cols",
        type=str,
        nargs="+",
        help="""If input file is too large, subsetting to read in only required columns can be helpful. 
                        Provide these column names.""",
    )
    parser.add_argument("--log", help="""Report stdout and stderr to this file.""")
    parser.add_argument(
        "--outputfile",
        required=True,
        help="""output file name. Will contain an extra column named 'SNP'.""",
    )
    args = parser.parse_args()
    return args


def peek(iterable):
    """Check if tabix output has some values. Return the third value that is the rsID.
    If the tabix output is empty, return NaN """
    try:
        value = next(iterable)[2]
    except StopIteration:
        return numpy.nan
    except TypeError:
        return numpy.nan
    return value


def query_function(chrom, pos, chromType, opened_tabix):
    """If some weird chromosome is not found in the vcf, tabix throws TabixError. Just
    return nan in such cases so rest of the positions can be queried"""
    function_dictionary = {"str": opened_tabix.query, "int": opened_tabix.queryi}
    try:
        value = function_dictionary[chromType](chrom, pos, pos)
    except tabix.TabixError:
        return numpy.nan
    except TypeError:
        # sys.stderr.write(f"query {chrom} {pos} {pos} with type {chromType} had a TypeError")
        return numpy.nan
    return value


def apply_on_chunk(
    chunk, opened_tabix, chrom_col, pos_col, subset_col=None, subset_val=None
):
    if subset_col is not None:
        chunk = chunk[chunk[subset_col] <= subset_val]
    if chunk.empty:
        return chunk
    chunk.dropna(inplace=True)
    chunk[chrom_col] = chunk[chrom_col].astype(str).replace("chr", "")
    chunk[chrom_col] = chunk[chrom_col].astype(str)

    # chromlist = chunk[chrom_col].drop_duplicates().tolist()

    chunk["SNP"] = chunk.apply(
        lambda x: peek(opened_tabix.query(x[chrom_col], x[pos_col], x[pos_col])), axis=1
    )

    return chunk


if __name__ == "__main__":

    args = getOpts()
    if args.log is not None:
        mylog = open(args.log, "w+")
        sys.stderr = sys.stdout = mylog

    if (args.header_pos is not None and args.header_names is not None) or (
        args.header_pos is None and args.header_names is None
    ):
        print("Provide either --header-pos or --header-names argument")
        sys.exit(1)

    d = pandas.read_csv(
        args.inputfile,
        sep="\t",
        header=args.header_pos,
        names=args.header_names,
        usecols=args.keep_cols,
        chunksize=args.chunksize,
    )

    tb = tabix.open(args.vcf)

    if args.chrom:
        chrom_col = args.chrom
    else:
        chrom_col = "chrom"

    if args.pos:
        pos_col = args.pos
    else:
        pos_col = "pos"

    l_apply = [
        apply_on_chunk(
            chunk,
            tb,
            chrom_col,
            pos_col,
            subset_col=args.subset_col,
            subset_val=args.subset_val,
        )
        for chunk in d
    ]

    out = pandas.concat(l_apply)

    if args.keep_best_in is not None:
        gene_col = args.keep_best_in[0]
        pval_col = args.keep_best_in[1]

        idx = out.groupby([gene_col])[pval_col].transform(min) == out[pval_col]
        out = out[idx]

        snpidx = out.groupby(["SNP"])[pval_col].transform(min) == out[pval_col]
        out = out[snpidx]
        out.rename(columns={pval_col: "P"}, inplace=True)

    out.to_csv(args.outputfile, index=False, sep="\t", na_rep="NA")

    if args.log is not None:
        mylog.close()
