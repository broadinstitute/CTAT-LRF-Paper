#!/usr/bin/env python

import sys, os, re
import gzip
import logging
import argparse
import pandas as pd

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():
    
    parser = argparse.ArgumentParser(description="__add_descr__", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--gff3_read_alignments", type=str, required=True, help="gff3 read alignments")
    parser.add_argument("--fusion_predictions", type=str, required=True, help="ctat-LRF-FI fusion predictions tsv")
    parser.add_argument("--output_prefix", required=True, help="output prefix")
    
    args = parser.parse_args()


    gff3_read_alignments_file = args.gff3_read_alignments
    fusion_predictions_tsv = args.fusion_predictions
    output_prefix = args.output_prefix
    

    # require SR and LR support:
    df = pd.read_csv(fusion_predictions_tsv, sep="\t")

    df = df.rename(columns={"#FusionName" : "FusionName",
                            "FFPM" : "SR_FFPM" } )

    df['num_SR'] = df['est_J'] + df['est_S']

    # require SR and LR support:
    df = df[ (df['LR_FFPM'] > 0) & (df['SR_FFPM'] > 0)]
    df = df.reset_index(drop=True)
    
    df = df[ ['FusionName', 'num_LR', 'num_SR', 'LeftLocalBreakpoint', 'LeftBreakpoint', 
          'RightLocalBreakpoint', 'RightBreakpoint', 'SpliceType', 'LR_FFPM', 'SR_FFPM', 'LR_accessions' ] ].copy()

    df['LR_accessions'] = df['LR_accessions'].apply(lambda x: x.split(","))
    df = df.explode('LR_accessions')

    
    # get list of the reads and the breakpoints of interest.
    long_read_to_local_brkpt_right = dict()
    long_read_to_local_brkpt_left = dict()
    for index,row in df.iterrows():
        fusion_name = row['FusionName']
        long_read = row['LR_accessions']
        long_read_fusion_token = "::".join([fusion_name, long_read])
        long_read_to_local_brkpt_right[long_read_fusion_token] = row['RightLocalBreakpoint']
        long_read_to_local_brkpt_left[long_read_fusion_token] = row['LeftLocalBreakpoint']

    # parse alignments, get coords for downstream alignment segments.
    alignments_df = pd.read_csv(gff3_read_alignments_file, sep="\t", header=None,
                            names=['FusionName', 'aligner', 'feat_type', 'lend', 'rend', 'per_id', 'strand', 'nothing', 'info'])

    alignments_df['read_name'] = alignments_df.apply(get_read_name, axis=1)
            
    alignments_df['long_read_fusion_token'] = alignments_df.apply(lambda x: x["FusionName"] + "::" + x['read_name'], axis=1)
    
    # restrict to alignments of interest, having long and short read support at breakpoints:
    alignments_df = alignments_df[ alignments_df['long_read_fusion_token'].isin(long_read_to_local_brkpt_right.keys()) ]

    fusion_3prime_lengths =  alignments_df.groupby('long_read_fusion_token').apply(sum_3prime_lens, long_read_to_local_brkpt_right).reset_index(drop=True)

    # merge fusion info with breakpoint 3' read length info
    df['long_read_fusion_token'] = df.apply(lambda x: x["FusionName"] + "::" + x['LR_accessions'], axis=1)

    fusion_3prime_lengths = df.merge(fusion_3prime_lengths, how='left', on='long_read_fusion_token')
    fusion_3prime_lengths.to_csv(output_prefix + ".per_read.tsv", sep="\t", index=False)
    
    # get median length for read alignment beyond breakpoint
    fusion_3prime_median_lengths = fusion_3prime_lengths \
      .groupby(['FusionName', 'LeftLocalBreakpoint', 'RightLocalBreakpoint']) \
      .agg({'threePrimeBrkLen': 'median',
            'align_len' : 'median'}) \
      .reset_index(drop=False)
    
    # incorporate the LR and SR support info:
    fusion_3prime_median_lengths = fusion_3prime_median_lengths.merge(df[['FusionName','LeftLocalBreakpoint', 'RightLocalBreakpoint','num_LR', 'num_SR', 'LR_FFPM', 'SR_FFPM']],
                                   on=['FusionName', 'LeftLocalBreakpoint', 'RightLocalBreakpoint']).drop_duplicates()

    fusion_3prime_median_lengths['SR/LR'] = fusion_3prime_median_lengths.apply(lambda x: (x['SR_FFPM']/x['LR_FFPM']), axis=1)


    fusion_3prime_median_lengths.to_csv(output_prefix + ".tsv", sep="\t", index=False)
    
    
    sys.exit(0)


def get_read_name(row):
    m  = re.search("Target=(\\S+)", row['info'])
    read_name = m.group(1)
    return read_name


            

# get length of 3' read alignment passed breakpoint.
def sum_3prime_lens(grp, long_read_to_local_brkpt_right):
    first_row = grp.head(1)
    #print(first_row)
    long_read_fusion_token = first_row['long_read_fusion_token'].values[0]
    
    local_breakpoint = long_read_to_local_brkpt_right[long_read_fusion_token]
    #print(f"local_brkpt: {local_breakpoint}")
    
    # might not be within snap distance...
    if not local_breakpoint in grp['lend'].values:
        print("Warning, missing brkpt {} for {}\n{}".format(local_breakpoint, long_read_fusion_token, grp), sys.stderr)
        return None
    
    
    threePrimeBrkLen = 0
    align_len = 0
    for row in grp.itertuples():
        #print(row)
        seglen = row.rend - row.lend + 1

        align_len += seglen
        
        if row.lend >= local_breakpoint:
            threePrimeBrkLen += seglen

    #print(str(threePrimeBrkLen))

    d = { 'long_read_fusion_token' : [long_read_fusion_token],
          'align_len' : [align_len],
          'threePrimeBrkLen' : [threePrimeBrkLen]
        }
    
    ret_df = pd.DataFrame.from_dict(d)

    return ret_df



if  __name__=='__main__':
    main()
