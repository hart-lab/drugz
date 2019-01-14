#!/usr/bin/env python

VERSION = "1.1.0.2"
BUILD   = 111

#---------------------------------
# DRUGZ:  Identify drug-gene interactions in paired sample genomic perturbation screens
# Special thanks to Matej Usaj
# Last modified 13 Jun 2018
# Free to modify and redistribute with attribtuion
#---------------------------------


# ------------------------------------
# python modules
# ------------------------------------
import sys

import six

import numpy as np
import pandas as pd
import scipy.stats as stats

pd.options.mode.chained_assignment = None  # default='warn'

# ------------------------------------
# constants
norm_value  = 1e7
min_reads_thresh = 1
#half_window_size = 500
# ------------------------------------

def drugz(readfile, drugz_outfile, control_samples, drug_samples, fc_outfile=None,
          remove_genes=None, pseudocount=5, minObs=1, half_window_size = 500, index_column=0, verbose=False):
    # parameter "nonessfile" removed
    def log_(msg):
        if verbose:
            six.print_(msg, file=sys.stderr)

    num_replicates = len(control_samples)

    log_('Control samples:  ' + str(control_samples))
    log_('Treated samples:  ' + str(drug_samples))

    #
    #read sgRNA reads counts file
    #
    reads = pd.read_table(readfile, index_col=index_column)

    # remove control genes
    # e.g. TKOv1 genes ['chr10Promiscuous','chr10Rand','chr10','EGFP','LacZ','luciferase']
    # TKOv3: 'EGFP','LacZ','luciferase'
    if ( remove_genes ):
        reads = reads.loc[~reads[ reads.columns.values[0] ].isin(remove_genes),:]

    numGuides, numSamples = reads.shape
    #

    #sample_sum = reads.ix[:,range(1,numSamples)].sum(0)
    #
    #normalize to norm_value reads
    #
    log_('Normalizing read counts')
    normed = norm_value * reads[control_samples+drug_samples] / reads[control_samples+drug_samples].sum().as_matrix()


    ##
    #Caculate fold change with normalized reads + pseudocount
    # maintain raw read counts for future filtering
    ##
    log_('Processing data')
    fc = pd.DataFrame(index=reads.index.values)
    fc['GENE'] = reads[ reads.columns.values[0] ]      # first column of input file MUST be gene name!

    for k in range(len(control_samples)):
        log_('Calculating raw fold change for replicate {0}'.format(k+1))
        fc[control_samples[k]] = reads[ control_samples[k] ]
        fc[drug_samples[k]] = reads[ drug_samples[k] ]
        fc['fc_{0}'.format(k)] = np.log2(( normed[ drug_samples[k] ] + pseudocount ) / ( normed[ control_samples[k] ]+ pseudocount))
        #
        # sort guides by readcount, descending:
        #
        fc.sort_values(control_samples[k], ascending=False, inplace=True)
        #
        # add columns for eb_mean, eb_std
        #
        #eb_mean_samplid = 'eb_mean_{0}'.format(k)
        eb_std_samplid  = 'eb_std_{0}'.format(k)
        #fc[eb_mean_samplid] = np.zeros(numGuides)
        fc[eb_std_samplid] = np.zeros(numGuides)
        #
        # get mean, std of fold changes based on 800 nearest fc
        #
        log_('Caculating smoothed Epirical Bayes estimates of stdev for replicate {0}'.format(k+1))
        #
        # initialize element at index 250
        #
        # do not mean-center. fc of 0 should be z=score of 0.
        #ebmean = fc.iloc[0:500]['fc_{0}'.format(k)].mean()
        #fc[eb_mean_samplid][0:250] = ebmean

        # parameter "half_window_size" added
        half_window_size = half_window_size
        ebstd  = fc.iloc[0:half_window_size*2]['fc_{0}'.format(k)].std()
        fc[eb_std_samplid][0:half_window_size]  = ebstd
        #
        # from 250..(end-250), calculate mean/std, update if >= previous (monotone smoothing)
        #
        for i in range(half_window_size, numGuides-half_window_size+25, 25):
            #every 25th guide, calculate stdev. binning/smoothing approach.
            ebstd  = fc.iloc[i-half_window_size:i+half_window_size]['fc_{0}'.format(k)].std()
            if (ebstd >= fc[eb_std_samplid][i-1]):
                fc[eb_std_samplid][i:i+25] = ebstd              #set new std in whole step size (25)
            else:
                fc[eb_std_samplid][i:i+25] = fc.iloc[i-1][eb_std_samplid]
        #
        # set ebmean, ebstd for bottom half-window set of guides
        #
        #log_('Smoothing estimated std for replicate {0}'.format(k+1))
        #fc[eb_mean_samplid][numGuides-250:] = fc.iloc[numGuides-251][eb_mean_samplid]
        fc[eb_std_samplid][numGuides-half_window_size:] = fc.iloc[numGuides-(half_window_size+1)][eb_std_samplid]
        #
        # calc z score of guide
        #
        log_('Caculating Zscores for replicate {0}'.format(k+1))
        #fc['Zlog_fc_{0}'.format(k)] = (fc['fc_{0}'.format(k)] - fc[eb_mean_samplid]) / fc[eb_std_samplid]
        fc['Zlog_fc_{0}'.format(k)] = fc['fc_{0}'.format(k)] / fc[eb_std_samplid]

    ##
    # write fc file as intermediate output
    ##
    if ( fc_outfile ):
    	fc.to_csv( fc_outfile, sep='\t', float_format='%4.3f')

    ##
    # sum guide-level zscores to gene-level drugz scores. keep track of how many elements (fold change observations) were summed.
    ##

    log_('Combining drugZ scores')

    # get unique list of genes in the data set
    usedColumns = ['Zlog_fc_{0}'.format(i) for i in range(num_replicates)]
    drugz = fc.groupby('GENE')[usedColumns].apply(lambda x: pd.Series([np.nansum(x.values), np.count_nonzero(x.values)]))
    drugz.columns = ['sumZ', 'numObs']
    #
    #
    log_('Writing output file')
    #
    # calculate normZ, pvals (from normal dist), and fdrs (by benjamini & hochberg).
    #
    drugz_minobs = drugz.ix[drugz.numObs>=minObs,:]
    numGenes, numCols = drugz_minobs.shape
    ##
    # normalize sumZ by number of observations
    ##
    drugz_minobs['normZ'] = drugz_minobs['sumZ'] / np.sqrt( drugz_minobs['numObs'])
    #
    # renormalize to fit uniform distribution of null p-vals
    #
    drugz_minobs['normZ'] = stats.zscore(drugz_minobs['normZ'])

    #
    # rank by normZ (ascending) to identify synthetic interactions
    #
    drugz_minobs.sort_values('normZ', ascending=True, inplace=True)
    drugz_minobs['pval_synth'] = stats.norm.sf( drugz_minobs['normZ'] * -1)
    drugz_minobs['rank_synth'] = np.arange(1,numGenes +1)
    drugz_minobs['fdr_synth'] = drugz_minobs['pval_synth']*numGenes/drugz_minobs['rank_synth']
    #
    # rerank by normZ (descending) to identify suppressor interactions
    #
    drugz_minobs = drugz_minobs.sort_values('normZ', ascending=False)
    drugz_minobs['pval_supp'] = stats.norm.sf( drugz_minobs['normZ'])
    drugz_minobs['rank_supp'] = np.arange(1,numGenes +1)
    drugz_minobs['fdr_supp']  = drugz_minobs['pval_supp'] * numGenes / drugz_minobs['rank_supp']
    drugz_minobs = drugz_minobs.sort_values('normZ', ascending=True)
    #
    # write output file
    #
    fout = drugz_outfile
    if not hasattr(fout, 'write'):
        fout = open(fout, 'w')
    fout.write('GENE')
    cols = drugz_minobs.columns.values
    for c in cols:
        fout.write('\t' + c)
    fout.write('\n')

    for i in drugz_minobs.index.values:
        #fout.write(i + '\t')
        fout.write( '{0:s}\t{1:3.2f}\t{2:d}\t{3:4.2f}\t{4:.3g}\t{5:d}\t{6:.3g}\t{7:.3g}\t{8:d}\t{9:.3g}\n'.format( \
            i, \
            drugz_minobs.loc[i,'sumZ'], \
            int(drugz_minobs.loc[i,'numObs']), \
            drugz_minobs.loc[i,'normZ'], \
            drugz_minobs.loc[i,'pval_synth'], \
            int(drugz_minobs.loc[i,'rank_synth']), \
            drugz_minobs.loc[i,'fdr_synth'], \
            drugz_minobs.loc[i,'pval_supp'], \
            int(drugz_minobs.loc[i,'rank_supp']), \
            drugz_minobs.loc[i,'fdr_supp'] ) )

    fout.close()


def main():
    import argparse

    ''' Parse arguments. '''
    p = argparse.ArgumentParser(description='DrugZ for chemogenetic interaction screens',epilog='dependencies: pylab, pandas')
    p._optionals.title = "Options"
    p.add_argument("-i", dest="infile", type=argparse.FileType('r'), metavar="sgRNA_count.txt", help="sgRNA readcount file", default=sys.stdin)
    p.add_argument("-o", dest="drugz", type=argparse.FileType('w'), metavar="drugz-output.txt", help="drugz output file", default=sys.stdout)
    p.add_argument("-f", dest="fc_outfile", type=argparse.FileType('w'), metavar="drugz-foldchange.txt", help="drugz normalized foldchange file (optional")
    #p.add_argument("-n", dest="ness", type=argparse.FileType('r'), metavar="NEG.txt", required=True, help="Non-essential gene list")
    p.add_argument("-c", dest="control_samples", metavar="control samples", required=True, help="control samples, comma delimited")
    p.add_argument("-x", dest="drug_samples", metavar="drug samples", required=True, help="treatment samples, comma delimited")
    p.add_argument("-r", dest="remove_genes", metavar="remove genes", help="genes to remove, comma delimited", default='')
    p.add_argument("-p", dest="pseudocount", type=int, metavar="pseudocount", help="pseudocount (default=5)", default=5)
    p.add_argument("-I", dest="index_column", type=int, help="Index column in the input file (default=0; GENE_CLONE column)", default=0)
    p.add_argument("--minobs", dest="minObs", type=int,metavar="minObs", help="min number of obs (default=1)", default=1)
    p.add_argument("--half_window_size", dest="half_window_size", type=int,metavar="half_window_size", help="width of variance-estimation window", default=500)
    p.add_argument("-q", dest="quiet", action='store_true', default=False, help='Be quiet, do not print log messages')

    args = p.parse_args()

    control_samples = args.control_samples.split(',')
    drug_samples = args.drug_samples.split(',')
    remove_genes = args.remove_genes.split(',')

    if len(control_samples) != len(drug_samples):
        p.error("Must have the same number of control and drug samples")

    #drugz(args.infile, args.ness, args.drugz, control_samples, drug_samples,
    #      remove_genes, args.pseudocount, args.minObs, args.index_column, not args.quiet)
    drugz(args.infile, args.drugz, control_samples, drug_samples,
    	args.fc_outfile, remove_genes, args.pseudocount, args.minObs, args.half_window_size, args.index_column, not args.quiet)


if __name__=="__main__":
    main()
