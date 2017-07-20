#!/usr/bin/env python
#Last-modified: 12 Jun 2017

#         Module/Scripts Description
# 
# Copyright (c) 2017 Traver Hart <traver@hart-lab.org>, Gang Wang <iskysinger@gmail.com>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  Bioinformatics
# @version: 0.1.1
# @authors:  Traver Hart, Gang Wang
# @contact: traver@hart-lab.org, iskysinger@gmail.com

# ------------------------------------
# python modules
# ------------------------------------
import sys

from matplotlib.mlab import find
import six

import numpy as np
import pandas as pd
import scipy.stats as stats


# ------------------------------------
# constants
norm_value  = 10000000.
qmin = 0.025
qmax = 0.975
min_reads_thresh = 0

# ------------------------------------

def drugz(readfile, nonessfile, drugz_outfile, control_samples, drug_samples, 
          remove_genes=None, pseudocount=5, minObs=6, index_column=0, verbose=False):
    num_replicates = len(control_samples)
    
    verbose and six.print_('Control samples:  ' + str(control_samples), file=sys.stderr)
    verbose and six.print_('Treated samples:  ' + str(drug_samples), file=sys.stderr)
    
    #
    #read non essential genes
    #
    non=pd.read_table(nonessfile,index_col=0);
    #
    #read sgRNA reads counts file
    reads = pd.read_table(readfile, index_col=index_column)
    numGuides, numSamples = reads.shape
    #
        
    #sample_sum = reads.ix[:,range(1,numSamples)].sum(0)
    #
    #normalize to norm_value reads
    #
    verbose and six.print_('Normalizing read counts', file=sys.stderr)
    normed = norm_value * reads.ix[:,control_samples+drug_samples] / reads.ix[:,control_samples+drug_samples].sum().as_matrix()
    
    
    ##
    #Caculate fold change with normalized reads + pseudocount
    # maintain raw read counts for future filtering
    ##
    
    verbose and six.print_('Caculating fold change', file=sys.stderr)
    fc = pd.DataFrame(index=reads.index.values)
    fc['GENE'] = reads.ix[:,'GENE']      # first column of input file MUST be gene name!
    for k in range(len(control_samples)):    
        fc[control_samples[k]] = reads[ control_samples[k] ]
        fc[drug_samples[k]] = reads[ drug_samples[k] ]
        fc['fc_{0}'.format(k)] = np.log2(( normed[ drug_samples[k] ] + pseudocount ) / ( normed[ control_samples[k] ]+ pseudocount))   
    
    
    # remove control genes
    # e.g. TKOv1 genes ['chr10Promiscuous','chr10Rand','chr10','EGFP','LacZ','luciferase']
    
    if ( remove_genes ):
        fc = fc.ix[~fc.GENE.isin(remove_genes),:]
    
    ##
    #create drugz foldchange dataframe
    ##
    
    dz_fc = pd.DataFrame(index=fc.index.values)
    dz_fc['GENE']=fc.GENE
    
    #
    # find nonessential/control reference 
    #
    nonidx = find( np.in1d(dz_fc.GENE, non.index.values))
    
    #
    # get fold changes from specficied samples
    #
    for i in range(len(control_samples)):
        f=find(fc.ix[:,control_samples[i]] > min_reads_thresh)
        dz_fc['dz_fc_{0}'.format(i)]=fc.ix[f,'fc_{0}'.format(i)]
    numGuides, numSamples = dz_fc.shape
    
    #
    # calculate moderated zscores for each gRNA
    #
    
    verbose and six.print_('Caculating Zscores', file=sys.stderr)
    for i in range(1, numSamples):
        sample = dz_fc.columns.values[i]
        fin    = find( np.isfinite(dz_fc.ix[nonidx,sample]))
        quants = dz_fc.ix[nonidx[fin],sample].quantile([qmin,qmax])
        g = find( (dz_fc.ix[nonidx,sample] > quants[qmin]) & (dz_fc.ix[nonidx,sample] < quants[qmax]) & ( np.isfinite(dz_fc.ix[nonidx,sample]) ) )
        sigmag = dz_fc.ix[nonidx[g],sample].std()
        mug = dz_fc.ix[nonidx[g],sample].mean()
        zsample = 'Z_' + sample
        dz_fc[zsample] = (dz_fc.ix[:,sample] - mug) / sigmag
    
    #
    # combine to gene-level drugz scores
    #
    
    verbose and six.print_('Combining drugZ scores', file=sys.stderr)
    
    # get unique list of genes in the data set
    usedColumns = ['Z_dz_fc_{0}'.format(i) for i in range(num_replicates)]
    drugz = dz_fc.groupby('GENE')[usedColumns].apply(lambda x: pd.Series([np.nansum(x.values), np.isfinite(x.values).sum()]))
    drugz.columns = ['sumZ', 'numObs']
    drugz.loc[:,'normZ'] = drugz.sumZ / np.sqrt(drugz.numObs)
    
    #
    #
    verbose and six.print_('Writing output file', file=sys.stderr)
    #
    # calculate numObs, pvals (from normal dist), and fdrs (by benjamini & hochberg).
    #
    drugz_minobs = drugz.ix[drugz.numObs>=minObs,:]
    numGenes, numCols = drugz_minobs.shape
    drugz_minobs=drugz_minobs.sort_values('normZ', ascending=True)
    drugz_minobs.loc[:,'pval_synth'] = stats.norm.sf( drugz_minobs.loc[:,'normZ'] * -1)
    drugz_minobs.loc[:,'rank_synth'] = np.arange(1,numGenes +1)
    drugz_minobs.loc[:,'fdr_synth'] = drugz_minobs['pval_synth']*numGenes/drugz_minobs.loc[:,'rank_synth']
    drugz_minobs = drugz_minobs.sort_values('normZ', ascending=False)
    drugz_minobs.loc[:,'pval_supp'] = stats.norm.sf( drugz_minobs.loc[:,'normZ'])
    drugz_minobs.loc[:,'rank_supp'] = np.arange(1,numGenes +1)
    drugz_minobs.loc[:,'fdr_supp']  = drugz_minobs.loc[:,'pval_supp'] * numGenes / drugz_minobs.loc[:,'rank_supp']
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
    #    fout.write(i + '\t')
        fout.write( '{0:s}\t{1:3.2f}\t{2:d}\t{3:4.2f}\t{4:.3g}\t{5:d}\t{6:.3g}\t{7:.3g}\t{8:d}\t{9:.3g}\n'.format( \
            i, \
            drugz_minobs.ix[i,'sumZ'], \
            int(drugz_minobs.ix[i,'numObs']), \
            drugz_minobs.ix[i,'normZ'], \
            drugz_minobs.ix[i,'pval_synth'], \
            int(drugz_minobs.ix[i,'rank_synth']), \
            drugz_minobs.ix[i,'fdr_synth'], \
            drugz_minobs.ix[i,'pval_supp'], \
            int(drugz_minobs.ix[i,'rank_supp']), \
            drugz_minobs.ix[i,'fdr_supp'] ) )
    
    fout.close()


def main():
    import argparse
    
    ''' Parse arguments. '''
    p = argparse.ArgumentParser(description='DrugZ for chemogenetic interaction screens',epilog='dependencies: pylab, pandas')
    p._optionals.title = "Options"
    p.add_argument("-i", dest="infile", type=argparse.FileType('r'), metavar="sgRNA_count.txt", help="sgRNA readcount file", default=sys.stdin)
    p.add_argument("-o", dest="drugz", type=argparse.FileType('w'), metavar="drugz-output.txt", help="drugz output file", default=sys.stdout) 
    p.add_argument("-n", dest="ness", type=argparse.FileType('r'), metavar="NEG.txt", required=True, help="Non-essential gene list")
    p.add_argument("-c", dest="control_samples", metavar="control samples", required=True, help="control samples, comma delimited")
    p.add_argument("-x", dest="drug_samples", metavar="drug samples", required=True, help="treatment samples, comma delimited")
    p.add_argument("-r", dest="remove_genes", metavar="remove genes", help="genes to remove, comma delimited", default='')
    p.add_argument("-p", dest="pseudocount", type=int, metavar="pseudocount", help="pseudocount (default=5)", default=5)
    p.add_argument("-I", dest="index_column", type=int, help="Index column in the input file (default=0; GENE_CLONE column)", default=0)
    p.add_argument("--minobs", dest="minObs", type=int,metavar="minObs", help="min number of obs (default=6)", default=6)
    p.add_argument("-q", dest="quiet", action='store_true', default=False, help='Be quiet, do not print log messages')

    args = p.parse_args()
    
    control_samples = args.control_samples.split(',')
    drug_samples = args.drug_samples.split(',')
    remove_genes = args.remove_genes.split(',')
    
    if len(control_samples) != len(drug_samples):
        p.error("Must have the same number of control and drug samples")
    
    drugz(args.infile, args.ness, args.drugz, control_samples, drug_samples, 
          remove_genes, args.pseudocount, args.minObs, args.index_column, not args.quiet)

if __name__=="__main__":
    main()
