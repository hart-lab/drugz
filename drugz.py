# !/usr/bin/env python

import argparse
import csv
import logging as log
import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import sys

VERSION = "1.1.0.2"
BUILD = 111

# ---------------------------------
# DRUGZ:  Identify drug-gene interactions in paired sample genomic perturbation screens
# Special thanks to Matej Usaj
# Last modified 13 Jun 2018
# Free to modify and redistribute with attribtuion
# ---------------------------------

##########################
#   GLOBAL VARIABLES     #
##########################

pd.options.mode.chained_assignment = None  # default='warn'

# ------------------------------------
# constants
NORM_VALUE = 1e7
min_reads_thresh = 1


# half_window_size = 500
# ------------------------------------
def get_args():
    parser = argparse.ArgumentParser(description='DrugZ for chemogenetic interaction screens',
                                     epilog='dependencies: pylab, pandas')

    parser._optionals.title = "Options"

    parser.add_argument("-i", dest="infile", type=argparse.FileType('r'), metavar="sgRNA_count.txt",
                        help="sgRNA readcount file", default=sys.stdin)
    parser.add_argument("-o", dest="drugz_output_file", type=argparse.FileType('w'), metavar="drugz-output.txt",
                        help="drugz output file", default=sys.stdout)
    parser.add_argument("-f", dest="fc_outfile", type=argparse.FileType('w'), metavar="drugz-foldchange.txt",
                        help="drugz normalized foldchange file (optional")
    # parser.add_argument("-n", dest="ness", type=argparse.FileType('r'),
    #                     metavar="NEG.txt", required=True, help="Non-essential gene list")
    parser.add_argument("-c", dest="control_samples", metavar="control samples", required=True,
                        help="control samples, comma delimited")
    parser.add_argument("-x", dest="drug_samples", metavar="drug samples", required=True,
                        help="treatment samples, comma delimited")
    parser.add_argument("-r", dest="remove_genes", metavar="remove genes", help="genes to remove, comma delimited")
    parser.add_argument("-p", dest="pseudocount", type=int, metavar="pseudocount", help="pseudocount (default=5)",
                        default=5)
    parser.add_argument("-I", dest="index_column", type=int,
                        help="Index column in the input file (default=0; GENE_CLONE column)", default=0)
    parser.add_argument("--minobs", dest="minObs", type=int, metavar="minObs", help="min number of obs (default=1)",
                        default=1)
    parser.add_argument("--half_window_size", dest="half_window_size", type=int, metavar="half_window_size",
                        help="width of variance-estimation window", default=500)
    parser.add_argument("-q", dest="quiet", action='store_true', default=False,
                        help='Be quiet, do not print log messages')
    return parser.parse_args()


def sort_by_sample_column(reads, fold_change, control_sample):
    """
    Add the current control sample column to our fold change dataframe (from the reads dataframe) and then sort the fold
    change by the current control sample in descending order.

    :param reads: read count dataframe object
    :param fold_change: fold change dataframe object
    :param control_sample: the column name (str) of the control sample for the current comparison.
    :return: fold change dataframe updated with the control sample column and ordered by this sample
    (largest to smallest)
    """
    fold_change[control_sample] = reads[control_sample]
    fold_change.sort_values(control_sample, ascending=False, inplace=True)
    return fold_change


def load_reads(filepath, index_column, genes_to_remove):
    """
    Load the input file and remove any genes provided through the remove_genes input.

    :param filepath: The path to the file to be loaded.
    :param index_column: The column to use as an index (not this will not be included in the resulting dataframe,
    default value is 0 e.g. the 1st column).
    :param genes_to_remove: A string of comma separated gene names.
    :return: reads: a dataframe containing the read counts to be used in the analysis.
    """

    reads = pd.read_csv(filepath, index_col=index_column, delimiter='\t')

    if genes_to_remove is not None:

        genes_to_remove = genes_to_remove.split(',')
        gene_column = reads.columns.values[0]

        cut_down_reads = reads.loc[~reads[gene_column].isin(genes_to_remove), :]

        return cut_down_reads

    else:
        return reads


def init_foldchange(reads):
    """
    Return a dataframe object containing the gene names for the analysis and with the correct index values.
    :param reads: A dataframe of the read count input file
    :return: A One column dataframe containing the gene names to be used in the analysis.
    """
    fold_change = pd.DataFrame(index=reads.index.values)
    fold_change['GENE'] = reads[reads.columns.values[0]]  # first column of input file MUST be gene name!
    # TODO: Should probably actually verify the above assumption prior to running program
    return fold_change


def normalise_readcounts(reads, treatment, control):
    """
    Normalise input read counts using the global variable NORM VALUE

    :param reads: Dataframe containing reads counts for each guide.
    :param treatment: List of columns names for the samples in the treatment group
    :param control: List of column names for the samples in the control group
    :return: normalised_counts: A dataframe of normalised read counts.
    """
    reads_to_normalise = reads[control + treatment]

    normalised_counts = (NORM_VALUE * reads_to_normalise) / reads_to_normalise.sum().values
    return normalised_counts


def empirical_bayes(fold_change, half_window_size, empiric_bayes_id, no_of_guides, replicate_id):
    """
    Calculate the variation present in log ratio of treatment to control read counts for bins of the data, smoothing
    the variation during this process to ensure that variation only ever increases or remains the same as the estimate
    for the previous bin.

    The estimate of variation is first initialise for 2x the range of the half_window_size. The default for this value
    if none provided via the argument is 500). Thus, for the first 1000 values we calculate the estimate of variation
    (std dev) and the set the 1st bin (0 to half_the_window_size - and so 0 to 500 by default) equal to this.

    For each subsequent bin between  j (the half-window-size) and n -  (where n is the number of rows) we then calculate
    the estimate of std dev for this bin and if the estimate is greater than the previous bin store the variance as
    this, otherwise we store it as the estimate of the previous bin. This acts to smooth the variance, so that varience
    estimates only ever tick upwards or remain flat between each bin.

    The final bin (n-j : n) is then set to the variance of the previous bin (e.g. the last estimate to be calculated).

    :param fold_change: A dataframe containing log ratio of treatment read counts to control read counts
    :param half_window_size: An integer value equal to the size of the first bin and half the size of the inital sample
    (window) to estimate StDev. By default this is 500.
    :param empiric_bayes_id: A string value that is the id for the current comparison. This acts as a column name under
    which to store the variance estimates.
    :param no_of_guides: An integer value euqal to the number of rows of the data frame.
    :param replicate_id: The column name under which the log ratio values are stored for the current comparison
    :return: fold_change: updated instance of the input dataframe with the zeros in the column specified by the
    empiric_bayes_id replace with smooothed (read doc string above) per bin estimate of StDev
    """

    # Calculate the standard deviation of the log ratio of the sample based on the a range 2x the defined window size
    std_dev = fold_change.iloc[0: half_window_size * 2][replicate_id].std()
    fold_change[empiric_bayes_id][0: half_window_size] = std_dev

    # do not mean-center. fc of 0 should be z=score of 0.
    # initialize element at index equal to half windows size (j), default = 500
    # so from j..(end-j), calculate mean/std, update if >= previous (this performs monotone smoothing)
    # Thus we generate an iterable of range(j, n-j, 25) where j = half_window_size and n is the number of rows
    for i in range(half_window_size, no_of_guides - half_window_size + 25, 25):
        # every 25th guide, calculate stdev. binning/smoothing approach.
        std_dev = fold_change.iloc[i - half_window_size:i + half_window_size][replicate_id].std()

        # If the current variation is greater than the previous bin then set variation equal to this
        if std_dev >= fold_change[empiric_bayes_id][i - 1]:
            fold_change[empiric_bayes_id][i:i + 25] = std_dev  # set new std in whole step size (25)

        # Otherwise, set it equal to the variation of the previous bin
        # In this way the variation estimate for each bin can only ever increase or stay the same as the previous.
        else:
            fold_change[empiric_bayes_id][i:i + 25] = fold_change.iloc[i - 1][empiric_bayes_id]

    # Get the variation estimate for the final bin and set the remaining values in the empirical bayes column equal to
    # this estimate.
    results = fold_change.iloc[no_of_guides - (half_window_size + 1)][empiric_bayes_id]
    fold_change[empiric_bayes_id][no_of_guides - half_window_size:] = results

    return fold_change


def calculate_fold_change(reads, normalised_counts, fold_change, control_samples, treatment_samples, pseudocount,
                          half_window_size, logger, replicate):
    # Generate a unique pair of ids for each replicate
    replicate_id = 'comparison_{replicate}'.format(replicate=replicate)
    zscore_id = 'zscore_id_' + replicate_id
    empirical_bayes_comparison_id = 'eb_std_{replicate}'.format(replicate=replicate)
    one_based_idx = replicate + 1

    # Get the control and treatment sample ids for this comparison:
    control_sample = control_samples[replicate]
    treatment_sample = treatment_samples[replicate]

    # Add the control sample column to the fold change dataframe and sort by this column
    fold_change = sort_by_sample_column(reads=reads, fold_change=fold_change, control_sample=control_sample)

    # Extract the row and columns numbers of the reads dataframe
    no_of_guides, no_of_samples = reads.shape

    logger.debug('Calculating raw fold change for replicate: {replicate}'.format(replicate=one_based_idx))

    # Calculate the log ratio of treatment normalised read counts to control
    # Psuedocount adds a constant value to all read counts for both treatment and control, default = 5
    fold_change[replicate_id] = np.log2(
        (normalised_counts[treatment_sample] + pseudocount) / (normalised_counts[control_sample] + pseudocount))

    # Make a column of zeros equal to the row no. of the dataframe to store the EB estimates in
    fold_change[empirical_bayes_comparison_id] = np.zeros(no_of_guides)

    logger.debug('Calculating smoothed Empirical Bayes estimates of stdev for replicate:'
                 '{replicate}'.format(replicate=replicate + 1))

    # Calculate a smoothed per bin estimate of the variance of the log ratio between the treatment and control sample
    # for their the normalised read counts for the current comparison
    fold_change = empirical_bayes(fold_change=fold_change, half_window_size=half_window_size, replicate_id=replicate_id,
                                  empiric_bayes_id=empirical_bayes_comparison_id, no_of_guides=no_of_guides)

    logger.debug('Calculating Zscores for replicate {replicate}'.format(replicate=one_based_idx))

    # calculate the z score of guides by dividing the log ratios values by the variance estimates
    # Remember that the values have been ordered by the control so the lower estimates of variance (if varience
    # isn't invariant) will correlate to the largest values of the control.
    fold_change[zscore_id] = fold_change[replicate_id] / fold_change[empirical_bayes_comparison_id]

    return fold_change, zscore_id


def calculate_drugz_score(fold_change, min_observations, columns):
    """
    Calculate per gene statistics for the zscores aggregated across all comparisons.

    The summed zscores and the number of observations for each gene are first aggregated. These zscores are then
    normalised and pvalue estimates (assuming guassian distribution), rank position, and FDR are calculated.

    The statistics are calculated twice first (with normalised zscores ranked smallest to largest) to identify synthetic
    interactions and then (with normalised zscores now ranked largest to smallest) to identify suppressor interactions.

    :param fold_change: Data frame containing calculated zscores per comparison
    :param min_observations: An integer value to act as a threshold for the minimum number observations to be included
    in the analysis. The default value is 1.
    :param columns: The column names (list) of the columns containing the per comparison zscores.
    :return: per_gene_results: a dataframe of summary statistics for each gene (see above doc string for more info).
    """
    ##
    # sum guide-level zscores to gene-level drugz scores. keep track of how many elements (fold change observations)
    # were summed.

    # This is going to produce a per gene summation of the zscores for each comparison. Missing values are converted
    # to zeros. The output is stored in a column (to be named sumZ).
    # The second column will be a per gene count of all non_zero observations for that gene
    per_gene_scores = fold_change.groupby('GENE')[columns].apply(lambda x: pd.Series([np.nansum(x.values),
                                                                                      np.count_nonzero(x.values)]))
    # Title these two columns.
    per_gene_scores.columns = ['sumZ', 'numObs']

    # Get a dataframe of values for genes where the number of observations is greater than the minimum threshold
    # By default the threshold is 1.
    per_gene_results = per_gene_scores.loc[per_gene_scores.numObs >= min_observations, :]

    # Update the column number and row number estimates for this new dataframe.
    no_of_genes, no_of_columns = per_gene_results.shape

    # normalize the sumZ values by number of observations, then renormalize these values to fit to fit uniform
    # distribution of null p-vals
    normalised_z_scores = stats.zscore(per_gene_results['sumZ'] / np.sqrt(per_gene_results['numObs']))
    per_gene_results['normZ'] = normalised_z_scores

    # Sort the data frame by normZ (ascending) to highlight synthetic interactions
    per_gene_results.sort_values('normZ', ascending=True, inplace=True)

    # Calculate pvals (from normal dist), and fdrs (by benjamini & hochberg).
    per_gene_results['pval_synth'] = stats.norm.sf(per_gene_results['normZ'] * -1)
    per_gene_results['rank_synth'] = np.arange(1, no_of_genes + 1)
    per_gene_results['fdr_synth'] = per_gene_results['pval_synth'] * no_of_genes / per_gene_results['rank_synth']

    # resort by normZ (descending) and recalculate above values to identify suppressor interactions
    per_gene_results = per_gene_results.sort_values('normZ', ascending=False)
    per_gene_results['pval_supp'] = stats.norm.sf(per_gene_results['normZ'])
    per_gene_results['rank_supp'] = np.arange(1, no_of_genes + 1)
    per_gene_results['fdr_supp'] = per_gene_results['pval_supp'] * no_of_genes / per_gene_results['rank_supp']
    per_gene_results = per_gene_results.sort_values('normZ', ascending=True)
    return per_gene_results


def run_drugz_analysis(args, logger):
    """
    Executes the drug z analysis, writing the output to a tab separated file.

    :param args: argument parser object
    :param logger: console logger object
    :return: None
    """

    # Check input control and treatment samples are of equal length and convert them into lists.
    control_samples, treatment_samples = handle_sample_input(args, logger)

    # Load the read counts as a dataframe object
    reads = load_reads(filepath=args.infile, index_column=args.index_column, genes_to_remove=args.remove_genes)

    # Normalise the read counts.
    logger.debug('Normalizing read counts using normalisation value: {value}'.format(value=NORM_VALUE))
    normalised_counts = normalise_readcounts(reads=reads, treatment=treatment_samples, control=control_samples)

    logger.debug('Initialising fold change data frame. This will be used to store the calculated fold changes')

    # Instantiate the fold change dataframe (starts as just a column of gene names) - this is where we will store our
    # results of the intermediate stages of the analysis.
    fold_change = init_foldchange(reads=reads)

    # Get the number of comparisons we will be iterating through.
    num_replicates = len(control_samples)

    # Ids for each comparison
    zscore_ids = list()

    # For each comparison calculate the zscore for read counts of that pair of samples
    # Each comparison consists read counts for a single drug sample and a single control sample
    for i in range(num_replicates):
        fold_change, zscore_id = calculate_fold_change(reads=reads, normalised_counts=normalised_counts,
                                                       fold_change=fold_change, control_samples=control_samples,
                                                       treatment_samples=treatment_samples,
                                                       pseudocount=args.pseudocount,
                                                       half_window_size=args.half_window_size, logger=logger,
                                                       replicate=i)
        zscore_ids.append(zscore_id)

    # Here we have the option, if we've chosen to, to write out the intermediate file, after we've calculated the
    # zscores for all of the comparisons.
    if args.fc_outfile:
        with args.fc_outfile as fold_change_file:
            fold_change.to_csv(fold_change_file, sep='\t', float_format='%4.3f')

    # Calculate the per gene aggregated summary statistics for the zscores calculated in the previous comparison.
    logger.debug('Combining drugZ scores')
    drug_z_results = calculate_drugz_score(fold_change=fold_change, min_observations=args.minObs, columns=zscore_ids)

    write_drugz_output(file_handler=args.drugz_output_file, output=drug_z_results)


def write_drugz_output(file_handler, output):
    """
    Write the ouput of the drugz analysis.

    :param file_handler: A file handler into which to write the output of the current drugz analysis.
    :param output: A dataframe containing the results of the analysis
    :return: None
    """

    with file_handler as output_file:
        # Define column names
        header = ['GENE', 'sumZ', 'numObs', 'normZ', 'pval_synth', 'rank_synth', 'fdr_synth', 'pval_supp',
                  'rank_supp', 'fdr_supp']

        writer = csv.DictWriter(output_file, delimiter='\t', lineterminator="\n", fieldnames=header)
        writer.writeheader()

        for gene in output.index.values:
            # For each gene make a dictionary of value to write
            row = {column_name: output.loc[gene, column_name] for column_name in header[1:]}

            # Format the values of this row to match the original drugz formatting
            row = format_row(row)

            # Add the gene entry
            row[header[0]] = gene

            # Write the row
            writer.writerow(row)


def format_row(row):
    """
    Ensure the format of each row matches the original drugz output, by converting all values contained in the dct
    either to integers or formatted strings.

    :param row: A dictionary containing the raw values for each of the columns specified in the header in
    write_drugz_output(), except 'GENE'.

    :return: row: A dictionary with formatted values.
    """

    # This is not best practices as we are having to specify each column individually, so any modification to column
    # names is going to require changes in multiple places in the code. Maybe just output raw values?
    row['sumZ'] = '{value:3.2f}'.format(value=row['sumZ'])
    row['numObs'] = int(row['numObs'])
    row['normZ'] = '{value:3.2f}'.format(value=row['normZ'])
    row['pval_synth'] = '{value:.3g}'.format(value=row['pval_synth'])
    row['rank_synth'] = int(row['rank_synth'])
    row['fdr_synth'] = '{value:.3g}'.format(value=row['fdr_synth'])
    row['pval_supp'] = '{value:.3g}'.format(value=row['pval_supp'])
    row['rank_supp'] = int(row['rank_supp'])
    row['fdr_supp'] = '{value:.3g}'.format(value=row['fdr_supp'])

    return row


def init_logger(name, quiet):
    """
    Instantiate a logging object and set the logging mode, either verbose or quiet.

    :param name: The name of logger (str)
    :param quiet: Boolean saying whether or not to run in quiet mode (default is True), this is reverse of how it is
    normally done e.g. it usually is by default quiet and the argument passed would be -v or whether to be verbose or
    not.

    :return: logger: a logger object that will output to the console/terminal (stdout).
    """

    if quiet:
        level = log.INFO
    else:
        level = log.DEBUG

    # Format how the output is presented to the screen.
    console_logger_format = log.Formatter(fmt='%(asctime)s %(message)s', datefmt='%H:%M:%S')

    # Set up a handler for console log messages
    console_handler = log.StreamHandler()
    console_handler.setFormatter(console_logger_format)

    # Set up a logger object and add the handler to the logging object.
    logger = log.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(console_handler)

    return logger


def handle_sample_input(args, logger):
    """
    Verify the control and treatment samples are of equal length, returning these samples in list format.

    Note if these two lists are of unequal length, then the program will exit

    :param args: argument parser object
    :param logger: console logger object
    :return:
    control_samples: a list of control samples
    treatment_samples: a list of treatment samples
    """
    control_samples = args.control_samples.split(',')
    treatment_samples = args.drug_samples.split(',')

    if len(control_samples) != len(treatment_samples):
        raise Exception("Must have the same number of control and drug samples.")
    else:
        logger.debug('Control samples:  {control_samples}'.format(control_samples=control_samples))
        logger.debug('Treated samples:  {drug_samples}'.format(drug_samples=treatment_samples))
        return control_samples, treatment_samples


def main():
    args = get_args()
    name = (os.path.basename(__file__)).split('.')[0]

    logger = init_logger(name=name, quiet=args.quiet)

    run_drugz_analysis(args=args, logger=logger)


if __name__ == "__main__":
    main()
