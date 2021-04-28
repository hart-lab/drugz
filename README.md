# drugz
DrugZ software from the Hart Lab
  
DrugZ detects synergistic and suppressor drug-gene interactions in CRISPR screens.  

```
usage: drugz.py [-h] [-i sgRNA_count.txt] [-o drugz-output.txt]  
                [-f drugz-foldchange.txt] -c control samples -x drug samples  
                [-r remove genes] [-p pseudocount] [-I INDEX_COLUMN]  
                [--minobs minObs] [--half_window_size half_window_size] [-q]  
  
-i      	Readcount file, tab-delimited text (input)  
-o      	DrugZ results file, tab-delimited text (output)  
-f      	DrugZ Z-transformed fold change file (optional)  
-c      	Control samples: comma-delimited list of column headers in readcount file  
-x      	Treated samples: comma-delimited list of column headers in readcount file  
-r      	Comma-delimited list of genes to remove before analysis  
-p      	Pseudocount to add to all readcounts; prevents log(0) problems (default=5) 
-I      	Index column (default=0)  
--minobs   	Ignore genes with fewer observations ( gRNA/gene x replicates) (default=1) 
--half_window_size  Size of the first bin and half the size of the inital sample
    (window) to estimate std (default=500) 
-unpaired Unpaired approach: compares mean(treated samples) to mean(control samples) (default=False)
```
  
The input file should be a tab-delimited file with the following format:

```
sgRNA	Gene	T0	T15_A_control	T15_B_control	T15_C_control	T15_A_olaparib	T15_B_olaparib	T15_C_olaparib
A1BG_CACCTTCGAGCTGCTGCGCG	A1BG	313	235	47	337	428	115	340
A1BG_AAGAGCGCCTCGGTCCCAGC	A1BG	99	8	1	13	26	5	28
A1BG_TGGACTTCCAGCTACGGCGC	A1BG	650	336	74	185	392	193	304
A1BG_CACTGGCGCCATCGAGAGCC	A1BG	718	192	34	296	178	69	185
A1BG_GCTCGGGCTTGTCCACAGGA	A1BG	180	230	29	122	394	148	364
A1BG_CAAGAGAAAGACCACGAGCA	A1BG	428	300	158	294	366	184	489
A1CF_CGTGGCTATTTGGCATACAC	A1CF	677	452	74	423	585	446	434
A1CF_GGTATACTCTCCTTGCAGCA	A1CF	138	69	43	109	96	184	127
A1CF_GACATGGTATTGCAGTAGAC	A1CF	396	183	38	106	193	120	198
(etc)
```

Critically, the "gene" column must be the first non-index column in the file, and the column headers are used on the command line. For example, to execute DrugZ analyzing just the A and B replicates of this file, the command line would be:

```
drugz.py -i [input_file] -o drugz-output.txt -c T15_A_control,T15_B_control -x T15_A_olaparib,T15_B_olaparib
```

To save the intermediate gRNA-level raw and normalized fold changes for other analyses, add the -f flag:

```
drugz.py -i [input_file] -o drugz-output.txt -f drugz-foldchange.txt -c T15_A_control,T15_B_control -x T15_A_olaparib,T15_B_olaparib
```
To run drugZ for an unpaired approach, add the -unpaired flag:

```
drugz.py -i [input_file] -o drugz-output.txt -f drugz-foldchange.txt -c T15_A_control,T15_B_control -x T15_A_olaparib,T15_B_olaparib -unpaired
```

To run drugZ analysis in a jupyter notebook, and save the output as variable:

```
# define the Arguments class (more convinient since iPython doesn't recognize argparse arguments)
# these are user-specified arguments

# infile = input readcounts matrix
# drugz_out_file = name of a file in which you will write the drugz results
# control_samples = the names of control samples (included in column names)
# drug_samples = the names of drug-treated samples (included in column names)
# unpaire = unpaired approach - compares mean(treated samples) to mean(control samples) 
# pseudocount = counts added to the observed readscounts, default = 5
# half_window_size = size of the first bin and half the size of the inital sample (window) to estimate std, default = 500 (for whole genome screens)

class Args:
    infile = "./sgRNA_count.txt"
    drugz_output_file = "./drugz_results.txt"
    fc_outfile = "./fc_results.txt"
    control_samples = 'T15_A_control,T15_B_control,T15_C_control'
    drug_samples = 'T15_A_olaparib,T15_B_olaparib,T15_C_olaparib'
    remove_genes = 'LacZ,luciferase,EGFR'
    unpaired = False
    pseudocount = 5
    half_window_size = 5 # 5 because of the size of test data set          (sgRNA_count.txt = 9 guides (i.e. rows))
    
drugz_results = dz.drugZ_analysis(Args())
```

For more option check drugZ_in_jupyter_notebook_tutorial.html
