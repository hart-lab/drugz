# drugz
DrugZ software  
  
DrugZ detects synergistic and suppressor drug-gene interactions in CRISPR screens.  
  
usage: drugz.py [-h] [-i sgRNA_count.txt] [-o drugz-output.txt]  
                [-f drugz-foldchange.txt] -c control samples -x drug samples  
                [-r remove genes] [-p pseudocount] [-I INDEX_COLUMN]  
                [--minobs minObs] [-q]  
  
-i      	Readcount file, tab-delimited text (input)  
-o      	DrugZ results file, tab-delimited text (output)  
-f      	DrugZ Z-transformed fold change file (optional)  
-c      	Control samples: comma-delimited list of column headers in readcount file  
-x      	Treated samples: comma-delimited list of column headers in readcount file  
-r      	(optional) comma-delimited list of genes to remove before analysis  
-p      	(default=5) pseudocount to add to all readcounts; prevents log(0) problems  
-I      	Index column (default=0)  
--minobs   	(default=6) Ignore genes with fewer observations ( gRNA/gene x replicates)  

