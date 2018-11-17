This file gives instructions about how to perform gene set analysis using magma. It was written by Dr Sarah Morgan (17/11/2018) and is based on code and directions kindly provided by [Dr Nicholas Clifton](https://www.cardiff.ac.uk/people/view/105079-clifton-nicholas).

# Downloads:

First, download the summary statistics from [Pardinas et al](https://doi.org/10.1038/s41588-018-0059-2) - http://walters.psycm.cf.ac.uk/.

Then download the magma software from: https://ctg.cncr.nl/software/magma. I found it helpful to use the static linking version.

You also need to download the gene locations file and the reference data from the same page. The build needs to match the build in which the summary stats were created. I used build 37 to match the summary stats from the Pardinas paper. We used European reference data because the brain images were collected in Europe.

# Preparing the input files:

The code can be run from the command line. Before running the code there are a few preparation steps:

The file `clozuk_pgc2.meta.sumstats.txt' have an "IMPUTE2" format (RSid:BP:A1:A2), whereas you just need the RSids (to match the g1000 file). In linux, to replace ids with RSids only:
```
sed 's/:[0-9A-Za-z:]*//g' clozuk_pgc2.meta.sumstats.txt > My_sumstats.txt
```

Then remove any lines which don't begin with 'rs':
```
awk '/^rs/' My_sumstats.txt > My_sumstats_no_rs.txt
```

Then put the header back in again:
```
echo -e "SNP\tFreq.A1\tCHR\tBP\tA1\tA2\tOR\tSE\tP" | cat - My_sumstats_no_rs.txt > My_sumstats_no_rs_head.txt
```

For the analysis you also need a .snp.loc file, which can be obtained by selecting columns 1, 3 and 4 from summarystats (SNP id, CHR and BP). We don't want the header here. The command to use is:
```
awk '{ print $1, $3, $4 }' My_sumstats_no_rs.txt > PGC_summarystats.snp.loc
```

Finally, prepare a file which contains your genesets. In the code below, `genesets.txt' is a file with two rows- the first row contains our first geneset (genes where Z>3) and the second row contains our second geneset (genes where Z<-3). Genes are given by their entrez IDs and should be tab-delimited. The first column should give the geneset names. You can use however many genesets you like.


# Run the analysis:

Now you're ready to run the analysis! There are three steps, as follows. Note that these commands require the inputs given above and you will have to change the paths etc to make sure it can find your input files.

Step 1:
```
magma --annotate window=35,10 --snp-loc ../clozuk_pgc2/PGC_summarystats.snp.loc --gene-loc ../NCBI37/NCBI37.3.gene.loc --out SUM_STATS.1
```

Step 2:
```
magma --bfile g1000/g1000_eur --pval ../clozuk_pgc2/My_sumstats_no_rs_head.txt use=SNP,P N=105318 --gene-annot SUM_STATS.1.genes.annot --out SUM_STATS.2
```

Step 3:
```
magma --gene-results SUM_STATS.2.genes.raw --set-annot genesets.txt --settings gene-include=maxintensity_bground.txt --out GSA_OUT
```
