## Key Learning Outcomes

After completing this practical the trainee should be able to:

-   Run the `IntOGen` analysis software on cohort mutation data.

-   Have gained experience of the structure of the analysis output files
    in order to identify potential driver genes.

-   Have gained overview knowledge of different methods for
    identification of genes important in cancers.

***
##Resources You’ll be Using

### Tools Used

IntOGen mutations platform:  
https://www.intogen.org/search

### Sources of Data

TCGA melanoma somatic SNV data from 338 tumour samples:  
https://tcga-data.nci.nih.gov/tcga/

### Useful Links

Mutation Annotation Format (MAF) specification:  
https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification

IntOGen installation instructions:
https://bitbucket.org/intogen/intogen-pipeline/overview

***
## Author Information

*Primary Author(s):*  
Ann-Marie Patch, QIMR Berghofer ann-marie.patch@qimrberghofer.edu.au  
Erdahl Teber, CMRI eteber@cmri.org.au  

*Contributor(s):*  
Scott Wood scott.wood@qimrberghofer.edu.au  

***
## Introduction

Cancer driver genes are commonly described as genes that when mutated
directly affect the potential of a cell to become cancerous. They are
important to a tumour cell as they confer a growth or survival advantage
over the normal surrounding cells. The mutations in these driver genes
are then clonally selected for as the population of tumour cells
increases. We think of the key genes driving tumour initiation
(development), progression, metastases, resistance and survival. Driver
gene mutations are often described as "early" events because they were
key in turning a normally functioning and regulated cell into a
dysregulated one. The logical assumption is that these key mutations
will be present in all tumour cells in a patient’s sample; although
sometime this is not true.

There are two major research goals that underline the need to identify
driver genes:

-   By identifying the early changes that take place researchers might
    be able to find a treatment to stop the root cause of why cells
    become malignant.

-   By identifying groups of patients with the same genes mutated then
    we can develop therapies that will work for all of them.

When we sequence tumour samples we tend to use samples that come from
fully developed cancers that can carry hundreds to thousands of
mutations in genes and many more outside of genes. The accumulation of
these passenger mutations in cancer cells can happen because often the
repair mechanisms or damage sensing processes are amongst the first
pathways to become disrupted accelerating the mutational rate. Mutations
that occur in genes after the cell has become cancerous may still affect
the growth rate, invasiveness and even the response to chemotherapy but
may not be present in all cells of a tumour. These genes may be drivers
of chemo-resistance or metastasis and are equally good targets for
therapies.

IntOGen-mutations is a platform that aims to identify driver mutations
using two methodologies from cancer cohort mutation data: the first
identifies mutations that are most likely to have a functional impact
combined with identifying genes that are frequently mutated; and the
second, genes that harbour clustered mutations. These measures are all
indicators of positive selection that occurs in cancer evolution and may
help the identification of driver genes.

***
## Analysing cancer cohort data with IntOGen

IntOGen-mutations is available as a web based service that can allow
users to run their analysis on the host’s servers or it can be downloaded
and run on a local server.

For the purposes of the course we will be using a local version of
`IntOGen` so that we don’t encounter any issues sharing resources.

- To begin open a terminal and navigate to the directory `somatic/intogen`.

    ```bash
    cd ~/somatic/intogen
    ```

In this directory you will find a Mutation Annotation Format (MAF) file
containing a cut down version of the somatic variant calls identified from
melanoma samples investigated as part of the TCGA cancer genomics projects.
You can see what files are in the directory by typing `ls`, look inside the
file using `less TCGA_Melanoma_SMgene.maf` and close the file and return to
the command line by typing `q`.

- Run the `IntOGen` analysis by typing

    ```bash
    intogen -i TCGA_Melanoma_slimSMgene.maf -o TCGA_Mela_out
    ```

The TCGA melanoma maf used in this practical has been modified from the
original to reduce processing time and only contains data for the top
680 mutated genes.

The tool will take around 10 minutes to run and the progress will be
indicated by the logging lines printed to the terminal. Once complete
the output can be explored.

Whilst the tool is running we can explore the options we have used to
run `IntOGen`.

- To get a list of `IntOGen` options open up a new terminal

  ```bash
  intogen --help
  ```   

This command will list the running options that you can alter as command line
inputs or in a configuration file. We are using the default options for this run
so we didn't have to supply a configuration file and we only used `-i` to set
the input and `-o` to control the mane of the output directory.

- To look at the default options open up the configuration file by typing  

```bash
less ~./intogen/task.conf  
less ~./intogen/system.conf
```


It is important to set the correct genome assembly in the `task.conf` to match
the one that you used as your reference when the variant were called. In our
`task.conf` this should be `hg19`.

***
## Exploring the output of IntOGen

When you run your data over the web on the remote site there is a browse
facility that allows you to explore your data using the web version of
the database. Running `IntOGen` locally provides the same tabular
information but in a flat file format.

There should be 14 files generated from a successful run of this version of `IntOGen`:

	gene.tsv
	gene.oncodriveclust
	pathway.recurrences
	gene.oncodrivefm
	sample_gene.impact
	gene.recurrences		     			     	
	sample_variant+transcript.impact
	summary.tsv
	transcript.recurrences
	TCGA_Melanoma_slimSMgene.smconfig
	oncodrivefm-pathways-MA_SCORE.tsv
	oncodrivefm-pathways-PPH2_SCORE.tsv  
	oncodrivefm-pathways-SIFT_SCORE.tsv  
	pathway.oncodrivefm


View these files by using `ls` as below.

  ```bash  
  ls ~/somatic/intogen/TCGA_Mela_out/project/TCGA_Melanoma_slimSMgene/
  ```

This practical will concentrate on the identification of driver genes so
we will look at the main output concerning genes.
The `gene.tsv` is the main gene centric output summary table.

- Open up the `gene.tsv` file in `LibreOffice` by double clicking on the icon on your
desktop.

- Select the file tab and click on open.

- Navigate to the results directory
`~/somatic/intogen/TCGA_Mela_out/project/TCGA_Melanoma_slimSMgene/`

- Double click on `gene.tsv`.

- In the pop-up box under the `Separator options` ensure only the tab box is
checked and click `OK`.

This file contains the overall summary results for the `IntOGen` pipeline
presented by gene and reports Q values (i.e. multiple testing corrected
P values) for the mutation frequency and cluster modules.

Significantly mutated genes from the cohort data are identified using both
the `OncodriveFM` and `OncodriveClUST` modules of `IntOGen`. The `OncodriveFM`
module detects genes that have accumulated mutations with a high functional
impact. It uses annotations from the Ensembl variant effect predictor (VEP, V.70)
that includes SIFT and Polyphen2 and precomputed MutationAssessor functional
impacts. It calculates a P value per gene from the number of mutations detected
across all possible coding bases of a gene with a positive weighting for mutations
with a high functional impact. The `OncodriveCLUST` module detects genes that
have more variants than would be expected across the cohort that alter the
same region of the gene.


The file is sorted to bring the most significantly altered genes to the top. The key
columns that help you identify the significantly mutated genes are the 3rd and 4th
(C and D) that indicate which of the modules identified a significant result and
the Q-values for the modules that are in 21st and 23rd (U and W)

The top twelve genes have significant Q-values for both modules and include BRAF,
NRAS and TP53. The next 35 are significant by only one of the modules.

All of these have small Q-values which means they are all
significantly mutated genes in this TCGA Melanoma cohort of 338
patients.

- Now look at their sample frequency count (column 9 `MUTS_CS_SAMPLES`) these
are the number of samples that contain at least one mutation in the gene.

!!! note "Question"
    a)   Which significantly mutated gene has mutations in the most samples?

	b)   Which gene/genes have the lowest Q-value from OncodriveFM and OncodriveCLUST?

    c)   Why don’t the genes with the lowest Q values also have the highest sample frequency value?

    !!! success ""
        ??? "**Answer**"
            a)  BRAF has 175 out of 327 cases with a mutation.

            b)  TP53 or PTEN have the lowest OncodriveFM Q-values and NRAS
            has the lowest Q-value for OncodriveCLUST.

            c)  The P value calculation takes into account the length of the
                coding sequence of the gene, the mutation rate of the nucleotides
                contained within it and for OncodriveFM the functional consequences
                of those changes. Therefore a small gene with a small number of deleterious
                mutations may have a lower P value and also Q value than a large
                gene with a high mutation frequency.


The results for the assessment of clustered mutations in genes carried
out by the `OncodriveCLUST` module of `IntOGen` are shown as  amino
acid residue positions of the encoded protein.

The three known oncogenes BRAF, NRAS and IDH1 have very low CLUST_QVALUEs
indicating that the mutations in these genes are highly clustered. The
`CLUST_COORDS` column reports that there are 160 samples with mutations
between the amino acid positions 594-601 of BRAF; 84 samples with mutations
at amino acid position 61 of NRAS; and 15 sample with mutations at amino
acid position 132 of IDH1.

!!! note "Question"
    Why are the oncogenes more likely to have clustered mutations and the
    tumour suppressor genes less likely?

    !!! success ""
        ??? "**Answer**"
            Gain of function mutations are required to activate oncogenes and so
            only key residues in the protein will result in activation. Tumour
            suppressors are frequently affected by loss of function mutations and
            deletions. A truncating mutation or frameshift indel can occur in any
            exon, except the last one, and have the same deleterious functional
            result.

<br>
The other files in the output support the information in this sheet.

The ` sample_variant+transcript.impact` file includes a summary of all mutations
found in each of the genes and protein coding transcripts of those genes for all
samples identified that have that mutation. It also reports the variant impact scores
from SIFT, PolyPhen2, MutationAssesor, reporting also impact categories of which
there are four; high, medium, low and none.

- Open up the ` sample_variant+transcript.impact` file and explore the data.

!!! note "Question"
    Can you find out what the nucleotide change details for the most common
    BRAF mutation that results in V600E amino acid change in the cohort?
	Sort the data by GENE, then TRANSCRIPT and then PROTEIN_POS to make this easier. The gene ID for BRAF is ` ENSG00000157764`.

    !!! success ""
        ??? "**Answer**"
            It is an A>T at position chr7:140453136 identified in 127 samples.

***
## References

Gunes et al. Nat. Methods 2010
:   http://www.nature.com/nmeth/journal/v7/n2/pdf/nmeth0210-92.pdf

Gonzalez-Perez et al. Nat. Methods 2013
:   http://www.nature.com/nmeth/journal/v10/n11/pdf/nmeth.2642.pdf
