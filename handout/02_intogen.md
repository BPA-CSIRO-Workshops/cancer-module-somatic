Key Learning Outcomes
---------------------

After completing this practical the trainee should be able to:

-   Run the IntOGen analysis software on cohort mutation data.

-   Have gained experience of the structure of the analysis output files
    in order to identify potential driver genes.

-   Have gained overview knowledge of different methods for
    identification of genes important in cancers.

Resources You’ll be Using
-------------------------

### Tools Used

IntOGen mutations platform
:   \
    <https://www.intogen.org/search>

### Sources of Data

TCGA melanoma somatic SNV data from 338 tumour samples
:   \
    <https://tcga-data.nci.nih.gov/tcga/>

### Useful Links

Mutation Annotation Format (MAF) specification
:   \
    <https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification>

Introduction
------------

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

Analysing cancer cohort data with IntOGen
-----------------------------------------

IntOGen-mutations is a web platform that can allow users to run their
analysis on the host’s servers or it can be downloaded and run without
any limits on the number of analyses on a local server.

For the purposes of the course we will be using a local version of
IntOGen so that we don’t encounter any issues sharing resources.

To start IntOGen open a terminal and navigate to the
somatic/intogen\_mutations\_analysis directory.

    cd ~/somatic/intogen_mutations_analysis

start IntOGen by typing

    ./run web

A browser window should open up with the IntOGen home page displayed.
IntoGen requires you to log in to Mozilla persona.

On the top right hand side of the screen there will be a "Sign in"
button.

Click the `Sign in` button.

If you can access your email from the web you can create yourself an
account. You do need to verify your email address by logging in to your
webmail.

Alternatively, you can log in using these dummy account details.

Email address:

amp.genomics@gmail.com.

Password if required:

letmepass.

Click on the `Analysis` button on the top tool bar. Click on the
`Analyse your data` button in the `Cohort analysis` panel.

A form will open up, click on the `Load` button and navigate to the
TCGA\_Melanoma\_SMgene.maf file located in the directory path
`somatic/intogen/ TCGA_Melanoma_SMgene.maf`.

Ensure the Genome assembly is set to `hg19(GRCh37)`.

The click the `Start analysis` button.

The TCGA melanoma maf used in this practical has been modified from the
original to reduce processing time and only contains data for the top
680 mutated genes.

The tool will take around 10 minutes to run and the progress will be
indicated by the green bar. You may need to refresh the browser if it
looks like there has been no progress for a while.

Once complete the output can be downloaded.

Click on the download link to the left of the green bar or go through
the results button at the top of the page.

Select to `save` the file which will end up in the download directory
` /Downloads`.

Click on the blue download arrow at the top right hand side of the
browser window.

Click on the file called `TCGA_Melanoma_slimSMgene.zip`.

Highlight the file name in the pop-up window and click on `Extract`.

Select the `/trainee/somatic/intogen/` directory and click on the
`Extract` button at the bottom of the window.

Click to show files once extracted.

Double click on the directory called `TCGA_Melanoma_slimSMgene`.

Exploring the output of IntOGen
-------------------------------

When you run your data over the web on the remote site there is a browse
facility that allows you to explore your data using the web version of
the database. Running IntOGen locally provides the same tabular
information but in a flat file format.

There are 7 files generated in a successful run of IntOGen:

-   genes.tsv - this is the main output summary table.

-   variant\_genes.tsv - describes the variants identified in the genes

-   variant\_samples.tsv - lists the variants together with the samples
    Ids

-   pathways.tsv - lists the Kegg pathway ID for perturbed pathways

-   consequences.tsv - lists consequences of the variants in all known
    transcripts

-   project.tsv - one line to summarise the input data

-   fimpact.gitools.tdm - the output of the functional impact tool

The scope of this practical will concentrate on identification of driver
genes so we will look at the main output concerning genes. OncodriveFM
does however also calculates a functional impact bias of high impact
mutations in annotated Kegg pathways
<http://www.genome.jp/kegg/pathway.html> as shown in the `pathways.tsv`
file but will not be used in this course.

Open up the `genes.tsv` in the spreadsheet software by double clicking
on the file.

This file contains the overall summary results for the IntOGen pipeline
presented by gene and reports Q values (i.e. multiple testing corrected
P values) for the mutation frequency and cluster modules.

Significantly mutated genes from the cohort data are identified using
the OncodriveFM module of IntOGen. This tool detects genes that have
accumulated mutations with a high functional impact. It uses annotations
from the Ensembl variant effect predictor (VEP, V.70) that includes SIFT
and Polyphen2 and precomputed MutationAssessor functional impacts. It
calculates a P value per gene from the number of mutations detected
across all possible coding bases of a gene with a positive weighting for
mutations with a high functional impact.

Sort the data in this file by two levels starting with the 5th column,
the FM\_QVALUE from smallest-to-biggest and then by the 12th column the
CLUST\_QVALUE also from smallest-to-biggest.

The top nine genes with the smallest Q values should be TP53, PTEN,
PPP6C, CDKN2A, BRAF, NRAS, ARID2, TTN, IDH1.

All of these have very small FM\_Qvalues which means they are all
significantly mutated genes in this TCGA Melanoma cohort of 338
patients.

Now look at their sample frequency (column 9 SAMPLE\_FREQ) values these
are the number of samples that contain at least one mutation in the
gene.

-   ​a) Which gene has mutations in the most samples?

-   ​b) Which gene had the lowest FM Q value?

-   ​c) Why don’t the genes with the lowest Q values also have the
    highest sample frequency value?

-   ​a) TTN has the highest number of samples with mutations. There are
    265 out of 327 samples with mutations in TTN.

-   ​b) TP53 or PTEN

-   ​c) The P value calculation takes into account the length of the
    coding sequence of the gene, the mutation rate of the nucleotides
    contained within it and the functional consequences of those
    changes. Therefore a small gene with a small number of deleterious
    mutations may have a lower P value and also Q value than a large
    gene with a high mutation frequency.

The results for the assessment of clustered mutations in genes carried
out by the OncodriveCLUST module of IntOGen are shown in the 10-14th
columns. The Q value is in column 12, CLUST\_QVALUE which indicates if
there is a significant grouping of mutations identified. The positions
of the mutations are now described in terms of the amino acid residue in
the encoded protein.

The three oncogenes BRAF, NRAS and IDH1 have very low CLUST\_QVALUEs
\<0.01 indicating that the mutations in these genes are highly
clustered. The CLUST\_COORDS column reports that there are 158 samples
with mutations between the amino acid positions 594-601 of BRAF; 84
samples with mutations at amino acid position 61 of NRAS; and 15 sample
with mutations at amino acid position 132 of IDH1.

Why are the oncogenes more likely to have clustered mutations and the
tumour suppressor genes less likely?

Gain of function mutations are required to activate oncogenes and so
only key residues in the protein will result in activation. Tumour
suppressors are frequently affected by loss of function mutations and
deletions. A truncating mutation or frameshift indel can occur in any
exon, except the last one, and have the same deleterious functional
result.

The INTOGEN\_DRIVER column indicates if this gene is a known cancer
driver gene with 1 for yes and 0 for no so it is promising to see the
majority of our top hit genes are known drivers.

The XREFS column indicates a mutation match in external databases. This
means the position and the base change has been seen before. For the
known Driver genes there are many COSMIC IDs indicating these mutations
have been recorded in the Catalogue of Somatic Mutations in Cancer
database <http://cancer.sanger.ac.uk/cosmic>. Remember the majority of
data in these databases have come from large scale cancer sequencing
projects carried out by TCGA and ICGC associated groups.

Have a look at the values in the XREFS field for the gene TTN.

There are a large number of COSMIC entries for mutations in this gene
but also a large number of ESP ID numbers that refer to the National
Heart, Lung, and Blood Institute (NHLBI), Exome Sequencing Project (ESP)
<http://evs.gs.washington.edu/EVS/>.

This data comes from a diverse population that typically don’t have
cancer but focus on patients with heart, lung and blood disorders.

Can we draw any conclusions from the high number of ESP mutations and
dbSNP references for TTN?

Functional studies would be the only way to prove conclusively if TTN
mutations were cancer driver mutations. The high number of samples in
the ESP cohort that have mutations in common with the cancer cohorts
could indicate that TTN mutation may not be unique to the cancer
setting. This conclusion seems to be backed up by the high number of
mutations also with dbSNP ids that indicates the potential for some of
these mutations to be present in the general population.

Amongst this high number of mutations there may be groups of patients
where the mutations are purely passenger and other groups where the
mutations could contribute to the tumour.

From this data there is no way to tell.

The other files in the output support the information in this sheet.

The `variant_genes.tsv` file includes a summary of all mutation found in
each of the genes with a count of the number of samples identified to
have that mutation. It also reports the variant impact score and
category of which there are four; high, medium, low and none.

Open up the `variant_genes.tsv` file and explore the data.

Can you find out what the nucleotide change details for the most common
BRAF mutation that results in V600E amino acid change in the cohort?
Sort the data on the gene symbol column to make this easier.

It is an A\>T at position chr7:140453136 identified in 127 samples.

The `variant_samples.tsv` worksheet allows you to find out the sample
identification numbers for the samples with the mutation that you are
interested in.

Open up the `variant_samples.tsv` file and explore the data.

Can you use the position that you found in the question above to locate
the sample IDs for the 3 cases with a V600R BRAF mutation that is caused
by a nucleotide change of AC/CT starting at chr7:140453136?

TCGA-D3-A1QA-06A-11D-A196-08, TCGA-EB-A3XC-01A-11D-A23B-08,
TCGA-ER-A193-06A-12D-A197-08

References
----------

Gunes et al. Nat. Methods 2010
:   \
    <http://www.nature.com/nmeth/journal/v7/n2/pdf/nmeth0210-92.pdf>

Gonzalez-Perez et al. Nat. Methods 2013
:   \
    <http://www.nature.com/nmeth/journal/v10/n11/pdf/nmeth.2642.pdf>

