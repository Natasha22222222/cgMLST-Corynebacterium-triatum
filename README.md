# cgMLST-Corynebacterium-triatum
# cgMLST to *Corynebacterium striatum*

The objective of this repository is to describe how we created a cgMLST for *C. striatum*. This scheme was created with the ChewBBACA pipeline (link below). Input files used to create and validate the scheme as well as output files from each step are available. 

Additionally, the steps to run the scheme and/or download the *C. striatum* selected target genes that compose the cgMLST scheme are also available in this repository. (Steps 7, 8 and 9)

## chewBBACA

To download ChewBBACA access the link:

* https://github.com/B-UMMI/chewBBACA_tutorial

## Workflow used to create the scheme

* Step 1: Scheme Creation
* Step 2: Allele calling
* Step 3: Scheme Validation (Allele call)
* Step 4: Extracting the Matrix target genes
* Step 5: Minimum Spanning Tree (MST)
* Step 6: Graphical evaluation of the scheme
* Step 7: How to identify the allelic profile of the genomes of interest using this cgMLST database for *C. striatum*
* Step 8: Analyze the results
* Step 9: Allele profile view by minimum spanning tree (MST)
 
## Softwares and Downloads (Main dependencies)

* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* Prodigal 2.6.0 https://github.com/hyattpd/prodigal/releases/ or above

## Other dependencies (for schema evaluation only):

* ClustalW2 http://www.clustal.org/download/current/
* MAFFT https://mafft.cbrc.jp/alignment/software/linux.html
* DATAMASH https://www.gnu.org/software/datamash/
* MLST https://github.com/sanger-pathogens/mlst_check

## Step 1: Schema Creation

## Selection of publicly available genomes for schema creation

As of August 29, 2021, 271 publicly available genome sequences were available at the NCBI (National Center for Biotechnology Information) genome sequence repository (https://www.ncbi.nlm.nih.gov/assembly).
A list of 271 publicly available genomes used to create this protocol can be found in the file 271-publicly-available-genome.txt.
C. striatum WP1a (RefSeq assembly accession: GCF_004138065.1) was used by Prodigal algorithm as reference to recognize coding sequences (CDs). Prodigal generated the GCF_004138065.1.trn file at this step.

## Step 1: Definition of CDs sequences  

## Command: 

```bash
# create schema
chewBBACA.py CreateSchema -i publicly-available-genome --cpu 46 -o schema_seed --ptf GCF_004138065.1.trn
```

The above command uses 46 CPU and creates a preliminary scheme (wgMLST) in the schema_seed folder using the trained product file GCF_004138065.1.trn that was generated using the reference genome WP1a (GCF_004138065.1) and 271 publicly available genome sequences. The wgMLST scheme generated contained 4934 *loci* based on the 271 publicly available genomes. 

**Note:** publicly-available-genome: Folder containing the list 271 publicly available genomes that created the schema.

## Step 2: Allele calling

In this step the allele calling is performed using the resulting set of *loci* determined in step 1.

## Command: 

```bash
# run allelecall
chewBBACA.py AlleleCall -i publicly-available-genome -g schema_seed/ -o results_cg --cpu 46 --ptf GCF_004138065.1.trn
```

The allele calling used the default BLAST Score Ratio (BSR) threshold of 0.6.

## Step 2.1: Paralog detection

In this step genes considered paralogous from result of the allelecall (see above) are removed

## Command: 

```bash
# run remove genes
chewBBACA.py RemoveGenes -i results_cg/results_20221004T122439/results_alleles.tsv -g results_cg/results_20221004T122439/RepeatedLoci.txt -o results_cg/alleleCallMatrix_cg
```
In this step, 66 genes were identified as possible paralogs and were removed from further analysis. The list RepeatedLoci.txt can be found at: ```results_cg/results_20221004T122439/RepeatedLoci.txt```

The output file can be found at: ```analysis_cg/alleleCallMatrix_cg.tsv```

## Step 2.2: Genome Quality Control

In this step we define a *Threshold* for the scheme that limits the loss of *loci* targets defined in the previous steps per genome and excludes genomes considered to be of low quality due to significant *loci* absence. 

With this analysis we define the percentage of *loci* that will constitute the scheme based on how many targets we want to keep in this phase. For example,  **100%, 99.5%, 99% and 95%** of the *loci* may present in the set of high quality genomes. This is one of the main steps in defining the cgMLST scheme targets.

## Command:

```bash
# run test genome quality
chewBBACA.py TestGenomeQuality -i alleleCallMatrix_cg.tsv -n 13 -t 300 -s 5
```

In the *Threshold* 15 a set of 1917 *loci* were found to be present in the 95% analyzed genomes. The output file can be found in the folder: ```analysis_cg/GenomeQualityPlot.html```. The list with the genes present in 95% of the genomes at the chosen *Threshold* can be retrieved in the folder ```analysis_cg/Genes_95%.txt```. 

In this stage we chose the *loci* present in 95% (*t0.95*) of the 271 publicly available genomes and the *Threshold* 15 that limited the loss of the *loci* in the genomes. In this *Threshold* (15) 38 genomes were removed due to loss of *loci* . The list of genomes removed at each *Threshold* can be retrieved in the folder ```analysis_cg/removedGenomes.txt```. From this list we created another (GenomeRemoved15thr.txt) with only the genomes removed at *Threshold* (15). The list of genomes removed at *Threshold 15* can be retrieved in the folder: ```analysis_cg/GenomeRemoved15thr.txt```

## Command:

```bash
# run ExtractCgMLST
chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_15 --t 0.95 --g GenomeRemoved15thr.txt
```

This script selects all * loci * present in the selected * Threshold *. The value * t * is associated with the percentage of * loci * we want to be present in the set of genomes, for example: * t1.0 * selects all * loci * present in the * Threshold * chosen in all genomes ie those present in 100% of genomes at that * Threshold *. Subsequently a cgMLST_15 folder receives the result of the allelic profile for each of the 1917 candidate * loci * (allelic profile matrix). The file in this folder (cgMLST.tsv) contains the allelic profile of 1917 selected * loci * and will be used to create the core gene list. In addition, another 3 output files are created at this point: *cgMLSTschema.txt; mdata_stats.tsv and Presence_Absence.tsv* and they can be found in the folder ```cgMLST_15/```.

## Step 2.3: Creating the core gene list

This command selects all target genes from the "cgMLST.tsv" spreadsheet.

```bash
# 10 list
head -n 1 cgMLST.tsv > Genes_95%_Core_15.txt
```
This step generated the file *Genes_95%_Core_15.txt* that can be retrieved in the folder ```results_cg/cgMLST_15/Genes_95%_Core_15.txt```. This list needs to be transposed so that each core gene name is reported in a single line:

## Command:

```bash
# transpose table
datamash -W transpose < Genes_95%_Core_15txt > Genes_Core_Al.txt 
```

This step generated the file > Genes_Core_Al.txt

You can see the list file with 1917 target genes at ```results_cg/cgMLST_15/Genes_Core_Al.txt``` and for the subsequent steps we added the full path to each locus fasta file.

This list ```results_cg/cgMLST_15/Genes_Core_Al.txt``` was then modified so that each name was preceeded by *schema_seed*:

## Command:

```bash
tail -n+1 Genes_Core_Al.txt | cut -f2 | perl -pe 's: :\n:g' | sort -Vu | awk '{print("/"$1)}' > list_genes_core.txt
```
This modified list can be found: ```results_cg/cgMLST_15/list_genes_core.txt```.

## Step 3: Scheme Validation (Allele calling)

For the validation step we selected 28 *C. striatum* strains from patients with hospital-acquired infection, 2 *C. striatum* strains from patients with community-acquired infections and 1 epidemiologically unrelated outgroup *C. striatum* strains. The list of all the Validation genomes used can be found in the file ```Genomes_Validation.txt```.
We this set of genomes (31 validation genomes) we repeated the allele call step using only the 1917 candidate target genes.

## Command:

```bash
chewBBACA.py AlleleCall -i ../Genomes_Validation/ --gl list_genes_core.txt -o ../results_all --cpu 46  -g ../schema_seed/schema_seed –ptf GCF_004138065.1.trn
```

The folder **Genomes_Validation** contains the 31 validation drafts genomes used for validation of the scheme.

The folder with the output file can be found at: ```results_all/ results_20221004T143059/```. This folder contains 5 files: "logging_info.txt; RepeatedLoci.txt; results_alleles.tsv; results_contigsInfo.tsv and results_statistics.tsv".

The ```results_all/results_20221004T143059/results_alleles.tsv``` file contains the allelic profile of the 31 validation genomes.

## Step 3.1: Concatenate the allelic profiles

The purpose of this step is to concatenate the matrix of the *loci* that defined the scheme and the matrix of the *loci* from the validation genomes. Thus, to concatenate the allelic profile matrix obtained from the creation of the scheme ```cgMLST_15/cgMLST.tsv``` with the matrix obtained for the validation genomes ```results_all/results_20221004T143059/results_alleles.tsv```.  The following command was used:

## Command:

```bash
# create header
head -n 1 results_cg/cgMLST_15/cgMLST.tsv > cgMLST_all.tsv
```

## Command:

```bash
# concatenate
grep -v ^FILE results_cg/cgMLST_15/cgMLST.tsv results_all/ results_20221004T143059/results_alleles.tsv >> cgMLST_all.tsv
```
The cgMLST_all.tsv file can be found in the folder: ```analysis_all/cgMLST_all.tsv```. This file (cgMLST_all.tsv) contains the allelic profile of the 264 genomes.

## Step 3.2: Evaluation of genome quality

After concatenation, we used the *TestGenomeQuality* to assess the impact of each validation genome on candidate *loci* in order to exclude low quality validation genomes. In this step you may need to define a new *Threshold*, as well as a new value of the parameter *p*, because *loci* that remain after the filters are the ones that will constituted the final scheme.

## Command:

```bash
 chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5
```
The folder with the output file can be found at: ```analysis_all/GenomeQualityPlot.html```. The list with the *loci* present in 95% of the genomes at the chosen *Threshold* can be retrieved in the folder ```analysis_all/Genes_95%.txt```. The list of genomes removed at each *Threshold* can be retrieved in the folder ```analysis_all/removedGenomes.txt```.

In order to exclude validation genomes that have left the scheme it is necessary to follow the steps described in **Step 2.2**

## Step 4: Extracting the Matrix loci

At this step we chose *loci* present in 99% (*t0.99*) of the validation genomes and the *Threshold* 15 to limit the loss of the *loci* in the genomes. In *Threshold* 15 a set of 1795 *loci* were found to be present in 99% the validation genomes. 

From the original "removedGenomes.txt" file that can be retrieved in the ```analysis_all/removedGenomes.txt``` folder we created another (removedGenomes15thr.txt) file with only the genomes removed at *Threshold* (15). The list of genomes removed at *Threshold 15* can be retrieved in the folder: ```analysis_all/removedGenomes15thr.txt``` 

Using *Threshold* (15) only 1 draft genomes were removed due to absence of *loci* targets.

## Command:

```bash
chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_15  --t 0.99 --g removedGenomes15thr.txt 
```

This script selects *loci* and genomes that remained in the *Threshold* 15 and excludes the validation genomes and *loci* that were excluded with this *Threshold*.

The folder with the output file can be found at: ```cgMLST_15```. This folder contains four files "cgMLST.tsv; cgMLSTschema.txt; mdata_stats.tsv and Presence_Absence.tsv".

The cgMLST targets can be found at: ```cgMLST_15/cgMLSTschema.txt``` It contains the list of 1795 genes in the core genome defined as targets for this scheme.
## Step 5: Minimum Spanning Tree

For the visualization of results, minimum spanning trees were buitl. Based on the allelic profiles obtained by the cgMLST scheme for each of the 263 genomes minimum spanning trees (MST) were constructed using the software GrapeTree (version 1.5.0) (https://github.com/achtman-lab/GrapeTree/releases) with parameters implemented in MSTree v2 ignoring missing values for the entire strain collection. The ```cgMLST_15/cgMLST.tsv ``` file contains the allelic profile of the 263 genomes typed by cgMLST.

## Step 6: Graphical evaluation of the scheme

To assess the variability of the gene targets of cgMLST as well explore and evaluate the type and extent of allelic variation detected in each of the chosen *loci*. We run this script and graphically visualize the data in a series of html files.

## Command:

```bash
chewBBACA.py SchemaEvaluator -i schema_seed/schema_seed -o rms/RmS.html  --cpu 46
```
## Step 7: How to identify the allelic profile of your *C. striatum* genomes of interest using this cgMLST scheme?

In the steps below we describe all the necessary steps to identify the allelic profile of the genomes of interest using this cgMLST scheme.

## Step 7.1: A master directory to run ChewBBACA

A master directory ```Analyze_genomes/``` containing the folders that are needed to run ChewBBACA with the 1795 cgMLST target genes was created as an example. This folder contains a directory called ```example_genomes/``` representing the folder of the genomes to be typed. The second folder is the ```example_schema_seed/``` representing the ```schema_seed/``` folder that must be downloaded. “gene_targets.txt” files containing the 1795 target genes of cgMLST and the file " GCF_004138065.trn" which is the prodigal's trained file to recognize CDs.

## Step 7.2: Install ChewBBACA and all its dependencies;

To download ChewBBACA and all its dependencies see the topic: Softwares and Downloads (Main dependencies); Other dependencies (for scheme evaluation only):

## Step 7.3: Download Scheme

To have access to the cgMLST scheme for * C. striatum * it is necessary to download the ```schema_seed/``` folder. 

**Note:** The ```schema_seed/``` folder was created in Step 1 where we identified all CDs of the 271 genomes creating the scheme.

## Step 7.4: Download the list of target genes

The list of cgMLST target genes can be obtained from the ```Analyze_genomes/gene_targets.txt```folder.

## Step 7.5: Download the Prodigal trained file

It is recommended to upload the trained folder (GCF_004138065.1.trn) because it is the reference folder for Prodigal to recognize the coding sequences (CDs). Prodigal was trained with the reference genome of *C. striatum* (WP1a) to recognize CDs.

## Step 7.6: Genomes of interest

The genomes of interest must be in a specific folder, for example the ```genomes/``` folder and must be in the same directory that contains the ```schema_seed/``` folder.

After this part, the next step is to run the command:

## Command: 

```bash
chewBBACA.py AlleleCall -i genomes -g gene_targets.txt -o results --cpu 46 --ptf GCF_004138065.1.trn
chewBBACA.py AlleleCall -i ./genomes --gl list_genes_core.txt -o results --cpu 46  -g ./schema_seed/schema_seed –ptf GCF_004138065.1.trn
```
**Note**:The folder “genomes" represents the folder with the genomes to be typed.

**Note**: The list "gene_targets.txt" containing the 1795 cgMLST target genes.
This command will release the output in the folder: ```results/```

## Step 8: Analyze the results

The allelic profile of the typed genomes will be in the folder: ```results/results_alleles.tsv``` which is the output of the file released by the script above. Other outputs will be in the folder ```results/```, as an example: RepeatedLoci.txt; logging_info.txt; results_contigsInfo.tsv and results_statistics.tsv.

**Note**: With the ```results/``` folder, analyze the file ```results/results_statistics.tsv.```It contains the number of genes that were found with 100% identity in the analyzed genomes (EXC - alleles which have exact matches (100% DNA identity) with previously identified alleles). **With the proposed scheme we consider that an isolate is considered typed when at least 95% of the targets are found in its genome**.

## Step 9: Allele profile view by minimum spanning tree (MST)

To view the allelic profile data of the typed genomes you need to access the output of the script present in the folder ```results/results_alleles.tsv```. It is possible to do so in two ways: the first is to download the free software GrapeTree (version1.5.0) with parameters implemented in MSTree v2 ignoring missing values for the entire strain collection available in (https://github.com/achtman-lab/GrapeTree/releases) or through the software Phyloviz, online.
