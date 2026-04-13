This repository contains code and data for:

**James ME, Melo MC, Roda F, Bernal-Franco D, Wilkinson MJ, Walter GM, Liu H, Engelstädter J, and Ortiz-Barrientos D. (2025) Ecological and mutation-order speciation in *Senecio*. Molecular Ecology [volume.page] [doi]**

## Study Overview

We demonstrate that ecological divergence and mutation-order processes can occur simultaneously as part of a continuum of selective pressures driving speciation in the *Senecio lautus* species complex.

## DNA Extraction and Sequencing

See [CTAB_protocol](laboratory/CTAB_protocol.docx) for the CTAB DNA extraction protocol used to extract DNA from leaf samples.

We undertook targeted re-sequencing of nuclear genomic regions (primarily intronic and intergenic regions) in the Access Array system. See [AccessArray_PrimerBinding](laboratory/AccessArray_PrimerBinding.docx) for details of the four-primer PCR process. Samples were sequenced on two lanes of the Roche GS FLX Titanium platform.

## Bioinformatics and Neutrality Tests

We used the standalone version of TagCleaner [tagcleaner.pl](scripts/tagcleaner.pl)
to trim forward and reverse barcodes. For each sample, the following code was run

```
tagcleaner.pl -fastq data.fastq -tag3 sequence1 -tag5 sequence2 -out trimmed -out_format 3 -64 -mm5 3 -mm3 5 -cont -info -log -trim_within 86 -split 3

```
, where `-fastq data.fastq` is the fastq file per individual, `-out trimmed` is the output filename prefix, `-out_format 3` sets output format to FASTQ, `-tag3 sequence1` defines the tag sequence expected at the 3'-end of reads (in this case the reverse complement of CS2), `-tag5 sequence2` defines the tag sequence expected at the 5'-end of reads (in this case the library specific key, barcode and CS1), `-64` tells the program to assume 64-bit, `-mm5 3` allows up to 3 mismatches when matching the 5'-end tag sequence, `-mm3 5` allows up to 5 mismatches when matching the 3'-end tag sequence, `-cont` enables continuous trimming of tag sequences, `-info` adds trimming information to the sequence header, `-log` creates a log file to track parameters and processing information, `-trim_within 86` searches for tag sequences within 86 bases from the sequence ends, and `-split 3` removes internal tag contaminations and splits fragment-to-fragment concatenations allowing up to 3 mismatches for internal tags.
A second round of tagcleaner was run
```
tagcleaner.pl -fastq data.fastq -tag3 sequence3 -tag5 sequence4 -out trimmed -out_format 3 -64 -mm5 3 -mm3 5 -cont -info -log -trim_within 86 -split 3
```
, this time removing the reverse complement of CS1 from the 3'-end of reads, and the 
in this case the library specific key, barcode and CS2 from the 5'-end of reads.

We then removed reads that had > 40% low quality bases (<Q20), > 2% Ns, or were < 50bp. 

We used the standalone version of ‘PRGmatic v1.6’ [PRGmatic.pl](scripts/PRGmatic.pl) to align reads into loci and construct haplotypes for each individual, using default parameters and an overlap of 100bp for read alignment.

The standalone version of ‘BLAT v35‘ [blat](scripts/blat) to map haplotypes from PRGmatic to each expected amplicon

```
./blat amplicons.txt fasta.trimmed -out=blast8 output.blast
```
, where ‘amplicons.txt’ is the list of expected amplicons, ‘fasta.trimmed’ is the query file containing the trimmed sequences, ‘-out=blast8’ specifies the output format as BLAST tabular format, and ‘output.blast’ is the output filename.

We then used ‘MUSCLE v3.8‘ [muscle](scripts/muscle) to align haplotypes for each amplicon
‘’’
./muscle -in locus1.aln.fna -out locus1.muscled.fasta
‘’’
, where ‘-in locus1.aln.fna’ specifies the input file containing sequences to be aligned, and ‘-out locus1.muscled.fasta’ is the output file for the multiple sequence alignment.

We only considered a locus for further if haplotypes were assigned to at least six individuals per population, excluding *S. madagascariensis* due to its smaller sample size. Following this criterion, we identified 26 loci. See [All_loci.fasta](data/All_loci.fasta) for a fasta file for all 26 loci (2 haplotypes per individual). The unaligned sequences (Ns and gaps removed) have been deposited on NCBI [Accession number to be added upon release]. 

We next assessed deviations from neutrality for each locus within each population using a combination of the following tests: Hudson-Kreitman-Aguad (HKA), Tajima’s D, Fu and Li’s D*, and Fu and Li’s F*. Analysis was performed in the GUI interface of ‘DNAsp v5‘, with default parameters. A locus was retained for phylogenetic and population structure analyses if it was considered neutral via the HKA test for all populations. In some cases, we were unable to calculate HKA due to the lack of sequencing data in the outgroup (*S. madagascariensis*). For these loci we assessed neutrality using a combination of the other three tests, and considered a locus neutral if, for all populations, it passed at least two tests. From this we obtained 13 neutral loci across all populations.

## Phylogenetic Analysis

To examine the phylogenetic independence of the *S. lautus* population pairs we conducted a Bayesian phylogenetic analysis using the 13 neutral loci in the GUI interface of ‘*BEAST v1.7.5‘ with default parameters, running the analysis with a chain length of 300,000,000 under a strict molecular clock with a Yules speciation process for species tree estimation. Input files were first created in the GUI interface of ‘Beauti v1.7.5‘, and ‘jModelTest v2.1.1‘ was used to determine the HKY model was the best fit for our loci. The GUI interface of ‘TreeAnnotator v1.7.5‘ was used to generate the maximum clade credibility tree with a burn-in of 10,000 steps, and was visualized in ‘FigTree v1.4.4‘ using *S. madagascariensis* as the outgroup. 

## Population Structure Analysis

To assess the major population structure in the data we identified the most likely number of genetic clusters across all populations using the GUI interface of ‘STRUCTURE v2.3.4’. We tested K values from 1 to 16, running 20 iterations per K with a burn-in of 100,000 and MCMC run length of 100,000. We confirmed convergence of model parameters through visual inspection of MCMC summary statistics.

## Intrinsic Reproductive Isolation 

We assessed intrinsic reproductive isolation by undertaking crosses both within and between ecotypes as well as within and between the two clades of the phylogeny, for two barriers: F1 seed set and F1 viability (germination). See [IntrinsicRIdata](data/IntrinsicRIdata) for the data file, and [IntrinsicAnalysis.R](scripts/IntrinsicAnalysis.R) for all analysis code.

## Extrinsic Reproductive Isolation

We assessed extrinsic reproductive isolation by reanalyzing data from a previous field transplant experiment by Walter et al. (2016). See [Transplant.csv](data/Transplant.csv) for the data file, and [ExtrinsicAnalysis.R](scripts/ExtrinsicAnalysis.R) for all analysis code.

## Environmental and Phenotypic Analyses

We asked whether there was a relationship between environmental distance and levels of reproductive isolation, as well as phenotypic distance and levels of reproductive isolation across the parapatric pairs. See [CrossesDH-ReciprocalsAveraged.txt](data/CrossesDH-ReciprocalsAveraged.txt) for the data file, and [PhenotypeEnvironment.R](scripts/PhenotypeEnvironment.R) for all analysis code.

## Theoretical Framework

To further understand how various forms of natural selection drive patterns of reproductive isolation in Senecio, we created a mathematical framework to explore how both ecological and mutation-order processes can jointly underlie patterns of reproductive isolation.

Our framework extends the classic Unckless-Orr model by introducing an environmental similarity parameter that captures how selection pressures differ between populations adapting to distinct environments. The framework also explores polygenic architectures, showing how genetic complexity amplifies the effects of environmental divergence on speciation probability. See [MathematicalFramework.R](scripts/MathematicalFramework.R)
