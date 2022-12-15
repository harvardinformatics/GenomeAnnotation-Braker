# BRAKER
[BRAKER](https://github.com/Gaius-Augustus/BRAKER) is a gene prediction tool that combines [GenemarkET](http://exon.gatech.edu/GeneMark/) and [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus), and that uses extrinsic evidence in the form homology information derived from protein or RNA-seq data to automatically train these tools along with directly inferred homology-based predictions of gene models. As with our implementation of other tools capable of using both protein and RNA-seq extrinsic evidence, we generate annotations that either only use protein evidence, or use both protein and RNA-seq evidence. Previously, Braker had options for simultaneously providing both types of evidence in a single analysis run, but this approach was classified as "experimental." Current best practice is to perform separate Braker analyses for each type of evidence, then merge the outputs with [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA). We have identified cases where TSEBRA erroneously discards large numbers of gene models because they purportedly are missing CDS features (when those features are present in the Braker gtf output files), and are waiting for fixes to the code before finalizing our annotations that combine RNA-seq and protein evidence.

## UTRs
When RNA-seq data are provided to BRAKER, it  has an option for predicting UTR portions of gene models. However, t is still viewed by its developers as "experimental". Furthermore, in initial testing, we discovered that, paradoxically, when running BRAKEr with UTR prediction, for some species it would predict far fewer CDS transcripts, leading to a non-trivial reduction in BUSCO recovery. Thus, in our evaluation, we do not implement UTR prediction.

## Singularity container
To aid users of Harvard's Cannon HPC cluster, and because successfully building Braker with all the required dependencies can be challenging, we run Braker using a Singularity image.

To build the Singularity image:

1. Clone this git repository
2. Register for GeneMark-ES/ET/EP (LINUX 64 kernel 3.10 - 5) at http://exon.gatech.edu/GeneMark/license_download.cgi 

    *NOTE: the license key expires after 400 days*
2. Download the genemark tarball (gmes_linux_64_4.tar.gz) and 64-bit license key (gm_key_64.gz) into the same directory as the same directory as the Singularity definition file (singularity.def) in this repository:

    curl -LO http://topaz.gatech.edu/GeneMark/tmp/GMtool_RANDOM_STRING/gmes_linux_64_4.tar..gz

    curl -LO http://topaz.gatech.edu/GeneMark/tmp/GMtool_RANDOM_STRING/gm_key_64.gz
3. Build the Singularity image:

   sudo singularity build Braker.sif singularity.def

## Pre-processing input files
In general, when Braker is provided fasta files that have more than one space-separated field in the header, it issues warning saying it may not work correctly in downstream steps. Therefore, we "clean" the headers of all genome and protein fasta files prior to use,simply by removing everything after the 1st field (which is usually the sequence name). In cases where a fasta file requires concatenation of more than one header field to provide the correct sequence name, custom tweaking of the below one-off awk command will be necessary:
```bash
awk '{print $1}' $my_input_fasta > headcleaned_${my_input_fasta}
```

## Protein-only
Running Braker with extrinsic protein evidence only requires you provide the genome sequence, and the protein (amino acid) fasta files. The developers of BRAKER indicate that it works best when multiple variants of individual proteins are provided, and recommend using [OrthoDB](https://www.orthodb.org/) proteins. Drop-down menus at this site all for easy extraction of taxonomically relevant protein fasta files. Guidelines for preparing the protein database (that BRAKER uses) from the protein fasta are found at the [ProtHint](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) github repository.

Once the database is prepared, assuming your Singularity container was called Braker.sif, one can run Braker as follows:
```
#!/bin/sh
# speciesname == name with no spaces
# genome fasta == full path to genome fasta
# orthodb == full path to protein database fasta

speciesname=$1
genomefasta=$2
orthodb=$3

singularity exec --cleanenv Braker.sif cp -Rs /opt/augustus/config augustus_config
singularity exec --no-home \
                 --home /opt/gm_key \
                 --cleanenv \
                 --env AUGUSTUS_CONFIG_PATH=${PWD}/augustus_config \
                 Braker.sif  braker.pl --cores 48 \
                 --species=${speciesname} --genome=${genomefasta} \
                 --prot_seq=${orthodb} \
                 --softmasking 
```

## RNA-seq-only
Braker provides a variety of options for how to integrate RNA-seq based evidence. The simplest way to do this is to supply an intron splice "hints" file in gff3 format as an input along with the genome fasta. However, to have gene predictions include UTR intervals, one must supply a bam file of RNA-seq spliced alignments to the genome (e.g. aligned with STAR or HISAT2). For our study, because we are comparing annotations to RNA-seq assemblers (Scallop and StringTie) that assemble UTRs and CDS, we include UTR predictions. Thus, we provide both an external hints gff (used for predicting features other than UTR) and a bam file representing the merge of all sample-level bam files used for the species of interest. To maximize sensitivity, our merged bam file includes sample-level bam files generated by both HISAT2 and STAR. The hints file required by Braker is an "intron hints" gff file generated with tools provided as part of [Augustus](https://github.com/Gaius-Augustus/Augustus). A detailed description of how we generate those hints can be found on this repositories [Comparative Augustus](https://github.com/harvardinformatics/GenomeAnnotation/tree/master/ComparativeAugustus) sub-page. 

Once this hints file is generated, on can simply execute Braker, supplying the genome fasta, the hints gff, and the merged bam file. Note, the merged bam file is the unfiltered bam file, not the filtered bam file produced during the course of generating the hints gff. Again, assuming the Singularity image is Braker.sif, BRAKER can be run as follows:
```
#!/bin/sh
# speciesname == name with no spaces
# genomefasta == full path to genome fasta
# bam == bam file of RNA-seq reads splice-aligned to genome with either Hisat2 or STAR

singularity exec --cleanenv Braker.sif cp -Rs /opt/augustus/config augustus_config
singularity exec --no-home \
                 --home /opt/gm_key \
                 --cleanenv \
                 --env AUGUSTUS_CONFIG_PATH=${PWD}/augustus_config\
                 Braker.sif braker.pl --cores 24 \
                 --bam=${bam} \
                 --species=${speciesname} --genome=${genomefasta} \
                 --softmasking 
```

### General cautionary notes
Particular features of braker outputs may lead to idisyncratic and undesireable behavior for downstream tools utilizing the braker annotation. A few that we have detected are as follows:
* Braker annotations come from AUGUSTUS and GeneMark. For single-exon gene, GeneMark does not appear to generate features of the type "exon", only "transcript" and "CDS" features. As a result, tool that extract nucleotide sequences based upon exon features (and their parent "transcript" features) will fail to report these. RSEM falls into this category.
