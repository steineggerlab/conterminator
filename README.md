# Conterminator
[![Build Status](https://dev.azure.com/themartinsteinegger/conterminator/_apis/build/status/martin-steinegger.conterminator?branchName=master)](https://dev.azure.com/themartinsteinegger/conterminator/_build/latest?definitionId=2&branchName=master)
## Detection of contamination in nucleotide and protein sequence sets
Conterminator is an efficient method for detecting incorrectly labeled sequences across kingdoms by an exhaustive all-against-all sequence comparison.
It is a free open-source GPLv3-licensed software for Linux and macOS, and is developed on top of modules provided by [MMseqs2](https://github.com/soedinglab/MMseqs2).


<p align="center"><img src="https://raw.githubusercontent.com/martin-steinegger/conterminator/master/.github/marv6.png" height="256" /></p>

[Terminating contamination: large-scale search identifies more than 2,000,000 contaminated entries in GenBank. Genome Biology, doi: 10.1186/s13059-020-02023-1 (2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02023-1)


# Install 
Conterminator requires a 64-bit Linux system (check that `uname -a | grep x86_64 | wc -l` is greater than `0`) with at least a SSE4.1 instruction set (check that `cat /proc/cpuinfo | grep sse4_1 | wc -l` yields an output greater than `0`).
   
    # SSE4.1
    wget https://mmseqs.com/conterminator/conterminator-linux-sse41.tar.gz; tar xvfz conterminator-linux-sse41.tar.gz; export PATH=$(pwd)/conterminator/:$PATH
    # AVX2
    wget https://mmseqs.com/conterminator/conterminator-linux-avx2.tar.gz; tar xvfz conterminator-linux-avx2.tar.gz; export PATH=$(pwd)/conterminator/:$PATH
    # conda
    conda install -c bioconda conterminator
    
# Getting started

Conterminator computes ungapped local alignments of all provided sequences and reports contamination across user-specifed specified taxa; by default this is done at the kingdom level.

Conterminator requires two input files: (1) a FASTA file containing all sequences (`example/dna.fna`/`example/prots.faa`) and (2) a `mappingFile` (`example/dna.mapping `/`examples/prots.mapping`), which maps FASTA identfiers to NCBI taxon identfiers. The program produces two output files with prefix (`${RESULT_PREFIX}`). More details are described in the Results section below.

To process nucleotide sequences use the following command:
    
    conterminator dna example/dna.fna example/dna.mapping ${RESULT_PREFIX} tmp     
    
Protein sequences can be processed as following:        

    conterminator protein example/prots.faa example/prots.mapping ${RESULT_PREFIX} tmp  
    
## Mapping file 

Conterminator needs a mapping file, which assigns each fasta identifier to a taxonomical identifier. The `mapping` file consists of two tab-delimited columns, (1) fasta identifier and (2) [NCBI taxonomy identifier] (taxonomy ID) (https://www.ncbi.nlm.nih.gov/taxonomy). 
By default, Conterminator takes the text up to the first blank space as the fasta identifier. However, with GenBank, Tremble, Swissprot, Conterminator extracts out only the unique identifier mapped to the taxonomy ID.

Example for detecting contamination in the NT database:

    blastdbcmd -db nt -entry all > nt.fna
    blastdbcmd -db nt -entry all -outfmt "%a %T" > nt.fna.taxidmapping
    conterminator dna nt.fna nt.fna.taxidmapping nt.result tmp
    
## Result

Conterminator produces two result files (1) `${RESULT_PREFIX}_conterm_prediction` and (2) `${RESULT_PREFIX}_all`.
The `{RESULT_PREFIX}_conterm_prediction` text file contains the predicted contamination. The file is TSV-seperated containing the following columns: 

```
1.) Numeric identifier
2.) Contaminated identifier
3.) Kingdom (default: 0: Bacteria&Archaea, 1: Fungi, 2: Metazoa, 3: Viridiplantae, 4: Other Eukaryotes)
4.) Species name
5.) Alignment start
6.) Alignment end
7.) Corrected contig length (length between flanking Ns)
8.) Identifier of the longest contaminating sequence
9.) Kingdom of the longest contaminating sequence
10.) Species name of the longest contaminating sequence
11.) Length of the longest contaminating sequence
12.) Count how often sequences from the contaminating kingdom align
```

Be aware that the result file may contain any contaminated identifier multiple times if multiple alignments were detected.  

The `{RESULT_PREFIX}_all` contains the information for all alignments that were used to predict contamination. 

```
1.) Numeric identifier
2.) Sequence identifier
3.) Alignment start
4.) Alignment end
5.) Corrected contig length (length between flanking Ns)
6.) Total sequence length
7.) Kingdom (default: 0: Bacteria&Archaea, 1: Fungi, 2: Metazoa, 3: Viridiplantae, 4: Other Eukaryotes)
8.) Species name 
```

## Important Parameters
### `--kingdom`

This parameter specifies which taxons across which contaminations should be considered.
Each taxon definition is seperated by a `,` e.g. to search for contamination between bacteria and human use `--kingdom 2,9606`. 
It is also possible to use more advanced expressions for contamination rules, through the following operators:

    ! NEGATION 
    || OR  
    && AND 

The default rule is as follows:

    2||2157,4751,33208,33090,2759&&!4751&&!33208&&!33090   
    
This searches for contamination between the following taxa:

    2||2157  # Bacteria OR Archaea 
    4751     # Fungi
    33208    # Metazoa
    33090    # Viridiplantae  
    2759&&!4751&&!33208&&!33090 # Eukaryota without Fungi Metazoa and Viridiplantae

# [OPTIONAL] Install by Compilation 
Users can install Conterminator using the commands specified above. However, Conterminator can also be installed by compiling directly from the source code using the following commands. 

    git clone --recursive https://github.com/martin-steinegger/conterminator 
    mkdir conterminator/build && cd conterminator/build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j 4
    make install
    export PATH=$(pwd)/bin/:$PATH 

