# Conterminator
## Detection of contamination in nucleotide and protein sequence sets
Conterminator is an efficient method to detect incorrectly labeled sequences across kingdoms by an exhaustive all-against-all sequence comparison.
It is free open-source GPLv3-licensed software for Linux and macOS, and is developed on top of modules provided by [MMseqs2](https://github.com/soedinglab/MMseqs2).

<p align="center"><img src="https://raw.githubusercontent.com/martin-steinegger/conterminator/master/.github/marv6.png" height="256" /></p>

 
# Installation
Conterminator can be installed by compiling from source. It requires a 64-bit system (check with uname -a | grep x86_64) with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on macOS).

    git clone --recursive https://github.com/martin-steinegger/conterminator 
    mkdir conterminator/build && cd conterminator/build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j 4
    make install
    export PATH=$(pwd)/bin/:$PATH 
    
# Getting started

Conterminator computes ungapped local alignments of all sequence and reports contamination across user-specifed specified taxa, by default this is done at kingdom level.

To process nucleotide sequences use the following command:
    
    conterminator dna sequence.fasta result tmp --tax-mapping-file mappingFile

Conterminator requires a `--tax-mapping-file` file, which maps FASTA identfiers to NCBI taxon identfiers.
    
Protein sequences can be processed as following:        

    conterminator protein sequence.fasta result tmp --tax-mapping-file mappingFile 

## Important Parameters
### `--taxon-list`

This parameters controls across which ranks contaminations should be considered. 
Each taxon definition is seperated by a `,` e.g. to search for contamination between bacteria and human use `--taxon-list 2,9606`. 
It is also possible to use more advanced expressions for contamination rules, through the following operators:

    ! NEGATION 
    | OR  
    & AND 

The default rule is as follows:

    (2|2157),4751,33208,33090,(2759&!4751&!33208&!33090)   
    
This searches for contamination between the following taxa:

    (2|2157) # Bacteria OR Archaea 
    4751     # Fungi
    33208    # Metazoa
    33090    # Viridiplantae  
    (2759&!4751&!33208&!33090) # Eukaryota without Fungi Metazoa and Viridiplantae

# Publication

Poster: [Terminating contamination—Large-scale search identifies >2,000,000 contaminated entries in GenBank](https://f1000research.com/posters/8-1861)