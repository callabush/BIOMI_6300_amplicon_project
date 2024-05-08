Analysis of the gut microbiome of P. sexstriatus.

Includes fastqc reports and a multiqc report. 

The following are summaries of .Rmd and .html files with the same names. 

01_DADA2 follows the DADA2 pipeline (https://benjjneb.github.io/dada2/ and Callahan 2016, Nature). 

02_PreProcessing takes objects generated through the DADA2 pipeline and creates an object for phyloseq (https://joey711.github.io/phyloseq/  McMurdie and Holmes 2013, PLOS ONE).

03a_Phylogenetic_Tree creates a phylogenetic tree with objects froms the DADA2 pipeline. This tree is used in downstream analysis for UniFrac calculations (Lozupone, ISMEJ). Alignment is done with MAFFT (https://mafft.cbrc.jp/alignment/software/ Katoh 2013, Mol. Biol. Evol.) and tree is generated with FastTree2 (http://www.microbesonline.org/fasttree/ Price 2010, PLOS ONE).

03b_Phylogenetic_Tree roots tree in ggtree (Gu 2017, Methods in Ecology and Evolution). 

04_Biodiversity measures alpha diversity measures, such as rarefraction using iNEXT (Hsieh and Chao 2024 R package, Chao et al 2014 in Ecological Monographs_

05_Community_Analysis looks at beta diveristy, comparing gut section IV to gut section V. Includes Sorenson and Bray Curtis dissimilarity metrics. As well as UniFrac. Graphs these using PCoA and NMDS. Statistical analysis is done using the adonis2 script in the vegan R package (https://search.r-project.org/CRAN/refmans/vegan/html/adonis.html). Statistical analysis done is Permanova and BetaDispR. 
