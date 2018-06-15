# BDTT

This repository contains R codes that can be used to run BDTT on a community ecology dataset. BDTT decomposes community dissimilarities by time-slicing the phylogenetic tree of all members (e.g. 16S rDNA reads) present in all communities. It allows the user to detect at which time or phylogenetic scale an environmental factor of interest shapes the most community compositions. It has recently been applied to analyse the composition of gut bacterial communities across mammals, to disentangle the effects of host phylogeny and host diet.

In this repository, we provide two files:

BDTT_functions.R, which contains functions that are needed to run BDTT on NON ULTRAMETRIC trees
BDTT_functions_UltrametricTree.R, which contains functions that are needed to run BDTT on ULTRAMETRIC trees
BDTT_Example_scripts_UltrametricTree.R, which contains descriptions of how BDTT works, along with illustrative examples on NON ULTRAMETRIC trees
If you have any questions, please contact either Mathieu Groussin (mgroussi@gmail.com) or myself(flo.mazel@gmail.com)
