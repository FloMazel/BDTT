# Beta-diversity through time (BDTT) 

This repository contains R codes that can be used to run BDTT on a community ecology dataset. BDTT decomposes community dissimilarities by time-slicing the phylogenetic tree of all members (e.g. 16S rDNA reads) present in all communities. It allows the user to detect at which time or phylogenetic scale an environmental factor of interest shapes the most community compositions. It has recently been applied to analyse the composition of mammalian communites around the world and the composition of gut bacterial communities within mammalian guts

The main function, "BDTT" has three main arguments: 

-- *similarity_slices*: a numerical vector giving the h heights at which the tree should be sliced

-- *tree*: a phylogeny in the "phylo" ape format, with tip name matching the names in the OTU matrix (the rownames of the sampleOTUs argument)

-- *sampleOTUs*: a matrix of n OTUS (row) by m samples (columns)

and returns an array of beta diversity values (h*b*m*m): m is the number of samples, h is the number of slices, and b is the number of beta-diversity metrics computed (3 so far: BrayCurtis, Jaccard and the true turnover component of Jaccard).   


Groussin, M., Mazel, F., Sanders, J. G., Smillie, C. S., Lavergne, S., Thuiller, W., & Alm, E. J. (2017). Unraveling the processes shaping mammalian gut microbiomes over evolutionary time. Nature communications, 8, 14319.

Mazel, F., Wüest, R. O., Lessard, J. P., Renaud, J., Ficetola, G. F., Lavergne, S., & Thuiller, W. (2017). Global patterns of β‐diversity along the phylogenetic time‐scale: The role of climate and plate tectonics. Global Ecology and Biogeography, 26(10), 1211-1221.

It has been written by myself with inputs from M. Groussin and J. Roy. In this repository, we provide a file that contains [functions](https://github.com/FloMazel/BDTT/blob/master/BDTT_functions.R) needed to run BDTT  and [example and tutorials](https://github.com/FloMazel/BDTT/blob/master/Tutorial_Examples.pdf). If you have any questions, please contact either Mathieu Groussin (mgroussi@gmail.com) or myself (flo.mazel@gmail.com)




