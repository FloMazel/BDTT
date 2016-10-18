
# Function to get Info on Edges
#---------------------------------
getTreeInfo=function(tr2)
{
  NH=node.depth.edgelength(tr2)
  treeH=max(NH)
  edgeInfo=cbind(tr2$edge,treeH-NH[tr2$edge[,1]],treeH-NH[tr2$edge[,2]])
  row.names(edgeInfo)=as.character(1:dim(edgeInfo)[1])
  return(edgeInfo)
}

# Function to get all branches and their descendant node at a given slice
#------------------------------------------------------------------------
GetBranchNode=function(slice,edgeInfo) {
  a=(edgeInfo[,3]>=slice)*(edgeInfo[,4]<=slice)
  return(edgeInfo[a==1,])
}

# Functions to fill the branch*sites matrices 
# (sites can be environmental sites or hosts for microbial communities)
#----------------------------------------------------------------------
vectorPres=function(node,tree)
{
  pres=rep(0,length(tree$tip.label))
  names(pres)=tree$tip.label
  pres[clade.members(node,tree,tip.label=T)]=1 
  return(pres)
}

vectorPresBigMat=function(nodes,tree,mat)
{
  for (i in nodes)
  {
    #print(100*match(i,nodes)/length(nodes))
    dd=clade.members(i,tree,tip.label=T)
    mat[as.character(i),dd]=1 
  }  
  
  return(mat)
}

# Function 'GetBranchOcc': computes a Branch*Sites matrix and directly save it in a file (not in Renv)
# Branch*sites matrices are also frequently named 'OTU tables' in the field of microbiology
#--------------------------------------------------------------------------------------------------------

#   slice: the age of the desired slice
#   tree: the community phylogenetic tree (must be ultrametric)
#   sitesp: the site * species matrice. Species names must match the tip names in the phylogenetic tree
#            In microbiology, sitesp is equivalent to the OTU table determining the distribution of unique 16S sequences across samples (100% similarity OTUs) 
#   pathtoSaveBranchOcc: directory where Branch*sites matrices are stored
#   bigmatrix=F: if site*species is a VERY big matrix, it is highly recommended to set bigmatrix=T (use of the bigmemory package)

GetBranchOcc=function(slice,tree,sitesp,pathtoSaveBranchOcc,bigmatrix=F)
{
  edgeInfo=getTreeInfo(tree)
  
  if (slice==0)
  {OTUmat=sitesp}
  else 
  {
    brCorres=GetBranchNode(slice=slice,edgeInfo=edgeInfo)
    NodeDes=brCorres[,2]    
    #print("getting OTUS host matrix")
    OTUnames=as.character(brCorres[,2])
    if (bigmatrix==F) {branchPAtotal=matrix(0,ncol=length(tree$tip.label),nrow=length(NodeDes),dimnames=list(OTUnames,tree$tip.label))}
    else if (bigmatrix==T) {branchPAtotal=filebacked.big.matrix(ncol=length(tree$tip.label),nrow=length(NodeDes),dimnames=list(OTUnames,tree$tip.label),backingfile = paste("branchPAtotal",slice,sep=""),backingpath =paste(pathtoSaveBranchOcc), descriptorfile = paste("branchPAtotalsauv",slice,sep=""))}   
      
      mat1=vectorPresBigMat(node=NodeDes,tree=tree,mat=branchPAtotal)
      mat1=mat1[,colnames(sitesp)]
      mat2=as.matrix(mat1)
      mat2=t(mat2)
      occM=as.matrix(sitesp)
      OTUmat=(occM)%*%(mat2)
  } 
      save(OTUmat,file=paste(pathtoSaveBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
}    

# Function to compute raw Jaccard and Sorensen Beta diversities 
#----------------------------------------------------------
getBeta=function(mat,ab=F)
{
  if (ab==T)
  {   
    h=bray.part(data.frame(mat))
    bctu=as.matrix(h[[1]])
    bcne=as.matrix(h[[2]])
    bc=as.matrix(h[[3]])
    res=abind(bctu,bcne,bc,along=0)
    dimnames(res)[[1]]=c("bctu","bcne","bc")
  }
  if (ab==F)
  {
    h=beta.pair(data.frame(mat), index.family="jaccard")
    hh=beta.pair(data.frame(mat), index.family="sorensen")
    jtu=as.matrix(h[[1]])
    jne=as.matrix(h[[2]])
    jac=as.matrix(h[[3]])
    stu=as.matrix(hh[[1]])
    sne=as.matrix(hh[[2]])
    sor=as.matrix(hh[[3]])  
    res=abind(jtu,jne,jac,stu,sne,sor,along=0)
    dimnames(res)[[1]]=c("jtu","jne","jac","stu","sne","sor")  
  }
  return(res)
}

# Function 'GetBetaDiv' to compute Beta-diversities from site*species matrices and to directly save the matrix of beta-diversities 
#---------------------------------------------------------------------------------------------------------------------------------

# INPUT VARIABLES
#   slice: the age of the desired slice
#   pathtoGetBranchOcc: the input file where the branch*sites matrix is saved (result of the function GetBranchOcc)
#   pathtoSaveBeta: the output file where the matrices of beta diversities are saved (several beta-diversity metrics are used)
#   bigmatrix=F: if you used bigmatrix=T to create the branch*sites matrix with 'GetBranchOcc' (OTU table), use bigmatrix=T again.

# OUTPUT
# save Beta diversity matrices as a 3D array in 'pathtoSaveBeta' 
# Array = array of sites*sites*beta diversity metrics 

#beta diversity metrics 

#"jtu" : True Turnover component of Jaccard (Presence/Absence)
#"jne" : Nestedness component of Jaccard (Presence/Absence)
#"jac" : Jaccard (Presence/Absence)
#"stu" : True Turnover component of Sorensen (Presence/Absence)
#"sne" : Nestedness component of Sorensen (Presence/Absence)
#"sor" : Sorensen (Presence/Absence)
#"bctu" : True Turnover component of Bray-Curtis (Abundance version of Sorensen)
#"bcne" : Nestedness component of Bray-Curtis (Abundance version of Sorensen)
#"bc" : Bray-Curtis (Abundance version of Sorensen)

GetBetaDiv=function(slice,pathtoGetBranchOcc,pathtoSaveBeta)
{
  load(paste(pathtoGetBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
  OTUmatPA=OTUmat
  OTUmatPA[OTUmatPA>0]=1
  Be1=getBeta(OTUmatPA,ab=F)  
  Be=getBeta(OTUmat,ab=T)
  Betaa=abind(Be1,Be,along=1)
  save(Betaa,file=paste(pathtoSaveBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
}   


# FUNCTION 'GetCorrelations': computes correlations between beta-diversities and environmental distances at the desired slice
#----------------------------------------------------------------------------------------------------------------------------

# INPUT VARIABLES
#   slice: the age of the desired slice
#   indice: the betadiv metric you want:
#		"jtu": True Turnover component of Jaccard (Presence/Absence)
#		"jne": Nestedness component of Jaccard (Presence/Absence)
#		"jac": Jaccard (Presence/Absence)
#		"stu": True Turnover component of Sorensen (Presence/Absence)
#		"sne": Nestedness component of Sorensen (Presence/Absence)
#		"sor": Sorensen (Presence/Absence)
#		"bctu": True Turnover component of Bray-Curtis (Abundance version of Sorensen)
#		"bcne": Nestedness component of Bray-Curtis (Abundance version of Sorensen)
#		"bc": Bray-Curtis (Abundance version of Sorensen)
#   pathtoGetBeta : the file where Beta-Diversity matrices were previously stored (output of the 'GetBetaDiv' function)
#   EnvDist: matrix of environmental distances. (e.g. geographic or climatic distances between sites, dietary or phylogenetic distances between hosts, etc).
#   TypeofMantel: type of correlation used to run the Mantel test (either "Spearman" or "Pearson")
#   nperm: number of permutations to perform to compute a p-value

# OUTPUT
# A 2*n matrix with R2 coefficients and their associated p-values


GetCorrelations=function(slice,indice="sor",pathtoGetBeta="",EnvDist,TypeofMantel="Spearman",nperm=1000)
{
  load(file=paste(pathtoGetBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
  #get same sites
  site=intersect(dimnames(Betaa)[[2]],colnames(EnvDist))
  Betaa=Betaa[indice,site,site]
  EnvDist=EnvDist[site,site]
  if (TypeofMantel=="Spearman") {  multiMant_SE=MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm,mrank=T)}
  if (TypeofMantel=="Pearson") {  multiMant_SE=MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm)}
  return(multiMant_SE$r.squared)
}

