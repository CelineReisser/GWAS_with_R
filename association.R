
#################################################
# R code used for carrying association analysis #
#################################################
## July 2016, Celine Reisser
## Published in Reisser et al. 2016. Molecular Biology and Evolution.



## association.R
## Note: some subroutines of this script are copied or adapted from Darren Obbard's "PHASE.R", also used in Joron et al. 2011, Nature.
## To execute, first copy and paste this page into R, then type association() to call the function. Select your input fasta file (formatted as described below).
## It might take a few minutes to run, depending of the size of your input file. 



## R Libraries to install prior to use this code: "sequinr"



##########################    INPUT INFORMATION    #########################
##
## Input format: fasta alignment
## IUPAC standard code: "ATCGYRKMSW"
## Missing data coded as "N"
##
## Example of input for three individuals at three SNP loci:
##   >ind1
##   ATY
##   >ind2
##   AAY
##   >ind3
##   TAN
##
##  Do number your SNPs in your reference VCF file (from 1 to length(input)), in order to re-identify them in the output ("site" column).
##  
##  Organise your individuals by phenotype: all individuals from phenotype 1, then inds from phenotype 2...
##  In the third paragraph of the code, change the value for "nNMP<-" and "nMP<-" according to your dataset (those name correspond to my phenotypes, NMP and MP, you can
##  modify the names, but be careful to do so throughout the rest of the code).



##########################    OUTPUT INFORMATION    #########################
## 
## Once finished the message "analysis completed, results written to file YOUR_INPUT_NAME.fasta_chisquare.out".
## It should contain the following information: site (locus ID number you have added to your reference VCF, see above); di-/multi (if your locus was diallelic or multiallelic);
## hyp (the allele corresponding to your hypothesis, here the most common allele); Chi2 (value of Chi2 statistics); P (corresponding P-value).
##
## site	di-/multiall	hyp	Chi2			P	
## 1	diallelic	g	0.00698597000483793	0.99651310837294	
## 2	diallelic	g	0.091387599469495	0.955334441069135	
## 3	diallelic	t	0.091387599469495	0.955334441069135	
## 4	diallelic	a	2.97826886015311	0.225567815781874	
## 5	diallelic	c	3.39279321355224	0.183342992073151	
## 6	diallelic	t	0.816850208892629	0.664696253421223	




########################  CODE TO CREATE THE "association" FUNCTION   ########################  

require("seqinr")

association<-function(infile.name){
  infile.name<-file.choose()
  indata<-read.fasta(infile.name)
  outfile.name<-paste(infile.name,"_chisquare.out",sep="")
  
  #check the format
  for (i in 1:length(indata)){
    if(is.SeqFastadna(indata[[i]])==FALSE)print("Incorrect format")
  }
  seqLength<-getLength(indata[[1]])  #saves doing it repeatedly
  nInds<-length(indata)
  nNMP<-XX # MODIFY according to your own dataset
  nMP<-XX # MODIFY according to your own dataset
  
  
  #  read the sequences into an array
  indata_matrix<-matrix(data = NA, nrow = (nInds), ncol = seqLength, byrow = FALSE,dimnames = NULL)
  for (i in 1:nInds) {
    indata_matrix[i,]<-getSequence(indata[[i]],as.string=FALSE)
  }
  
  # check which sites are variable 
  j<-1
  polymorphic_sites<-0  # note that this may screw up is the file has no variable sites at all, but with SNP data this shouldn't happen
  for (i in 1:seqLength) {
    if(sum(is.element(c("a","g","c","t","y","r","k","m","s","w"),indata_matrix[,i]))>1){
      polymorphic_sites[j]<-i
      j<-j+1
    }
  }
  polys<-length(polymorphic_sites)
  
  #make matrix with polymorphic sites only
  poly_matrix<-matrix(data = NA, nrow = nInds, ncol = polys, byrow = FALSE, dimnames = NULL)
  for (i in 1:nInds){
    poly_matrix[i,]<-indata_matrix[i,polymorphic_sites]
  }
  
  
  #make an allele matrix by arbitrarily splitting heterozygotes (arbitrary phasing, obtaining genotypes in two rows)
  #note: the IUPAC code is used here, so genotypes are encoded using W,S,K...
  allele_matrix<-matrix(data = NA, nrow = (2*nInds), ncol = polys, byrow = FALSE,dimnames = NULL)
  for (i in 1:nInds){
    i1<-((2*i)-1)
    i2<-(2*i)
    allele_matrix[i1,]<-poly_matrix[i,]
    allele_matrix[i2,]<-poly_matrix[i,]
    
    for (j in 1:polys) {
      if(allele_matrix[i1,j]=="y"){
        allele_matrix[i1,j]<-"c"
        allele_matrix[i2,j]<-"t"
        next  
      }
      if(allele_matrix[i1,j]=="r"){
        allele_matrix[i1,j]<-"a"
        allele_matrix[i2,j]<-"g"
        next
      }
      if(allele_matrix[i1,j]=="m"){
        allele_matrix[i1,j]<-"c"
        allele_matrix[i2,j]<-"a"
        next	
      }
      if(allele_matrix[i1,j]=="s"){
        allele_matrix[i1,j]<-"c"
        allele_matrix[i2,j]<-"g"
        next
      }			
      if(allele_matrix[i1,j]=="k"){
        allele_matrix[i1,j]<-"g"
        allele_matrix[i2,j]<-"t"
        next
      }
      if(allele_matrix[i1,j]=="w"){
        allele_matrix[i1,j]<-"a"
        allele_matrix[i2,j]<-"t"
        next
      }
    }	
  }
  
  #determine multiallelic sites, and only keep bi-allelic sites (if SNP call done by Stacks, only bi-allelic sites should be present)
  multi<-0
  for (j in 1:polys){
    if(sum(is.element(c("a","g","c","t"),allele_matrix[,j]))>2){
      multi[j]<-"multiallelic:_check_hypothesis!"
    }else{
      multi[j]<-"diallelic"
    }
  }
  
  #determine allele frequencies in both phenotypes
  allfreq.nmp<-matrix(data = NA, nrow = 4, ncol = polys, byrow = FALSE,dimnames = NULL)
  allfreq.mp<-matrix(data = NA, nrow = 4, ncol = polys, byrow = FALSE,dimnames = NULL)
  allfreq.all<-matrix(data = NA, nrow = 4, ncol = polys, byrow = FALSE,dimnames = NULL)
  for (j in 1:polys){ #iterate through all loci
    a<-0; c<-0; g<-0; t<-0
    for (i in 1:(2*nNMP)){# iterate through all individuals in NMP (here NMP were the first 17 individuals)
      if(allele_matrix[i,j]=="a"){
        a<-a+1
        next
      }
      if(allele_matrix[i,j]=="c"){
        c<-c+1
        next
      }
      if(allele_matrix[i,j]=="g"){
        g<-g+1
        next
      }
      if(allele_matrix[i,j]=="t"){
        t<-t+1
        next
      }
    }
    allfreq.nmp[1,j]<-a
    allfreq.nmp[2,j]<-c
    allfreq.nmp[3,j]<-g
    allfreq.nmp[4,j]<-t
    
    a<-0; c<-0; g<-0; t<-0
    for (i in (2*nNMP+1):(2*nInds)){ # iterate through all individuals in MP (ordered after the NMP, hence 2*NMP+1, since genotypes are in two rows)
      if(allele_matrix[i,j]=="a"){
        a<-a+1
        next
      }
      if(allele_matrix[i,j]=="c"){
        c<-c+1
        next
      }
      if(allele_matrix[i,j]=="g"){
        g<-g+1
        next
      }
      if(allele_matrix[i,j]=="t"){
        t<-t+1
        next
      }
    }
    allfreq.mp[1,j]<-a
    allfreq.mp[2,j]<-c
    allfreq.mp[3,j]<-g
    allfreq.mp[4,j]<-t
    allfreq.all[1,j]<-allfreq.nmp[1,j]+allfreq.mp[1,j]
    allfreq.all[2,j]<-allfreq.nmp[2,j]+allfreq.mp[2,j]
    allfreq.all[3,j]<-allfreq.nmp[3,j]+allfreq.mp[3,j]
    allfreq.all[4,j]<-allfreq.nmp[4,j]+allfreq.mp[4,j]
  }
  
  
  #determine the number of valid NMP and MP genotypes at each site
  nnNMP<-0
  nnMP<-0
  for (j in 1:polys){
    nnNMP[j]<-(sum(allfreq.nmp[,j]))/2		
    nnMP[j]<-(sum(allfreq.mp[,j]))/2		
  }
  
  #determine hypothesis to be tested: the most common allele is found at the homozygous state in MP but at the heterozygous state in NMP
  hypothesis<-0
  first<-0
  second<-0
  allfreq.first.all<-0
  allfreq.first.nmp<-0
  allfreq.first.mp<-0
  allfreq.sec.all<-0
  allfreq.sec.nmp<-0
  allfreq.sec.mp<-0
  alleles<-c("a","c","g","t")
  for (j in 1:polys){
    helpv.all<-allfreq.all[,j]/(2*(nnNMP[j]+nnMP[j]))  #definition of the most common allele in the whole pop
    helpv.nmp<-allfreq.nmp[,j]/(2*nnNMP[j])
    helpv.mp<-allfreq.mp[,j]/(2*nnMP[j])
    helpv2<-rank(helpv.all, na.last=FALSE, ties.method= "first") #ordering allele freq in the whole pop from smallest to largest
    helpv3<-c(0,0,0,0)
    for (i in 1:4){
      if(helpv2[i]==4){
        helpv3[1]<-i    #if the rank is 4, it is the most common allele. place the alleles
        next
      }
      if(helpv2[i]==3){
        helpv3[2]<-i
        next
      }
      if(helpv2[i]==2){
        helpv3[3]<-i
        next
      }
      if(helpv2[i]==1){
        helpv3[4]<-i
        next
      }
    }
    first[j]<-alleles[helpv3[1]]						#the most common allele in the population
    allfreq.first.all[j]<-helpv.all[helpv3[1]]  			#allele frequency of the most common allele
    second[j]<-alleles[helpv3[2]]						#the second most common allele
    allfreq.sec.all[j]<-helpv.all[helpv3[2]]				#allele frequency of the second most common allele
    hyp<-helpv3[1]
    hypothesis[j]<-alleles[hyp]						#the most common allele in population
  }
  
  
  
  #Chi-square test
  P<-0; n1<-0; n2<-0; n3<-0; n4<-0; exp1<-0; exp2<-0; exp3<-0; exp4<-0; nbNMP<-0; nbMP<-0
  Chi1<-0; Chi2<-0; Chi1.n1<-0; Chi1.n2<-0; Chi1.n3<-0; Chi1.n4<-0;
  for (j in 1:polys){
    NMP.A1A1<-0;NMP.A1A2<-0;NMP.A2A2<-0;MP.A1A1<-0;MP.A1A2<-0;MP.A2A2<-0
    for (i in 1:nNMP){
      if(poly_matrix[i,j]==hypothesis[j]){
        NMP.A1A1<-NMP.A1A1+1
        next
      }else{	
        if(is.element(poly_matrix[i,j],c("a","g","c","t"))){
          NMP.A2A2<-NMP.A2A2+1 	
          next
        }else{
          if(is.element(poly_matrix[i,j],c("y","r","k","m","s","w"))){
            NMP.A1A2<-NMP.A1A2+1 		  
          }
        }
      }	
    }		
    for (i in (nNMP+1):nInds){
      if(poly_matrix[i,j]==hypothesis[j]){
        MP.A1A1<-MP.A1A1+1 
        next
      }else{	
        if(is.element(poly_matrix[i,j],c("a","g","c","t"))){
          MP.A2A2<-MP.A2A2+1 	
          next
        }else{
          if(is.element(poly_matrix[i,j],c("y","r","k","m","s","w"))){
            MP.A1A2<-MP.A1A2+1 		  
          }
        }
      }	
    }		
    n1[j]<-NMP.A1A2
    n2[j]<-NMP.A1A1+NMP.A2A2
    n3[j]<-MP.A1A1
    n4[j]<-MP.A1A2+MP.A2A2
    nbNMP[j]<-n1[j]+n2[j]
    nbMP[j]<-n3[j]+n4[j]
    exp1[j]<-nbNMP[j]*2*allfreq.first.all[j]*(1-allfreq.first.all[j])
    exp2[j]<-nbNMP[j]-exp1[j]
    exp3[j]<-nbMP[j]*allfreq.first.all[j]*allfreq.first.all[j]
    exp4[j]<-nbMP[j]-exp3[j]
    
    Chi1.n1[j]<-((n1[j]-exp1[j])^2)/exp1[j]
    Chi1.n2[j]<-((n2[j]-exp2[j])^2)/exp2[j]
    Chi1.n3[j]<-((n3[j]-exp3[j])^2)/exp3[j]
    Chi1.n4[j]<-((n4[j]-exp4[j])^2)/exp4[j]
    
    Chi1[j]<-Chi1.n1[j]+Chi1.n2[j]+Chi1.n3[j]+Chi1.n4[j]
    
    if((n1[j]-exp1[j])>0){
      if((n3[j]-exp3[j])>0){
        Chi2[j]<-Chi1[j]
        next
      }else{
        Chi2[j]<-Chi1.n1[j]+Chi1.n2[j]
      }
      next
    }else{
      if((n3[j]-exp3[j])>0) {
        Chi2[j]<-Chi1.n3[j]+Chi1.n4[j]
        next
      }else{
        Chi2[j]<-0
      }
    }
  }
  for (j in 1:polys){
    P[j]<-pchisq(Chi2[j],2,lower.tail=FALSE)
  }
  
  
  #write the outfile
  outmatrix<-matrix(data = NA, nrow = polys+1, ncol = 5, byrow = FALSE,dimnames = NULL)
  firstrow<-c("site","di-/multiall","hyp","Chi2","P")
  for (j in 1:5){
    outmatrix[1,j]<-firstrow[j]
  }
  for (i in 2:(polys+1)){
    outmatrix[i,1]<-polymorphic_sites[i-1]
    outmatrix[i,2]<-multi[i-1]
    outmatrix[i,3]<-hypothesis[i-1]
    outmatrix[i,4]<-Chi2[i-1]
    outmatrix[i,5]<-P[i-1]
  }
  
  cat(outmatrix[1,],"\n",file=outfile.name,sep="\t",append=FALSE)
  for (i in 2:(polys+1)){
    cat(outmatrix[i,],"\n",file=outfile.name,sep="\t",append=TRUE)
  }
  cat("analysis completed, results written to file ", paste(outfile.name),"\n",file ="")
}
