
#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
ProteoformMap_4 = function(Peptide_FC_P_Selected,Peptide_ProteinSeq_2){
  data_peptide<-Peptide_FC_P_Selected
  #data_protein_seq<-read.xlsx(address_2,rowNames = F,colNames = T)
  data_protein_seq<-Peptide_ProteinSeq_2
  for(n in 1: nrow(data_peptide)){# begin 0
    
    #n=1
    data_peptide[n,"Protein.Id"]
    data_protein_seq_test<-subset(data_protein_seq,UNIPROTKB==data_peptide[n,"Protein.Id"])
    if(!is.na(data_protein_seq_test[1,"SEQUENCE"])&&nchar(data_protein_seq_test[1,"SEQUENCE"])>10){# begin1
      a=data_protein_seq_test[1,"SEQUENCE"]
      data_peptide[n,"PeptideLocation"]=indexOf(a,data_peptide[n,"PeptideSeq"])
      data_peptide[n,"PeptideLen"]=nchar(data_peptide[n,"PeptideSeq"])
      data_peptide[n,"ProteinLen"]=nchar(data_protein_seq_test[1,"SEQUENCE"])
      data_peptide[n,"RTPosition"]=indexOf(a,data_peptide[n,"PeptideSeq"])/nchar(data_protein_seq_test[1,"SEQUENCE"])
      #print(data_peptide[n,1])
      #print(n)
      
    }# final1
  }# final 0
  #write.csv(data_peptide,file=address_4,row.names=FALSE)
  
  # ########analysis the RTPosition and the Fold use the data_peptide from above section
  for(n in 1: nrow(data_protein_seq)){# begin 0 use the data_protein_seq as the list name This is the for cycle for the C,N,M adugement
    
    
    test2<-subset(data_peptide,Protein.Id==data_protein_seq[n,"UNIPROTKB"])
    # 1 C case
    if((test2[1,"Fold"]<0)&&(test2[nrow(test2),"Fold"]>0)){
      test2_small<-subset(test2,Fold==-2)
      test2_large<-subset(test2,Fold==1)
      if(min(test2_large$RTPosition)>max(test2_small$RTPosition)){
        for(n in 1: nrow(data_peptide)){
          if(data_peptide[n,"Protein.Id"]==test2[1,1]){
            data_peptide[n,"Proteoform"]="C"
          }# all change as C
          
        }# This is for to change the all peptide select as C
        
        
        
      }# this is the test small and test large compare
      
    }# this is the for the list in proteinseq list all name
    #2 N case
    if((test2[1,"Fold"]>0)&&(test2[nrow(test2),"Fold"]<0)){
      test2_small<-subset(test2,Fold==-2)
      test2_large<-subset(test2,Fold==1)
      if(min(test2_small$RTPosition)>max(test2_large$RTPosition)){
        for(n in 1: nrow(data_peptide)){# all change as C
          if(data_peptide[n,"Protein.Id"]==test2[1,1]){
            data_peptide[n,"Proteoform"]="N"
          }# all change as C
          
        }# This is for to change the all peptide select as C
        
        
        
      }# this is the test small and test large compare
      
    }# this is the for the list in proteinseq list all name
    
    
    ##############
    
    # test use talbe for M case
    
    # 3 M case   method 3
    
    if((test2[1,"Fold"]<0)&&(test2[nrow(test2),"Fold"]<0)){
      proteoM="M"
      # first order the test2
      #TEST4<-test2[order(test2[,"RTPosition"]),test2[,"F"]),]
      test2_M_large=subset(test2,Fold>0)
      test2_M_small=subset(TEST4,Fold<0)
      test2_M_smallRTP<-min(test2_M_large$RTPosition)
      test2_M_largeRTP<-max(test2_M_large$RTPosition)
      for(n in 1: nrow(test2_M_small)){
        if((test2_M_small[n,"RTPosition"]>test2_M_smallRTP)&&(test2_M_small[n,"RTPosition"]<test2_M_largeRTP)){
          proteoM=0
          break
        }
        #print(proteoM)
        #print(n)
      }
      for(n in 1: nrow(data_peptide)){# all change as C
        if(data_peptide[n,"Protein.Id"]==test2[1,1]){
          data_peptide[n,"Proteoform"]="proteoM"
        }# all change as C
      }
      
      
      ############# 3 M case   method 3
      
    }# this is for for } not use for test
    
    
    test2_small<-subset(test2,Fold==-2)
    test2_large<-subset(test2,Fold==1)
    if(min(test2_small$RTPosition)>max(test2_large$RTPosition)){
      for(n in 1: nrow(data_peptide)){
        if(data_peptide[n,"Protein.Id"]==test2[1,1]){
          data_peptide[n,"Proteoform"]="N"
        }# all change as C
        
      }# This is for to change the all peptide select as C
      
      
      
    }# this is the test small and test large compare
    
  }# this is the for the list in proteinseq list all name This is the for cycle for the C,N,M adugement
  #above is judge N and C 
  ##################### 
  
  # output the preotorm not NA from data_peptide above
  
  Peptide_proteoform=data_peptide[complete.cases(data_peptide[,"Proteoform"]),]
  Peptide_proteoform_list=Peptide_proteoform[!duplicated(Peptide_proteoform$Protein.Id),]
  return(Peptide_proteoform)
}# the part 4 end