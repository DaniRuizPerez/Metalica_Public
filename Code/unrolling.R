

rm(list=ls())
library(scales)
library(stringr)
library(dplyr)
options("scipen"=100, "digits"=4)
library(sqldf)

source("../AuxFunctions.R")


getUnrolledTGMTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    genes =  unique(as.character(tinteraction[grep("g__.*",tinteraction)]))
    for (g in genes){
      # for each gene, get the metabolites it produces
      ginteraction = unique(as.character(pairs[grep(g,pairs[,1]),2]))
      metabolites =  unique(as.character(ginteraction[grep("m__.*",ginteraction)]))
      bootscoreTG = pairs[intersect(grep(t,pairs[,1]), grep(g,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (m in metabolites){
        minteraction = unique(as.character(pairs[grep(m,pairs[,1]),2]))
        taxainfluenced = unique(as.character(minteraction[grep("s__.*",minteraction)]))
        bootscoreGM = pairs[intersect(grep(g,pairs[,1]), grep(m,pairs[,2]))[1],3]
        for (ti in taxainfluenced){
          # print(paste(t,"->",g,"->",m,"->",ti))
          # write(paste(".*",t,".*",ti), outputFile,  append=TRUE)
          # write(paste(t,"->",g,"->",m,"->",ti), outputFile,  append=TRUE)
          
          #Get the bootscore of each
          bootscoreMTi = pairs[intersect(grep(m,pairs[,1]), grep(ti,pairs[,2]))[1],3]
          unrolledpairs = rbind(unrolledpairs,c(t,"",g,m,ti,bootscoreTG,bootscoreGM,bootscoreMTi))
          
        }
      }
    }
  }
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","","g","m","ti","bTG","bGM","bMTi")
  return(unrolledpairs)
}

getUnrolledTGTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    genes =  unique(as.character(tinteraction[grep("g__.*",tinteraction)]))
    for (g in genes){
      # for each gene, get the metabolites it produces
      ginteraction = unique(as.character(pairs[grep(g,pairs[,1]),2]))
      taxainfluenced =  unique(as.character(ginteraction[grep("s__.*",ginteraction)]))
      bootscoreTG = pairs[intersect(grep(t,pairs[,1]), grep(g,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (ti in taxainfluenced){
        # write(paste(t,"->",g,"->",ti), outputFile,  append=TRUE)
        # unrolledpairs = rbind(unrolledpairs,c(t,g,"",ti))
        
        #Get the bootscore of each
        bootscoreGTi = pairs[intersect(grep(g,pairs[,1]), grep(ti,pairs[,2]))[1],3]
        unrolledpairs = rbind(unrolledpairs,c(t,"",g,"",ti,bootscoreTG,"",bootscoreGTi))
        
      }
    }
  }
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","","g","m","ti","bTG","","bGTi")
  return(unrolledpairs)
}

getUnrolledTMTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    metabolites =  unique(as.character(tinteraction[grep("m__.*",tinteraction)]))
    for (m in metabolites){
      # for each gene, get the metabolites it produces
      minteraction = unique(as.character(pairs[grep(m,pairs[,1]),2]))
      taxainfluenced =  unique(as.character(minteraction[grep("s__.*",minteraction)]))
      bootscoreTM = pairs[intersect(grep(t,pairs[,1]), grep(m,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (ti in taxainfluenced){
        # write(paste(t,"->",m,"->",ti), outputFile,  append=TRUE)
        # unrolledpairs = rbind(unrolledpairs,c(t,"",m,ti))
        
        #Get the bootscore of each
        bootscoreMTi = pairs[intersect(grep(m,pairs[,1]), grep(ti,pairs[,2]))[1],3]
        unrolledpairs = rbind(unrolledpairs,c(t,"","",m,ti,"",bootscoreTM,bootscoreMTi))
      }
    }
  }
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","","g","m","ti","","bTM","bMTi")
  return(unrolledpairs)
}

getUnrolledTTTpairs <- function(pairs){
  unrolledpairs = c()
  #get the taxa ti+1
  taxa = unique(as.character(pairs[grep("s__.*",pairs[,1]),1]))
  for (t in taxa){
    #for each taxa, get the genes it expresses
    tinteraction = unique(as.character(pairs[grep(t,pairs[,1]),2]))
    middleTaxa =  unique(as.character(tinteraction[grep("s__.*",tinteraction)]))
    for (tx in middleTaxa){
      # for each gene, get the metabolites it produces
      middleSinteraction = unique(as.character(pairs[grep(tx,pairs[,1]),2]))
      finalTaxaInfluenced =  unique(as.character(middleSinteraction[grep("s__.*",middleSinteraction)]))
      bootscoreTti = pairs[intersect(grep(t,pairs[,1]), grep(tx,pairs[,2]))[1],3]
      #now try to see if any of these metabolites itneract with the taxa
      for (ti in finalTaxaInfluenced){
        #Get the bootscore of each
        bootscoreTTiTj = pairs[intersect(grep(tx,pairs[,1]), grep(ti,pairs[,2]))[1],3]
        unrolledpairs = rbind(unrolledpairs,c(t,tx,"","",ti,bootscoreTti,"",bootscoreTTiTj))
      }
    }
  }
  # unrolledpairs = rbind(unrolledpairs,c(t,g,m,ti,bootscoreTG,bootscoreGM,bootscoreMTi))
  
  if (length(unrolledpairs) > 0)
    colnames(unrolledpairs) = c("t","tx","g","m","ti","","bTM","bMTi")
  return(unrolledpairs)
}

intersectUnrolledAndOriginal <- function(unrolledpairs,pairsTT){
  #now try to do the itnereseciton between pairsTT and unrolledpairs
  intersected = c()
  
  if (NROW(pairsTT) == 0 || NROW(unrolledpairs) == 0){
    return(NULL)
  }
  
  for (i in 1:NROW(pairsTT)){
      indexOfverifiedinteraction = intersect(grep(pairsTT[i,1],unrolledpairs[,1]), grep(pairsTT[i,2],unrolledpairs[,5]))[1]
      if (!is.na(indexOfverifiedinteraction)){
        # print(unrolledpairs[indexOfverifiedinteraction,])
        intersected = rbind(intersected,c(unrolledpairs[indexOfverifiedinteraction,],pairsTT[i,3]))
      }
    # }
  }
  return(intersected)
}

getfinalTGTM <- function(intersectedAll,TGTchains,TMTchains){
  finalTGTM = c()
  for (i in 1:NROW(intersectedAll)){
    tgm = intersectedAll[i,]
    #Check if this interaction is also in TGT
    isInTGT = intersect(which(TGTchains[,1] %in% tgm[1]),intersect(which(TGTchains[,5] %in% tgm[5]),which(TGTchains[,3] %in% tgm[3])))
    #Check if this interaction is also in TMT
    isInTMT = intersect(which(TMTchains[,1] %in% tgm[1]),intersect(which(TMTchains[,5] %in% tgm[5]),which(TMTchains[,4] %in% tgm[4])))
    finalTGTM = rbind(finalTGTM,c(tgm,length(isInTGT)!=0,length(isInTMT)!=0))
  }
  

  
  #remove self loops
  selfLoopsToRemove = c()
  for (i in 1:NROW(finalTGTM)){
    if (finalTGTM[i,1] == finalTGTM[i,2] && finalTGTM[i,1] == finalTGTM[i,5]){
      selfLoopsToRemove = c(selfLoopsToRemove,i)
    }
  }
  if (length(selfLoopsToRemove) > 0){
    finalTGTM = finalTGTM[-selfLoopsToRemove,]
  }
  
  
  
  # Don't remove duplicates, since they are triangles, so independent interactions
  ## finalTGTM= removeDuplicates(finalTGTM)
  finalTGTM = unique(finalTGTM)
  
  if (NROW(finalTGTM) == 0){
    return (NULL)
  }
  
  #Add the overalls core as the multiplication of all bootscores
  aux = finalTGTM[,6:9]
  aux[aux==""] = 1
  overallScore = as.numeric(aux[,1])*as.numeric(aux[,2])*as.numeric(aux[,3])*as.numeric(aux[,4])
  if (length(overallScore) == 0 || is.null(overallScore) || is.na(overallScore)){
    overallScore = 0
  }
  finalTGTM = cbind(finalTGTM[,1:9],overallScore, finalTGTM[,10:13])

  colnames(finalTGTM) = c("T","Tx","G","M","Ti", "T->(G|Tx)_Bs","(T|G)->M_Bs","(M|G|Tx)->Ti_Bs" ,"T->Ti_Bs","overallScore","ChainInTTxTi","ChainInTGMTi","ChainInTGTi","ChainInTMTi")
  
  finalTGTM = finalTGTM[order(finalTGTM[,12],finalTGTM[,13],finalTGTM[,14], finalTGTM[,11],finalTGTM[,10],decreasing=T),]
  # finalTGTM = finalTGTM[duplicated(finalTGTM), ]
  return(finalTGTM)
}

removeDuplicates <- function(d){
  remove = c()
  for (i in 1:(NROW(d)-1)){
    for (j in (i+1):NROW(d)){
      if (all(d[i,-c(6,7,8,9)] == d[j,-c(6,7,8,9)])){
        remove = c(remove,i)
        break
      }
    }
  }
  if (length(remove) != 0){
    d = d[-remove,]
  }
  return(d)
}

getNumberOfInteractionTypes <- function(pairs){
  
  et = length(intersect(grep("(week |Disease |Age_).*",pairs[,1]),grep("s__.*",pairs[,2])))
  tg = length(intersect(grep("s__.*",pairs[,1]),grep("g__.*",pairs[,2])))
  gm = length(intersect(grep("g__.*",pairs[,1]),grep("m__.*",pairs[,2])))
  mt = length(intersect(grep("m__.*",pairs[,1]),grep("s__.*",pairs[,2])))
  tt = length(intersect(grep("s__.*",pairs[,1]),grep("s__.*",pairs[,2])))
  mm = length(intersect(grep("m__.*",pairs[,1]),grep("m__.*",pairs[,2])))
  gg = length(intersect(grep("g__.*",pairs[,1]),grep("g__.*",pairs[,2])))
  gt = length(intersect(grep("g__.*",pairs[,1]),grep("s__.*",pairs[,2])))
  
  stats = c(et,gt,tg,gm,mt,tt,mm,gg)
  names(stats) = c("#e->t","#g->t","#t->g","#g->m","#m->t","#t->t","#m->m","#g->g")
  
  return(stats)
  # return(paste("Number of t->g: ", tg, ", g->m: ", gm, ", m->t: ", mt, ", t->t: " , tt, ", m->m: ", mm, ", g->g: ", gg,sep=""))
}

validateChains <- function(folder,TGMName,TGName,TMName,TName,TTName,outputFile, method){
  
  pairsTGMT = readGenericNetwork(paste(folder,TGMName,sep="/"),method)
  pairsTT = readGenericNetwork(paste(folder,TName,sep="/"),method)
  pairsTGT = readGenericNetwork(paste(folder,TGName,sep="/"),method)
  pairsTMT = readGenericNetwork(paste(folder,TMName,sep="/"),method)
  pairsTTT = readGenericNetwork(paste(folder,TTName,sep="/"),method)
  
  interactionTypespairsTT = getNumberOfInteractionTypes(pairsTT)
  interactionTypespairsTGMT = getNumberOfInteractionTypes(pairsTGMT)
  interactionTypespairsTGT = getNumberOfInteractionTypes(pairsTGT)
  interactionTypespairsTMT = getNumberOfInteractionTypes(pairsTMT)
  interactionTypespairspairsTTT = getNumberOfInteractionTypes(pairsTTT)
  print("pairsTT")
  print(interactionTypespairsTT)
  print("pairsTGMT")
  print(interactionTypespairsTGMT)
  print("pairsTGT")
  print(interactionTypespairsTGT)
  print("pairsTMT")
  print(interactionTypespairsTMT)
  print("pairsTTT")
  print(interactionTypespairspairsTTT)
  
  TGMTchains = getUnrolledTGMTpairs(pairsTGMT)
  TGTchains = getUnrolledTGTpairs(pairsTGT)
  TMTchains = getUnrolledTMTpairs(pairsTMT)
  TTTchains = getUnrolledTTTpairs(pairsTTT)
  
  intersectedTGTM = intersectUnrolledAndOriginal(TGMTchains,pairsTT)
  intersectedTGT = intersectUnrolledAndOriginal(TGTchains,pairsTT)
  intersectedTMT = intersectUnrolledAndOriginal(TMTchains,pairsTT)
  intersectedTTT = intersectUnrolledAndOriginal(TTTchains,pairsTT)

  if (is.null(intersectedTGTM) && is.null(intersectedTGT) && is.null(intersectedTMT) &&is.null(intersectedTTT)){
    finalAll = NULL
  }else{
 
    intersectedAll= c()
    if (!is.null(intersectedTGTM)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTGTM,rep("FALSE",NROW(intersectedTGTM)),rep("TRUE",NROW(intersectedTGTM))))
    }
    if (!is.null(intersectedTGT)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTGT,rep("FALSE",NROW(intersectedTGT)),rep("FALSE",NROW(intersectedTGT))))
    }
    if (!is.null(intersectedTMT)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTMT,rep("FALSE",NROW(intersectedTMT)),rep("FALSE", NROW(intersectedTMT))))
    }
    if (!is.null(intersectedTTT)){
      intersectedAll = rbind(intersectedAll,cbind(intersectedTTT,rep("TRUE",NROW(intersectedTTT)),rep("FALSE", NROW(intersectedTTT))))
    }
      finalAll = getfinalTGTM(intersectedAll,TGTchains,TMTchains)
  }
  if (!is.null(finalAll)){
    write.table(finalAll, outputFile, sep=",", append=F, row.names = F,quote=F, col.names = T)
  }
  fractionUnrolledTGMT =if (is.null(intersectedTGTM)) 0 else (NROW(intersectedTGTM))/NROW(pairsTT)
  fractionUnrolledTGT = if (is.null(intersectedTGT)) 0 else (NROW(intersectedTGT))/NROW(pairsTT)
  fractionUnrolledTMT = if (is.null(intersectedTMT)) 0 else (NROW(intersectedTMT))/NROW(pairsTT)
  fractionUnrolledTTT = if (is.null(intersectedTTT)) 0 else (NROW(intersectedTTT))/NROW(pairsTTT)
  
  if (!is.null(finalAll)){
    overallScore = sum(as.numeric(finalAll[,10]))/nrow(finalAll)
  }else{
    overallScore = 0
  }
  
  stats = c( if (is.null(TTTchains)) 0 else NROW(TTTchains),
              if (is.null(TGTchains)) 0 else NROW(TGTchains),
              if (is.null(TMTchains)) 0 else NROW(TMTchains),
              if (is.null(TGMTchains)) 0 else NROW(TGMTchains),
              if (is.null(intersectedTTT)) 0 else NROW(intersectedTTT),
              if (is.null(intersectedTGTM)) 0 else NROW(intersectedTGTM),
              if (is.null(intersectedTGT)) 0 else NROW(intersectedTGT),
              if (is.null(intersectedTMT)) 0 else NROW(intersectedTMT),
              fractionUnrolledTTT,fractionUnrolledTGMT,fractionUnrolledTGT,fractionUnrolledTMT,overallScore)

  #TGTchains: Number of chains t1->g1->t2 in the TG network.
  #intersectedTGT: How many from  #TGTchains are actually an unrolling. This means that the interaction t1->g1->t2 appears in the form t1->t2 in the TT network. This is expanded into the individual interactions list that I generate in another file 
  #%UnrolledTGT: simply dividing  intersectedTGT /#TGTchains
    
  firstStatsNames = c("#TTTchains","#TGTchains","#TMTchains","#TGMTchains","intersectedTTT","intersectedTGTM","intersectedTGT","intersectedTMT","%UnrolledTTT","%UnrolledTGMT","%UnrolledTGT","%UnrolledTMT","overallScore")
  #Augment the stasts with the number of each interaction type for each network
  stats = c(interactionTypespairsTT,interactionTypespairsTGMT,interactionTypespairsTGT,interactionTypespairsTMT,interactionTypespairspairsTTT,stats)
  names(stats) = c(paste("TT:",names(interactionTypespairsTT)),paste("TGMT:",names(interactionTypespairsTGMT)),paste("TGT:",names(interactionTypespairsTGT)),
                   paste("TMT:",names(interactionTypespairsTMT)),paste("TTT:",names(interactionTypespairspairsTTT)),firstStatsNames)
  
  
  stats = c(TGMName,stats)
  names(stats)[1]  = "Name"
  
  return(list(stats,finalAll))
}
  
main <- function(){
  folder = "../networks/" 
  resultsFolder = "../unrolling/" 
  
  # matrix = "HEGTM_Skeleton.*"
  statsAll = c()
  unrollingsAll = c()
  srList = c("_sr7d")
  dataset = "filteredAllSubjects_ibd_"
  methods = c("PyCausal_","Tigramite_","DBN_")
  methods = c("DBN_")
  alignments = c("noalignment","alignment")
  # alignments = c("alignment")
  t = "t_"
  g = "g_"
  m = "m_"
  for (method in methods){
      if (method == "Tigramite_"){
      bs = "_bs0.1"
      testList = c("_test_ParCorr","_test_GPDC","_test_CMIknn")
      # testList = c("_test_CMIknn")
      aList = c("_a0.1","_a0.01","_a0.001","_a0.0001")
      scoreList = c("")
      nparentsList = c("")
      matrix = ".*"
      p = "_ptreshold0.05"
      nboots = ""
      tau = "_T1"
      sr = "_sr7d"
    }else if (method == "PyCausal_") {
      bs = ""
      scoreList = c("_scoreDiscreteMixedBicScore","_scoreConditionalGaussianBicScore","_scoreDegenerateGaussianBicScore","_scoreFisherZScore","_scoreMNLRBicScore","_scoreMVPBicScore","_scorePeterScore","_scoreSemBicScore","_scoreSemBicScoreDeterministic")
      testList = c("_testMultinomialLogisticRegressionWald","_testConditionalGaussianLRT", "DegenerateGaussianLRT","_testFisherZ","_testKci","_testMNLRLRT","_testPositiveCorr","_testSemBicTest")
      nparentsList = c("")
      
      scoreList = c("_scoreFisherZScore")
      testList = c("_testPositiveCorr")
      aList = c("_a0.1","_a0.01","_a0.001","_a0.0001")
      p = ""
      sr = "_sr7d"
      nboots = "_nboots10"
      matrix = ".*"
      tau = ""
    }else if (method == "DBN_") {
      # DBN
      testList = c("")
      aList = c("")
      scoreList = c("")
      bs = ""
      p = ""
      nparentsList = c("_nParents3","_nParents4","_nParents5","_nParents6")
      tau = ""
      test =""
      sr = "_sr7d"
      matrix =".*"
      nboots = "_nboots100"
    }else{
      # DBN
      nparents = c("")
      testList = c("")
      aList = c("")
      scoreList = c("")
      bs = ""
      p = ""
      tau = ""  
      nparentsList = c("")
      test =""
      sr = "_sr14d"
      matrix ="_dbnIntraBoot"
      nboots = ""
    }
    for (alignment in alignments){
      for (a in aList){
        for (test in testList){
          for (score in scoreList){
            for (sr in srList){
              for (nparents in nparentsList){
                TGMName = list.files(path = folder, pattern = paste(method,dataset,t,g,m,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                TGName = list.files(path = folder, pattern = paste(method,dataset,t,g,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                TMName = list.files(path = folder, pattern = paste(method,dataset,t,m,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                TName = list.files(path = folder, pattern = paste(method,dataset,t,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                TTName = list.files(path = folder, pattern = paste(method,dataset,t,t,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                
                
                # Tigramite_filtered_ibd_t_g_m_noalignment_sr7dTrim_bs0.1_a0.5_ptreshold0.05_T1_test_ParCorrHEGTM_Skeleton.matrix
                # print("Tigramite_filtered_ibd_t_g_m_noalignment_sr7d_bs0.1_a0.5_ptreshold0.05_T1_test_ParCorrHEGTM_Skeleton.matrix")
                # print(paste(method,dataset,t,g,m,alignment,sr,bs,a,nboots,p,tau,score,test,matrix,sep=""))
                print(paste(method,dataset,t,t,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                
                if (length(TGMName) == 0){
                  print("Aborting, TGMName is:")
                  print(TGMName)
                  next
                }
                if (length(TGName) == 0){
                  print("Aborting, TGName is:")
                  print(TGName)
                  next
                }
                if (length(TMName) == 0){
                  print("Aborting, TMName is:")
                  print(TMName)
                  next
                }
                if (length(TName) == 0){
                  print("Aborting, TName is:")
                  print(TName)
                  next
                }
                if (length(TGMName) > 1){
                  TGMName = TGMName[1]
                }
                if (length(TGName) > 1){
                  TGName = TGName[1]
                }
                if (length(TMName) > 1){
                  TMName = TMName[1]
                }
                if (length(TName) > 1){
                  TName = TName[1]
                }
                print(TGMName)
                print(TGName)
                print(TMName)
                print(TName)
                
                outputFile = paste(resultsFolder,method,dataset,t,g,m,alignment,sr,bs,a,nboots,nparents,p,tau,test,".csv",sep="")
                
                validationOutput = validateChains(folder = folder,TGMName = TGMName, TGName = TGName, TMName = TMName, TName= TName,TTName = TTName, method = method, outputFile = outputFile)

                stats = validationOutput[[1]]
                unrollings = validationOutput[[2]]
                unrollings = cbind(rep(paste(method,dataset,t,g,m,alignment,sr,bs,a,nboots,nparents,p,tau,test,sep=""),NROW(unrollings)),unrollings)
                # print(stats)
                if (nrow(unrollings) > 0){
                  unrollingsAll = rbind(unrollingsAll,unrollings)
                }
                  
                statsAll = rbind(statsAll,stats)
              }
            }
          }
        }
      }
    }
}
  #Output overall stats
  print(statsAll)
  write.table(statsAll, paste(resultsFolder,"/statsAllAugmented.csv",sep=""), sep=",", append=F, row.names = F,quote=F, col.names = T)
  
  #Output overall stats for each method
  overallStatsForEachMethod = c()
  statsAll = as.data.frame(statsAll)
  for (method in methods){
    for (alignment in alignments){
      foundIn = intersect(grep(method, statsAll$Name),grep(paste("_",alignment,sep=""), statsAll$Name))
      statsAsNumeric = mutate_all(statsAll[foundIn,2:(ncol(statsAll))], function(x) as.numeric(as.character(x)))
      overallStatsThisMethod = colSums(statsAsNumeric)/length(foundIn)
      overallStatsForEachMethod= rbind(overallStatsForEachMethod,c(method,alignment,overallStatsThisMethod))
    }
  }
  write.table(overallStatsForEachMethod, paste(resultsFolder,"/overallAVGStatsForEachMethod.csv",sep=""), sep=",", append=F, row.names = F,quote=F, col.names = T)
  

  
  #Compute most comon unrollings througout all methods, and all the unrollings
  colnames(unrollingsAll)[1] = c("NetworkName")
  unrollingsAll = unrollingsAll[order(unrollingsAll[,13],unrollingsAll[,14],unrollingsAll[,15], unrollingsAll[,11],unrollingsAll[,10],decreasing=T),]
  unrollingsAll = as.data.frame(unrollingsAll)
  countUnrollings = sqldf("select T, Tx, G, M,Ti, COUNT(*) AS NetworksThisUnrollingWasFoundIn  from unrollingsAll  GROUP BY T, G, M  ORDER BY NetworksThisUnrollingWasFoundIn desc") #HAVING T NOT LIKE '' AND G NOT LIKE '' AND M NOT LIKE '' AND Ti NOT LIKE ''
  combinedUnrollingsAndCount = sqldf("select a.*, NetworksThisUnrollingWasFoundIn  from unrollingsAll a join countUnrollings b on a.T = b.T and a.G=b.G and a.M=b.M and a.Ti=b.Ti and a.Tx=b.Tx" )
  duplicated = duplicated(combinedUnrollingsAndCount[,2:(NCOL(combinedUnrollingsAndCount))]) 
  combinedUnrollingsAndCountNoDuplicates = combinedUnrollingsAndCount[!duplicated,]
  write.table(combinedUnrollingsAndCountNoDuplicates, paste(resultsFolder,"/combinedUnrollingsAndCount.csv",sep=""), sep=",", append=F, row.names = F,quote=F, col.names = T)
  
}


main()

#See if the unrolled interactions get validated more frequently







