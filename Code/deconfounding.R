#This file calculates and outputs the confounders and confounders of confounders for each file

rm(list=ls())
source("./AuxFunctions.R")

library(hash)
library(sqldf)
library(dplyr)


folder =  "../networks/" 
resultsFolder =  "../deconfounding/" 
setwd(resultsFolder)

#deconfound
# if we want to deconfound pairsTT with pairsTMT, build a hash with pairsTMT where the key is the child node and the value the parent, then search 
# for every pair in pairstTT for both nodes in the hash
findCounfounders = function(basePairs, extraOmicPairs,method, baseName, extraOmicName){
        
        if(is.null(basePairs) || is.null(extraOmicPairs)){
                return (c())
        }
        
        ExtraOmic_HashFromChildToParent <- hash() 
        
        # build lookup structure
        for (i in 1:nrow(extraOmicPairs)){
                pair = extraOmicPairs[i,]
                
                if (!is.null(ExtraOmic_HashFromChildToParent[[pair[2]]])){
                        ExtraOmic_HashFromChildToParent[[pair[2]]] = c(paste(pair[1],pair[3],sep="&"),ExtraOmic_HashFromChildToParent[[pair[2]]])
                }else{
                        ExtraOmic_HashFromChildToParent[[pair[2]]] = c(paste(pair[1],pair[3],sep="&"))
                }
                
        }
        
        confoundingTable = c()
        for (i in 1:nrow(basePairs)){
                pair = basePairs[i,]
                
                #We ignore self lops
                if (pair[[1]] == pair[[2]]){
                        next
                }
                
                #We get all parents of the confounded parent
                parents1InExtraOmicA = ExtraOmic_HashFromChildToParent[[pair[1]]] 
                #We get all parents of the confounded child
                parents2InExtraOmicA = ExtraOmic_HashFromChildToParent[[pair[2]]] 
                
                #We remove the bootscore -- for now
                parents1InExtraOmic = c()
                for (p in parents1InExtraOmicA){
                        parents1InExtraOmic = c(parents1InExtraOmic,strsplit(p,"&")[[1]][1])
                }
                parents2InExtraOmic = c()
                for (p in parents2InExtraOmicA){
                        parents2InExtraOmic = c(parents2InExtraOmic,strsplit(p,"&")[[1]][1])
                }
                
                #we make sure the "confounded" interaction doesn't exist in the new network
                if (pair[1] %in% parents2InExtraOmic) {
                        next
                }

                #Get the nodes that are parents for both confounded nodes in the new extra omic
                intersection = intersect(parents1InExtraOmic,parents2InExtraOmic)
                if(length(intersection)>0){
                        #we found a confounder
                        for(confounder in intersection){
                                # I need to filter out confounders that were already in the first newtork
                                if (!confounder %in% basePairs[,1] && !confounder %in% basePairs[,2]){
                                        #Get the bootsore for the deconfounding interacitons
                                        bootscorePC1 = strsplit(parents1InExtraOmicA[grep(confounder,parents1InExtraOmicA)],"&")[[1]][2]
                                        bootscorePC2 = strsplit(parents2InExtraOmicA[grep(confounder,parents2InExtraOmicA)],"&")[[1]][2]
                                        
                                        
                                        confoundingTable = rbind(confoundingTable,c(confounder,pair[1],pair[2],bootscorePC1,bootscorePC2,pair[3],as.numeric(bootscorePC1)*as.numeric(bootscorePC2)*as.numeric(pair[3]),method,baseName,extraOmicName))
                                }
                        }
                }
                
        }
        if (!is.null(confoundingTable)){
                colnames(confoundingTable) = c("confounder","confoundedParent","confoundedChild","BootscoreConfounderP","BootscoreConfounderC","BootscoreConfounded","overalScore","method","sourceLevel","deconfoundedWith")
                confoundingTable =unique(confoundingTable)
                confoundingTable = confoundingTable[order(confoundingTable[,7],decreasing=T),]
        }else{
                confoundingTable = c()
        }
        
        return(confoundingTable)
}

confounderListToPairs = function(confounderList) {
        if (is.null(confounderList)){
                return (c())
        }
        pairsList = c()

        if (is.null(nrow(confounderList))){
                confounderList = as.data.frame(t(confounderList))
        }
        for (i in 1:nrow(confounderList)){
                pairsList = rbind(pairsList,c(confounderList[i,1],confounderList[i,2]))
                pairsList = rbind(pairsList,c(confounderList[i,1],confounderList[i,3]))
        }       
        return(pairsList)
}

deconfound = function(){
        
        overallStats = data.frame()
        overallStatsMethod = data.frame()
        overallStatsMethodAcum = data.frame()
        allInteractions = c()

        
        # matrix = "HEGTM_Skeleton.*"
        srList = c("_sr7d")
        dataset = "filteredAllSubjects_ibd_"
        dataset = "filteredAllSubjects_ibd_"
        alignments = c("noalignment","alignment")
        # alignments = c("noalignment")
        methods = c("PyCausal_","Tigramite_","DBN_")
        # methods = c("DBN_")
        t = "t_"
        g = "g_"
        m = "m_"
                for (method in methods){
                        if (method == "Tigramite_"){
                        bs = "_bs0.1"
                        testList = c("_test_ParCorr","_test_GPDC","_test_CMIknn")
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
                        # nparentsList = c("_nParents4","_nParents5","_nParents6")
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
                                                                GMName = list.files(path = folder, pattern = paste(method,dataset,g,m,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                                TGName = list.files(path = folder, pattern = paste(method,dataset,t,g,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                                TTGName = list.files(path = folder, pattern = paste(method,dataset,t,t,g,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                                TMName = list.files(path = folder, pattern = paste(method,dataset,t,m,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                                TTMName = list.files(path = folder, pattern = paste(method,dataset,t,t,m,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                                TName = list.files(path = folder, pattern = paste(method,dataset,t,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                                TTName = list.files(path = folder, pattern = paste(method,dataset,t,t,alignment,sr,bs,a,nboots,nparents,p,tau,score,test,matrix,sep=""))
                                                               
                                                                

                                                                if (length(TGMName) > 1){
                                                                        TGMName = TGMName[1]
                                                                }
                                                                if (length(TGName) > 1){
                                                                        TGName = TGName[1]
                                                                }
                                                                if (length(TTGName) > 1){
                                                                        TTGName = TTGName[1]
                                                                }
                                                                if (length(TMName) > 1){
                                                                        TMName = TMName[1]
                                                                }
                                                                if (length(TTMName) > 1){
                                                                        TTMName = TTMName[1]
                                                                }
                                                                if (length(TName) > 1){
                                                                        TName = TName[1]
                                                                }
                                                                print(TGMName)
                                                                print(TGName)
                                                                print(TMName)
                                                                print(TName)
                                                                
                                                                outputFile = paste(method,dataset,t,g,m,alignment,sr,bs,a,nboots,nparents,p,tau,test,".csv",sep="")
                                                
                                        
                                                                pairsT = readGenericNetwork(paste(folder,TName,sep="/"),method)
                                                                pairsTT = readGenericNetwork(paste(folder,TTName,sep="/"),method)
                                                                pairsTG = readGenericNetwork(paste(folder,TGName,sep="/"),method)
                                                                pairsTTG = readGenericNetwork(paste(folder,TTGName,sep="/"),method)
                                                                pairsTM = readGenericNetwork(paste(folder,TMName,sep="/"),method)
                                                                pairsTTM = readGenericNetwork(paste(folder,TTMName,sep="/"),method)
                                                                pairsTGM = readGenericNetwork(paste(folder,TGMName,sep="/"),method)
                                                                pairsGM = readGenericNetwork(paste(folder,GMName,sep="/"),method)
                                                                
                                                                #Deconfound base networks
                                                                confoundersT = findCounfounders(pairsT, pairsTT,method,"T","TT")
                                                                confoundersM = findCounfounders(pairsT, pairsTM,method,"T","TM")
                                                                confoundersG = findCounfounders(pairsT, pairsTG,method,"T","TG")

                                                                #Deconfound multimic networks
                                                                confoundersGM = findCounfounders(pairsGM, pairsTGM,method,"GM","TGM")
                                                                confoundersTM = findCounfounders(pairsTM, pairsTTM,method,"TM","TTM")
                                                                confoundersTTM = findCounfounders(pairsTT, pairsTTM,method,"TT","TTM")
                                                                confoundersTG = findCounfounders(pairsTG, pairsTTG,method,"TG","TTG")
                                                                confoundersTTG = findCounfounders(pairsTT, pairsTTG,method,"TT","TTG")

                                                                #Deconfound the deconfounder
                                                                confoundersOfConfoundersTTM = findCounfounders(confounderListToPairs(confoundersT), pairsTTM,method,"confoundersOfT","TTM")
                                                                confoundersOfConfoundersTGM = findCounfounders(confounderListToPairs(confoundersT), pairsTGM,method,"confoundersOfT","TGM")
                                                                confoundersOfConfoundersTT = findCounfounders(confounderListToPairs(confoundersT), pairsTT,method,"confoundersOfT","TT")
     
                                                                
                                                                # print(stats)
                                                                stats = rbind(confoundersT,confoundersM,confoundersG,confoundersTM,confoundersGM,confoundersTTM,confoundersTG,confoundersTTG,
                                                                              confoundersOfConfoundersTTM,confoundersOfConfoundersTGM,confoundersOfConfoundersTT)
                                                                if (is.null(stats)){
                                                                        next
                                                                }
                                                                
                                                                write.table(stats, paste(resultsFolder,outputFile, sep=""), sep=",", append=F, row.names = F,quote=F, col.names = T)
                                                                print(stats)
                                                                

                                                                # amountOfConfounding = c()
                                                                # colnames(amountOfConfounding) = c("Name","#EdgesInT","#EdgesInTT","#EdgesInGM","#EdgesInTG","#EdgesInTM","#EdgesInTGM","#EdgesInTTG","#EdgesInTTM","#EdgesInTTGM")
                                                                
                                                                
        
                                                                
                                                                #combine all the small outputs into one large matrix, together with the file that they
                                                                # came from, and then do a group by with a count
                                                                stats = cbind(rep(outputFile,nrow(stats)),stats)
                                                                allInteractions = rbind (allInteractions,stats)
                                                                
                                                                
                                                                #compute some stats about this network
                                                                combinedLevels = c()
                                                                for (i in 1:nrow(stats)){
                                                                        combinedLevels = c(combinedLevels,paste(stats[i,10],stats[i,11],sep="_dw_"))        
                                                                }
                                                                
                                                                dfs = as.data.frame(combinedLevels)
                                                                
                                                                statsOfNetwork = sqldf("select combinedLevels,count(*) from dfs group by combinedLevels")
                                                                colnames = statsOfNetwork[,1]
                                                                statsOfNetwork = statsOfNetwork[,-1]
                                                                statsOfNetwork = as.data.frame(t(statsOfNetwork))
                                                                colnames(statsOfNetwork) = colnames
                                                                # statsOfNetwork$NumberOfEdges = outputFile
                                                                statsOfNetwork$Method = method
                                                                statsOfNetwork$TotalEdgesAllNetworks = NROW(pairsT)+ NROW(pairsTT)+ NROW(pairsTG)+ NROW(pairsTTG)+ NROW(pairsTM)+ NROW(pairsTTM)+ NROW(pairsTGM)+ NROW(pairsGM)
                                   
                                                                
                                                                statsOfNetwork$Percentage_T_dw_TG_Deconfounding = ifelse(length(statsOfNetwork$T_dw_TG/NROW(pairsT)*100)!=0, statsOfNetwork$T_dw_TG/NROW(pairsT)*100, 0)
                                                                statsOfNetwork$Percentage_T_dw_TM_Deconfounding = ifelse(length(statsOfNetwork$T_dw_TM/NROW(pairsT)*100)!=0, statsOfNetwork$T_dw_TM/NROW(pairsT)*100, 0)
                                                                statsOfNetwork$Percentage_T_dw_TT_Deconfounding = ifelse(length(statsOfNetwork$T_dw_TT/NROW(pairsT)*100)!=0, statsOfNetwork$T_dw_TT/NROW(pairsT)*100, 0)
                                                                statsOfNetwork$Percentage_GM_dw_TGM_Deconfounding = ifelse(length(statsOfNetwork$GM_dw_TGM/NROW(pairsGM)*100)!=0, statsOfNetwork$GM_dw_TGM/NROW(pairsGM)*100, 0)
                                                                statsOfNetwork$Percentage_TG_dw_TTG_Deconfounding = ifelse(length(statsOfNetwork$TG_dw_TTG/NROW(pairsTG)*100)!=0, statsOfNetwork$TG_dw_TTG/NROW(pairsTG)*100, 0)
                                                                statsOfNetwork$Percentage_TM_dw_TTM_Deconfounding = ifelse(length(statsOfNetwork$TM_dw_TTM/NROW(pairsTM)*100)!=0, statsOfNetwork$TM_dw_TTM/NROW(pairsTM)*100, 0)
                                                                statsOfNetwork$Percentage_TT_dw_TTG_Deconfounding = ifelse(length(statsOfNetwork$TT_dw_TTG/NROW(pairsTT)*100)!=0, statsOfNetwork$TT_dw_TTG/NROW(pairsTT)*100, 0)
                                                                statsOfNetwork$Percentage_TT_dw_TTM_Deconfounding = ifelse(length(statsOfNetwork$TT_dw_TTM/NROW(pairsTT)*100)!=0, statsOfNetwork$TT_dw_TTM/NROW(pairsTT)*100, 0)
                                                            
                                                                
                                                                statsOfNetwork$edgesT = NROW(pairsT)
                                                                statsOfNetwork$edgesTT = NROW(pairsTT)
                                                                statsOfNetwork$edgesTG = NROW(pairsTG)
                                                                statsOfNetwork$edgesTM = NROW(pairsTM)
                                                                statsOfNetwork$edgesGM = NROW(pairsGM)
                                                                statsOfNetwork$edgesTTG = NROW(pairsTTG)
                                                                statsOfNetwork$edgesTTM = NROW(pairsTTM)
                                                                statsOfNetwork$edgesTGM = NROW(pairsTGM)
                                                                
                                                                statsOfNetwork$overallScore = mean(as.numeric(c(pairsT[,3],pairsTT[,3],pairsTG[,3],pairsTM[,3],pairsGM[,3],pairsTTG[,3],pairsTTM[,3],pairsTGM[,3])))
                                                                
                                                                statsOfNetwork$Name = outputFile
                                                                statsOfNetwork$Alignment = alignment
                                                                statsOfNetwork[is.na(statsOfNetwork)] = 0
                                                        
                                                                overallStats = dplyr::bind_rows(statsOfNetwork, overallStats)
                                                                #Also adds some percentages, not just absolute values. And add info to see if T deconfounds more than M or G
                                                                #unrolled T_dw_TT is calculated dividing by #T
                                                                print("a")
                                                                
                                                               
                                                                
                                                        }
                                                }
                                        }
                                }
                        }
                        closeAllConnections()
                }
                        
        }
        

        #compute overall stats for method and alignment
        properColumnNames = c(0,0,0,0,0,0,0,0)
        names(properColumnNames) = c("GM_dw_TGM","TG_dw_TTG","TM_dw_TTM","TT_dw_TTG","TT_dw_TTM","T_dw_TG","T_dw_TM","T_dw_TT")
        overallStats = dplyr::bind_rows(overallStats, properColumnNames) # make sure all the column names are here
        overallStats = overallStats[-nrow(overallStats),]
        overallStatsByMethod = sqldf("select Method,Alignment,avg(Percentage_T_dw_TG_Deconfounding),avg(Percentage_T_dw_TM_Deconfounding),avg(Percentage_T_dw_TT_Deconfounding),avg(Percentage_GM_dw_TGM_Deconfounding),avg(Percentage_TG_dw_TTG_Deconfounding),avg(Percentage_TM_dw_TTM_Deconfounding),avg(Percentage_TT_dw_TTG_Deconfounding),avg(Percentage_TT_dw_TTM_Deconfounding),avg(GM_dw_TGM),avg(TG_dw_TTG),avg(TM_dw_TTM),avg(TT_dw_TTG),avg(TT_dw_TTM),avg(T_dw_TG),avg(T_dw_TM),avg(T_dw_TT),avg(edgesT),avg(edgesTT),avg(edgesTG),avg(edgesTM),avg(edgesGM),avg(edgesTTG),avg(edgesTTM),avg(edgesTGM),avg(overallScore),sum(TotalEdgesAllNetworks),avg(TotalEdgesAllNetworks) from overallStats group by Method, Alignment")
        colnames(overallStatsByMethod) =   gsub("\\)","",gsub("avg\\(","",gsub("Percentage","%", (gsub("edges","#Edges",colnames(overallStatsByMethod))))))
        write.table(overallStatsByMethod, "overallDeconfoundingStatsByMethod.csv", sep=",", append=F, row.names = F,quote=F, col.names = T, na="0")

        colnames(overallStats) = gsub("Percentage","%", (gsub("edges","#Edges",colnames(overallStats))))
        
        overallStats = overallStats[,c("Method","Alignment","Name","%_T_dw_TG_Deconfounding" ,"%_T_dw_TM_Deconfounding", "%_T_dw_TT_Deconfounding", "%_GM_dw_TGM_Deconfounding",
                                       "%_TG_dw_TTG_Deconfounding", "%_TM_dw_TTM_Deconfounding", "%_TT_dw_TTG_Deconfounding", "%_TT_dw_TTM_Deconfounding",
                                       "GM_dw_TGM","TG_dw_TTG","TM_dw_TTM","TT_dw_TTG","TT_dw_TTM","T_dw_TG","T_dw_TM","T_dw_TT",      
                                       "#EdgesT","#EdgesTT","#EdgesTG","#EdgesTM","#EdgesGM","#EdgesTTG","#EdgesTTM","#EdgesTGM",
                                       "overallScore","TotalEdgesAllNetworks")]
        
        # overallStats = overallStats[,c(10,30,29,12:19,1:9,20:27, 11,28)]
        write.table(overallStats, "overallDeconfoundingStats.csv", sep=",", append=F, row.names = F,quote=F, col.names = T, na="0")
        
        return(allInteractions)
}

#Compute the overall stats by seeing what interactios appear more times
findMostComonConfounder = function (allInteractions){
        combineStatsHash <- hash() 
        
        for (i in 1:nrow(allInteractions)){
                thisRowConf = paste(allInteractions[i,2],"&",allInteractions[i,3],"&",allInteractions[i,4], sep="")
                thisRowIdent =  paste(allInteractions[i,1],"&",allInteractions[i,10],"&",allInteractions[i,11],"$",allInteractions[i,8], sep="")
                if(!is.null(combineStatsHash[[thisRowConf]])){
                        combineStatsHash[[thisRowConf]] = c(combineStatsHash[[thisRowConf]],thisRowIdent)
                }else{
                        combineStatsHash[[thisRowConf]] = c(thisRowIdent)
                }
        }
        
        combinedStats = c()
        for (k in keys(combineStatsHash)){
                splitKey = strsplit(k,"&")
                #Calculate the average overall score
                overallScore = 0
                for (instance in combineStatsHash[[k]]){
                        if (is.na(as.numeric(strsplit(instance,"\\$")[[1]][2]))){
                                next
                        }
                        overallScore = overallScore + as.numeric(strsplit(instance,"\\$")[[1]][2])
                }
                overallScore = overallScore/length(combineStatsHash[[k]])
               
                combinedStats = rbind(combinedStats,c(splitKey[[1]][1],splitKey[[1]][2],splitKey[[1]][3],overallScore,length(combineStatsHash[[k]]),paste(combineStatsHash[[k]],collapse="|")))
        }
        colnames(combinedStats) = c("confounder","confoundedParent","confoundedChild","overallScore","timesItAppears","WhereItAppears")
        combinedStats <-combinedStats[order(combinedStats[,5],combinedStats[,4],decreasing=T),]
        
        closeAllConnections()
        write.table(combinedStats,  paste(resultsFolder,"mostCommonConfounders.csv",sep=""), sep=",", append=F, row.names = F,quote=F, col.names = T)
        
}




allInteractions = deconfound()

findMostComonConfounder(allInteractions)


