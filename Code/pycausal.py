# -*- coding: utf-8 -*-
import pandas as pd
from pycausal import prior as p
from pycausal import search as s
from pycausal.pycausal import pycausal as pc
import networkx as nx
import javabridge
import os
import glob
import pandas as pd
import pydot
from IPython.display import SVG
import matplotlib.pyplot as plt
import re
from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import time
import re

'''

Receives as input multiple subjects temporal data and a matrix of allowed interactions
Creates proper time lag dataset removing the timepoints between subjects to avoid temporal inconsistencies and creates object with prior knowledge
Performs bootstrap (sampling with replacement" by returning all edges learnt in any repetition
Outputs a graphml with separated edge for every edge type, making sure they are in the correct order

module load jdk-11.0.1-gcc-8.2.0-ref6fpq
module load miniconda3-4.5.11-gcc-8.2.0-oqs2mbg
#conda create -n pycausal python=3.7
source activate pycausal 
conda install numpy
conda install pandas
conda install pydot
conda install networkx
conda install IPython
conda install matplotlib
conda install joblib
conda install cython
conda install -c pwwang glibc214
conda install -c conda-forge javabridge
'''


def startJVM():
    tetrad_libdir = os.path.join(os.getcwd(), 'src', 'pycausal', 'lib')
    for l in glob.glob(tetrad_libdir + os.sep + "*.jar"):
        javabridge.JARS.append(str(l))

    javabridge.start_vm(run_headless=True, max_heap_size='4000M')
    javabridge.attach()
    return javabridge

def readDataAllSubjects(filename):
    with open(filename) as f:
        data = f.readlines()
    prevSubject = -1
    acum = []
    indexOfSubjectChange = []
    var_names_ascii = ""
    for i,line in enumerate(data, start=0):
        line = line.split("\t")
        if (i == 0):
            var_names = line[2:]
            var_names_ascii = [name.replace(":","").replace("-","").replace(" ","_").replace("(","").replace(")","").replace("[","").replace("]","").replace("+","").replace("\n","") for name in var_names]
            continue
        line = [float(k) for k in line]
        if (i > 1 and prevSubject != line[0]):
            indexOfSubjectChange += [i-2]
        acum += [line[2:]]
        prevSubject = line[0]
    #We save the last subject
    data = pd.DataFrame(acum)
    data.columns = [name.split("$")[0] for name in var_names_ascii]
    return data,var_names_ascii,indexOfSubjectChange

def combineSubjectData(dataPerSubject):
    acum = pd.DataFrame()
    for k,v in dataPerSubject.items():
        if (len(acum) == 0):
            acum = v
        else:
            acum  = pd.concat([acum, v], axis=0)
    return acum

def readData(filename):
    with open(filename) as f:
        data = f.readlines()
    acum = []
    var_names = ""
    var_names_ascii = ""
    for i,line in enumerate(data, start=0):
        line = line.split("\t")
        if (i == 0):
            var_names = line[2:]
            var_names_ascii = [name.replace(":","").replace("-","").replace(" ","_").replace("(","").replace(")","").replace("[","").replace("]","").replace("+","").replace("\n","") for name in var_names]
            continue
        line = [float(k) for k in line]
        acum += [line[2:]]

    data = pd.DataFrame(acum)
    data.columns = [name.split("$")[0] for name in var_names_ascii]
    return data,var_names_ascii, var_names

def transformData(dframe,javabridge, indexOfSubjectChange):
    node_list = javabridge.JClassWrapper("java.util.ArrayList")()
    for col in dframe.columns:
        nodi = javabridge.JClassWrapper("edu.cmu.tetrad.data.ContinuousVariable")(col)
        node_list.add(nodi)

    dataBox = javabridge.JClassWrapper("edu.cmu.tetrad.data.DoubleDataBox")(len(dframe.index), dframe.columns.size)

    for col in range(0, dframe.columns.size):
        for row in dframe.index:
            value = javabridge.JClassWrapper("java.lang.Double")(dframe.iloc[row, col])
            dataBox.set(row, col, value)

    boxData = javabridge.JClassWrapper("edu.cmu.tetrad.data.BoxDataSet")(dataBox, node_list)

    tetradData = javabridge.static_call('edu/cmu/tetrad/search/TimeSeriesUtils', 'createLagData','(Ledu/cmu/tetrad/data/DataSet;I)Ledu/cmu/tetrad/data/DataSet;', boxData, numLags)

    #We remove the rows in the indices of subject change, to avoid creating a lag involving two subjects
    td = javabridge.JClassWrapper("edu.cmu.tetrad.data.BoxDataSet")(tetradData)
    rowsToRemove = javabridge.get_env().make_int_array(np.array(indexOfSubjectChange, np.int32))
    # rowsToRemove1 = javabridge.JClassWrapper("java.util.Arrays").copyOf(rowsToRemove, 5)
    td.removeRows(rowsToRemove)

    return tetradData

#Export network to graphml where each edge link is separated in a new edge with its own bootscore
def exportToGraphml(tetradGraph, graphName, var_names_ascii):

    def formatName(oldName):
        if (oldName[-2:] == ":1"):
            return(oldName[:-2]+"_ti")
        else:
            return (oldName+"_ti+1")

    text_str = javabridge.static_call('edu/cmu/tetrad/graph/GraphUtils','graphToText','(Ledu/cmu/tetrad/graph/Graph;)Ljava/lang/String;',tetradGraph)
    print(text_str)
    G = nx.MultiDiGraph()
    G.add_nodes_from([formatName(name.split("$")[0]) for name in var_names_ascii])

    for line in text_str.split("\n")[4:]:
        if line != "\r" and line != '':
            edgeList = re.findall(r"\[.+$", line[:-1])[0].split(";")
            for edgeIter in edgeList:
                edge, bootScore = edgeIter.split("]:")
                if (edge != "[no edge"):
                    source, link, target = edge[1:].split(" ")
                    #if it is encoded target <- source instead of source -> target, reverse it
                    if (target[-2:] == ":1" and source[-2:] != ":1"):
                        if (link[0] == "<"):
                            if  link[2] == ">":
                                link = "<-->"
                            else:
                                link = link[2]+"->"
                            source, target = target, source
                    G.add_edges_from([(formatName(source),formatName(target) ,{'link': link, 'bootscore': bootScore})])
    nx.write_graphml_xml(G, graphName)

#Export to graphml where  all the links of one edge are exported in the same edge
def exportToGraphmlOld(tetradGraph, graphName):
    print(tetradGraph)
    dot_str = javabridge.static_call('edu/cmu/tetrad/graph/GraphUtils','graphToDot','(Ledu/cmu/tetrad/graph/Graph;)Ljava/lang/String;',tetradGraph)
    graphs = pydot.graph_from_dot_data(dot_str)
    G = nx.nx_pydot.from_pydot(graphs[0])
    # nx.draw(G)
    # mAdj = nx.adjacency_matrix(G, nodelist=var_names_ascii)
    nx.draw_networkx(G,node_size =100,font_size=8 )
    plt.show()
    nx.write_graphml_xml(G, graphName)

#Import restriction matrix and generate prior knowledge object
def generateAllowedMatrix(filename,var_names_ascii,javabridge):
    #  0 := no inter-edges and no intra-edges
    #  1 := self-edges only
    #  2 := inter-edges only
    #  3 := intra-edges only
    #  4 := inter- and intra-edges

    prior = javabridge.JClassWrapper('edu.cmu.tetrad.data.Knowledge2')()

    with open(filename, 'r') as f:
        types = f.readline().replace("\n","").split('\t')
        m = [[int(num) for num in line.split('\t')] for line in f.readlines()]

        for name in var_names_ascii:
            prior.addToTier(0, name.split("$")[0]  +':1')
            prior.addToTier(1, name.split("$")[0])

        # .split("$")[0]
        for i, typei in enumerate(types):
            namesOfThisTypei = [name.split("$")[0] for name in list(filter(re.compile(".*\$" + typei).match, var_names_ascii))]
            for j, typej in enumerate(types):
                namesOfThisTypej = [name.split("$")[0] for name in list(filter(re.compile(".*\$"+typej).match, var_names_ascii))]
                for namei in namesOfThisTypei:
                    #Add self edges Unless they are forbidden
                    if (m[i][j] == 0):
                        prior.setForbidden(namei+':1', namei)
                    # else:
                        prior.setRequired(namei+':1', namei)

                    for namej in namesOfThisTypej:
                        #Forbid all incoming edges to the _ti tier
                        # prior.setForbidden(namei+':1', namej+':1')
                        # prior.setForbidden(namei, namej+':1')
                        # prior.setForbidden(namej + ':1', namei + ':1')
                        # prior.setForbidden(namej, namei + ':1')

                        #If it is inter edge, then forbid all intra edges
                        if (m[i][j] == 2 or m[i][j] == 1 or m[i][j] == 0):
                            prior.setForbidden(namei,namej)
                            prior.setForbidden(namei+ ':1', namej + ':1')
                        #Add it is intra edge, then forbid all inter edges
                        if (m[i][j] == 3 or m[i][j] == 1 or m[i][j] == 0):
                            prior.setForbidden(namei+':1',namej)

    return prior

#Run without bootstrap
def runAlgorithm(tetradData,alpha, javabridge, prior, indTestName, scoreName, maxIndegree):
    IndTest = javabridge.JClassWrapper('edu.cmu.tetrad.search.IndTestFisherZ')(tetradData, alpha)
    # IndTest = javabridge.JClassWrapper(indTestName)(tetradData, 0.01)
    score = javabridge.JClassWrapper('edu.cmu.tetrad.search.SemBicScore')(tetradData)
    # score = javabridge.JClassWrapper(scoreName)(tetradData)
    score.setPenaltyDiscount(penaltydiscount)
    tsgfci = javabridge.JClassWrapper('edu.cmu.tetrad.search.TsGFci')(IndTest, score)
    tsgfci.setMaxPathLength(maxPathLength)
    tsgfci.setCompleteRuleSetUsed(completeRuleSetUsed)
    tsgfci.setFaithfulnessAssumed(faithfulnessAssumed)
    tsgfci.setVerbose(verbose)
    tsgfci.setMaxIndegree(maxIndegree)
    tsgfci.setKnowledge(prior)

    tetradGraph = tsgfci.search()
    return tetradGraph

def runAlgorithmBootstrap(tetradData,alpha, javabridge, prior,numBootstrap, indTestName, scoreName, maxIndegree):
    #alpha is getting passed
    indTest = javabridge.JClassWrapper(indTestName)()
    score = javabridge.JClassWrapper(scoreName)()
    algorithm = javabridge.JClassWrapper('edu.cmu.tetrad.algcomparison.algorithm.oracle.pag.TsGfci')(indTest, score)

    parameters = javabridge.JClassWrapper('edu.cmu.tetrad.util.Parameters')()
    parameters.set('penaltyDiscount', penaltydiscount)
    parameters.set('maxPathLength', maxPathLength)
    parameters.set('alpha', alpha)
    parameters.set('completeRuleSetUsed', completeRuleSetUsed)
    parameters.set('faithfulnessAssumed', faithfulnessAssumed)
    parameters.set('verbose', verbose)
    parameters.set('percentResampleSize',percentResampleSize)
    parameters.set('fasRule', fasRule) #Use concurrent PC Stable
    parameters.set('maxIndegree', maxIndegree) #Max number of parents
    parameters.set('maxDegree', maxIndegree) #Max number of parents
    # parameters.set('saveLatentVars', True) #Max number of parents
    # parameters.set('numLatents', 5) #Max number of parents

    tsgfci = javabridge.JClassWrapper('edu.pitt.dbmi.algo.resampling.GeneralResamplingTest')(tetradData, algorithm, numBootstrap)
    edgeEnsemble = javabridge.get_static_field('edu/pitt/dbmi/algo/resampling/ResamplingEdgeEnsemble',ensembleMethod,'Ledu/pitt/dbmi/algo/resampling/ResamplingEdgeEnsemble;')

    tsgfci.setEdgeEnsemble(edgeEnsemble)
    tsgfci.setParameters(parameters)
    tsgfci.setVerbose(verbose)

    tsgfci.setKnowledge(prior)
    tetradGraph = tsgfci.search()
    return tetradGraph

def main(javabridge, filename, alpha, numBootstrap, indTestName, scoreName, maxIndegree):
    print("STARTING Executing " + filename + " for " + str(numBootstrap) + " repetitions. Alpha: "+str(alpha)+". scoreName: " + scoreName + " and indTestName: " + indTestName )

    data, var_names_ascii, indexOfSubjectChange = readDataAllSubjects(filename=filename+".tsv")

    if ("_t_g_m" in filename):
        allowedMatrix = "../data/HEGTM_Skeleton.matrix"
    elif ("_t_g" in filename):
        allowedMatrix = "../data/TG.matrix"
    elif ("t_m" in filename):
        allowedMatrix = "../data/TM.matrix"
    else:
        allowedMatrix = "../data/HEGTM.matrix"
    prior = generateAllowedMatrix(allowedMatrix, var_names_ascii,javabridge)
    tetradData = transformData(dframe=data,javabridge=javabridge, indexOfSubjectChange=indexOfSubjectChange)

    start_time = time.time()

    # tetradGraph = runAlgorithm(tetradData = tetradData, alpha=alpha, javabridge = javabridge, prior = prior, indTestName = indTestName, scoreName = scoreName,maxIndegree = maxIndegree)
    tetradGraph = runAlgorithmBootstrap(numBootstrap= numBootstrap,tetradData = tetradData, alpha=alpha, javabridge = javabridge, prior = prior, indTestName = indTestName, scoreName = scoreName, maxIndegree = maxIndegree)

    execTime = time.time() - start_time
    print("--- %s seconds ---" % (execTime))

    # + "_bs" + str(bootstrap_level)\
    outputName = "../networks/" + "PyCausal_" + filename.split("/")[len(filename.split("/")) - 1] + "_a" + str(alpha)\
                 + "_nboots" + str(numBootstrap) + "_score" + str(scoreName.split(".")[len(scoreName.split(".")) - 1]) +\
                 "_test" + indTestName.split(".")[len(indTestName.split(".")) - 1] +\
                 "_nParents" +str(maxIndegree) + "_matrix" + allowedMatrix.split("/")[len(allowedMatrix.split("/")) - 1] +\
                 "_s"+str(int(execTime))+".graphml"

    exportToGraphml(tetradGraph=tetradGraph, graphName = outputName, var_names_ascii = var_names_ascii)
    print("FINISHED Executing " + filename + " for " + str(numBootstrap) + " repetitions. ALpha: "+str(alpha)+" with allowed matrix: " + allowedMatrix + ". scoreName: " + scoreName + " and indTestName: " + indTestName )

    return javabridge






numLags=1
penaltydiscount = 4 # set to 2 if variable# <= 50 otherwise set it to 4
maxPathLength = -1 #The maximum length for any discriminating path. -1 if unlimited
completeRuleSetUsed = True
faithfulnessAssumed = False
verbose = False
ensembleMethod = 'Preserved' #Preserved (all edges), Highest (returns the edge orientation which the highest proportion of sample graphs returned), and Majority (at least 50 percent of the sample graphs agree on an edge orientation in order to return any edge at all)
percentResampleSize = 0.9
fasRule = 3  # Use concurrent PC Stable
maxIndegree = 3 #-1 is unlimited. But it looks like this parameter is somewhat ignored
# indTestName = 'edu.cmu.tetrad.search.IndTestFisherZ', scoreName = 'edu.cmu.tetrad.search.SemBicScore'
#https://rdrr.io/github/bd2kccd/r-causal/f/

javabridge = startJVM()

indTestNames = ["MultinomialLogisticRegressionWald","ConditionalGaussianLRT", "DegenerateGaussianLRT","FisherZ","Kci","MNLRLRT","PositiveCorr","SemBicTest"]
scoreNames = ["DiscreteMixedBicScore", "ConditionalGaussianBicScore","DegenerateGaussianBicScore","FisherZScore","MNLRBicScore","MVPBicScore","PeterScore","SemBicScore","SemBicScoreDeterministic"]

indTestNames = ["PositiveCorr"]
scoreNames = ["FisherZScore"]
# 0.01

alphaList = [0.0001,0.001]
alphaList = [0.001,0.01]
# alphaList = [0.01, 0.1,0.5]
dataFolder = "../data/"
filenames = []

for file in os.listdir(dataFolder):
    if file.endswith(".tsv"):
        filenames += [dataFolder + file.split(".tsv")[0]]

for alpha in alphaList:
    for test in indTestNames:
        for score in scoreNames:
            for file in filenames:
                javabridge = main(javabridge= javabridge, filename=dataFolder+file, alpha = alpha, numBootstrap = 10,  indTestName = 'edu.cmu.tetrad.algcomparison.independence.'+test,scoreName = 'edu.cmu.tetrad.algcomparison.score.'+score,maxIndegree = maxIndegree)

javabridge.detach()
javabridge.kill_vm()
