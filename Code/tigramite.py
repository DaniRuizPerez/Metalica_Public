'''

Because tigramite doesn't support runing multiple subjects at the same time I have to run one by one
because of this, bootstrapping doesn't make sense because I would have to select different subjects, but each subect's execution yields the same results
bootstrapping at the timepoint level doesn't make sense because subjects already have too few timepoints, so subsampling there is a bad idea
I'm removing the columns of a subject that are mostly zeroes.

'''


import copy
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import sklearn
from sklearn import preprocessing
import tigramite
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb, ParCorrBIC
from tigramite.models import LinearMediation, Prediction
from joblib import Parallel, delayed
import multiprocessing
import networkx as nx
import re
import time
import os


#import rpy2.robjects.packages as rpackages

# import R's utility package
#utils = rpackages.importr('utils')

# select a mirror for R packages
#utils.chooseCRANmirror(ind=1) # select the first mirror in the list
#utils.chooseBioCmirror(ind=1)

#utils.install_packages('remotes')
#remotes = rpackages.importr('remotes')
#remotes.install_github("ericstrobl/RCIT")
##utils.install_packages('devtools')
##utils.install_packages('cli')
##utils.install_packages('usethis')



'''
#conda create -n tigramite python=3.7
source activate tigramite 
conda install numpy
conda install scipy
conda install scikit-learn
conda install matplotlib
conda install networkx
conda install cython
conda install mpi4py
conda install rpy2
'''


# Setup analysis
np.random.seed(42)  # Fix random seed

def readData(filename):
    scale = False #Doesn't work very well
    scaler = None #Used to scale all data with the same tansformation
    with open(filename) as f:
        data = f.readlines()
    dataPerSubject = dict()
    prevSubject = -1
    acum = []
    var_names = ""
    for i,line in enumerate(data, start=0):
        line = line.split("\t")
        if (i == 0):
            var_names = line[2:]
            continue
        line = [float(k) for k in line]
        if (i > 1 and prevSubject != line[0]):
            if (scaler is None and scale):
                scaler = preprocessing.StandardScaler().fit( np.asarray(acum))
            if scale:
                subject = scaler.transform(np.asarray(acum))
            else:
                subject = (np.asarray(acum))
            dataPerSubject[prevSubject] = subject
            acum = [line[2:]]
        else:
            acum += [line[2:]]
        prevSubject = line[0]
    #We save the last subject
    if scale:
        subject = scaler.transform(np.asarray(acum))
    else:
        subject = (np.asarray(acum))
    dataPerSubject[prevSubject] = subject
    return dataPerSubject,var_names

# @deprecated
def get_sig_links_Ïgnore_constant(data,tau_max,alpha_level,allowedEdges,test):
    #remove constant columns.
    notConstantColumns = data.var(0) != 0
    indexOfFalse = np.asarray(np.where(notConstantColumns == False))[0]
    data = data[:,notConstantColumns]

    # PCMCI
    pcmci = PCMCI(dataframe=pp.DataFrame(data), cond_ind_test=test)
    # pcmci = PCMCI(dataframe=pp.DataFrame(data), cond_ind_test=CMIknn())
    results = pcmci.run_pcmci(selected_links=allowedEdges,tau_min=0,tau_max=tau_max, pc_alpha=alpha_level)
    # results = pcmci.run_pcmci(tau_min=0,tau_max=tau_max, pc_alpha=alpha_level)

    #Reshape the matrix to have zeroes where the attribute is constant
    for index in indexOfFalse:
        results['p_matrix'] = np.insert(results['p_matrix'],obj = index,values = 1, axis = 1 )
        results['p_matrix'] = np.insert(results['p_matrix'],obj = index,values = 1, axis = 0 )
        results['val_matrix'] = np.insert(results['val_matrix'],obj = index,values = 0, axis = 1 )
        results['val_matrix'] = np.insert(results['val_matrix'],obj = index,values = 0, axis = 0 )

    return  results['p_matrix'], results['val_matrix']

#This sets to a random vector the constant attributes
def runPCMCI(data,tau_max,alpha_level,allowedEdges,test):
    #remove constant columns.
    notConstantColumns = data.var(0) != 0
    notMostlyEmptyColumns = np.count_nonzero(data, axis=0)/data.shape[0] >= 0.2
    notConstantFinal = np.logical_and(notConstantColumns,notMostlyEmptyColumns)
    # print(notConstantFinal)
    indexOfFalse = np.asarray(np.where(notConstantFinal == False))[0]
    # print(indexOfFalse)

    for element in indexOfFalse:
        data[:,element] =  np.random.rand(1,len(data))/100000
    # PCMCI

    #Scale the data
    data  = preprocessing.scale(data)

    pcmci = PCMCI(dataframe=pp.DataFrame(data), cond_ind_test=test, verbosity= 0)
    # pcmci = PCMCI(dataframe=pp.DataFrame(data), cond_ind_test=CMIknn())
    # results = pcmci.run_pcmci(selected_links=allowedEdges,tau_min=0,tau_max=tau_max, fdr_method='fdr_bh', pc_alpha= alpha_level)
    results = pcmci.run_pcmciplus(selected_links=allowedEdges,tau_min=0,tau_max=tau_max, fdr_method='fdr_bh', pc_alpha= alpha_level)

    return  results['p_matrix'], results['val_matrix'], results['q_matrix'],

def generateGraph (var_names,ave_val_matrix,sig_links,tau_max,bootstrap_matrix,outputName):
    pastVariables = ([e +"_ti" for e in var_names])
    futureVariables = ([e +"_ti+1" for e in var_names])

    G = nx.MultiDiGraph()
    G.add_nodes_from(pastVariables)
    G.add_nodes_from(futureVariables)

    for i in range(len(pastVariables)):
        for j in range(len(futureVariables)):
            for t in range(tau_max+1):
                if (bootstrap_matrix[i,j,t]):
                    if (t == 0):
                        # G.add_edges_from([(pastVariables[i], pastVariables[j], {'weight': ave_val_matrix[i,j,t],'bootscore': sig_links[i,j,t]})])
                        G.add_edges_from([(futureVariables[i], futureVariables[j], {'weight': ave_val_matrix[i,j,t],'bootscore': sig_links[i,j,t]})])
                    elif(t==1):
                        G.add_edges_from([(pastVariables[i], futureVariables[j], {'weight': ave_val_matrix[i,j,t],'bootscore': sig_links[i,j,t]})])

    # for line in nx.generate_graphml(G):
    #     print(line)

    nx.write_graphml_xml(G, outputName)
    return G

def generateAllowedMatrix(filename,var_names):
    #  0 := no inter-edges and no intra-edges
    #  1 := self-edges only
    #  2 := inter-edges only
    #  3 := intra-edges only
    #  4 := inter- and intra-edges
    if ("_t_g_m" in filename):
        matrixname = "HEGTM_Skeleton.matrix"
    elif ("_t_g" in filename):
        matrixname = "TG.matrix"
    elif ("t_m" in filename):
        matrixname = "TM.matrix"
    else:
        matrixname = "HEGTM.matrix"
    allowedEdges = dict()
    with open(dataFolder + matrixname, 'r') as f:
        types = f.readline().replace("\n","").split('\t')
        m = [[int(num) for num in line.split('\t')] for line in f.readlines()]

        for i,varname in enumerate(var_names):
            allowedEdges[i] = set()

        for i, typei in enumerate(types):
            namesOfThisTypei = list(filter(re.compile(".*\$" + typei).match, var_names))
            for j, typej in enumerate(types):
                namesOfThisTypej = list(filter(re.compile(".*\$"+typej).match, var_names))
                for namei in namesOfThisTypei:
                    n1 = var_names.index(namei)
                    #Add self edges
                    if (m[i][j] > 0):
                        allowedEdges[n1].add((n1, -1))
                    for namej in namesOfThisTypej:
                        n2 = var_names.index(namej)
                        #Add inter edges
                        if (m[i][j] == 2 or m[i][j] == 4):
                            allowedEdges[n2].add((n1, -1))
                        #Add intra edges
                        if (m[i][j] == 3 or m[i][j] == 4):
                            allowedEdges[n2].add((n1, 0))
    for i,varname in enumerate(var_names):
        allowedEdges[i] = [(a,b) for a,b in allowedEdges[i]]
        aux = varname + "" + str([(var_names[a],b) for a,b in allowedEdges[i]])
        print(aux)
    return allowedEdges,matrixname
# Randomly select as many subjects as the number of subjects we have

def main (filenames, tau_max,alpha_levels,pvalue_threshold, bootstrap_level,testNames):
    parallel = True
    for filename in filenames:
        for alpha_level in alpha_levels:
            for testName in testNames:
                start_time = time.time()
                print("Processing: " + filename.split("/")[len(filename.split("/"))-1]+"_bs"+str(bootstrap_level)+"_a"+str(alpha_level)+\
                             "_ptreshold"+str(pvalue_threshold)+"_T"+str(tau_max)+"_test_"+testName)
                test = ParCorr()
                if testName =="ParCorr":
                    test = ParCorr()
                if testName == "ParCorrBIC":
                    test = ParCorrBIC()
                elif testName =="GPDC":
                    test = GPDC()
                elif testName =="CMIknn":
                    test = CMIknn()
                elif testName =="CMIsymb":
                    test = CMIsymb()

                dataPerSubject, var_names = readData(filename+".tsv")
                allowedEdges,matrixname = generateAllowedMatrix(filename,var_names)

                var_names = [s.split('$')[0] for s in var_names]

                p_matrix = np.ones((len(dataPerSubject.keys()), len(var_names), len(var_names), tau_max + 1))
                val_matrix = np.zeros((len(dataPerSubject.keys()), len(var_names), len(var_names), tau_max + 1))
                q_matrix = np.zeros((len(dataPerSubject.keys()), len(var_names), len(var_names), tau_max + 1))

                if (parallel):
                    result = Parallel(n_jobs=min(multiprocessing.cpu_count()-3,len(dataPerSubject.keys())))\
                        (delayed(runPCMCI)(dataPerSubject[subject],tau_max,alpha_level, allowedEdges,test) for subject in dataPerSubject.keys())
                else:
                    for subject in dataPerSubject.keys():
                        if (subject != 23):
                            continue
                        print("Subject: " + str(subject))
                        runPCMCI(dataPerSubject[subject], tau_max, alpha_level, allowedEdges, test)
                    continue

                for i in range(len(dataPerSubject.keys())):
                    p_matrix[i] = result[i][0]
                    val_matrix[i] = result[i][1]
                    q_matrix[i] = result[i][2]

                # Get true positive rate (=power) and false positive rate
                sig_links = (q_matrix <= pvalue_threshold).mean(axis=0)
                ave_val_matrix = val_matrix.mean(axis=0)

                bootstrap_matrix = (sig_links >= bootstrap_level)
                # tp.plot_graph(val_matrix=ave_val_matrix,link_matrix=bootstrap_matrix, var_names=var_names,arrow_linewidth=70)
                fig, ax = tp.plot_time_series_graph(val_matrix=ave_val_matrix,link_matrix=bootstrap_matrix,var_names=var_names, link_colorbar_label='MCI')
                fig.show()

                execTime = time.time() - start_time
                print("--- %s seconds ---" % (execTime))

                outputName = "../networks/"+"Tigramite_"+filename.split("/")[len(filename.split("/"))-1]+"_bs"+str(bootstrap_level)+"_a"+str(alpha_level)+\
                             "_ptreshold"+str(pvalue_threshold)+"_T"+str(tau_max)+"_test_"+testName+matrixname+"_s"+str(int(execTime))+".graphml"

                generateGraph (var_names,ave_val_matrix,sig_links,tau_max,bootstrap_matrix,outputName)





dataFolder = "../data/"
filenames = []
testNames = ["ParCorr","GPDC","CMIknn"]
alphaList = [0.0001,0.001,0.01,0.1]
for file in os.listdir(dataFolder):
    if file.endswith(".tsv"):
        filenames += [dataFolder + file.split(".tsv")[0]]
print(filenames)


main(filenames=filenames, tau_max=1, alpha_levels=alphaList, pvalue_threshold = 0.05, bootstrap_level = 0.1, testNames=testNames)

