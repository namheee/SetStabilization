from pyboolnet.state_transition_graphs import primes2stg
from pyboolnet.attractors import compute_attractors_tarjan
from pyboolnet.file_exchange import bnet2primes
import pyboolnet
import re
import pandas as pd
import numpy as np
from choonFunction import modeltext_transform
from sympy.logic.boolalg import to_dnf, simplify_logic
import itertools
import random
import collections
import copy


def computeAtt(primes):
    stg = primes2stg(primes, "synchronous")
    steady, cyclic = compute_attractors_tarjan(stg)
    return steady, cyclic

def cyclic2str(cyclicSet):
    cyclicList = []
    for cyc in cyclicSet:
        catt0 = pd.DataFrame([list(x) for x in cyc])
        catt = catt0.apply(lambda x: len(set(x.values.tolist())) == 1, axis=0).tolist()
        cyclicList.append(''.join([str(x) if i in np.where([x == True for x in catt])[0] else '*' for i,x in enumerate(catt0.iloc[0,:])]))
    return cyclicList

def readBNETfromPrimes(model_file):
    primes = bnet2primes(model_file)
    net = pyboolnet.file_exchange.primes2bnet(primes)
    net = re.sub("&", "and", net)
    net = re.sub("[|]", "or", net)
    net = re.sub("!", "not ",net)
    net = re.sub(",   ", " = ",net)
    net = re.sub("\n\n", "\n",net)

    net_splited = [(x.split(' ')[0],x) for x in net.split('\n')]
    net_ordered = []
    for nodex in primes.keys():
        for eq in net_splited:
            if eq[0] == nodex:
                net_ordered.append(eq[1])

    # print('\n'.join(net_ordered), end="\n\n")
    modeltext = '\n'.join(net_ordered)
    
    return modeltext


def returnINPUT(model_file, desired, undesired, phenotype, recurL=0):
    modeltext = readBNETfromPrimes(model_file)
    F_net = modeltext_transform(modeltext)
    # print('F_net:','\n', F_net)


    F_net_eq = {x.split(' = ')[0]:x.split(' = ')[1] for x in F_net.strip().split('\n')}

    V_idx = {k[:-1]: k for k in F_net_eq.keys()}
    if recurL != 0: V_l = recurL
    else: V_l = len(V_idx)
    
    V_pair = [(v,'x'+str(int(v[1:-1])+V_l)+'_') if int(v[1:-1])+V_l>= 10 else (v,'x0'+str(int(v[1:-1])+V_l)+'_') for v in V_idx.values()]
    
    # print('V_pair:', V_pair)

    negV = ['x'+str(int(v[1:-1])+V_l)+'_' if int(v[1:-1])+V_l>= 10 else 'x0'+str(int(v[1:-1])+V_l)+'_' for v in V_idx.values()]
    V_idx.update({'~'+k[:-1]:n for n,k in zip(negV, F_net_eq.keys())})
    # print('V_idx:', V_idx)

    def desired2Var(desired, F_net_eq, V_idx):
        desiredVariable = []
        cyclicnode = []
        for datt in desired:
            valueList = []
            for k, dv in zip(F_net_eq.keys(), list(datt)):
                if dv == '0': 
                    valueList.append(V_idx['~'+k[:-1]])
                elif dv == '*':
                    cyclicnode.append(V_idx['~'+k[:-1]])
                    cyclicnode.append(V_idx[k[:-1]])
                    # continue
                else:
                    valueList.append(V_idx[k[:-1]])
            desiredVariable.append(valueList)
        return desiredVariable, cyclicnode

    desiredVariable, cyclicnode = desired2Var(desired, F_net_eq, V_idx)
    undesiredVariable, undesired_cyclicnode = desired2Var(undesired, F_net_eq, V_idx)
    # print(desiredVariable, cyclicnode)
    
    phenotype = [V_idx[x] for x in phenotype]
    cphenotype = [v[1] for v in V_pair if v[0] in phenotype] + [v[0] for v in V_pair if v[1] in phenotype]
    # print(phenotype, cphenotype)

    F_net_simple = ''
    for k,v in F_net_eq.items():    
        expression_dnf = to_dnf(v, simplify=True, force=True)
        F_net_simple += k + ' = ' + ' | '.join(str(expression_dnf).split("|")) + '\n'

    negF_net = ''
    for k,v in F_net_eq.items():    
        expression = " ~ ( " + v + ")"
        expression_dnf = to_dnf(expression, simplify=True, force=True)
        negF_net += V_idx['~'+k[:-1]] + ' = ' + ' | '.join(str(expression_dnf).split("|")) + '\n'

    G_net = F_net_simple + negF_net

    for k,v in V_idx.items():
        G_net = G_net.replace(k, v[:-1])
     
    # print('G_net:','\n', G_net)
    desiredSet = (desiredVariable, cyclicnode)
    undesiredSet = (undesiredVariable, undesired_cyclicnode)
    return G_net, V_pair, V_idx, desiredSet, undesiredSet, phenotype, cphenotype


def checkMutExclusive(equation, V_pair, checkMore):
    vSet = set([x.strip() for x in itertools.chain(*[x.split('&') for x in equation])])  
    checked = True
    
    for v in V_pair:
        if (v[0] in vSet) & (v[1] in vSet):
            chekced = False
            break
        else: continue

    if len(vSet & set(checkMore)) != 0:
        checked = False

    return checked


def ppCombination(G_net, V_pair, layered, checkMore):
    G_net_eq = {x.split(' = ')[0]:x.split(' = ')[1] for x in G_net.strip().split('\n')}

    ppCombiList = []
    for p in layered:
        ppList = []
        for pp in G_net_eq[p].split(' | '):
            pp = pp.replace('(','')
            pp = pp.replace(')','')
            pp = pp.strip()
            ppList.append(pp)
        ppCombiList.append(ppList)


    ppCombiList_updated = []
    for equation in list(itertools.product(*ppCombiList)):
        if checkMutExclusive(equation, V_pair, checkMore):
            ppCombiList_updated.append(equation)

    return layered, ppCombiList_updated

def makeReducedG(G_net, S, Spp):
    G_net_line = G_net.splitlines()
    S_pp = {x:y for x,y in zip(S, Spp)}
    G_net_new = ""
    for line in G_net_line:
        logicList = line.split("=")
        logicOutput = logicList[0].strip()
        logicInput = logicList[1].strip()
        
        if logicOutput in S: # 
            G_net_new += logicOutput + " = " + S_pp[logicOutput] + '\n'
        else:
            G_net_new += logicOutput + " = " + logicInput + '\n'        

    return G_net_new 

def makeReducedForm(G_net):
    # without ''
    G_net_line = G_net.splitlines()
    G_net_new = ""
    for line in G_net_line:
        logicList = line.split("=")
        logicOutput = logicList[0].strip()
        logicInput = logicList[1].strip()
        
        if logicInput != '': # 
            G_net_new += line + '\n'

    return G_net_new


def eliminateComplementS(G_net, complementS):
    G_net_line = G_net.splitlines()

    G_net_new = ""
    for line in G_net_line:
        logicList = line.split("=")
        logicOutput = logicList[0].strip()
        logicInput = logicList[1].strip()

        if logicOutput in complementS: # 
            G_net_new += logicOutput + " = " + '' + '\n'
        else:
            remainedInput = []
            for pp in logicInput.split(' | '):
                pp = pp.strip()
                pp = pp.replace('(','')
                pp = pp.replace(')','')
                pp = pp.split('&')
                pp = set([x.strip() for x in pp])

                if (set(pp) & set(complementS)) != set():
                    continue
                else:
                    remainedInput.append(' & '.join(pp))
            # print(remainedInput)

            G_net_new += logicOutput + " = " + ' | '.join(remainedInput) + '\n'        

    return G_net_new


def recursiveCanalization(G_new, cS):
    # delete complement variables
    
    G_eliminated = eliminateComplementS(G_new, cS)
    # print(G_eliminated)
    canalizedV = [x.split(' = ')[0] for x in G_eliminated.splitlines() if x.split(' = ')[1].strip() == '']
    i = 1
    while 1:
        G_eliminated = eliminateComplementS(G_eliminated, canalizedV)
        # print('r',G_eliminated)
        canalizedV_new = [x.split(' = ')[0] for x in G_eliminated.splitlines() if x.split(' = ')[1].strip() == '']
        if set(canalizedV)==set(canalizedV_new):
            break
        else:
            i += 1
            canalizedV = canalizedV_new

        print('recursiveCanalization:',i)

    return G_eliminated

# integrity condition
def checkIN(nodeNum, desiredVariable, S):
    didxList = []
    for didx, datt in enumerate(desiredVariable):
        if set(S)-set(datt) == set(): # point attractor
            didxList.append(didx)

        else:
            continue            
    return didxList
            

import networkx as nx
from collections import defaultdict

def reducedwithS(net, S, depth = 5):
    netline = makeReducedForm(net).strip()
    net_eq = {x.split(' = ')[0]:x.split(' = ')[1] for x in netline.split('\n')}

    for _ in range(depth):
        net_eq_new = copy.deepcopy(net_eq)
        deleteVariable = set(net_eq.keys())-set(S)

        for k, vline in net_eq.items():    
            if k not in deleteVariable:
                for x in deleteVariable:
                    vline = vline.replace(x, net_eq[x])
                    net_eq_new[k] = vline
                    net_eq_new[x] = ''
                
        if np.any([x in ''.join(net_eq_new.values()) for x in deleteVariable]):
            continue
        else:
            break
    new_net_line = '\n'.join([k+' = ' +v for k,v in net_eq_new.items() if v != ''])
    return new_net_line

def defineUD(desiredSet, undesiredSet):
    D_state = desiredSet.T.loc[(desiredSet.sum(axis=0)/(desiredSet.shape[0]) == 0) | (desiredSet.sum(axis=0)/(desiredSet.shape[0]) == 1),0].to_dict()
    # U_state = undesiredSet.T.loc[(undesiredSet.sum(axis=0)/(undesiredSet.shape[0]) == 0) | (undesiredSet.sum(axis=0)/(undesiredSet.shape[0]) == 1),0].to_dict()
    # D_list = []
    # for k,v in D_state.items():
    #     if (k<=8) & (v==0): D_list.append('~x0'+str(k+1))
    #     elif (k<=8) & (v==1): D_list.append('x0'+str(k+1))
    #     elif (k>8) & (v==0): D_list.append('~x'+str(k+1))
    #     else: D_list.append('x'+str(k+1))
    D_list = []
    for k,v in D_state.items():
        if (v==0) : D_list.append('~'+k)
        else: D_list.append(k)
            
    U_lists = [] # individual 
    for u_idx in range(undesiredSet.shape[0]): 
        u_states = {x:int(y) for x,y in undesiredSet.iloc[u_idx,:].to_dict().items() if y != '*'}
        U_list = []
        for k,v in u_states.items():
            if (v==0) : U_list.append('~'+k)
            else: U_list.append(k)
            # if (k<=8) & (v==0): U_list.append('~x0'+str(k+1))
            # elif (k<=8) & (v==1): U_list.append('x0'+str(k+1))
            # elif (k>8) & (v==0): U_list.append('~x'+str(k+1))
            # else: U_list.append('x'+str(k+1)) 
        U_lists.append(U_list)

            
    return D_list, U_lists

def defineUD_individual(desiredSet, undesiredSet):
    #D_state = desiredSet.T.loc[(desiredSet.sum(axis=0)/(desiredSet.shape[0]) == 0) | (desiredSet.sum(axis=0)/(desiredSet.shape[0]) == 1),0].to_dict()
    #U_state = undesiredSet.T.loc[(undesiredSet.sum(axis=0)/(undesiredSet.shape[0]) == 0) | (undesiredSet.sum(axis=0)/(undesiredSet.shape[0]) == 1),0].to_dict()
    D_lists = [] # individual 
    for d_idx in range(desiredSet.shape[0]): 
        d_states = {x:int(y) for x,y in desiredSet.iloc[d_idx,:].to_dict().items() if y != '*'}
        D_list = []
        for k,v in d_states.items():
            if (v==0) : D_list.append('~'+k)
            else: D_list.append(k)
            # if (k<=8) & (v==0): D_list.append('~x0'+str(k+1))
            # elif (k<=8) & (v==1): D_list.append('x0'+str(k+1))
            # elif (k>8) & (v==0): D_list.append('~x'+str(k+1))
            # else: D_list.append('x'+str(k+1)) 
        D_lists.append(D_list)
    D_list = list(set(itertools.chain(*D_lists)))
    
    U_lists = [] # individual 
    for u_idx in range(undesiredSet.shape[0]): 
        u_states = {x:int(y) for x,y in undesiredSet.iloc[u_idx,:].to_dict().items() if y != '*'}
        U_list = []
        for k,v in u_states.items():
            if (v==0) : U_list.append('~'+k)
            else: U_list.append(k)
            # if (k<=8) & (v==0): U_list.append('~x0'+str(k+1))
            # elif (k<=8) & (v==1): U_list.append('x0'+str(k+1))
            # elif (k>8) & (v==0): U_list.append('~x'+str(k+1))
            # else: U_list.append('x'+str(k+1)) 
        U_lists.append(U_list)

            
    return D_list, U_lists


def makeDict(x):
    perturb = defaultdict()
    for p in x:
        if '~' not in p: perturb[p] = True
        else: perturb[p[1:]] = False
    return perturb


def transformedTarget(addTarget, V_pair, desiredVariable):
    addTarget2 = []
    for addt in addTarget:
        # changed the variables corresponding to the desired values
        for desiredV in desiredVariable:
            addTarget2.append([list(set([x for x in V_pair if add_ in x][0]) & set(desiredV))[0] for add_ in addt if len(set([x for x in V_pair if add_ in x][0]) & set(desiredV)) > 0 ])    
    addTarget2 = list(set([tuple(x) for x in addTarget2 if len(x) == len(addt)]))

    
    return addTarget2



