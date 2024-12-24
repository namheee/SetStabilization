from choonFunction import mFVSs, modeltext_transform, canalizing
from setFunction import *

from scipy.spatial import distance
import numpy as np
from collections import Counter
from pyboolnet.digraphs import _primes2signed_digraph
from sklearn.cluster import KMeans
import itertools
import copy
import pickle



def model2attBSIM(model_file, simulationPickle, n_clusters=2, pThres=2):
    primes = bnet2primes(model_file)
    with open(simulationPickle,'rb') as f:
        nodeList_ = pickle.load(f)
        str_steady = pickle.load(f)
    
    # simulation matching
    str_steady = set(str_steady)
    print('Matched?',list(primes.keys()) == nodeList_)
    if len(str_steady)<2: return [], [], [], [], [], []
    
    # compute similarity using hamming distance (K-means clustering)
    similarity = np.array([[1-distance.hamming(list(w1),list(w2)) for w1 in str_steady] for w2 in str_steady])
    kmeans = KMeans(n_clusters=n_clusters, random_state=123, n_init="auto").fit(similarity)
    
    print(model_file)
    print('Cluster#', len(np.unique(kmeans.labels_)), kmeans.labels_)
    labels = dict(Counter(kmeans.labels_))
    
    max_labels = [k for k,v in labels.items() if v == max(labels.values())]

    desired_ = [x for i,x in enumerate(str_steady) if i in np.where(kmeans.labels_ == max_labels[0])[0]]
    undesired_ = list(set(str_steady) - set(desired_))

    desiredSet = pd.DataFrame([list(x) for x in desired_], columns = list(primes.keys()))
    desiredSet = desiredSet.loc[:,~np.any(desiredSet == '*', axis=0)].astype('int')
    undesiredSet = pd.DataFrame([list(x) for x in undesired_], columns = list(primes.keys()))
 
    D_list, U_list = defineUD(desiredSet, undesiredSet)
    desired = list(set(desired_))
    undesired = list(set(str_steady) - set(desired))

    if pThres == 100 : pThres = len(U_list)*2/3
    # if a set of undesired attractor has three attractors, it means that two out of three choose different one from desired one
    else: pThres = pThres
    phenotypeNode = [k for k,v in Counter(itertools.chain(*[set(D_list)-set(u_list) for u_list in U_list])).items() if v >= pThres][:2]
    return str_steady, D_list, U_list, desired, undesired, phenotypeNode
    
def makefTarget_updated(fTarget, model_file, V_idx, V_pair, desired, undesired, desiredVariable, undesiredVariable, S2, remainedThres = 5):
    renamed = {v:k for k,v in V_idx.items()}
    fTarget2 = []
    remained = defaultdict(list)
    for candTarget in fTarget:
        modeltext = readBNETfromPrimes(model_file)
        F_net = modeltext_transform(modeltext)
    
        F_net_remained = canalizing(F_net, makeDict([renamed[x]+'_' for x in candTarget]))
        # G_ = _primes2signed_digraph(primes)
        # scc_num = len(list(nx.strongly_connected_components(G_)))
        
        if len(F_net_remained.split('=')) <= remainedThres : # applying FVS to remained network
            # print('F_net_remained:','\n', F_net_remained)
            addTargeto = mFVSs(F_net_remained)
            addTargetd = transformedTarget(addTargeto, V_pair, desiredVariable)

            # check if controlloing FVS nodes avoid any undesired state
            addTarget = []
            for add in addTargetd:
                uidx = checkIN(len(V_pair), undesiredVariable , (S2 | set(add)))
                didx = checkIN(len(V_pair), desiredVariable, (S2 | set(add)))
                isNpoint = np.all(['*' in x for i,x in enumerate(desired) if i in didx]) # is any point att?
                if (len(uidx) == 0 ) & (~isNpoint):
                    addTarget.append(add)
            
         
            if (len(addTargeto)>0) & (len(addTarget)==0):
                addTarget = [[]] # at least one value of addTargets is * ,i.e., there is no correct answer to drvie desired state
                fTarget2 += []
            else:
                fTarget2 += [list(set(itertools.chain(*x))) for x in (itertools.product(addTarget, [candTarget]))]
                
        else: # iteratively apply the proposed  algorithm
            
            remainedName = (model_file[:-5] + ':' + ':'.join(candTarget)+'.bnet').replace('/algorithm/network/','/algorithm/tmp/')
            F_net_remained_ = re.sub(" = ", ", ", F_net_remained)
            F_net_remained_ = re.sub("~", "!", F_net_remained_)
            F_net_remained_ = re.sub("_", "", F_net_remained_)
            with open(remainedName,'w') as f: f.write(F_net_remained_)
                
            primes = bnet2primes(remainedName)
            newNode = [x+'_' for x in primes.keys()]
            
            newD = []
            for atts in desired:
                newd_ = ''.join([y for i,y in enumerate(list(atts)) if i in np.where([x[0] in newNode for x in V_pair])[0]])
                newD.append(newd_)
    
            newU = []
            candTargetidx = np.where([len(set(x)&set(candTarget)) == 1 for x in V_pair])[0]
            controlled = makeDict([renamed[x]+'_' for x in candTarget])
            for atts in undesired:
                check1 = [y != '*' for i,y in enumerate(list(atts)) if i in candTargetidx] # if fixed nodes are '*', this att. cannot be arrived.
                check2 = []
                for ki,vi in makeDict([renamed[x]+'_' for x in candTarget]).items():
                    kidx = (np.where([len(set(x)&set([ki])) == 1 for x in V_pair])[0])
                    check2.append(atts[kidx[0]] == str(int(vi)))
                if np.all(check1) & np.all(check2):
                    newu_ = ''.join([y for i,y in enumerate(list(atts)) if i in np.where([x[0] in newNode for x in V_pair])[0]])
                    newU.append(newu_)
    
            desiredSet = pd.DataFrame([list(x) for x in newD], columns = list(primes.keys()))
            desiredSet = desiredSet.loc[:,~np.any(desiredSet == '*', axis=0)].astype('int')
            undesiredSet = pd.DataFrame([list(x) for x in newU], columns = list(primes.keys()))
            D_list, U_list = defineUD_individual(desiredSet, undesiredSet)
            # D_list, U_list = defineUD(desiredSet, undesiredSet)
            
            pThres = 1 
            newP = [k for k,v in Counter(itertools.chain(*[set(D_list)-set(u_list) for u_list in U_list])).items() if v >= pThres][:1]
    
            # print(newD, newU)
            remained[':'.join(candTarget)].append((remainedName, newNode, set(newD), set(newU), newP))
    
            
    return fTarget2, remained

def main4_updated(model_file, desired, undesired, phenotypeNode, oriN=0, remainedThres=5):

    if oriN == 0 : recurL = 0;
    else: recurL = oriN
    
    G_net, V_pair, V_idx, desiredSet, undesiredSet, phenotype, cphenotype = returnINPUT(model_file, desired, undesired, phenotypeNode, recurL)
    nodeNum = len(V_pair)
    desiredVariable, cyclicnode = desiredSet
    undesiredVariable, undesired_cyclicnode = undesiredSet
    
    pp_seq, ppCombiList_updated = ppCombination(G_net, V_pair, phenotype, cphenotype)
    
    
    allresult = defaultdict(dict)
    allremained = defaultdict(dict)
    for cidx, ppcombi in enumerate(ppCombiList_updated):
    
        # initial step
        S0 = set(phenotype)
        S1 = set([y.strip() for y in itertools.chain(*[x.split('&') for x in ppcombi])])
        S = S0 | S1
        pairedS = set(itertools.chain(*[x for x in V_pair if (x[0] in S) or (x[1] in S) ]))
        cS = pairedS-S
    
        G_new = makeReducedG(G_net, pp_seq, ppcombi)
        desiredIdx = checkIN(nodeNum, desiredVariable, S) # make R
    
        if len(desiredIdx) > 0: # start
            # print('n',G_new)
            G_eliminated = recursiveCanalization(G_new, cS)
            # print('e',G_eliminated)
            result = defaultdict()
            if (len(S) == nodeNum):
                desiredIdx = checkIN(nodeNum, desiredVariable, S)
                undesiredIdx = checkIN(nodeNum, undesiredVariable, S)
                # print('END1', desiredIdx, undesiredIdx)
                controlTarget = set(S) - set([x.split(' = ')[0] for x in makeReducedForm(G_eliminated).splitlines()])
                setTarget = mFVSs(reducedwithS(G_eliminated, S))
                fTarget = [list(set(itertools.chain(*x))) for x in (itertools.product(setTarget, list([controlTarget])))]
                result['pp2_*'] = (0, desiredIdx, undesiredIdx, controlTarget, fTarget, tuple(S))

    
            else:
                def updatedS_G(G_eliminated, V_pair, S, cS):
                    pp_seqS, ppCombiS_updated = ppCombination(makeReducedForm(G_eliminated), V_pair, list(S), cS)
    
                    # ============================================================================================
                    S_updated = []
                    for ppidx, ppC in enumerate(ppCombiS_updated):
                        S_new = set(S) | set(itertools.chain(*[[y.strip() for y in x.split('&')] for x in ppC]))
                        lenR = len(checkIN(nodeNum, desiredVariable, S_new))
                        S_updated.append([ppidx, lenR, len(S_new)])
    
                    df = pd.DataFrame(S_updated, columns = ['ppidx','lenR','lenS'])
                    df1 = df.loc[(df.lenR == max(df.lenR)),]
                    df2 = df1.loc[(df1.lenS == min(df1.lenS)),]
                    # =============================================================================================
    
                    G_o = copy.deepcopy(G_eliminated)
                    S_o = copy.deepcopy(S)
                    pp_combIdx = df2.ppidx.tolist()
                    return G_o, S_o, pp_combIdx, pp_seqS, ppCombiS_updated
    
    
    
                def iteratewithS_updated(ppidx, ppCombiS_updated, pp_seq, V_pair, G_eliminated, S):
                    G_new = makeReducedG(makeReducedForm(G_eliminated), pp_seq, ppCombiS_updated[ppidx]) 
    
                    S_new = set(S) | set(itertools.chain(*[[y.strip() for y in x.split('&')] for x in ppCombiS_updated[ppidx]]))
                    pairedS = set(itertools.chain(*[x for x in V_pair if (x[0] in S_new) or (x[1] in S_new) ]))
                    cS_new = pairedS-S_new
    
                    G_eliminated_new = recursiveCanalization(G_new, cS_new)
    
                    return S_new, G_eliminated_new            
    
                # select the largest lenR & the smallest lenS
                G_o, S_o, pp_combIdx, pp_seqS, ppCombiS_updated = updatedS_G(G_eliminated, V_pair, S, cS)
    
                for usidx in pp_combIdx:
                    j = 1; S = S_o; G_eliminated = G_o
                    check_idx = usidx
                    ppCombiS_updated2, pp_seqS2 = ppCombiS_updated, pp_seqS
                    print('S',S)
                    while 1:
                        S2, G_eliminated2 = iteratewithS_updated(check_idx, ppCombiS_updated2, pp_seqS2, V_pair, G_eliminated, S)
                        desiredIdx = checkIN(nodeNum, desiredVariable, S2)
                        undesiredIdx = checkIN(nodeNum, undesiredVariable, S2)   
                        if len(desiredIdx) == 0 : break
                        print(desiredIdx, undesiredIdx)
                        print('S',S2)
                        if (S == S2) | (len(S2) == nodeNum):        
                            print('END2', desiredIdx, undesiredIdx)
                            if len(S2) != nodeNum: j -= 1
    
                            controlTarget = set(S2) - set([x.split(' = ')[0] for x in makeReducedForm(G_eliminated2).splitlines()])                        
                            setTarget = mFVSs(reducedwithS(G_eliminated2, S2))     
                            fTarget = [list(set(itertools.chain(*x))) for x in (itertools.product(setTarget, list([controlTarget])))]                           
                            
                            if len(undesiredIdx) != 0:
                                while 1:
                                    
                                    fTarget2, remained = makefTarget_updated(fTarget, model_file, V_idx, V_pair, 
                                                                             desired, undesired, desiredVariable, 
                                                                             undesiredVariable, S2, remainedThres)
                                
                                    if len(remained) > 0 : allremained.update(remained); break
                                    elif len(fTarget2) == 0 : remainedThres = 1
                                    else: 
                                        result['pp3_'+str(usidx)] = (j, desiredIdx, undesiredIdx, controlTarget, fTarget2, tuple(S2))
                                        break
                                
                            else: result['pp2_'+str(usidx)] = (j, desiredIdx, undesiredIdx, controlTarget, fTarget, tuple(S2))
                            break
                        else: # S update
                            pairedS2 = set(itertools.chain(*[x for x in V_pair if (x[0] in S2) or (x[1] in S2) ]))
                            cS2 = pairedS2-S2
                            G_eliminated, S, pp_combIdx2, pp_seqS2, ppCombiS_updated2 = updatedS_G(G_eliminated2, V_pair, S2, cS2)
                            check_idx = pp_combIdx2[0] # randomly select one pp (it does not guarantee the least)
                            j += 1
                    print(j)
            allresult['pp1_'+str(cidx)] = result
        else: continue
    return V_idx, allresult, allremained

def dict2df(renamed, allresult, plusCtrl):
    rList = []
    for x,y in allresult.items():
        for value in list(y.values()):
            solution = [tuple([renamed[y] for y in x]) for x in value[4]]
            solution = tuple(itertools.product(solution, [[renamed[y] for y in plusCtrl]]))
            solution = [tuple(set(itertools.chain(*x))) for x in solution]
            rList.append((len(value[1]), len(value[2]), len(value[3]), tuple(solution) , tuple(plusCtrl)))
    rdf = pd.DataFrame(rList, columns = ['num_D','num_U','ctrl','solution','plusCtrl'])  
    return (rdf)