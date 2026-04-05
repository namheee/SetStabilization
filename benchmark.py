import os
import pickle
import numpy as np
import itertools
import pandas as pd
import cabean
import time
import pystablemotifs as sm
import biobalm
import pystablemotifs.export as ex
import networkx as nx
from pyboolnet.file_exchange import bnet2primes
import pyboolnet
import itertools
import cana
import re
    
def pert_dict2str(pert_i, annotType = 'int'):
    if annotType == 'str':      
        return ('&'.join(sorted(['~'+x if y == 'false' else x for x,y in pert_i.items()])))
    else:
        return ('&'.join(sorted(['~'+x if y == 0 else x for x,y in pert_i.items()])))


def minBF_LDOI(primes, target_state, G_iter = 2000):
    # minimal (logical) drivers using brute-force search ============================
    start = time.time()
    min_bf = sm.drivers.minimal_drivers(target_state, primes)
    t_bf = time.time()-start

    # LDOI ===========================================================================
    start = time.time()
    # Search for drivers of target in primes using the method of Yang et al. 2018.
    target_ldoi = sm.drivers.GRASP(target_state, primes, GRASP_iterations = G_iter)   
    min_ = [len(x.keys()) for x in target_ldoi]
    target_ldoi = [x for x in target_ldoi if len(x.keys())==min(min_)]
    t_ldoi = time.time()-start 

    return (tuple([pert_dict2str(x) for x in min_bf]), t_bf), (tuple([pert_dict2str(x) for x in target_ldoi]), t_ldoi)

def SM(primes, target_state):
    # appproximated Stable motif ======================================================
    start = time.time()
    max_simulate_size = 20
    ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)    
    # Controlling Attractors
    sm_result = ar.succession_diagram.reprogram_to_trap_spaces(logically_fixed=target_state,
                                                   target_method='history',
                                                   driver_method='internal')
    
    t_sm = time.time()-start    
    return (tuple([pert_dict2str(x) for x in sm_result]), t_sm)


def balm(model_file, target_state):
    start = time.time()
    primes = bnet2primes(model_file)
    net = pyboolnet.file_exchange.primes2bnet(primes)
    sd = biobalm.SuccessionDiagram.from_rules(net)
    interventions = biobalm.control.succession_control(sd, target_state)
    try:
        intervention = interventions[0]
        ctrl_ = []
        for int1 in intervention.control:
            tmp1 = [pert_dict2str(x) for x in int1]
            min1 = np.min([len(x.keys()) for x in int1])
            tmp2 = [x for x in tmp1 if len(x.split('&'))==min1]
            ctrl_.append(tmp2)
        ctrl_targets = ['&'.join(x) for x in itertools.product(*ctrl_)]
        ftime = time.time() - start
    except IndexError:
        ctrl_targets = []
        ftime = 0.0
    return (tuple(ctrl_targets), ftime)


def mFVSs(modeltext):
    modeltext = modeltext.replace("=", "*=")
    net = cana.boolean_network.BooleanNetwork.from_string_boolean(modeltext)

    # Mapping nodes
    mappind_dic = {}
    for node in net.nodes:
        mappind_dic[node.id] = node.name

    # FVSs
    FVS_bruteforce = net.feedback_vertex_set_driver_nodes(graph='structural', method='bruteforce', max_search=10, keep_self_loops=True)  # brutuforce

    FVS_list_list = []
    for FVS in FVS_bruteforce:
        FVS_list = []
        for node in FVS:
            FVS_list.append(mappind_dic[node])
        FVS_list_list.append(FVS_list)

    return FVS_list_list

def compare_mFVS(model_file, target):
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

    modeltext = '\n'.join(net_ordered)
    nodeList = list(primes.keys())  
    # =================================
    # FVS
    # =================================
    start = time.time()
    fvsresult = tuple(mFVSs(modeltext))
    fvsresult = ['&'.join(['~'+x if target[nodeList.index(x)] == '0' else x for x in r_id]) for r_id in fvsresult]
    ftime = time.time() - start
    return (tuple(fvsresult), ftime)