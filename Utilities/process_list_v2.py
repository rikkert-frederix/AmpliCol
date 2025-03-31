#!/usr/bin/env python

import itertools
import copy
import re
import argparse
from collections import Counter

def ValidColorOrd(proc,perm):
    # Check if 'perm' is a valid color order for the process 'proc'
    found=[False,False,False,False] # quarks, antiquarks,singlets,gluons
    for i in range(len(perm)):
        if proc[perm[i]] in quarks:
            if found[0] or found[3]:
                validColorOrd=False
                return
            if found[1] and perm[i] < perm[0]: # for 2-quark line process, the removes one of the two orders
                validColorOrd=False
                return
            found[0]=True
            found[1]=False
            found[2]=False
            found[3]=False
        if proc[perm[i]] in antiquarks:
            if found[1] or found[2] or not(found[0]) :
                validColorOrd=False
                return
            found[0]=False
            found[1]=True
            found[2]=False
            found[3]=False
        if proc[perm[i]] in singlets:
            if found[3] or found[0] or not(found[1]):
                validColorOrd=False
                return
            found[2]=True
        if proc[perm[i]] in gluons:
            if found[1] or found[2]:
                validColorOrd=False
                return
            if i == 0 and perm[i] != 0: # for all-gluon procs, this removes cyclic permutations
                validColorOrd=False
                return
            found[3]=True
    return True

def UniqueColorOrd(proc,perm):
    # Check if 'perm' is a color order in canonical order for the process 'proc'.

    # Start at the first incoming particle. From there, the indentical
    # final state particles should each come in *increasing* order.
    zero=perm.index(0)
    perm_mapped=perm[zero:]+perm[:zero]
    for particle in all_coloured:
        particle_positions=[]
        for i,part in enumerate(proc[2:]):
            if part==particle:
                particle_positions.append(i+2)
        previous_position=0        
        for i in particle_positions:
            if perm_mapped.index(i) < previous_position:
                return False
            else:
                previous_position=perm_mapped.index(i)
    return True

def OrderProcPerm(proc,perm):
    zero=perm.index(0)
    perm_mapped=list(perm[zero:]+perm[:zero])
    proc_mapped=proc[zero:]+proc[:zero]
    elements_to_order=[]
    for i in perm_mapped:
        if proc[i] in massless_QCD and i > 1:
            elements_to_order.append(i)
    indices=[perm_mapped.index(x) for x in elements_to_order]
    sorted_elements=sorted(elements_to_order)
    for i, val in zip(indices,sorted_elements):
        perm_mapped[i]=val
    perm_ordered=tuple(perm_mapped[len(perm)-zero:]+perm_mapped[:len(perm)-zero])
    
    proc_ordered=[None]*len(proc)
    for i in range(len(perm_ordered)):
        proc_ordered[perm_ordered[i]]=proc[perm[i]]

    return tuple(proc_ordered),perm_ordered

def AddProcPermToPhaseSpaceOrder(proc,perm,phase_space_orders):
    zero=perm.index(0)
    perm_mapped=tuple(perm[zero:]+perm[:zero])
    if perm_mapped in phase_space_orders:
        phase_space_orders[perm_mapped].append((proc,perm))
    else:
        phase_space_orders[perm_mapped]=[(proc,perm)]


def ParseCollision(input_string):
    parts=input_string.split(">")
    if len(parts) != 2:
        raise ValueError("Invalid collision format. Expected 'p p > ...'.")
    initial_state=parts[0].strip().split()
    final_state=parts[1].strip().split()
    jet_match=re.match(r"(\d+)j",final_state[-1]) if final_state else None
    jet_count=int(jet_match.group(1)) if jet_match else 0
    rest=final_state[:-1] if jet_match else final_state
    return {"initial_state":initial_state,"jet_count":jet_count,"rest":rest}

def count_matching_elements(main_list,check_list):
    main_counts=Counter(main_list)
    count=sum(main_counts[item] for item in check_list if item in main_counts)
    return count

def ValidProc(proc):
    nq=count_matching_elements(proc,quarks)
    naq=count_matching_elements(proc,antiquarks)
    if nq > 2 : return False
    if naq > 2 : return False
    if nq != naq : return False
    # check flavour changing currents:
    for q in quarks:
        if count_matching_elements(proc,[q]) != count_matching_elements(proc,[q+'bar']) : return False
    # need at least one quark line if there are colour singlets
    if nq == 0 and count_matching_elements(proc,singlets) > 0 : return False
    return True

def GenerateAllUniqueProcs(process):
    procs=[[]]
    for part in process['initial_state']:
        if part != 'p':
            raise ValueError("Initial state should be a proton ('p').")
    for part in range(process['jet_count']+2):
        procs_new=[]
        for proc in procs:
            for p in jet:
                if part > 3 and p not in gluons: continue
                if not procs_new:
                    procs_new=[sorted(proc+[p])]
                else:
                    procs_new.append(sorted(proc+[p]))
        procs=copy.deepcopy(procs_new)
    for part in process['rest']:
        for proc in procs:
            proc.append(part)
    unique_procs=[]
    for proc in procs:
        if (ValidProc(proc)):
            unique_procs.append(tuple(proc))
    return set(unique_procs)

def GenerateAllProcs(unique_procs):
    procs=set()
    for proc in unique_procs:
        for i,j in itertools.combinations(range(len(proc)),2):
            if proc[i] in proton and proc[j] in proton:
                pair1=[proc[i],proc[j]]
                pair2=[proc[j],proc[i]]
                remaining=[e for k,e in enumerate(proc) if k not in (i,j)]
                procs.add(tuple(pair1+remaining))
                procs.add(tuple(pair2+remaining))
    return procs

    

quarks=['d','u','s','c','b','t']
antiquarks=['dbar','ubar','sbar','cbar','bbar','tbar']
#flavour_scheme=['d','u','s','c','b'] # all the massless quarks
flavour_scheme=['d'] # all the massless quarks
singlets=['a','z','w+','w-','e+','e-','mu+','mu-','ta+','ta-','ve','ve~','vm','vm~','vt','vt~','h']
gluons=['g']
all_coloured=quarks+antiquarks+gluons
massless_QCD=flavour_scheme+[q+'bar' for q in flavour_scheme]+gluons
proton=massless_QCD
jet=massless_QCD

if jet != proton:
    raise ValueError("definition of 'jet' and 'proton' should be the same")

parser=argparse.ArgumentParser(description="Generate the full list of processes, ordered by phase-space order")
parser.add_argument("process_string",type=str,help="Process to consider (e.g., 'p p > w+ z 4j')")
args=parser.parse_args()
                            
process=ParseCollision(args.process_string)

all_unique_procs=GenerateAllUniqueProcs(process)
all_procs=GenerateAllProcs(all_unique_procs)

#print(len(all_unique_procs),all_unique_procs)
#print(len(all_procs),all_procs)

phase_space_orders={}

for proc in all_procs:
    all_possible_color_ord=list(itertools.permutations([i for i in range(len(proc))]))
    valid_color_ord=[]
    for perm in all_possible_color_ord:
        if ValidColorOrd(proc,perm):
            valid_color_ord.append(perm)
    unique_color_ord=[]
    for perm in valid_color_ord:
        if UniqueColorOrd(proc,perm):
            unique_color_ord.append(perm)
    for perm in unique_color_ord:
        ordered_proc,ordered_perm=OrderProcPerm(proc,perm)
        AddProcPermToPhaseSpaceOrder(ordered_proc,ordered_perm,phase_space_orders)
        
print(len(phase_space_orders.keys()))
for key in phase_space_orders.keys():
    print(key,':',len(phase_space_orders[key]))
#    print(key,':',phase_space_orders[key])
    
