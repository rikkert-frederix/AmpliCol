#!/usr/bin/env python

# multi-jet


import itertools
import math
import copy



def create_all_procs_from_base_procs(proc):
    # given a flavour configuration, list all subprocesses that can be
    # obtained from that flavour configuration. Effectively, this
    # means looping over all possible pairs of initial state
    # particles.
    all_procs = []
    # Iterate over all permutations of the first two elements
    for perm in itertools.permutations(proc, 2):
        if perm[0] not in proton or perm[1] not in proton : continue
        # Find the remaining elements
        remaining = list(proc)
        remaining.remove(perm[0])
        remaining.remove(perm[1])
        
        # Sort remaining elements since final state is unordered
        remaining_sorted = sorted(remaining)
        
        # Combine the permuted first two elements with the sorted remaining
        combined = list(perm)+remaining_sorted

        if combined not in all_procs :
            all_procs.append(combined)
    return all_procs

def generate_unique_final_state_permutations(arr):
    all_procs=[]
    for perm in itertools.permutations(proc[2:]):
        combined=[proc[0]]+[proc[1]]+list(perm)
        if combined not in all_procs :
            all_procs.append(combined)
    return all_procs


def order_permutation_with_psorder(order,perm):
    return [perm[i-1] for i in order]

def valid_perm(perm):
    # check that all anti-quarks are followed by a quark
    n = len(perm)
    for i in range(n):
        if perm[i] in antiquarks :  # If the current element is anti-quark
            next_index = (i + 1) % n  # Cyclic next element
            if (perm[next_index] not in quarks ):  # next element must be quark
                return -1
    for i in range(n):
        if perm[i] in quarks:
            return i
    return 0

def shift_order(shift,order):
    order_shifted=order[shift:]+order[:shift]
    return order_shifted

    
def convert_to_pdg(perm):
    return " ".join([str(pdgs[antipart[perm[0]]]),str(pdgs[antipart[perm[1]]])]+[str(pdgs[i]) for i in perm[2:]])



nfinal=4
quarks=['d','u','s','c','b','t']
antiquarks=['dbar','ubar','sbar','cbar','bbar','tbar']
flavour_scheme=['d','u','s','c','b'] # all the massless quarks
#flavour_scheme=['d','u'] # all the massless quarks
proton=['g','d','u','s','c','b','dbar','ubar','sbar','cbar','bbar'] # all partons that can be element of proton
pdgs={'g':'21','d':'1','u':'2','s':'3','c':'4','b':'5','t':'6','dbar':'-1','ubar':'-2','sbar':'-3','cbar':'-4','bbar':'-5','tbar':'-6','a':'22'}
antipart={'g':'g','d':'dbar','u':'ubar','s':'sbar','c':'cbar','b':'bbar','t':'tbar','dbar':'d','ubar':'u','sbar':'s','cbar':'c','bbar':'b','tbar':'t','a':'a'}
all_part=['g','d','u','s','c','b','t','dbar','ubar','sbar','cbar','bbar','tbar','a']

# color-singlets
color_singlets=[]
#color_singlets=['a']
#color_singlets=['a','a']


# all-gluon process
#base_procs=[[]]
base_procs=[]
## one-quark-line process
#if nfinal+2 >= len(color_singlets)+2 : 
#    for q in flavour_scheme:
#        base_procs.append([q,q+'bar'])
# two-quark-line process
if nfinal+2 >= len(color_singlets)+4 : 
    # different-flavour
    for i,q in enumerate(flavour_scheme[:-1]):
        for qp in flavour_scheme[i+1:]:
            base_procs.append([q,qp,q+'bar',qp+'bar'])
    # same-flavour (put after different flavour)
    for q in flavour_scheme:
        base_procs.append([q,q,q+'bar',q+'bar'])

# add the gluons
for proc in base_procs:
    while len(proc) < nfinal+2 -len(color_singlets):
        proc.append('g')

unique_procs=copy.deepcopy(base_procs)
# add the 2qq_df processes with the two incoming particles interchanged
i=0
while i < len(unique_procs):
    proc = unique_procs[i]
    if proc[0] in quarks and proc[1] in quarks and proc[0] != proc[1]:
        swapped_proc=proc[:]
        swapped_proc[0],swapped_proc[1]=swapped_proc[1],swapped_proc[0]
        unique_procs.insert(i+1,swapped_proc)
        i+=1
    i+=1
        
    
# add the colour singlets
for proc in unique_procs:
    for s in color_singlets:
        proc.append(s)

        
all_procs=[]
for proc in base_procs:
    all_procs+=create_all_procs_from_base_procs(proc)

    
# Get symmetry factor
iden_fac={}
for i,proc in enumerate(all_procs):
    i_fac=1
    for p in all_part:
        i_fac*=max(1,math.factorial(proc[2:].count(p)))
    iden_fac[i]=i_fac

    
    
order=[1]+[i for i in range(3,nfinal+3)]
phase_space_order=[]
for i in range(1,nfinal+2):
    phase_space_order.append(order[:i]+[2]+order[i:])


towrite=[[] for _ in phase_space_order]
    
for pso,psorder in enumerate(phase_space_order):
    # The phase-space order fixes where particles 1 and 2 go in the
    # order. For each process, take all possible permutations of the
    # other particles (i.e., all the final state ones), and
    #
    # 1. Check that it is a unique configurations
    #
    # 2. Check that it is compatible with the phase-space-order (i.e.,
    # if there's a quark line, the gluons are between the quark and
    # anti-quark (up to cyclic permutations) etc.
    for i,proc in enumerate(all_procs):
        unique_permutations=generate_unique_final_state_permutations(proc)
        for perm in unique_permutations:
            perm_ordered=order_permutation_with_psorder(psorder,perm)
            shift=valid_perm(perm_ordered)
            if shift != -1:
                order_shifted=shift_order(shift,psorder)
                towrite[pso].append(convert_to_pdg(perm)+"   "+" ".join([str(i) for i in order_shifted])+"   "+str(iden_fac[i]))
#                towrite.append([perm,order_shifted,proc])
                

                




with open('processes.txt','w') as f:
    f.write(str(nfinal+2)+' '+str(len(unique_procs))+'\n')
    for proc in unique_procs:
        f.write(' '.join([pdgs[ele] for ele in proc])+'\n')
    f.write('\n')
    f.write('\n')
    f.write(str(len(towrite))+'\n')
    f.write('\n')
    for pso in towrite:
        f.write(str(len(pso))+'\n')
        for ele in pso:
            f.write(ele+'\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')


    
#    valid_perms=generate_permutations(proc)
#    for pso,psorder in enumerate(phase_space_order):
#        for valid_perm in valid_perms:
#            if valid_perm[0]==proc[0] and valid_perm[pso+1]==proc[1]:
#                to_write[pso].append([psorder,valid_perm,proc])
#print(to_write[0])
