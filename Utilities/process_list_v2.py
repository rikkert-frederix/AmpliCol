#!/usr/bin/env python

import itertools
import copy
import re
import argparse
from collections import Counter
import multiprocessing
import math

# Global sets (make then 'frozenset' so that they are immutable):
quarks=frozenset({'d','u','s','c','b','t'})
antiquarks=frozenset({'dbar','ubar','sbar','cbar','bbar','tbar'})
singlets=frozenset({'a','z','w+','w-','e+','e-','mu+','mu-','ta+','ta-','ve','ve~','vm','vm~','vt','vt~','h'})
gluons=frozenset({'g'})
#flavour_scheme=frozenset({'d','u','s','c','b'}) # all the massless quarks
flavour_scheme=frozenset({'d','u'}) # all the massless quarks
all_coloured=quarks | antiquarks | gluons
massless_QCD=flavour_scheme | frozenset([q+'bar' for q in flavour_scheme]) | gluons
proton=massless_QCD
jet=massless_QCD
if jet != proton:
    raise ValueError("definition of 'jet' and 'proton' should be the same")
pdgs={'g':'21','d':'1','u':'2','s':'3','c':'4','b':'5','t':'6','dbar':'-1','ubar':'-2','sbar':'-3','cbar':'-4','bbar':'-5','tbar':'-6','a':'22','z':'23','w+':'24','w-':'-24','e+':'-11','e-':'11','mu+':'-13','mu-':'13','ta+':'-15','ta-':'15','ve':'12','ve~':'-12','vm':'14','vm~':'-14','vt':'16','vt~':'-16','h':'25'}
anti_particle={'g':'g','d':'dbar','u':'ubar','s':'sbar','c':'cbar','b':'bbar','t':'tbar','dbar':'d','ubar':'u','sbar':'s','cbar':'c','bbar':'b','tbar':'t','a':'a','z':'z','w+':'w-','w-':'w+','e+':'e-','e-':'e+','mu+':'mu-','mu-':'mu+','ta+':'ta-','ta-':'ta+','ve':'ve~','ve~':'ve','vm':'vm~','vm~':'vm','vt':'vt~','vt~':'vt','h':'h'}


def ProcessProcess(proc):
    """Function to process each 'proc' in parallel"""
    phase_space_orders_local = {}  # Local dictionary to avoid race conditions

    all_possible_color_ord = list(itertools.permutations(range(len(proc))))
    valid_color_ord = [perm for perm in all_possible_color_ord if ValidColorOrd(proc, perm)]
    unique_color_ord = [perm for perm in valid_color_ord if UniqueColorOrd(proc, perm)]

    for perm in unique_color_ord:
        ordered_proc, ordered_perm = OrderProcPerm(proc, perm)
        zero = ordered_perm.index(0)
        perm_mapped = tuple(ordered_perm[zero:] + ordered_perm[:zero])

        if perm_mapped in phase_space_orders_local:
            phase_space_orders_local[perm_mapped].append((ordered_proc, ordered_perm, []))
        else:
            phase_space_orders_local[perm_mapped] = [(ordered_proc, ordered_perm, [])]

    return phase_space_orders_local


def ValidColorOrd(proc,perm):
    # Check if 'perm' is a valid color order for the process 'proc'
    found_quark = found_antiquark = found_singlet = found_gluon = False
    for idx in perm:
        particle = proc[idx]
        if particle in quarks:
            if found_quark or found_gluon:
                return False
            if found_antiquark and idx < perm[0]:  # Two-quark line process, reject one ordering
                return False
            found_quark = True
            found_antiquark = found_singlet = found_gluon = False
        elif particle in antiquarks:
            if found_antiquark or found_singlet or not found_quark:
                return False
            found_antiquark = True
            found_quark = found_singlet = found_gluon = False
        elif particle in gluons:
            if found_antiquark or found_singlet:
                return False
            found_gluon = True
        else:  # Assuming the rest are singlets
            if found_quark or found_gluon or not found_antiquark :
                return False
            found_singlet = True
    if found_gluon: # Only for all-gluon process found_gluon is True here. Use it to remove cyclic permutations
        if perm[0] != 0 :
            return False
    return True

def UniqueColorOrd(proc,perm):
    # Check if 'perm' is a color order in canonical order for the process 'proc'.
    #
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
    if nq > 2 : return False    # at most two quarks
    if naq > 2 : return False   # at most two anti-quarks
    if nq != naq : return False # same number of quarks and anti-quarks
    # remove flavour changing currents:
    for q in quarks:
        if count_matching_elements(proc,[q]) != count_matching_elements(proc,[q+'bar']) : return False
    # need at least one quark line if there are colour singlets:
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
        procs=procs_new.copy()
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

def CombineResults(results):
    phase_space_orders = {}
    for result in results:
        for key, value in result.items():
            if key in phase_space_orders:
                phase_space_orders[key].extend(value)
            else:
                phase_space_orders[key] = value
    return phase_space_orders
                
def ParseArgument():
    parser=argparse.ArgumentParser(description="Generate the full list of processes, ordered by phase-space order")
    parser.add_argument("process_string",type=str,help="Process to consider (e.g., 'p p > w+ z 4j')")
    args=parser.parse_args()
    return ParseCollision(args.process_string)

def IdenticalParticleSymmetryFactor(proc):
    i_fac=1
    for p in all_coloured:
        i_fac*=max(1,math.factorial(proc[2:].count(p)))
    return i_fac


def MultiChannelPartners(proc,perm,k,l):
    all_possible_perms={perm}
    colour_singlets_in_proc=[perm.index(i) for i,p in enumerate(proc) if p in singlets]
    anti_quarks_in_proc=tuple([perm.index(i) for i,p in enumerate(proc) if p in antiquarks])
    if colour_singlets_in_proc:
        all_singlet_orders=tuple(itertools.permutations(colour_singlets_in_proc))
    else:
        all_singlet_orders=()
    if len(anti_quarks_in_proc) == 1:
        for s in all_singlet_orders:
            order=[]
            for i in range(len(perm)):
                if i == anti_quarks_in_proc[0]:
                    order.extend([perm[p] for p in list(anti_quarks_in_proc+s)])
                elif i not in colour_singlets_in_proc:
                    order.append(perm[i])
            all_possible_perms.add(tuple(order))
    elif len(anti_quarks_in_proc) == 2:
        for j in range(len(all_singlet_orders)+1):
            for s in all_singlet_orders:
                order=[]
                for i in range(len(perm)):
                    if i == anti_quarks_in_proc[0]:
                        order.extend([perm[p] for p in list((anti_quarks_in_proc[0],)+s[:j])])
                    elif i == anti_quarks_in_proc[1]:
                        order.extend([perm[p] for p in list((anti_quarks_in_proc[1],)+s[j:])])
                    elif i not in colour_singlets_in_proc:
                        order.append(perm[i])
                all_possible_perms.add(tuple(order))
    mt=[]
    for o in all_possible_perms:
        found=False
        for i,key in enumerate(all_keys_sorted):
            for (process,order,multichannel) in phase_space_orders[key]:
                if process == proc and order==o:
                    if found:
                        print('FOUND DOUBLE')
                    else:
                        Found=True
                        mt.append(i)
#                        print('found',process,proc,order,o,i)
        if not Found:
            print('NOT FOUND')
    phase_space_orders[k][l]=(proc,perm,tuple(sorted(mt)))
#    print(proc,':',all_singlet_orders,':',anti_quarks_in_proc,':',all_possible_perms)

def DetermineMultiChannelPartnersAndSymmetryFactor():
    for j,key in enumerate(all_keys_sorted):
        for i,(process,order,multichannel) in enumerate(phase_space_orders[key]):
            MultiChannelPartners(process,order,key,i)
    for key in all_keys_sorted:
        for i,(process,order,multichannel) in enumerate(phase_space_orders[key]):
            phase_space_orders[key][i]=(process,order,multichannel,IdenticalParticleSymmetryFactor(process))


if __name__ == "__main__":    
    process=ParseArgument()
    all_unique_procs=GenerateAllUniqueProcs(process)
    all_procs=GenerateAllProcs(all_unique_procs)
    
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(ProcessProcess, all_procs)  # Parallelize across procs
    phase_space_orders=CombineResults(results)
    all_keys_sorted=sorted(phase_space_orders.keys())

    DetermineMultiChannelPartnersAndSymmetryFactor()

    for key in all_keys_sorted:
#        print(key,':',len(phase_space_orders[key]))
        print(key,':',phase_space_orders[key])
        
