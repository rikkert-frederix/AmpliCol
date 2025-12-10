#!/usr/bin/env python3.8

import itertools
import copy
import re
import argparse
from collections import Counter
import multiprocessing
import math



# Global sets (make then 'frozenset' so that they are immutable):
quarks=frozenset({'d','u','s','c','b','t'})
antiquarks=frozenset({'d~','u~','s~','c~','b~','t~'})
singlets=frozenset({'a','z','w+','w-','e+','e-','mu+','mu-','ta+','ta-','ve','ve~','vm','vm~','vt','vt~','h'})
gluons=frozenset({'g'})
all_coloured=quarks | antiquarks | gluons
flavour_scheme=frozenset({'d','u','s','c','b'}) # all the massless quarks
massless_QCD=flavour_scheme | frozenset([q+'~' for q in flavour_scheme]) | gluons
proton=massless_QCD
jet=massless_QCD
if jet != proton:
    raise ValueError("definition of 'jet' and 'proton' should be the same")
pdgs={'g':'21','d':'1','u':'2','s':'3','c':'4','b':'5','t':'6','d~':'-1','u~':'-2','s~':'-3','c~':'-4','b~':'-5','t~':'-6','a':'22','z':'23','w+':'24','w-':'-24','e+':'-11','e-':'11','mu+':'-13','mu-':'13','ta+':'-15','ta-':'15','ve':'12','ve~':'-12','vm':'14','vm~':'-14','vt':'16','vt~':'-16','h':'25'}
anti_particle={'g':'g','d':'d~','u':'u~','s':'s~','c':'c~','b':'b~','t':'t~','d~':'d','u~':'u','s~':'s','c~':'c','b~':'b','t~':'t','a':'a','z':'z','w+':'w-','w-':'w+','e+':'e-','e-':'e+','mu+':'mu-','mu-':'mu+','ta+':'ta-','ta-':'ta+','ve':'ve~','ve~':'ve','vm':'vm~','vm~':'vm','vt':'vt~','vt~':'vt','h':'h'}
#sort_particles={'g':0,'d':1,'u':2,'s':3,'c':4,'b':5,'t':6,'d~':7,'u~':8,'s~':9,'c~':10,'b~':11,'t~':12,'a':99,'z':99,'w+':99,'w-':99,'e+':99,'e-':99,'mu+':99,'mu-':99,'ta+':99,'ta-':99,'ve':99,'ve~':99,'vm':99,'vm~':99,'vt':99,'vt~':99,'h':99}
sort_particles={'g':13,'d':1,'u':2,'s':3,'c':4,'b':5,'t':6,'d~':7,'u~':8,'s~':9,'c~':10,'b~':11,'t~':12,'a':80,'z':81,'w+':82,'w-':83,'e+':84,'e-':85,'mu+':86,'mu-':87,'ta+':88,'ta-':89,'ve':90,'ve~':91,'vm':92,'vm~':93,'vt':94,'vt~':95,'h':96}
charges3={'g':0,'d':-1,'u':2,'s':-1,'c':2,'b':-1,'t':2,'d~':1,'u~':-2,'s~':1,'c~':-2,'b~':1,'t~':-2,'a':0,'z':0,'w+':3,'w-':-3,'e+':3,'e-':-3,'mu+':3,'mu-':-3,'ta+':3,'ta-':-3,'ve':0,'ve~':0,'vm':0,'vm~':0,'vt':0,'vt~':0,'h':0}
family={'g':0,'d':1,'u':1,'s':11,'c':11,'b':21,'t':21,'d~':-1,'u~':-1,'s~':-11,'c~':-11,'b~':-21,'t~':-21,'a':0,'z':0,'w+':0,'w-':0,'e+':-31,'e-':31,'mu+':-41,'mu-':41,'ta+':-51,'ta-':51,'ve':31,'ve~':-31,'vm':41,'vm~':-41,'vt':51,'vt~':-51,'h':0}

options = {}
process_order_to_index = {}

def SwitchFlavourScheme(FS):
    # Overwrite the relevant global variables so that we switch to the
    # 'FS' flavour-scheme process definition.
    global flavour_scheme
    global massless_QCD
    global proton
    global jet

    if FS==1:
        flavour_scheme=frozenset({'d'}) # all the massless quarks
    elif FS==2:
        flavour_scheme=frozenset({'d','u'}) # all the massless quarks
    elif FS==3:
        flavour_scheme=frozenset({'d','u','s'}) # all the massless quarks
    elif FS==4:
        flavour_scheme=frozenset({'d','u','s','c'}) # all the massless quarks
    elif FS==5:
        flavour_scheme=frozenset({'d','u','s','c','b'}) # all the massless quarks
    elif FS==6:
        flavour_scheme=frozenset({'d','u','s','c','b','t'}) # all the massless quarks
    else:
        print("ERROR: unknown flavour scheme",FS)
        quit()
    massless_QCD=flavour_scheme | frozenset([q+'~' for q in flavour_scheme]) | gluons
    proton=massless_QCD
    jet=massless_QCD

def ProcessProcess(proc):
    # For the process 'proc', consider all the possible colour
    # orderings compatible with that input process. Then, rearrange
    # them into a phase-space ordering dictionary.
    phase_space_orders_local = {}

    # All n! possible colour orderings:
    all_possible_color_ord = list(itertools.permutations(range(len(proc))))
    # Restrict them to the ones compatible with the process under consideration:
    valid_color_ord = [perm for perm in all_possible_color_ord if ValidColorOrd(proc, perm)]
    # In case we have identical final state particles, there can be
    # multiple colour orderings that are identical. Filter those out:
    unique_color_ord = [perm for perm in valid_color_ord if UniqueColorOrd(proc, perm)]

    # Group the possible colour orderings into their corresponding
    # phase-space orderings. The OrderProcPerm() function rearranges
    # the massless_QCD final state particles to reduce the number of
    # required phase-space orderings. 
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
    # Check if 'perm' is a valid color order for the process
    # 'proc'. Valid orders are of the form:
    #
    # [q,g,...,g,qbar,s,...,s,q,g,...g,qbar,s,...,s]
    #
    # Note that we keep all possible orderings of the colour singlets;
    # we will remove that overcounting by considering those processes
    # to be multi-channel partners. This allows us to map all the
    # peaks from e.g., close-to-collinear photons to any quark or
    # anti-quark.
    found_quark = found_antiquark = found_singlet = found_gluon = found_first=False
    for idx in perm:
        if idx == 0: found_first=True
        particle = proc[idx]
        if particle in quarks:
            # A quark should not come directly after another quark or
            # a gluon in the colour ordering:
            if found_quark or found_gluon:
                return False
            # In case of a multi-quark line process, we need to reject
            # some orderings---interchanging 'q,...,qbar' lines in the
            # ordering, results, in fact, not in a different
            # ordering---it does not matter which ones we remove, as
            # long as we keep one possible one. We take the one where
            # the colour-order is such that they are increasing in the
            # labels for the quarks (current idx is larger than the
            # idx from previous quark in the ordering).
            #
            # ALTERNATIVELY: we could keep all and take care of it
            # through multi-channeling (just as we do for the
            # colour-singlet orderings).
            # UPDATE: do this below when particle==antiquark
#            if found_antiquark and idx < quark_idx: 
#                return False
            found_quark = True
            found_antiquark = found_singlet = found_gluon = False
            quark_idx=idx
        elif particle in antiquarks:
            # An antiquark should not come directly after another
            # anti-quark or a colour singlet. Moreover, there should
            # have been already a quark earlier in the colour order:
            if found_antiquark or found_singlet or not found_quark:
                return False
            # UPDATE: Only remove duplicates from cyclic ordering. For
            # three quark lines this will give two times as many dual
            # amplitudes as we need; these will be taken care of
            # through multi-channeling
            if not found_first: return False
            found_antiquark = True
            found_quark = found_singlet = found_gluon = False
        elif particle in gluons:
            # Gluons should not come directly after antiquarks or
            # singlets:
            if found_antiquark or found_singlet:
                return False
            found_gluon = True
        else:  # Assuming the rest are singlets.  Singlets should come
               # directly after other singlets or antiquarks (i.e.,
               # not directly after quarks or gluons).
            if found_quark or found_gluon or not found_antiquark :
                return False
            found_singlet = True
    if found_gluon: # Only for all-gluon process found_gluon is True here. Use it to remove cyclic permutations
        if perm[0] != 0 :
            return False
    return True

def UniqueColorOrd(proc,perm):
    # Check if 'perm' is a color order in canonical order for the
    # process 'proc' in case there are identical final state
    # particles.
    #
    # Start at the first incoming particle. From there, the identical
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

def second_tuple_index(tup):
    return tup[1]

def OrderProcPerm(proc,perm):
    # In order to reduce the number of needed phase-space orderings,
    # as well as having the largest chance of finding identical matrix
    # elements, we bring the final state massless_QCD particles (and
    # corresponding colour orderings) into a canonical order.
    #
    # The canonical order is obtained by rearranging the massless_QCD
    # particles such that their corresponding colour ordering is
    # strictly increasing, when one starts counting as if the first
    # incoming particle would be at the first position in the order
    # (i.e., as if cyclicly permuting the order such that particle 1
    # is at position 1).
    #
    # Start by bringing the colour-order ('perm') to the canonical
    # ordering.
    zero=perm.index(0)
    perm_mapped=list(perm[zero:]+perm[:zero]) # cyclicly permute
    elements_to_order=[]
    for i in perm_mapped:
        if proc[i] in massless_QCD and i > 1:
            elements_to_order.append(i)
    indices=[perm_mapped.index(x) for x in elements_to_order]
    sorted_elements=sorted(elements_to_order)
    for i, val in zip(indices,sorted_elements):
        perm_mapped[i]=val
    perm_ordered=perm_mapped[len(perm)-zero:]+perm_mapped[:len(perm)-zero] # undo the cyclic permutation
    # Rearrange the process following the rearrangement of the colour ordering
    proc_ordered=[None]*len(proc)
    for i in range(len(perm_ordered)):
        proc_ordered[perm_ordered[i]]=proc[perm[i]]
    return tuple(proc_ordered),tuple(perm_ordered)


def ParseCollision(input_string):
    # Parse the process string, and cross the initial state particles.
    input_string=input_string.replace('bar','~')
    input_string=input_string.replace(' j',' 1j')
    parts=input_string.split(">")
    if len(parts) != 2:
        raise ValueError("Invalid collision format. Expected 'p p > ...'.")
    initial_state=parts[0].strip().split()
    for i,p in enumerate(initial_state):
        if p != 'p': 
            initial_state[i]=anti_particle[p]
    final_state=parts[1].strip().split()
    jet_match=re.match(r"(\d+)j",final_state[-1]) if final_state else None
    jet_count=int(jet_match.group(1)) if jet_match else 0
    rest=final_state[:-1] if jet_match else final_state
    lep=[]
    if (options['include_resonance']):
        i=0
        for k in rest:
            if (abs(int(pdgs[k])) >= 11 and abs(int(pdgs[k])) <= 16): 
                lep.append(k)
        for l in lep:
            rest.remove(l)
        if (sum(charges3[l] for l in lep) == 0): rest.append('z')
        if (sum(charges3[l] for l in lep) < 0): rest.append('w-')
        if (sum(charges3[l] for l in lep) > 0): rest.append('w+')
    return {"initial_state":initial_state,"jet_count":jet_count,"rest":rest,"lep":lep}

def count_matching_elements(main_list,check_list):
    # Counts how many elements of main_list are there in check_list
    # e.g., main_list=[a,b,b,b,c,c,d] ; check_list=[b,d,e] ; --> count_matching_elements=4
    main_counts=Counter(main_list)
    count=sum(main_counts[item] for item in check_list if item in main_counts)
    return count

def ValidProc(proc):
    # Checks if 'proc' is a valid process based on simple QFT rules
    # (such as charge conservation) and some restrictions imposed by
    # this code (e.g., maximum number of quark lines).
    nq=count_matching_elements(proc,quarks)
    naq=count_matching_elements(proc,antiquarks)
    if (options["include_3qqbar_processes"]) :
        if nq > 3 : return False    # at most three quarks
        if naq > 3 : return False   # at most three anti-quarks
    else :
        if nq > 2 : return False    # at most two quarks
        if naq > 2 : return False   # at most two anti-quarks
#    if nq < 2 : return False    # at least two quarks
#    if naq < 2 : return False   # at least two anti-quarks
    if nq != naq : return False # same number of quarks and anti-quarks
    # check charge conservation:
    if sum([charges3[x] for x in proc]) != 0 : return False
    # remove flavour changing currents:
    if (not options["include_cc_processes"]) :
        if sum([family[x] for x in proc]) != 0 : return False
        if 'w+' not in proc and 'w-' not in proc:
            for q in quarks:
                if count_matching_elements(proc,[q]) != count_matching_elements(proc,[q+'~']) : return False
    # need at least one quark line if there are colour singlets:
    if nq == 0 and count_matching_elements(proc,singlets) > 0 : return False
    return True

def CompatibleUniqueProc(process,proc):
    # Determines what are the required particles that need to be there
    # in the process (skipping 'p' and 'j'), and checks if all these
    # particles are actually part of 'proc'. Ignores if particles are
    # initial or final state.
    mandatory=[]
    proc_local=proc.copy()
    for part in process['initial_state']:
        if part != 'p':
            mandatory.append(part)
    mandatory.extend(process['rest'])
    try:
        for p in mandatory:
            proc_local.remove(p)
        return True
    except:
        return False
    
def CompatibleProc(process,proc):
    # Similar to CompatibleUniqueProc, but checks for both initial and
    # final state separately.
    proc_local=list(proc[2:])
    for i in [0,1]:
        if process['initial_state'][i] != 'p':
            if proc[i] != process['initial_state'][i]:
                return False
    try:
        for p in process['rest']:
            proc_local.remove(p)
    except:
        return False
    return True
    
def GenerateAllUniqueProcs(process):
    # Returns a set of unique_procs with the 'massless_QCD' particles
    # in a canonical (sorted) order. Even though the particles in
    # these processes are compatible with the input process, they are
    # not necessarily in the same order. The way this function works
    # is not really optimised: it first counts how many particles are
    # in the definition of 'massless_QCD', and then simply generates all 
    procs=[[]]
    for part in process['initial_state']:
        if part != 'p' and part not in massless_QCD:
            raise ValueError("Initial state should be a proton ('p').")
    jp=0
    for part in process['rest']:
        if part not in massless_QCD : continue
        jp=jp+1
    # The following for-loop will generate all processes of length
    # 'part' that contain all possible massless_QCD particles. Hence,
    # procs will have a size n^part, where n is the number of
    # particles in massless_QCD (by default 11=10quarks+1gluon). (This
    # is slightly reduced when part is large, since only up to 2 (or
    # 3) qqbar pairs are considered). To check that these are valid
    # processes (e.g., equal number of quarks and anti-quarks) is done
    # later in this function. 
    for part in range(process['jet_count']+2+jp): # number of jets + two incoming + other massless_QCD
        procs_new=[]
        for proc in procs:
            for p in massless_QCD:
                if (not options["include_3qqbar_processes"]) and part > 3 and p not in gluons: continue
                if part > 5 and p not in gluons: continue
                if not procs_new:
                    procs_new=[sorted(proc+[p])]
                else:
                    procs_new.append(sorted(proc+[p]))
        procs=procs_new.copy()
    # Add the non-massless_QCD particles in the process to all the
    # procs.
    for part in process['rest']:
        if part in massless_QCD : continue
        for proc in procs:
            proc.append(part)
    unique_procs=[]
    # Only at this stage check if they are valid (e.g., equal number
    # of quarks and anti-quarks) and compatible with the input
    # process.
    for proc in procs:
        if ValidProc(proc) and CompatibleUniqueProc(process,proc):
            unique_procs.append(tuple(proc))
    return set(unique_procs)

def GenerateAllProcs(unique_procs,process):
    # Determine which particles will become the actual initial state
    # ones. Do this, for each of the unique_procs, by considering all
    # possible pairs of two particles. Put them at then beginning of
    # the process (in both orders), and if this results in a process
    # compatible with the into process, its a valid process to
    # consider and add it to the 'procs'. Note that the 'procs' is a
    # set. Hence, the same process won't appear twice, so we do not
    # need to worry about identical particles.
    procs=set()
    for proc in unique_procs:
        # pick among all particles in the process two that become the two incoming ones.
        for i,j in itertools.combinations(range(len(proc)),2):
            if proc[i] in jet and proc[j] in jet:
                pair1=[proc[i],proc[j]]
                pair2=[proc[j],proc[i]]
                remaining=[e for k,e in enumerate(proc) if k not in (i,j)]
                proc1=tuple(pair1+remaining) # the two incoming + all others
                if CompatibleProc(process,proc1):
                    procs.add(proc1)
                proc2=tuple(pair2+remaining) # the two incoming (in other order) + all others
                if CompatibleProc(process,proc2):
                    procs.add(proc2)
    return procs

def CombineResults(results):
    # Results is an array of phase-space orderings
    # dictionaries. Concatenate/combine this into one big phase-space
    # orderings dictionary containing everything.
    phase_space_orders = {}
    for result in results:
        for key, value in result.items():
            if key in phase_space_orders:
                phase_space_orders[key].extend(value)
            else:
                phase_space_orders[key] = value
    return phase_space_orders
                
def ParseArgument():
    parser = argparse.ArgumentParser(description="Generate the full list of processes, ordered by phase-space order")
    parser.add_argument("process_string", type=str, help="Process to consider (e.g., 'p p > w+ z 4j')")
    parser.add_argument("-FS", "--flavour_scheme", type=int, choices=range(1, 6), metavar="[1-5]",
                        help="Switch to N-flavor scheme (NFS), where N is 1-5 (default=5)")
    parser.add_argument("-3", "--include_3qqbar", action='store_true', help="Include processes with up to 3 quark lines")
    parser.add_argument("-s", "--serial", action='store_true', help="Do not use multi-processes (parallel execution). Useful for debugging.")
    parser.add_argument("-cc", "--include_cc", action='store_true', help="Include flavour-changing processes")
    parser.add_argument("-res", "--resonance", action='store_true', help="Treat the lepton pair as a resonance")
    args=parser.parse_args()
    if (args.flavour_scheme):
        SwitchFlavourScheme(args.flavour_scheme)
    if args.include_3qqbar:
        options["include_3qqbar_processes"] = True
    else:
        options["include_3qqbar_processes"] = False
    if args.include_cc:
        options["include_cc_processes"] = True
    else:
        options["include_cc_processes"] = False
    if args.resonance:
        options["include_resonance"] = True
    else:
        options["include_resonance"] = False
    if args.serial:
        options["serial"] = True
    else:
        options["serial"] = False
    return ParseCollision(args.process_string)

def IdenticalParticleSymmetryFactor(proc):
    # Determine the symmetry factor that takes into account that we
    # are considering only a subset of all possible LC colour
    # orderings (based on the symmetry of the phase-space); this
    # factor is exactly identical to the usual factor coming from the
    # identical final state particles for all the colour particles:
    i_fac=1.0
    for p in all_coloured:
        i_fac*=max(1,math.factorial(proc[2:].count(p)))
    return i_fac


def build_process_index():
    # Dictionary to go from process + color order directly to the
    # phase-space order.
    for i, key in enumerate(all_keys_sorted):
        for (process, order, _) in phase_space_orders[key]:
            process_order_to_index[(process, order)] = i

def MultiChannelPartners(proc, perm, k, l):
    # Determine the multi-channel partners for the process 'proc' with
    # "colour-ordering" 'perm' (that belong to phase-space order 'k',
    # with process index 'l'; these latter two are only used to
    # overwrite the current proc+perm element in the
    # phase_space_orders dictionary with the one that includes the
    # multi-channel partners determined in this function). Note that
    # in the 'perm' also contains the "colour-ordering" of the
    # colour-singlet particles.
    all_possible_perms = []
    singlet_indices = [perm.index(i) for i, p in enumerate(proc) if p in singlets]
    anti_quark_indices = tuple([perm.index(i) for i, p in enumerate(proc) if p in antiquarks])
    # Precompute singlet permutations
    if len(singlet_indices) > 1:
        singlet_perms = tuple(itertools.permutations(singlet_indices))
    elif len(singlet_indices) == 1:
        singlet_perms = (tuple(singlet_indices),)
    else:
        singlet_perms = ()

    # Build all possible permutations of the colour-ordering (and
    # therefore have different phase-space orderings) but that have
    # the same matrix elements. These will be the multi-channel
    # partner processes. Currently, this is only based on the order of
    # the colour-singlets and how they are distributed among the
    # colour-strings. For processes with three quark lines, each
    # contribution is counted for twice. Hence, 'iden' should be set to
    # two, so that their contribution is halved.
    if len(anti_quark_indices) == 3:
        iden=2.
    else:
        iden=1.
    if not singlet_perms:
        all_possible_perms.append((perm,proc))
    else:
        if len(anti_quark_indices) == 1:
            # For a single quark line, the only thing we need to consider
            # is all the permutations of the colour singlets.
            for s in singlet_perms:
                order = []
                for i in range(len(perm)):
                    if i == anti_quark_indices[0]:
                        # Add all singlets in one go after the anti-quark
                        order.extend(perm[p] for p in (anti_quark_indices + s))
                    elif i not in singlet_indices:
                        # Add QCD particles one at the time
                        order.append(perm[i])
                all_possible_perms.append((tuple(order),tuple(proc)))
        elif len(anti_quark_indices) == 2:
            # For two quark lines, we need to consider both all the
            # permutations of the colour singlets AND how they are
            # distributed between the two quark-antiquark colour
            # groupings.
            for j in range(len(singlet_perms) + 1): # 'j' is the number of singlets attached to the 1st quark line
                for s in singlet_perms:
                    order = []
                    for i in range(len(perm)):
                        if i == anti_quark_indices[0]:
                            # Add j singlets after first anti-quark:
                            order.extend(perm[p] for p in ((anti_quark_indices[0],) + s[:j]))
                        elif i == anti_quark_indices[1]:
                            # Add rest of singlets after second anti-quark:
                            order.extend(perm[p] for p in ((anti_quark_indices[1],) + s[j:]))
                        elif i not in singlet_indices:
                            # Add the QCD particles one at the time
                            order.append(perm[i])
                    all_possible_perms.append((tuple(order),tuple(proc)))
        elif len(anti_quark_indices) == 3:
            for j1 in range(len(singlet_perms) + 1): # 'j1' is the number of singlets attached to the 1st quark line
                for j2 in range(len(singlet_perms)-j1 + 1): # 'j2' is the number of singlets attached to the 2nd quark line
                    for s in singlet_perms:
                        order = []
                        for i in range(len(perm)):
                            if i == anti_quark_indices[0]:
                                # Add j1 singlets after first anti-quark:
                                order.extend(perm[p] for p in ((anti_quark_indices[0],) + s[:j1]))
                            elif i == anti_quark_indices[1]:
                                # Add j2 of singlets after second anti-quark:
                                order.extend(perm[p] for p in ((anti_quark_indices[1],) + s[j1:j1+j2]))
                            elif i == anti_quark_indices[2]:
                                # Add rest of singlets after third anti-quark:
                                order.extend(perm[p] for p in ((anti_quark_indices[2],) + s[j1+j2:]))
                            elif i not in singlet_indices:
                                # Add the QCD particles one at the time
                                order.append(perm[i])
                        all_possible_perms.append((tuple(order),tuple(proc)))

    # The possible permutations should be processes that are already
    # included into other phase-space orderings. Look-up to which
    # phase-space orders these permutations belong. These are the
    # multi-channel partners.
    mt = []
    for (o,p) in all_possible_perms:
        idx = process_order_to_index.get((p, o))
        if idx is not None:
            if idx in mt:
                print("ERROR: found double. Each permutation should be unique for the multi-channel partners",p,o)
                quit()
            else:
                mt.append(idx)
        else:
            print("ERROR: expected multi-channel partner not found among phase-space orderings",p,o)
            quit()
    # Overwrite the current proc+perm element with the one that also
    # includes the multi-channel partners:
    phase_space_orders[k][l] = (proc, perm, tuple(sorted(mt)),1/iden)

def DetermineMultiChannelPartnersAndSymmetryFactor():
    # Using what's in all_keys_sorted, determine if a process needs
    # multi-channeling, and what are their multi-channel partners. It
    # also adds the symmetry factor for each process determined by
    # considering that we are including only the non-unique colour
    # orderings (when considering that the phase-space for identical
    # particles is symmetric) and we need to compensate for this.
    #
    # Build a dictionary to get the phase-space order directly from
    # the process + color order so that we can quickly find in which
    # phase-space order we find the multi-channel partners:
    build_process_index()
    # Determine the multi-channel partners:
    for key in all_keys_sorted:
        for i,(process,order,multichannel) in enumerate(phase_space_orders[key]):
            MultiChannelPartners(process,order,key,i)
    # Add the identical particle symmetry factor:
    for key in all_keys_sorted:
        for i,(process,order,multichannel,iden) in enumerate(phase_space_orders[key]):
            phase_space_orders[key][i]=(process,order,multichannel,iden*IdenticalParticleSymmetryFactor(process))

def ConvertProcToString(proc):
    # Convert the process 'proc' into a string.
    process,order,multi_channel,iden=proc
    crossed=[pdgs[p] if i>1 else pdgs[anti_particle[p]] for i,p in enumerate(process)] # cross intial state
    line=str(len(multi_channel))
    line=line+'   '+' '.join([str(m+1) for m in multi_channel])
    line=line+'   '+' '.join(crossed)
    line=line+'   '+' '.join([str(o+1) for o in order])
    line=line+'   '+str(iden)
    return line

def sort_by_pdg_codes(process):
    # Give a label to each process that can be used to order the
    # process in a process_list. This is not really important, but
    # makes the processes.txt file look neater.
    nq=count_matching_elements(process,quarks)
    same_flavour=max([count_matching_elements(process,[q]) for q in quarks])
    # first sort by 'nq', then by same_flavour, then by (modified) PDG codes:
    return (nq,same_flavour,[sort_particles[p] for p in process])

def sort_by_pdg_codes2(proc):
    process=proc[0]
    return sort_by_pdg_codes(process)

def WriteAllProcsIntoList(lep):
    # Convert all the information stored in the phase_space_orders
    # (or, actually, all_keys_sorted) into strings that can be
    # directly written in the processes.txt file.
    towrite=[]
    towrite.append(str(len(all_keys_sorted))) # number of phase-space orderings to consider
    towrite.append('')
    all_keys_sorted_new=copy.copy(all_keys_sorted)
    phase_space_orders_new={}
    # re-shuffle the phase space orders if there is a resonance
    if (options['include_resonance']):
        for i,key in enumerate(all_keys_sorted):
            first = phase_space_orders[key][0]
            if 'z' in first[0]:
                j = key.index(first[0].index("z"))
                new = tuple([first[0].index("z"),first[0].index("z")+1])
                t = key[0:j] + new + key[j+1:]
                all_keys_sorted_new[i]=t
                phase_space_orders_new[t]=phase_space_orders[key]
    else:
        phase_space_orders_new=copy.copy(phase_space_orders)

    for i,key in enumerate(all_keys_sorted_new):
        towrite.append(str(i+1)+'   '+str(len(phase_space_orders_new[key]))+'   '+str(max(len(proc[2]) for proc in phase_space_orders_new[key]))+'   '+' '.join([str(k+1) for k in key]))
        # order the processes in the process_list, so that we get a neat processes.txt file:
        process_list=sorted(phase_space_orders_new[key],key=sort_by_pdg_codes2)
        for proc in process_list:
            if (options['include_resonance']):
                if 'z' in proc[0]:
                    ib = proc[0].index("z")
                    t = proc[0][:ib] + tuple(lep) + proc[0][ib+1:]
                    proc = (t,) + proc[1:]
                    j = proc[1].index(ib)
                    new = tuple([ib,ib+1])
                    t = proc[1][0:j] + new + proc[1][j+1:]
                    proc = (proc[0],) + (t,) + proc[2:]
            process_line=ConvertProcToString(proc)
            towrite.append(process_line)
        towrite.append('')
        towrite.append('')
        towrite.append('')
    return towrite

def Addqq_dfProcesses(sorted_procs):
    # add all the multiple qq diff flavour processes by flipping orders of the quarks and anti-quarks
    i=0
    while i < len(sorted_procs):
        proc = sorted_procs[i]
        nq=count_matching_elements(proc,quarks)
        if nq >= 2 :  # multi-quark-line process
            qs=proc[0:nq]      # assume that first nq elements are the quarks
            aqs=proc[nq:2*nq]  # then come all the anti-quarks
            for q_perm in itertools.permutations(qs):
                for a_perm in itertools.permutations(aqs):
                    swapped_proc=list(q_perm)+list(a_perm)+proc[2*nq:]
                    if swapped_proc not in sorted_procs:
                        sorted_procs.insert(i+1,swapped_proc)
                        i+=1
        i+=1
    return sorted_procs

def WriteUniqueProcsIntoList(procs,lep):
    # Sort the unique processes to make them ready to be written to the file
    sorted_procs=sorted([sorted(proc,key=lambda x: sort_particles[x]) for proc in procs],key=sort_by_pdg_codes)
    if (options['include_resonance']):
        for proc in sorted_procs:
            if 'z' in proc:
                proc[proc.index("z"):proc.index("z")+1] = lep
    # in case of different flavour multiple-quark line processes, add all the possible orders:
    sorted_procs=Addqq_dfProcesses(sorted_procs)
    try:
        line=[str(len(sorted_procs[0]))+' '+str(len(sorted_procs))]
    except IndexError:
        print("ERROR: no processes found. Try './process_list.py --help' to get more information on usage")
        quit()
    for proc in sorted_procs:
        line.append(' '.join(pdgs[p] for p in proc))
    line.append('')
    line.append('')
    line.append('')
    return line

def CheckConsistency():
    allprocs={}
    for key in all_keys_sorted:
        for (process,order,multichannel,iden) in phase_space_orders[key]:
            proc=[process[0],process[1]]
            proc.extend(sorted(process[2:],key=lambda x: sort_particles[x]))
            proc=tuple(proc)
            if proc in allprocs:
                allprocs[proc]=allprocs[proc]+iden/len(multichannel)
            else:
                allprocs[proc]=iden/len(multichannel)
    for proc in allprocs.keys():
        if abs(allprocs[proc]-ExpectedNumberOfDualAmplitudes(proc)) > 1e-5:
            print('ERROR: inconsistent number of dual amplitudes for process:',proc,'. Found:',allprocs[proc],'. Expected:',ExpectedNumberOfDualAmplitudes(proc))
            quit()

def ExpectedNumberOfDualAmplitudes(proc):
    nq=count_matching_elements(proc,quarks)
    if nq == 0 :
        ng=count_matching_elements(proc,gluons)
        return math.factorial(ng-1)
    elif nq == 1 :
        ng=count_matching_elements(proc,gluons)
        return math.factorial(ng)
    elif nq == 2 :
        # number of gluons times number of ways gluons can be divided
        # times two ways of connecting quarks with anti-quarks
        ng=count_matching_elements(proc,gluons)
        return math.factorial(ng)*(ng+1)*2
    elif nq == 3 :
        ng=count_matching_elements(proc,gluons)
        return math.factorial(ng)*((ng+2)*(ng+1)/2)*6
    else:
        print("ERROR: unknown number of quarks",nq)
        quit()
        
    
if __name__ == "__main__":
    # Parse the argument. Cross the initial state particles to the
    # final state to avoid confusion about initial state quarks (that
    # are treated as anti-quarks when considering
    # e.g. colour-ordering)
    process=ParseArgument()
    lep = process["lep"]
    all_unique_procs=GenerateAllUniqueProcs(process)
    all_procs=GenerateAllProcs(all_unique_procs,process)
    # At this stage, we have two sets:
    # 'all_unique_procs' contains all processes compatible with the
    #     input process, where the massless_QCD particles are put in
    #     canonical order. No distinction between incoming and
    #     outgoing particles is made.
    # 'all_procs' contains all processes compatible with the input
    #     process. It has all possibilities for this two incoming
    #     particles (even though they are not crossed to the initial
    #     state), and only one order for all the final state
    #     particles, (i.e, only one of 'u u~ > d d~ g ' and 'u u~ > d
    #     g d~', since they are the same process).
    # No knowledge on colour orderings or phase-space orderings have
    # been considered up to now.
    #
    # Consider all the possible colour-orderings and collect all the
    # info into the phase_space_orders dictionary.
    if not options["serial"]:
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(ProcessProcess, all_procs)  # Parallelize across procs
    else:
        results=[ProcessProcess(x) for x in all_procs]
    phase_space_orders=CombineResults(results)
    all_keys_sorted=sorted(phase_space_orders.keys())
    DetermineMultiChannelPartnersAndSymmetryFactor() # updates the phase_space_orders dictionary
    # Check the consistency of the generated processes
    CheckConsistency()
    # write to disk
    towriteunique=WriteUniqueProcsIntoList(all_unique_procs,lep)
    towriteallprocs=WriteAllProcsIntoList(lep) # puts the phase_space_orders dictionary in a writable list
    
    with open('processes.txt','w') as f:
        f.write('\n'.join(towriteunique))
        f.write('\n'.join(towriteallprocs))
        f.write('\n')
