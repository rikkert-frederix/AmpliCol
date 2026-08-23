#!/usr/bin/env python
"""Build the ``processes.txt`` input used by the Fortran event generator.

The script takes a high-level process string, for example ``p p > w+ 2j``,
and expands it into the concrete partonic subprocesses, leading-colour
orderings, phase-space orderings, multichannel partners, and symmetry factors
needed by ``amplicol_generate``.

Important internal representations:

* ``process`` is the parsed user request returned by ``ParseCollision``. It
  contains the requested incoming particles, the number of inclusive jets, and
  the non-jet final-state particles that must be present.
* ``proc`` is a concrete subprocess tuple. Its first two entries are the
  incoming particles, and the remaining entries are final-state particles.
* ``perm`` is a zero-based colour-order permutation of the entries in
  ``proc``. The file written to disk converts these indices to one-based
  Fortran-style indices.
* ``phase_space_orders`` maps an internal phase-map key to records
  ``(proc, perm, multichannel_partners, factor, P, partner_Ps)``.  ``P`` maps
  canonical phase-space labels to the fixed labels used by ``proc`` and
  ``perm``; ``partner_Ps`` contains the corresponding map for every listed
  multichannel partner.  This lets relabelling-related coefficients share a
  canonical adaptive grid without changing which fixed-label amplitude is
  evaluated.
"""

import itertools
import re
import argparse
from collections import Counter, defaultdict
import multiprocessing
import math

PROCESS_FILE_VERSION = 4


# Global sets (make then 'frozenset' so that they are immutable):
quarks=frozenset({'d','u','s','c','b','t'})
antiquarks=frozenset({'d~','u~','s~','c~','b~','t~'})
singlets=frozenset({'a','z','w+','w-','e+','e-','mu+','mu-','ta+','ta-','ve','ve~','vm','vm~','vt','vt~','h'})
gluons=frozenset({'g'})
all_coloured=quarks | antiquarks | gluons
flavour_scheme=frozenset({'d','u','s','c','b'}) # all the massless quarks
flavour_scheme_number=5
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
channel_resonances = {}
process_provenance = ""
_resonance_history_cache = {}


def _build_three_point_vertices(active_flavours):
    """Return the physical cubic vertices needed for resonance discovery.

    Particles are represented in the all-outgoing convention used by the
    process builder.  Auxiliary fields used to factorise four-point vertices
    are deliberately absent: they are not physical Breit--Wigner poles.
    """
    vertices = {(21, 21, 21), (24, -24, 22), (24, -24, 23),
                (24, -24, 25), (23, 23, 25), (25, 25, 25)}
    for flavour in range(1, 7):
        vertices.add((21, flavour, -flavour))
        vertices.add((22, flavour, -flavour))
        vertices.add((23, flavour, -flavour))
    # A quark has a Higgs Yukawa current precisely when it is outside the
    # active, exactly-massless flavour scheme.  Top is therefore always
    # included, while b/c/s/u enter successively in lower schemes.
    for flavour in range(active_flavours + 1, 7):
        vertices.add((25, flavour, -flavour))
    for down in (1, 3, 5):
        up = down + 1
        vertices.add((24, down, -up))
        vertices.add((-24, -down, up))
    for charged, neutrino in ((11, 12), (13, 14), (15, 16)):
        vertices.add((22, charged, -charged))
        vertices.add((23, charged, -charged))
        vertices.add((23, neutrino, -neutrino))
        vertices.add((24, charged, -neutrino))
        vertices.add((-24, -charged, neutrino))
    return tuple(vertices)


THREE_POINT_VERTICES = _build_three_point_vertices(flavour_scheme_number)


def _anti_pdg(pdg):
    """Return the antiparticle PDG for the physical fields used here."""
    return pdg if pdg in (21, 22, 23, 25) else -pdg


def _build_current_combinations(vertices):
    """Map two outgoing currents to every current allowed by a cubic vertex."""
    combinations = defaultdict(set)
    for vertex in vertices:
        for current_index in range(3):
            current = _anti_pdg(vertex[current_index])
            children = [vertex[index] for index in range(3)
                        if index != current_index]
            combinations[tuple(sorted(children))].add(current)
    return {key: frozenset(value) for key, value in combinations.items()}


CURRENT_COMBINATIONS = _build_current_combinations(THREE_POINT_VERTICES)

def SwitchFlavourScheme(FS):
    """Set which quark flavours are treated as massless proton/jet content.

    ``p`` and ``j`` are intentionally kept identical in this code: both expand
    to gluons plus all quarks and antiquarks in the chosen flavour scheme.
    """
    # Overwrite the relevant global variables so that we switch to the
    # 'FS' flavour-scheme process definition.
    global flavour_scheme
    global flavour_scheme_number
    global massless_QCD
    global proton
    global jet
    global THREE_POINT_VERTICES
    global CURRENT_COMBINATIONS

    ordered_flavours=('d','u','s','c','b')
    if FS < 1 or FS > len(ordered_flavours):
        raise ValueError(f"unknown flavour scheme {FS}; expected 1--5")
    flavour_scheme_number=FS
    flavour_scheme=frozenset(ordered_flavours[:FS])
    massless_QCD=flavour_scheme | frozenset([q+'~' for q in flavour_scheme]) | gluons
    proton=massless_QCD
    jet=massless_QCD
    THREE_POINT_VERTICES=_build_three_point_vertices(FS)
    CURRENT_COMBINATIONS=_build_current_combinations(THREE_POINT_VERTICES)
    _resonance_history_cache.clear()

def ProcessProcess(proc):
    """Return phase-space groups for all valid colour orderings of ``proc``.

    The result is local to one subprocess and is later merged with the results
    from all other subprocesses.  The internal dictionary key contains both the
    canonical phase-space order and the (pre-canonicalisation) open-string
    connection.  The latter is essential: final-state relabelling can otherwise
    put a leading connected flow and an auxiliary-U(1) flow on the same adaptive
    integration grid even though their normalisations and peak histories differ.
    """
    phase_space_orders_local = {}

    valid_color_ord = GenerateValidColorOrders(proc)
    # In case we have identical final state particles, there can be
    # multiple colour orderings that are identical. Filter those out:
    unique_color_ord = [perm for perm in valid_color_ord if UniqueColorOrd(proc, perm)]

    # Group the possible colour orderings into their corresponding
    # phase-space orderings. The OrderProcPerm() function rearranges
    # the massless_QCD final state particles to reduce the number of
    # required phase-space orderings. 
    for perm in unique_color_ord:
        orientation = ThreeQuarkLineOrientation(proc, perm)
        ordered_proc, ordered_perm = OrderProcPerm(proc, perm)
        # Express the connection in the same canonical external-leg labelling
        # as the phase-space map.  This lets genuinely like-shaped flavour
        # subprocesses share a grid without collapsing distinct three-line
        # endpoint permutations into one coarse cycle-count class.
        topology = ColourTopologyTag(ordered_proc, ordered_perm) + (orientation,)
        zero = ordered_perm.index(0)
        perm_mapped = tuple(ordered_perm[zero:] + ordered_perm[:zero])
        channel_key = (perm_mapped, topology)
        if channel_key in phase_space_orders_local:
            phase_space_orders_local[channel_key].append((ordered_proc, ordered_perm, []))
        else:
            phase_space_orders_local[channel_key] = [(ordered_proc, ordered_perm, [])]
    return phase_space_orders_local


def GenerateValidColorOrders(proc):
    """Generate valid orders without permuting singlets through invalid slots.

    A brute-force permutation of every external occurrence scales as ``n!``
    and becomes prohibitive as soon as explicit four- or six-lepton decays are
    kept in the process.  First construct the much smaller coloured skeleton,
    then use the existing singlet-placement routine to insert every singlet
    permutation after an antiquark.  The final validity check keeps this
    exactly equivalent to the historical enumeration.
    """
    coloured = tuple(index for index, particle in enumerate(proc)
                     if particle in all_coloured)
    singlet_labels = tuple(index for index, particle in enumerate(proc)
                           if particle not in all_coloured)
    valid = set()
    for coloured_order in itertools.permutations(coloured):
        if not ValidColorOrd(proc, coloured_order):
            continue
        if not singlet_labels:
            valid.add(coloured_order)
            continue
        antiquark_positions = [
            position for position, label in enumerate(coloured_order)
            if proc[label] in antiquarks
        ]
        if not antiquark_positions:
            continue
        insertion = antiquark_positions[-1] + 1
        seed = coloured_order[:insertion] + singlet_labels + \
            coloured_order[insertion:]
        for order, _ in SingletMultiChannelOrders(proc, seed):
            if ValidColorOrd(proc, order):
                valid.add(order)
    return tuple(sorted(valid))


def ValidColorOrd(proc,perm):
    """Return whether ``perm`` is an allowed leading-colour ordering.

    Valid colour strings have the schematic form
    ``q, g, ..., g, qbar, s, ..., s`` repeated for each quark line, where
    ``s`` denotes a colour-singlet particle. Colour singlet permutations are
    deliberately kept at this stage; later they become multichannel partners,
    which lets the phase-space generator cover singular regions associated
    with singlets attached to different colour lines.
    """
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
            # UPDATE: Only remove duplicates from cyclic ordering. For three
            # quark lines this retains both orders of the two unanchored open
            # strings. Later they become multichannel partners when they have
            # the same fixed external-leg labelling; otherwise their old
            # half-weight normalization is retained.
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
    """Remove duplicate colour orderings caused by identical final particles.

    Starting at the first incoming particle, identical coloured final-state
    particles must appear in increasing original-index order. This keeps one
    representative of each physically equivalent ordering.
    """
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
    """Canonicalize a subprocess and its colour order.

    For up to two quark lines, final-state massless QCD particles are relabelled
    so their colour-order positions increase.  Three-line phase maps instead
    retain every distinguishable external label: relabelling those endpoints
    moves the physical fermion-transfer pole and was found to destroy the grid
    efficiency.  Only genuinely identical final particles are canonicalised in
    that case.
    """
    # Start by bringing the colour order into the canonical frame.
    zero=perm.index(0)
    perm_mapped=list(perm[zero:]+perm[:zero]) # cyclicly permute
    if count_matching_elements(proc,quarks) == 3:
        for particle in all_coloured:
            particle_positions = [
                i for i, part in enumerate(proc) if i > 1 and part == particle
            ]
            positions_in_order = [
                i for i, idx in enumerate(perm_mapped)
                if idx in particle_positions
            ]
            for i, val in zip(positions_in_order, sorted(particle_positions)):
                perm_mapped[i] = val
        perm_ordered=perm_mapped[len(perm)-zero:]+perm_mapped[:len(perm)-zero]
        return tuple(proc),tuple(perm_ordered)

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
    """Parse a process string such as ``p p > w+ 2j``.

    Non-proton initial states are crossed to the final-state convention used
    by the colour-ordering logic. A final token like ``3j`` is interpreted as
    three inclusive massless-QCD jets.  Decay products remain explicit;
    resonance histories are discovered later for each concrete subprocess.
    """
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
    return {"initial_state":initial_state,"jet_count":jet_count,"rest":rest}

def count_matching_elements(main_list,check_list):
    """Count entries of ``main_list`` whose value is present in ``check_list``."""
    # Counts how many elements of main_list are there in check_list
    # e.g., main_list=[a,b,b,b,c,c,d] ; check_list=[b,d,e] ; --> count_matching_elements=4
    main_counts=Counter(main_list)
    count=sum(main_counts[item] for item in check_list if item in main_counts)
    return count


def _is_lepton_pdg(pdg):
    return 11 <= abs(pdg) <= 16


def _lepton_pair_currents(first_pdg, second_pdg):
    """Return physical vector currents that can produce two leptons.

    The cubic-vertex table already contains the flavour and charge rules.  In
    particular, charged same-flavour pairs admit both ``gamma*`` and ``Z``,
    neutrino-antineutrino pairs admit ``Z``, and charged-current partners of
    one family admit exactly one of ``W+`` or ``W-``.
    """
    currents = CURRENT_COMBINATIONS.get(
        tuple(sorted((first_pdg, second_pdg))), ()
    )
    return tuple(sorted(set(currents) & {22, 23, 24, -24}))


def _lepton_pairings(process, labels=None):
    """Enumerate all complete physical pairings of lepton occurrences.

    Each returned pairing is a tuple of ``(labels, currents)`` records.  The
    occurrence labels make crossed pairings of identical four-lepton states
    distinct even when their particle names are the same.
    """
    external_pdgs = tuple(int(pdgs[particle]) for particle in process)
    if labels is None:
        labels = tuple(index for index, pdg in enumerate(external_pdgs)
                       if _is_lepton_pdg(pdg))
    else:
        labels = tuple(labels)
    if not labels:
        return ((),)
    if len(labels) % 2:
        return ()

    pairings = []

    def pair_remaining(remaining, prefix):
        if not remaining:
            pairings.append(tuple(sorted(prefix, key=lambda item: item[0])))
            return
        first = remaining[0]
        for position in range(1, len(remaining)):
            second = remaining[position]
            currents = _lepton_pair_currents(
                external_pdgs[first], external_pdgs[second]
            )
            if not currents:
                continue
            pair = (tuple(sorted((first, second))), currents)
            pair_remaining(
                remaining[1:position] + remaining[position + 1:],
                prefix + (pair,),
            )

    pair_remaining(tuple(sorted(labels)), ())
    return tuple(sorted(set(pairings)))


def _current_builder(external_pdgs):
    """Create a memoised tree-current finder for occurrence-tagged legs."""
    cache = {}

    def currents(mask):
        if mask in cache:
            return cache[mask]
        if mask == 0:
            result = frozenset()
        elif mask & (mask - 1) == 0:
            index = mask.bit_length() - 1
            result = frozenset((external_pdgs[index],))
        else:
            found = set()
            anchor = mask & -mask
            subset = (mask - 1) & mask
            while subset:
                other = mask ^ subset
                if other and subset & anchor:
                    left = currents(subset)
                    right = currents(other)
                    for first in left:
                        for second in right:
                            found.update(CURRENT_COMBINATIONS.get(
                                tuple(sorted((first, second))), ()))
                subset = (subset - 1) & mask
            result = frozenset(found)
        cache[mask] = result
        return result

    return currents


def _has_higgs_vector_split(mask, currents):
    """Check that a Higgs current is specifically reachable through VV."""
    anchor = mask & -mask
    subset = (mask - 1) & mask
    while subset:
        other = mask ^ subset
        if other and subset & anchor:
            left = currents(subset)
            right = currents(other)
            if (23 in left and 23 in right) or \
                    (24 in left and -24 in right) or \
                    (-24 in left and 24 in right):
                return True
        subset = (subset - 1) & mask
    return False


def _compatible_resonances(first, second):
    """Return whether two descendant sets form a laminar family."""
    first_labels = frozenset(first[1])
    second_labels = frozenset(second[1])
    if first_labels == second_labels:
        return False
    overlap = first_labels & second_labels
    return not overlap or first_labels <= second_labels or \
        second_labels <= first_labels


def _canonical_history(entries):
    """Canonical inner-to-outer ordering for one resonance history."""
    return tuple(sorted(set(entries), key=lambda item: (
        len(item[1]), item[1], item[0]
    )))


def _history_has_physical_tree(process, history):
    """Check that all requested poles coexist in one cubic tree.

    Laminar descendant masks are necessary but not sufficient: for example,
    an ``H`` mask and two nested photon masks are geometrically compatible but
    this model has no tree-level ``H gamma gamma`` vertex. Required masks are
    therefore treated as indivisible currents while the same recursive current
    builder used for discovery closes the complete external tree.
    """
    if not history:
        return True
    external_pdgs = tuple(int(pdgs[particle]) for particle in process)
    full_mask = (1 << len(process)) - 1
    required = {}
    for resonance, labels in history:
        mask = sum(1 << label for label in labels)
        if mask in required and required[mask] != resonance:
            return False
        required[mask] = resonance

    cache = {}

    def constrained_currents(mask):
        if mask in cache:
            return cache[mask]
        if mask == 0:
            result = frozenset()
        elif mask & (mask - 1) == 0:
            index = mask.bit_length() - 1
            result = frozenset((external_pdgs[index],))
        else:
            found = set()
            anchor = mask & -mask
            subset = (mask - 1) & mask
            while subset:
                other = mask ^ subset
                if other and subset & anchor:
                    splits_required_current = any(
                        required_mask != mask and
                        (required_mask & mask) == required_mask and
                        required_mask & subset and required_mask & other
                        for required_mask in required
                    )
                    if not splits_required_current:
                        for first in constrained_currents(subset):
                            for second in constrained_currents(other):
                                found.update(CURRENT_COMBINATIONS.get(
                                    tuple(sorted((first, second))), ()
                                ))
                subset = (subset - 1) & mask
            result = frozenset(found)
        if mask in required:
            result = frozenset((required[mask],)) \
                if required[mask] in result else frozenset()
        cache[mask] = result
        return result

    anchor = full_mask & -full_mask
    subset = (full_mask - 1) & full_mask
    while subset:
        other = full_mask ^ subset
        if other and subset & anchor:
            if any(required_mask & subset and required_mask & other
                   for required_mask in required):
                subset = (subset - 1) & full_mask
                continue
            left = constrained_currents(subset)
            right = constrained_currents(other)
            if any(_anti_pdg(current) in right for current in left):
                return True
        subset = (subset - 1) & full_mask
    return False


def DiscoverResonanceHistories(process):
    """Discover physical tree-level resonance histories for ``process``.

    External occurrences, rather than particle names, label descendants.  A
    current is retained only when the complementary external legs can supply
    the antiparticle current, which excludes poles that cannot occur in a
    complete tree diagram for the concrete subprocess. Photon pair assignments
    use an ordinary pair-block density and therefore do not appear as mapped
    finite-width poles in the returned histories.
    """
    cache_key = (flavour_scheme_number, process)
    if cache_key in _resonance_history_cache:
        return _resonance_history_cache[cache_key]

    external_pdgs = tuple(int(pdgs[particle]) for particle in process)
    final_indices = tuple(range(2, len(process)))
    final_leptons = tuple(index for index in final_indices
                          if _is_lepton_pdg(external_pdgs[index]))
    has_top_decay_seed = any(
        external_pdgs[index] in (5, -5, 24, -24)
        for index in final_indices
    )
    if len(final_leptons) < 2 and not has_top_decay_seed:
        result = ((),)
        _resonance_history_cache[cache_key] = result
        return result

    currents = _current_builder(external_pdgs)
    full_mask = (1 << len(process)) - 1
    final_mask = full_mask ^ 3
    candidates = set()
    subset = final_mask
    while subset:
        if subset.bit_count() >= 2:
            labels = tuple(index for index in final_indices
                           if subset & (1 << index))
            values = currents(subset)
            complement_values = currents(full_mask ^ subset)
            all_leptons = all(_is_lepton_pdg(external_pdgs[index])
                              for index in labels)

            if all_leptons:
                # A photon assignment is meaningful only for a genuine
                # two-body charged-lepton current. Massive vector parents may
                # contain more leptons through a nested decay tree.
                if len(labels) == 2 and 22 in values and \
                        22 in complement_values:
                    candidates.add((22, labels))
                for resonance in (23, 24, -24):
                    if resonance in values and \
                            _anti_pdg(resonance) in complement_values:
                        candidates.add((resonance, labels))

            higgs_to_leptons = all_leptons and len(labels) >= 4
            if higgs_to_leptons and 25 in values and \
                    25 in complement_values and \
                    _has_higgs_vector_split(subset, currents):
                candidates.add((25, labels))

            for resonance, bottom, weak_boson in (
                    (6, 5, 24), (-6, -5, -24)):
                if resonance not in values or \
                        _anti_pdg(resonance) not in complement_values:
                    continue
                for bottom_label in labels:
                    if external_pdgs[bottom_label] != bottom:
                        continue
                    daughters = subset ^ (1 << bottom_label)
                    daughter_labels = tuple(
                        index for index in final_indices
                        if daughters & (1 << index)
                    )
                    explicit_w = len(daughter_labels) == 1 and \
                        external_pdgs[daughter_labels[0]] == weak_boson
                    leptonic_w = daughter_labels and all(
                        _is_lepton_pdg(external_pdgs[index])
                        for index in daughter_labels
                    )
                    if (explicit_w or leptonic_w) and \
                            weak_boson in currents(daughters):
                        candidates.add((resonance, labels))
                        break
        subset = (subset - 1) & final_mask

    candidates = tuple(sorted(candidates, key=lambda item: (
        len(item[1]), item[1], item[0]
    )))
    histories = {()}
    candidate_set = set(candidates)
    pair_candidates = {
        candidate for candidate in candidates
        if len(candidate[1]) == 2 and
        all(label in final_leptons for label in candidate[1]) and
        candidate[0] in (22, 23, 24, -24)
    }

    # A lepton-bearing density starts from one *complete* perfect matching.
    # Choosing one vector current for every pair gives gamma/Z alternatives
    # for charged neutral currents and the competing WW/ZZ histories in four-
    # lepton states. A gamma* assignment is represented by the ordinary
    # pair-block density: unlike Z/W it has no finite-width pole to serialize.
    pair_histories = set()
    if not final_leptons:
        pair_histories.add(())
    else:
        for pairing in _lepton_pairings(process, final_leptons):
            choices = []
            for labels, allowed_currents in pairing:
                available = tuple(
                    (resonance, labels) for resonance in allowed_currents
                    if (resonance, labels) in candidate_set
                )
                if not available:
                    break
                choices.append(available)
            else:
                for combination in itertools.product(*choices):
                    pair_histories.add(_canonical_history(combination))

    histories.update(pair_histories)

    # Higgs, top, and broader vector poles can be layered on a complete pair
    # history. This retains H->VV and t->bW mappings while ensuring that no
    # added density leaves another external lepton unpaired.
    structural = tuple(candidate for candidate in candidates
                       if candidate not in pair_candidates)

    def add_structural_combinations(seed, prefix=(), start=0):
        combined_so_far = seed + prefix
        for index in range(start, len(structural)):
            candidate = structural[index]
            if not all(_compatible_resonances(candidate, existing)
                       for existing in combined_so_far):
                continue
            combined = prefix + (candidate,)
            histories.add(_canonical_history(seed + combined))
            add_structural_combinations(seed, combined, index + 1)

    for pair_history in pair_histories:
        add_structural_combinations(pair_history)

    histories = {
        history for history in histories
        if _history_has_physical_tree(process, history)
    }
    # Do not turn the massive companion of a gamma* pair into a partial
    # resonance history. The complete gamma-bearing assignment is represented
    # by its ordinary pair-block density; fully massive assignments retain all
    # of their Z/W poles.
    histories = {
        () if any(candidate[0] == 22 for candidate in history) else history
        for history in histories
    }

    result = tuple(sorted(histories, key=lambda history: (
        len(history), tuple((len(item[1]), item[1], item[0])
                            for item in history)
    )))
    _resonance_history_cache[cache_key] = result
    return result

def ValidProc(proc):
    """Apply process-level physics and code-support constraints."""
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
    # Every external lepton line must close through a physical neutral- or
    # charged-vector current. This rejects odd and flavour-mismatched lepton
    # collections before factorial colour-order generation starts.
    lepton_labels = tuple(index for index, particle in enumerate(proc)
                          if _is_lepton_pdg(int(pdgs[particle])))
    if lepton_labels and not _lepton_pairings(proc, lepton_labels):
        return False
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
    """Check a canonical particle multiset against the user request.

    This ignores the distinction between initial and final state and only
    checks that all mandatory non-inclusive particles are present.
    """
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
    """Check a concrete subprocess against the requested initial/final state."""
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
    """Generate canonical particle multisets compatible with the user request.

    The returned tuples do not yet decide which particles are incoming. They
    contain all requested particles plus enough massless-QCD particles to
    account for the two beams and requested inclusive jets.  Generate those
    particles as multisets directly: enumerating every ordered choice and
    sorting it afterwards repeats each six-parton multiset up to ``6!`` times.
    """
    for part in process['initial_state']:
        if part != 'p' and part not in massless_QCD:
            raise ValueError("Initial state should be a proton ('p').")

    requested_qcd = sum(part in massless_QCD for part in process['rest'])
    qcd_slots = process['jet_count'] + 2 + requested_qcd
    fixed_particles = [
        part for part in process['rest'] if part not in massless_QCD
    ]

    # ValidProc requires the complete process to contain equally many quarks
    # and antiquarks, with at most two (or three) lines.  Enumerate only count
    # combinations that can meet that condition.  The fixed-particle offsets
    # matter for processes containing explicit massive quarks, such as t W.
    fixed_quarks = count_matching_elements(fixed_particles, quarks)
    fixed_antiquarks = count_matching_elements(fixed_particles, antiquarks)
    max_quark_lines = 3 if options["include_3qqbar_processes"] else 2
    massless_quarks = sorted(quarks & massless_QCD)
    massless_antiquarks = sorted(antiquarks & massless_QCD)

    unique_procs = set()
    first_line_count = max(fixed_quarks, fixed_antiquarks)
    for line_count in range(first_line_count, max_quark_lines + 1):
        nquarks = line_count - fixed_quarks
        nantiquarks = line_count - fixed_antiquarks
        ngluons = qcd_slots - nquarks - nantiquarks
        if ngluons < 0:
            continue
        quark_choices = itertools.combinations_with_replacement(
            massless_quarks, nquarks
        )
        for quark_choice in quark_choices:
            antiquark_choices = itertools.combinations_with_replacement(
                massless_antiquarks, nantiquarks
            )
            for antiquark_choice in antiquark_choices:
                proc = sorted(
                    quark_choice + antiquark_choice + ('g',) * ngluons
                )
                proc.extend(fixed_particles)
                if ValidProc(proc) and CompatibleUniqueProc(process, proc):
                    unique_procs.add(tuple(proc))
    return unique_procs

def GenerateAllProcs(unique_procs,process):
    """Expand canonical multisets into concrete incoming/outgoing processes.

    For each unique particle multiset, every valid pair of QCD particles is
    tried as the two incoming partons, in both beam orders. The set return
    value removes duplicates from identical particles automatically.
    """
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
    """Merge per-subprocess phase-space dictionaries."""
    phase_space_orders = {}
    for result in results:
        for key, value in result.items():
            if key in phase_space_orders:
                phase_space_orders[key].extend(value)
            else:
                phase_space_orders[key] = value
    return phase_space_orders
                
def ParseArgument():
    """Parse command-line options and return the normalized process request."""
    global process_provenance
    parser = argparse.ArgumentParser(description="Generate the full list of processes, ordered by phase-space order")
    parser.add_argument("process_string", type=str, help="Process to consider (e.g., 'p p > w+ z 4j', (including quotation marks))")
    parser.add_argument("-FS", "--flavour_scheme", type=int, choices=range(1, 6), metavar="[1-5]",
                        help="Switch to N-flavor scheme (NFS), where N is 1-5 (default=5)")
    parser.add_argument("-3", "--include_3qqbar", action='store_true', help="Include processes with up to three quark lines (instead of just two).")
    parser.add_argument("-s", "--serial", action='store_true', help="Do not use multi-processes (parallel execution). Useful for debugging.")
    parser.add_argument("-cc", "--include_cc", action='store_true', help="Include flavour-changing processes")
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
    if args.serial:
        options["serial"] = True
    else:
        options["serial"] = False
    options["flavour_scheme"] = args.flavour_scheme or 5
    process_provenance = args.process_string
    return ParseCollision(args.process_string)

def IdenticalParticleSymmetryFactor(proc):
    """Return the identical-particle factor for coloured final states."""
    # Since equivalent leading-colour orderings are removed using phase-space
    # symmetry, the compensating factor is the usual final-state identical
    # particle symmetry factor, restricted here to coloured particles.
    i_fac=1.0
    for p in all_coloured:
        i_fac*=max(1,math.factorial(proc[2:].count(p)))
    return i_fac


def PhaseSpaceOrderFromKey(key):
    """Return the external-leg permutation stored in an internal channel key.

    Distinct matrix-element coefficients can require independent adaptive grids
    even when they use the same phase-space permutation.  After multichannel
    finalisation such duplicate channels are keyed by ``(order, topology)``;
    the historical plain tuple remains accepted for compatibility.
    """
    if len(key) >= 2 and isinstance(key[0], tuple):
        return key[0]
    return tuple(key)


def ColourTopologyFromKey(key):
    """Return the immutable colour-connection tag in an internal key."""
    if len(key) >= 2 and isinstance(key[0], tuple):
        return key[1]
    # Backward-compatible fallback for callers constructing historical keys.
    return (0, (), 0, 0)


def build_process_index():
    """Index amplitude targets by every phase-space channel containing them."""
    process_order_to_index.clear()
    indices = defaultdict(list)
    for channel, key in enumerate(all_keys_sorted):
        for row in phase_space_orders[key]:
            process, order = row[:2]
            if channel not in indices[(process, order)]:
                indices[(process, order)].append(channel)
    process_order_to_index.update(
        (target, tuple(channels)) for target, channels in indices.items()
    )


def QuarkLineBlocks(proc, perm, expected_lines=None):
    """Split a cyclic colour order into fixed-label open-string blocks.

    Each returned block starts with a quark and contains its ordered gluons,
    closing antiquark, and any following colour singlets.  The block containing
    external leg zero is placed first.  No external labels are canonicalised.
    """
    quark_positions = [i for i, idx in enumerate(perm) if proc[idx] in quarks]
    if expected_lines is not None and len(quark_positions) != expected_lines:
        raise ValueError(
            f"expected {expected_lines} open colour strings, found "
            f"{len(quark_positions)}"
        )
    if not quark_positions:
        return ()

    zero_position = perm.index(0)
    starts_before_zero = [i for i in quark_positions if i <= zero_position]
    anchor_start = starts_before_zero[-1] if starts_before_zero \
        else quark_positions[-1]
    anchored_order = perm[anchor_start:] + perm[:anchor_start]

    blocks = []
    for idx in anchored_order:
        if proc[idx] in quarks:
            blocks.append([])
        if not blocks:
            raise ValueError("colour order does not start on an open string")
        blocks[-1].append(idx)
    if expected_lines is not None and len(blocks) != expected_lines:
        raise ValueError("could not split colour order into open strings")
    if 0 not in blocks[0]:
        raise ValueError("could not anchor colour strings on external leg zero")
    for block in blocks:
        if sum(proc[idx] in antiquarks for idx in block) != 1:
            raise ValueError("an open colour string must contain one antiquark")
    return tuple(tuple(block) for block in blocks)


def _cycle_count(permutation):
    """Return the number of cycles of a permutation in one-line notation."""
    seen = set()
    cycles = 0
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycles += 1
        current = start
        while current not in seen:
            seen.add(current)
            current = permutation[current]
    return cycles


def ColourTopologyTag(proc, perm):
    """Describe the labelled open-string connection before canonicalisation.

    The endpoint permutation distinguishes colour coefficients that may share a
    numerical phase-space order after :func:`OrderProcPerm`.  When a unique
    flavour-preserving fermion-line matching exists, the final entry is the
    cycle count relative to that matching.  It is ``-1`` for identical-flavour
    or flavour-changing cases where external PDGs do not identify a unique Wick
    contraction; those cases are augmented conservatively later.
    """
    nlines = sum(p in quarks for p in proc)
    # Zero- and one-line coefficients have no disconnected colour-flow
    # topology.  Keep their historical shared phase-space grouping, including
    # mixed inclusive subprocess groups.
    if nlines <= 1:
        return (0, (), 0)

    blocks = QuarkLineBlocks(proc, perm, nlines)
    quark_legs = sorted(block[0] for block in blocks)
    antiquark_legs = sorted(
        next(idx for idx in block if proc[idx] in antiquarks)
        for block in blocks
    )
    endpoint_by_quark = {}
    for block in blocks:
        endpoint_by_quark[block[0]] = next(
            idx for idx in block if proc[idx] in antiquarks
        )
    endpoint_permutation = tuple(
        antiquark_legs.index(endpoint_by_quark[q]) for q in quark_legs
    )

    # Enumerating at most 3! matchings is cheap and avoids choosing an arbitrary
    # Wick contraction for repeated flavours.
    physical_matchings = []
    for matching in itertools.permutations(range(nlines)):
        if all(
            proc[antiquark_legs[matching[i]]] == anti_particle[proc[q]]
            for i, q in enumerate(quark_legs)
        ):
            physical_matchings.append(matching)
    if len(physical_matchings) == 1:
        matching = physical_matchings[0]
        inverse_matching = [0] * nlines
        for q_index, antiquark_index in enumerate(matching):
            inverse_matching[antiquark_index] = q_index
        relative_permutation = tuple(
            inverse_matching[endpoint] for endpoint in endpoint_permutation
        )
        relative_cycles = _cycle_count(relative_permutation)
        connection = relative_permutation
    else:
        relative_cycles = -1
        connection = endpoint_permutation

    return (nlines, connection, relative_cycles)


def CanonicalPhaseMap(proc, perm):
    """Canonicalize only a phase header, retaining the amplitude labels.

    The returned permutation maps canonical (``base``) external labels to the
    fixed labels used by ``proc`` and ``perm``::

        target_label = permutation[base_label]

    Initial-state and colour-singlet labels are fixed.  Consequently a point
    generated with the canonical header can be evaluated for the fixed-label
    target by assigning ``p_target[P[base]] = p_base[base]``.  The returned
    ``base_proc`` and ``base_perm`` are used only to derive a canonical colour
    topology; they are not substituted for the target amplitude record.
    """
    zero = perm.index(0)
    natural_header = tuple(perm[zero:] + perm[:zero])
    final_qcd_labels = tuple(sorted(
        i for i, particle in enumerate(proc)
        if i > 1 and particle in massless_QCD
    ))
    qcd_positions = [
        position for position, label in enumerate(natural_header)
        if label in final_qcd_labels
    ]
    if len(qcd_positions) != len(final_qcd_labels):
        raise ValueError("could not locate every final-state QCD label")

    canonical_header = list(natural_header)
    base_to_target = list(range(len(proc)))
    for position, base_label in zip(qcd_positions, final_qcd_labels):
        base_to_target[base_label] = natural_header[position]
        canonical_header[position] = base_label

    if sorted(base_to_target) != list(range(len(proc))):
        raise ValueError("phase-map label transformation is not a permutation")
    target_to_base = [None] * len(proc)
    for base_label, target_label in enumerate(base_to_target):
        target_to_base[target_label] = base_label

    base_proc = tuple(proc[base_to_target[label]] for label in range(len(proc)))
    base_perm = tuple(target_to_base[label] for label in perm)
    canonical_header = tuple(canonical_header)
    base_to_target = tuple(base_to_target)

    if tuple(base_to_target[label] for label in canonical_header) != \
            natural_header:
        raise ValueError("canonical phase header does not reconstruct its target")
    base_zero = base_perm.index(0)
    if tuple(base_perm[base_zero:] + base_perm[:base_zero]) != canonical_header:
        raise ValueError("canonical colour order and phase header disagree")
    return canonical_header, base_to_target, base_proc, base_perm


def AdaptiveTopologyClass(topology):
    """Return the compact adaptive class for a canonical colour topology.

    The historical compact classes are retained through two quark lines.
    Three-line coefficients use the complete endpoint permutation: the three
    transpositions and two oriented three-cycles put fermion-transfer poles on
    different labelled momenta and do not train one grid efficiently.  Actual
    flavour names are omitted, so kinematically equivalent flavour variants
    continue to share a grid.
    """
    nlines, connection, cycles = topology[:3]
    if nlines == 3:
        if cycles < 0:
            return ("mixed", 3, connection)
        return ("flow", 3, connection)
    if nlines < 2 or cycles == 1:
        return ("leading", 0)
    if cycles < 0:
        return ("mixed", nlines)
    return ("disconnected", cycles - 1)


def ThreeQuarkLineBlocks(proc, perm):
    """Split ``perm`` into three strings, anchored on the one containing zero."""
    return QuarkLineBlocks(proc, perm, 3)


def ThreeQuarkLineOrientation(proc, perm):
    """Label the two cyclic representations of three open colour strings."""
    if sum(p in quarks for p in proc) != 3:
        return 0
    blocks = ThreeQuarkLineBlocks(proc, perm)
    return int(blocks[1][0] > blocks[2][0])


def SwapThreeQuarkLineOrder(proc, perm):
    """Swap the two unanchored open colour strings in a three-line order.

    Cyclic invariance fixes the string containing external leg zero.  The two
    remaining strings are the redundant representations retained by
    :func:`ValidColorOrd`.  Return the canonical process/order pair used by the
    phase-space-group lookup.
    """
    blocks = ThreeQuarkLineBlocks(proc, perm)

    swapped_order = tuple(blocks[0] + blocks[2] + blocks[1])
    return OrderProcPerm(proc, swapped_order)


def ThreeQuarkLineBlockSignature(proc, perm):
    """Return the exact coloured-leg strings represented by ``perm``.

    The ordering of the three independent strings is irrelevant, but the leg
    coloured-leg labels *within* every string are not.  Colour singlets are
    deliberately omitted because their placements are already canonicalized by
    the amplitude reader.  In particular, canonicalizing two identical final
    coloured particles can move their momenta between strings.  Such a
    relabelled ordering is not the same colour coefficient point by point and
    must not be used as a multichannel partner without an explicit momentum
    permutation.
    """
    blocks = ThreeQuarkLineBlocks(proc, perm)
    return tuple(
        sorted(tuple(idx for idx in block if proc[idx] in all_coloured)
               for block in blocks)
    )


def SingletMultiChannelOrders(proc, perm):
    """Return equivalent phase-space orders with physical lepton blocks.

    Non-leptonic singlets retain the historical permutation and colour-line
    distribution. Leptons are different: every occurrence stays next to its
    partner from a physical ``gamma*/Z/W`` current. Four- and six-lepton
    states enumerate all complete pairings, but never split a pair over two
    colour lines or construct an ordering with an unpaired lepton.
    """
    singlet_labels = tuple(index for index, particle in enumerate(proc)
                           if particle in singlets)
    if not singlet_labels:
        return ((perm, tuple(proc)),)

    lepton_labels = tuple(
        index for index in singlet_labels
        if _is_lepton_pdg(int(pdgs[proc[index]]))
    )
    pairings = _lepton_pairings(proc, lepton_labels)
    if lepton_labels and not pairings:
        return ()

    other_blocks = tuple((index,) for index in singlet_labels
                         if index not in lepton_labels)
    anti_positions = tuple(
        perm.index(index) for index, particle in enumerate(proc)
        if particle in antiquarks
    )
    if not anti_positions:
        return ()
    singlet_positions = frozenset(perm.index(index)
                                  for index in singlet_labels)
    all_possible_perms = set()

    for pairing in pairings:
        lepton_blocks = tuple(labels for labels, _ in pairing)
        blocks = lepton_blocks + other_blocks
        for ordered_blocks in itertools.permutations(blocks):
            # A weak composition assigns consecutive blocks to each open
            # colour string. Blocks, rather than individual leptons, are the
            # unit of distribution.
            separators_iter = itertools.combinations_with_replacement(
                range(len(blocks) + 1), max(len(anti_positions) - 1, 0)
            )
            for separators in separators_iter:
                boundaries = (0,) + separators + (len(blocks),)
                blocks_by_antiquark = {
                    position: ordered_blocks[
                        boundaries[line]:boundaries[line + 1]
                    ]
                    for line, position in enumerate(anti_positions)
                }
                order = []
                for position, label in enumerate(perm):
                    if position in anti_positions:
                        order.append(label)
                        for block in blocks_by_antiquark[position]:
                            order.extend(block)
                    elif position not in singlet_positions:
                        order.append(label)
                all_possible_perms.add((tuple(order), tuple(proc)))

    return tuple(sorted(all_possible_perms))


def BuildTopologyAwareMultiChannels():
    """Build compact grids without relabelling fixed amplitude targets.

    For three quark lines, the phase header alone is canonicalized and its
    base-to-target leg permutation is retained on the row.  Singlet placements
    and exact fixed-label swaps of the two unanchored colour strings are joined
    into one multichannel component.  Swaps that additionally relabel identical
    final particles share an adaptive grid but remain different MIS targets. If
    two maps in one component canonicalize to the same header and topology, a
    small component-local source rank keeps both maps present.
    """
    global phase_space_orders, all_keys_sorted

    nodes = []
    node_index = defaultdict(list)
    exact_target_index = defaultdict(list)
    for channel, key in enumerate(all_keys_sorted):
        for row_index, row in enumerate(phase_space_orders[key]):
            process, order = row[:2]
            raw_topology = ColourTopologyFromKey(key)
            nlines = raw_topology[0]
            if nlines == 3:
                current_map, permutation, base_process, base_order = \
                    CanonicalPhaseMap(process, order)
                canonical_topology = ColourTopologyTag(
                    base_process, base_order
                )
            else:
                current_map = PhaseSpaceOrderFromKey(key)
                permutation = tuple(range(len(process)))
                canonical_topology = raw_topology[:3]
            base_key = (
                current_map,
                AdaptiveTopologyClass(canonical_topology),
            )
            node = len(nodes)
            nodes.append({
                "source_key": key,
                "row_index": row_index,
                "old_channel": channel,
                "process": process,
                "order": order,
                "raw_topology": raw_topology,
                "nlines": nlines,
                "base_key": base_key,
                "permutation": permutation,
            })
            node_index[(process, order, raw_topology)].append(node)
            exact_target_index[(process, order)].append(node)

    parent = list(range(len(nodes)))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(first, second):
        first = find(first)
        second = find(second)
        if first != second:
            parent[second] = first

    singlet_channel_count = [1] * len(nodes)
    for node, record in enumerate(nodes):
        process = record["process"]
        order = record["order"]
        topology = record["raw_topology"]
        singlet_channels = set()
        singlet_orders = SingletMultiChannelOrders(process, order)
        for partner_order, partner_process in singlet_orders:
            if record["nlines"] == 3:
                matches = exact_target_index.get(
                    (partner_process, partner_order), ()
                )
            else:
                matches = node_index.get(
                    (partner_process, partner_order, topology), ()
                )
            if len(matches) != 1:
                raise ValueError(
                    "expected one singlet multichannel target, found "
                    f"{len(matches)}: {partner_process} {partner_order}"
                )
            partner = matches[0]
            union(node, partner)
            singlet_channels.add(nodes[partner]["old_channel"])
        if len(singlet_channels) != len(singlet_orders):
            raise ValueError(
                "singlet permutations unexpectedly share a phase-space group: "
                f"{process} {order}"
            )
        singlet_channel_count[node] = len(singlet_channels)

        if record["nlines"] == 3:
            blocks = ThreeQuarkLineBlocks(process, order)
            raw_swapped_order = tuple(blocks[0] + blocks[2] + blocks[1])
            swapped_process, swapped_order = OrderProcPerm(
                process, raw_swapped_order
            )
            # Only a swap that leaves every fixed external label unchanged is
            # a pointwise identity and hence a valid MIS partner.  If
            # canonicalizing identical final particles changes the order, the
            # two rows obey A_2(p)=A_1(Pp), not A_2(p)=A_1(p).  They may still
            # share the canonical adaptive grid through ``base_key``, but must
            # retain separate components and their original one-half weights.
            if swapped_process == process and \
                    swapped_order == raw_swapped_order:
                matches = exact_target_index.get(
                    (swapped_process, swapped_order), ()
                )
                if len(matches) != 1:
                    raise ValueError(
                        "expected one exact three-string partner target, found "
                        f"{len(matches)}: {swapped_process} {swapped_order}"
                    )
                union(node, matches[0])

    components = defaultdict(list)
    for node in range(len(nodes)):
        components.setdefault(find(node), []).append(node)

    component_maps = []
    all_channel_keys = set()
    for component in components.values():
        old_weight = math.fsum(
            (0.5 if nodes[node]["nlines"] == 3 else 1.0)
            / singlet_channel_count[node]
            for node in component
        )

        # Rank only maps of this coefficient that would otherwise collapse to
        # one grid.  Ranks are local to a component, so unrelated flavour
        # variants retain the maximum useful grid sharing.
        buckets = defaultdict(list)
        for node in component:
            buckets[nodes[node]["base_key"]].append(node)

        maps = {}
        for base_key, bucket in buckets.items():
            ordered_bucket = sorted(
                bucket,
                key=lambda item: (
                    nodes[item]["permutation"],
                    nodes[item]["order"],
                    nodes[item]["process"],
                ),
            )
            for source_rank, node in enumerate(ordered_bucket):
                record = nodes[node]
                channel_key = base_key + (source_rank,)
                if channel_key in maps:
                    raise ValueError(
                        "two coefficient maps share a ranked channel key"
                    )
                maps[channel_key] = (
                    record["process"],
                    record["order"],
                    record["permutation"],
                )
                all_channel_keys.add(channel_key)
        component_maps.append((maps, old_weight))

    all_keys_sorted = sorted(all_channel_keys)
    channel_number = {key: i for i, key in enumerate(all_keys_sorted)}
    rebuilt = {key: [] for key in all_keys_sorted}
    for maps, old_weight in component_maps:
        partner_keys = tuple(sorted(maps, key=channel_number.__getitem__))
        partners = tuple(channel_number[key] for key in partner_keys)
        partner_permutations = tuple(maps[key][2] for key in partner_keys)
        for key, (process, order, permutation) in maps.items():
            factor = old_weight * IdenticalParticleSymmetryFactor(process)
            rebuilt[key].append((
                process,
                order,
                partners,
                factor,
                permutation,
                partner_permutations,
            ))

    phase_space_orders = rebuilt
    build_process_index()


def DetermineMultiChannelPartnersAndSymmetryFactor():
    """Finalize each subprocess record with channels and symmetry factors."""
    BuildTopologyAwareMultiChannels()
    AddResonanceMultiChannels()


def _history_in_density_labels(history, base_to_target):
    """Translate target occurrence labels into one density's base labels."""
    target_to_base = [None] * len(base_to_target)
    for base_label, target_label in enumerate(base_to_target):
        target_to_base[target_label] = base_label
    translated = []
    for resonance, labels in history:
        translated.append((resonance, tuple(sorted(
            target_to_base[label] for label in labels
        ))))
    return _canonical_history(translated)


def AddResonanceMultiChannels():
    """Add automatic resonance densities to every existing MIS component.

    The matrix-element process and colour order are never changed.  Each
    existing density is cloned for the physical histories supported by its
    target subprocess; the stored external-leg map translates occurrence
    labels into that density's canonical labelling.  Every clone remains a
    full-support map and all clones are linked as multichannel partners.
    """
    global phase_space_orders, all_keys_sorted, channel_resonances

    old_keys = tuple(all_keys_sorted)
    old_rows = tuple(tuple(phase_space_orders[key]) for key in old_keys)
    signatures_by_channel = [set() for _ in old_keys]
    for channel, rows in enumerate(old_rows):
        for row in rows:
            process = row[0]
            permutation = row[4]
            signatures_by_channel[channel].add(())
            representative = min(row[2])
            if channel != representative:
                continue
            for history in DiscoverResonanceHistories(process):
                if not history:
                    continue
                signatures_by_channel[channel].add(
                    _history_in_density_labels(history, permutation)
                )

    registry = {}
    new_keys = []
    for channel, old_key in enumerate(old_keys):
        signatures = sorted(signatures_by_channel[channel], key=lambda history: (
            len(history), tuple((len(item[1]), item[1], item[0])
                                for item in history)
        ))
        for history in signatures:
            new_key = old_key + (history,)
            registry[(channel, history)] = new_key
            new_keys.append(new_key)

    channel_number = {key: index for index, key in enumerate(new_keys)}
    rebuilt = {key: [] for key in new_keys}
    for source_channel, rows in enumerate(old_rows):
        for row in rows:
            process, order, old_partners, factor, permutation, \
                old_partner_permutations = row
            physical_histories = DiscoverResonanceHistories(process)
            partner_records = []
            # Keep every original QCD/singlet map as a nonresonant density.
            for partner, partner_permutation in zip(
                    old_partners, old_partner_permutations):
                key = registry[(partner, ())]
                partner_records.append((channel_number[key], partner_permutation))

            # A recursive resonance density does not need to be repeated for
            # every singlet placement of the same coefficient.  Attach one
            # copy of each physical history to the component's deterministic
            # representative; all rows still see it in the same MIS sum.
            representative = min(old_partners)
            representative_index = old_partners.index(representative)
            representative_permutation = old_partner_permutations[
                representative_index
            ]
            for history in physical_histories:
                if not history:
                    continue
                signature = _history_in_density_labels(
                    history, representative_permutation
                )
                try:
                    key = registry[(representative, signature)]
                except KeyError as error:
                    raise ValueError(
                        "resonance history is missing from its multichannel "
                        f"representative: channel={representative} "
                        f"history={signature}"
                    ) from error
                partner_records.append((
                    channel_number[key], representative_permutation
                ))
            partner_records.sort(key=lambda item: (item[0], item[1]))
            partners = tuple(item[0] for item in partner_records)
            partner_permutations = tuple(item[1] for item in partner_records)

            source_histories = [()]
            if source_channel == representative:
                source_histories.extend(
                    history for history in physical_histories if history
                )
            for history in source_histories:
                signature = _history_in_density_labels(history, permutation)
                key = registry[(source_channel, signature)]
                rebuilt[key].append((
                    process,
                    order,
                    partners,
                    factor,
                    permutation,
                    partner_permutations,
                ))

    all_keys_sorted = new_keys
    phase_space_orders = rebuilt
    channel_resonances = {
        key: key[-1] for key in all_keys_sorted
    }
    build_process_index()


def ConvertProcToString(proc):
    """Convert one subprocess record into the ``processes.txt`` row format."""
    process, order, multi_channel, iden = proc[:4]
    if len(proc) >= 5:
        permutation = proc[4]
    else:
        # Preserve the historical helper API for callers that construct a
        # four-field row directly.
        permutation = tuple(range(len(process)))
    if len(proc) >= 6:
        partner_permutations = proc[5]
    else:
        partner_permutations = tuple(
            permutation for _ in multi_channel
        )
    if len(permutation) != len(process):
        raise ValueError("row leg permutation has the wrong size")
    if sorted(permutation) != list(range(len(process))):
        raise ValueError("row leg map is not a permutation")
    if len(partner_permutations) != len(multi_channel):
        raise ValueError("partner permutations are not aligned with channels")
    if any(len(item) != len(process) for item in partner_permutations):
        raise ValueError("partner leg permutation has the wrong size")
    crossed=[pdgs[p] if i>1 else pdgs[anti_particle[p]] for i,p in enumerate(process)] # cross intial state
    line=str(len(multi_channel))
    line=line+'   '+' '.join([str(m+1) for m in multi_channel])
    line=line+'   '+' '.join(crossed)
    line=line+'   '+' '.join([str(o+1) for o in order])
    line=line+'   '+str(iden)
    line=line+'   '+' '.join(str(label+1) for label in permutation)
    line=line+'   '+' '.join(
        str(label+1)
        for partner_permutation in partner_permutations
        for label in partner_permutation
    )
    return line

def sort_by_pdg_codes(process):
    """Return a stable sort key for readable, deterministic output."""
    nq=count_matching_elements(process,quarks)
    same_flavour=max([count_matching_elements(process,[q]) for q in quarks])
    # first sort by 'nq', then by same_flavour, then by (modified) PDG codes:
    return (nq,same_flavour,[sort_particles[p] for p in process])

def sort_by_pdg_codes2(proc):
    process=proc[0]
    return sort_by_pdg_codes(process)

def WriteAllProcsIntoList():
    """Serialize the phase-space groups and subprocess rows.

    This is the second block in ``processes.txt``. Each group starts with the
    group id, number of subprocess rows, maximum number of multichannel
    partners, and phase-space order. The following rows are emitted by
    ``ConvertProcToString``.
    """
    towrite=[]
    towrite.append(str(len(all_keys_sorted))) # number of phase-space orderings to consider
    towrite.append('')
    output_channels=[]
    for internal_key in all_keys_sorted:
        key = PhaseSpaceOrderFromKey(internal_key)
        rows = phase_space_orders[internal_key]
        resonances = channel_resonances.get(internal_key, ())
        output_channels.append((key, resonances, rows))

    for i,(key,resonances,rows) in enumerate(output_channels):
        towrite.append(str(i+1)+'   '+str(len(rows))+'   '+str(max(len(proc[2]) for proc in rows))+'   '+' '.join([str(k+1) for k in key])+'   '+str(len(resonances)))
        for resonance, labels in resonances:
            towrite.append(str(resonance)+'   '+str(len(labels))+'   '+' '.join(str(label+1) for label in labels))
        # order the processes in the process_list, so that we get a neat processes.txt file:
        process_list=sorted(rows,key=sort_by_pdg_codes2)
        for proc in process_list:
            process_line=ConvertProcToString(proc)
            towrite.append(process_line)
        towrite.append('')
        towrite.append('')
        towrite.append('')
    return towrite

def Addqq_dfProcesses(sorted_procs):
    """Include all quark-line permutations for different-flavour processes."""
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

def WriteUniqueProcsIntoList(procs):
    """Serialize the unique-process block at the top of ``processes.txt``."""
    sorted_procs=sorted([sorted(proc,key=lambda x: sort_particles[x]) for proc in procs],key=sort_by_pdg_codes)
    # in case of different flavour multiple-quark line processes, add all the possible orders:
    sorted_procs=Addqq_dfProcesses(sorted_procs)
    try:
        line=[str(len(sorted_procs[0]))+' '+str(len(sorted_procs))+' '+str(PROCESS_FILE_VERSION)]
    except IndexError:
        print("ERROR: no processes found. Try './process_list.py --help' to get more information on usage")
        quit()
    line.append('# process: '+process_provenance)
    line.append('# options: flavour_scheme='+str(flavour_scheme_number)+
                ' include_3qqbar='+str(options['include_3qqbar_processes']).lower()+
                ' include_cc='+str(options['include_cc_processes']).lower()+
                ' resonance_discovery=automatic')
    for proc in sorted_procs:
        line.append(' '.join(pdgs[p] for p in proc))
    line.append('')
    line.append('')
    line.append('')
    return line

def CheckConsistency():
    """Verify that the generated rows contain the expected dual amplitudes."""
    allprocs={}
    for key in all_keys_sorted:
        for row in phase_space_orders[key]:
            process, order, multichannel, iden = row[:4]
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
    """Return the expected number of leading-colour dual amplitudes."""
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
        with multiprocessing.Pool(
                processes=multiprocessing.cpu_count(),
                initializer=SwitchFlavourScheme,
                initargs=(flavour_scheme_number,)) as pool:
            results = pool.map(ProcessProcess, all_procs)  # Parallelize across procs
    else:
        results=[ProcessProcess(x) for x in all_procs]
    phase_space_orders=CombineResults(results)
    all_keys_sorted=sorted(phase_space_orders.keys())
    DetermineMultiChannelPartnersAndSymmetryFactor() # updates the phase_space_orders dictionary
    # Check the consistency of the generated processes
    CheckConsistency()
    # write to disk
    towriteunique=WriteUniqueProcsIntoList(all_unique_procs)
    towriteallprocs=WriteAllProcsIntoList() # puts the phase_space_orders dictionary in a writable list
    
    with open('processes.txt','w') as f:
        f.write('\n'.join(towriteunique))
        f.write('\n'.join(towriteallprocs))
        f.write('\n')
