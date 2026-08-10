********************************
***        AmpliCol          ***
********************************
*   R. Frederix and T. Vitos   *
* Please cite arXiv:2601.19483 *
********************************

There are three steps in running the code: 
1. process generation, 
2. event generation,
3. reweighting.


1. Process generation

The goal is to create a list of subprocesses, integration channels, and
multiplication factors that the integrator code can use to generate
events. The code is a python script, 'process_list.py', that takes a
process as an argument, and creates a 'processes.txt' file with the
required information. Usage of the code is as follows:

*****************************
$ ./process_list.py -h
usage: process_list.py [-h] [-FS [1-5]] [-3] [-s] [-cc] [-res] process_string

Generate the full list of processes, ordered by phase-space order

positional arguments:
  process_string        Process to consider (e.g., 'p p > w+ z 4j', (including quotation marks))

options:
  -h, --help            show this help message and exit
  -FS [1-5], --flavour_scheme [1-5]
                        Switch to N-flavor scheme (NFS), where N is 1-5 (default=5)
  -3, --include_3qqbar  Include processes with up to 3 quark lines
  -s, --serial          Do not use multi-processes (parallel execution). Useful for debugging.
  -cc, --include_cc     Include flavour-changing processes
  -res, --resonance     Treat the lepton pair as a resonance
*****************************

Note that, you should NOT use 'p p > j j j' as process_string for
three jet production. The correct syntax is 'p p > 3j'. By default it
will use a 5 Flavour scheme, (both proton 'p' and jet 'j' will include
the 5 lightest quarks (and the gluon)). Furthermore, flavour
configurations with three or more quark lines are not included by
default, but can be included by providing the '--include_3qqbar'
option. If these are included, the performance of the event generation
is significantly reduced. Furthermore, these contributions will only
be calculated at leading-ccolour accuracy (i.e., the reweighting code,
see below, will skip these contributions).



2. Event generation

Using the 'processes.txt' file, the 'amplicol_generate' code will
generate events at leading-colour accuracy for that process. The
Amplicol code is a fortran code that can be compiled with

make amplicol_generate

It requires 'LHAPDF' as an external dependency, for which it is
expected that the correct link commands are provided through
`lhapdf-config --ldflags`, see the makefile lines 98-99.

There is currently no "run_card" or "process_card" to set input
parameters. All parameteres are hard-coded in the fortran code
---mostly in 'common.f03' (e.g., collision energy, PDFset to use,
generation cuts, and renormalisation/factorisation scale choice are
all in common.f03). Changing any of these parameters requires
recompiling amplicol_generate.

The possible arguments when running the event generation code are:

*****************************
./amplicol_generate -h

Usage: amplicol_generate <arguments>'. Possible arguments are

  --help,           -h      : Show this message.
  --process=[X],    -p=[X]  : Born process specified in file [X] (default is './processes.txt').
  --real-process=[X]         : Real-emission process file; integrates B + R-sum(D) in one run.
  --dim-reg=[hv|fdh]         : Dimensional scheme for integrated dipoles (default: hv).
  --nevents=[X],    -n=[X]  : Number of unweighted events to generate (default is 10000).
  --phasespace=[X], -ps=[X] : Phase-space parametrisation to use -- 1=gen23 (default), 2=HAAG, 3=pT-based, 4=t-channel.
  --seed=[X],       -s=[X]  : The random number seed to use (default is read from randinit file).
  --itmax=[X],      -i=[X]  : The maximum number of iterations to use (detault is 128).
  --library=[X],    -l=[X]  : To create or use a library for the amplitudes, set [X] to 'create' or 'use', respectively. (To use a library, re-compile code with 'make amplicol_generate_library' after a library has been created). Default is 'none'.
  --tag=[X],        -t=[X]  : Event file (and log file) names will be prepended with a tag '[X]_'.
  --me_test=[X],    -mt=[X] : Perform ME-level test against MG with [X] points tested (single PS kinematics)
  --limit_test,     -lt     : Test soft and collinear Catani-Seymour limits and exit.
  --tail-replay=[FILE]      : Re-evaluate the saved maximum-weight real-subtraction point and exit.
*****************************

Most important are the '--nevents=X' to set the number of events to
generate, '--seed=X' to set the random seed.

The limit diagnostic can be run with '--limit_test'. It tests each channel at
100 generated phase-space points. A limit passes at a given channel and limit
when fewer than 5 of the tested points fail, and the screen reports the
failure fraction separately for every integral (`iint`). At an individual
point, a limit passes when three consecutive finite
deformation points have a positive matrix-element/dipole ratio within 1% of
one and vary by less than 2%. If roundoff spoils this window, a two-point
convergent tail or a local minimum within 1%, with both neighbours within 5%,
is accepted. Soft non-gluons are checked for bounded `lambda * matrix element`,
which permits integrable `1/lambda` growth. Collinear pairs without a matching
Catani-Seymour dipole are checked for bounded matrix elements. In both cases,
three consecutive finite samples must not grow by more than 5% as the limit is
approached. These checks are included in the reported failure fractions and
are labelled '(integrable)' or '(finite)' in the summary.
Soft limits for massive particles, and collinear limits involving at least one
massive particle, are skipped.
The base point for every test is regenerated until it has positive phase-space
weight and passes the built-in cuts before the limiting deformation starts.
For failed limits, the complete deformation table is written to
'Outputs/limit_test_failures.log' (or the corresponding tagged filename).

For real Catani-Seymour dipoles, the optional '--alpha=X' argument restricts the
dipole phase space. A single value applies to all topologies; alternatively use
four comma-separated values in the order FF, FI, IF, II. Values must satisfy
0 < X <= 1. The default is one. In a combined Born-plus-real run, the matching topology-specific
endpoint corrections are applied to the integrated terms as well.

The optional '--real-process=FILE' mode loads the Born processes from
'--process=FILE' and real-emission processes from the additional file into one
integration. The Born channels and the local real contribution `R - sum(D)`
are sampled by the same multichannel integrator, even when their phase spaces
have different dimensions. The local contribution is split point by point into
two exact, disjoint strata. The regular stratum contains points for which the
real event and every alpha-active mapped dipole have the same cut decision. The
migration stratum contains points for which at least one active mapped decision
differs from the real decision. Each real subprocess consequently has separate
regular and migration estimators and point quotas. Within a physical channel,
all regular leaves share one adaptive grid and all migration leaves share a
second, independent grid. Their sum is identically the unsplit local
contribution. The final report prints both strata and checks this closure.

For jet-cut migrations, the diagnostics record the signed distance of the real
and mapped configurations from the relevant jet-pT threshold: the required
Nth-hardest accepted jet pT minus `pTj_min`. When the real and mapped pT
decisions cross that threshold, the smallest absolute distance involved is
reported in the tail record; migrations caused only by another cut use -1.
This derived distance is not an extra integration coordinate; the separate
migration grid allows the existing phase-space coordinates to adapt to that
boundary without biasing the integral.

Each Born subprocess supplies an endpoint integral and one P+K convolution
integral for each incoming beam, all in that same run. It requires
'--accuracy=X' and does not produce unweighted events. The final report lists
Born, the total/regular/migration R-sum(D), the three Laurent coefficients of I,
P, K, and the finite subtotal B+(R-D)+I0+P+K. The virtual term is not included.
Accuracy-only allocation reserves 25% of every post-pilot iteration for uniform
exploration across channel/integral leaves; the remaining points follow the
measured variances. A leaf quota may grow by at most a factor of 16 per
iteration. These safeguards affect sampling only and never clip or modify an
integration weight.

The requested-accuracy stop also requires the largest single migration point's
current variance proxy to be at most 20% of the total migration variance proxy.
Set another fraction with '--migration-tail-fraction=X'; zero disables this
additional convergence gate. A failed gate forces another grid refinement and
doubles the next global point budget, subject to the per-leaf growth cap above.
The user-specified iteration cap is never overridden, and a warning is printed
if it is reached while the gate is still unsatisfied.

Combined runs retain separate top-eight lists for the signed real-subtraction
residual and the individual R/D counterevent envelope in
`Outputs/<tag>_tail_diagnostics.log`. Each record contains the random
coordinates, momenta, real matrix element, individual dipoles, alpha cut
variables, mapped-cut decisions, PDF/phase-space factors, and its contribution
to the second moment. Compact fixtures for the global component and residual
maxima are written to `Outputs/<tag>_tail_replay.dat` and
`Outputs/<tag>_tail_residual_replay.dat`. Version-3 fixtures retain the sampled
regular/migration leaf; version-2 fixtures from unsplit runs remain accepted.
Re-run with the original process, phase-space and alpha arguments plus
`--tail-replay=FILE` to reproduce either point; each fixture retains its
historical multichannel sampling factors, and the executable verifies the
signed weight before exiting.

Integrated histories are constructed exclusively from the local
dipoles and are checked against the supplied Born processes and leading-colour
orders. The supported dimensional schemes are HV (default) and FDH, selected
with '--dim-reg=hv' or '--dim-reg=fdh'.
Process files generated by `process_list.py` record their active flavour
scheme (`-FS 1` through `-FS 5`) in a header. A combined run requires this
header in both files and requires the two schemes to agree; regenerate legacy
files with `process_list.py`. The P/K and massive IF kernels use that declared
scheme. Use an LHAPDF set with the same flavour convention. Without LHAPDF,
combined runs are supported only for the default five-flavour scheme.
Endpoint I histories still come only from local real dipoles, so deliberately
omitted higher-quark-line subprocesses are not restored by the integrated
terms.
The integrated kernels require a massless unresolved parton. In addition to
the fully massless dipoles, they support a massive final-state quark emitter
for Q -> Qg in FF and FI topologies, a massive final-state spectator in FF and
IF topologies, and an FF dipole in which both the quark emitter and spectator
are massive. Massive unresolved partons and massive initial-state legs are not
supported; a combined run stops during initialisation if such a history is
encountered. The massive finite terms include the topology-specific alpha
restriction and HV/FDH scheme terms. The current executable chooses equal
renormalisation and factorisation scales, but the integrated equations retain
their distinct mu_R and mu_F dependence (including the massive-IF finite
scale-ratio term), so a later interface can vary them independently.
Identical-flavour copies are expanded explicitly, and reduced histories are
canonicalised onto the supplied Born leg and colour-order conventions before
duplicate insertions are removed. At initialisation, purely massless
physical-parent histories are checked to reconstruct the universal quark and
gluon poles for every light-flavour sector present in the supplied
real-process file.
Initial-state flavour-changing kernels are normalised to the spin- and
colour-averaged reduced Born contribution. In particular, a real incoming
quark that reduces to a gluon includes the quark/gluon averaging ratio before
the ordered-history colour weight is applied; auxiliary U(1)-parent beam
convolutions carry one additional factor 1/Nc.
 Its real-emission measurement clusters massless QCD partons with inclusive
kT and radius `DRjj_min`, requiring at least one fewer accepted jets than the
number of real final-state jet legs. Each mapped dipole is cut independently
at parton level: every reduced final-state jet must pass `pTj_min`, `etaj_max`,
and `DRjj_min`.

For the number of events, it is best to require somewhere between
100000-1000000 events per run. Requiring too few, and then combining
many separate runs together might undersample some phase-space regions
(there is no 'grid-pack' mode), while requesting too many events in
one go might require too much memory, since a large fraction of the
generated events will be kept in memory. Indeed, in general, the code
has not been optimised to be memory efficient---one might need several
GB's of RAM (or even more) for process with ~6 final state particles.

Running the code will generate an event file, 'events.lhe', with
events at leading-colour accuracy in the ./Outputs/ directory. An
extensive logfile is provided in the same location. Using the
'--tag=X' option will prepent the event file name and log file name
with the tag 'X_'.

Finally, it might be useful to use the '--library=X' option. This can
speed up the event generation at the cost of a longer compilation
time. Running the amplicol_generate code with '--library=create' will
not generate events, but rather the code will create a bunch of
fortran source files (in the ./Library/ directory) with all the matrix
elements relevant for this processes written to disk. These can then
be compiled (with optimisation flags) into a set of process libraries
that can be linked to the amplicol_generate code, using 'make
amplicol_generate_library'; you might want to use '-jX', with X the
number of CPU cores available to speed up the compilation. The
resulting code can be significantly faster, but compilation time is
non-negligible, so only relevant when generating many events. The
usage in this case would be:

make amplicol_generate
./amplicol_generate --library=create --process=processes.txt
make amplicol_generate_library
./amplicol_generate --library=use --nevents=1000000 --seed=101

It might be needed to include the "-mcmodel=large" compilation flag
when the process library is large, see also lines 5-6 of the
makefile. Using this option requires a 'make clean' first, before
recompiling the amplicol_generate code. (This also means re-creating
the process library, since a make clean removes all of the process
library source code).


3. Reweighting

Once the leading-colour accurate events are generated (by default in
the ./Outputs/events.lhe file), these events can be reweighted to
full-colour accuracy by the amplicol_reweight code, compiled with

make amplicol_reweight

Running the code, requires the leading-colour accurate event file as
an argument. The output is an event file with the same name, apart
from '.rwgt' appended. This LHEF contains the events at full colour
accuracy. Hence, by default

./amplicol_reweight Outputs/events.lhe

will generate the full-colour accurate events in the file
Outputs/events.lhe.rwgt. Possible options are

*****************************
./amplicol_reweight -h

Usage: 'amplicol_reweight <event_file> <arguments>', where <event_file> is the leading colour event file to reweight to full colour.
The code creates an LHEF, '<event_file>.rwgt' containing the full colour events.
Possible arguments are

  --help,   -h      : Show this message.
  --unwgt           : Unweight the reweight events.
  --remove_comments : Remove comment lines in the final LHEF.
*****************************

By default, the code reweights all events and update the event weights
to include the full-colour information. The weights are such that they
average to the total cross section. Note that this means that it will
produce *weighted* events. These weights need to be kept and passed
through to the analysis code. If, instead the '--unwgt' option is
give, a (secondary) unweighting step is performed and the events are
given again unit weights. This will reduce the number of events in the
event file ---typically by 30-50% or so. In practice, it is generally
more efficient to use the weighted events (since the Kish effective
sample size ratio is typically well over 98%), except when
post-processing (parton showering, dectector simulation, etc.) is
extremely time-consuming.

By default, the produced LHEF also contains some lines starting with
'#' that contain some information for debugging. The writing of these
lines is skipped by adding the '--remove_comments' option to the
execution of the amplicol_reweight code.



NOTE: both the leading-colour as well as the reweighted events are not
completely randomised when they are written to file! (Some flavour
configurations are over-represented in the beginning of the file;
others towards the end of the file).  Hence, when only a subset of all
events are passed to the analysis, the events should first be
randomised to have an unbiased sample!
