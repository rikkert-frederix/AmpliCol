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
  --process=[X],    -p=[X]  : Process specified in file [X] (default is './processes.txt').
  --nevents=[X],    -n=[X]  : Number of unweighted events to generate (default is 10000).
  --phasespace=[X], -ps=[X] : Phase-space parametrisation to use -- 1=gen23 (default), 2=HAAG, 3=pT-based, 4=t-channel.
  --seed=[X],       -s=[X]  : The random number seed to use (default is read from randinit file).
  --itmax=[X],      -i=[X]  : The maximum number of iterations to use (detault is 128).
  --library=[X],    -l=[X]  : To create or use a library for the amplitudes, set [X] to 'create' or 'use', respectively. (To use a library, re-compile code with 'make amplicol_generate_library' after a library has been created). Default is 'none'.
  --tag=[X],        -t=[X]  : Event file (and log file) names will be prepended with a tag '[X]_'.
  --me_test=[X],    -mt=[X] : Perform ME-level test against MG with [X] points tested (single PS kinematics)
  --gauge=[X]               : Vector-current gauge: unitary (default) or fd.
*****************************

The optional FD gauge uses the five-component massive-vector
representation (four reduced-vector components plus the associated
Goldstone component) and the momentum-dependent gauge direction
n(q)=(sgn(q0),-q_vec/|q_vec|).  It follows the formulation in
arXiv:2003.03003, arXiv:2203.10440, and arXiv:2405.01256.  The default
'unitary' value preserves the pre-existing implementation.

Most important are the '--nevents=X' to set the number of events to
generate, '--seed=X' to set the random seed.

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

An amplitude library is gauge-specific.  Use the same '--gauge=X'
option when creating and using it; a mismatched or old library is
rejected and must be recreated.

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
  --gauge=[X]       : Vector-current gauge: unitary (default) or fd.
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


4. FD-gauge regression tests

The component checks for the five-component propagator, external
states, and auxiliary-particle decompositions are run with

make test_fd_gauge

The full matrix-element fixture is generated locally because its
golden file is large:

make update_matrix_goldens
make test_matrix_elements

The latter command checks both the default unitary gauge and the FD
gauge against the same golden.dat file. These gauge-equivalence tests
use zero particle widths: fixed widths violate Ward identities and
would introduce a gauge difference unrelated to the FD vertices being
tested.



NOTE: both the leading-colour as well as the reweighted events are not
completely randomised when they are written to file! (Some flavour
configurations are over-represented in the beginning of the file;
others towards the end of the file).  Hence, when only a subset of all
events are passed to the analysis, the events should first be
randomised to have an unbiased sample!
