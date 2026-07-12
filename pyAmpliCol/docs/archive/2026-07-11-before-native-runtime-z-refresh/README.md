# Pre-native-runtime Z-table refresh

This snapshot was taken after completing the first n=9 O1/ASM/O3 refresh and
before rebuilding Rusticol with the redundant momentum-copy and direct NumPy
output-transfer improvements. The n=1 through n=8 Z timings predate that
runtime build; n=9 was generated with the current shared-LC generator but uses
the previous `time-process` wall-timing path.

The live Z table is subsequently regenerated in increasing multiplicity. All
process artifacts remain under `docs/.z_performance_outputs` and are preserved
with timestamped names when replaced.
