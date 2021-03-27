
EXE=../../local/MFP-PHM.1d.gnu.DEBUG.MPI.CVODE.PARTICLES.HDF5.ex

INPUTS=../Collisions/Collisions.inputs

rm *.hdf5

valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --undef-value-errors=no --suppressions=suppress  ${EXE} ${INPUTS} 2>&1 | tee run_log.txt
