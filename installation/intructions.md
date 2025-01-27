
Mark Munyi - August 2024

This is the procedure for installling and compiling `Quantum Package 2`, on Nots cluster:
Start by cloning the clone with https://github.com/QuantumPackage/qp2.git.
`cd` into the qp2 directory, then run:

`module purge`
`module load gomkl/2023a Python Python-bundle-PyPI SciPy-bundle GMP HDF5 Ninja`
`nodule load ZeroMQ `

This loads Ninja, python among other related dependances like HDF5 for trexio.

You'll need some python packages:

`pip install --user resultsFile`
`pip install --user bse`

Then you'll have to build the remaining external dependencies:

`./configure -i trexio `
`./configure -i f77zmq`
`./configure -i ocaml`

finally, load the compiler. This will depend on the cluster you are using, i.e Nots Vs Notsx, but gfortran should generally work:

`./configure -c config/gfortran_mkl.cfg`

run module save to save these modules and ensure than you dont have to load them all once again.

then run 
`make`

to compile the code.

youll need to activate your python environment once in the qp2 directory with

`conda activate qp2.`

Incase you make any changes to the code, i.e editing the .cfg files, you'll need to run make again.

:::::::::::::::::::::::::::::::
