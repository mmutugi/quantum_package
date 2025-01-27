'''Instruction to load compile Quantum Package 2.2 on Expanse
    Mark Munyi
'''

1. Clone the repository: https://github.com/QuantumPackage/qp2.git. Check for the updated one.

2. `cd` into the `qp2` directory and run `./configure`
    This runs some compatibility and dependency checks and then returns a `configure` error that `Ninja` is not installed.

3. Run `./configure -i ninja`(ninja with small `n`) -> This takes a minute or so and complains that `ZeroMQ` is not installed.

4. Run `./configure -i zeromq` -> This takes ~3 mins and complains that `TREXIO (trexio | trexio-nohdf5) is not installed`.

5. Run `./configure -i trexio-nohdf5` -> This takes like 10 mins and then complains `Fortran binding of ZeroMQ (f77zmq) is not installed`

6. Run `./configure -i f77zmq` -> This also takes about 5 mins and the complains `OCaml (ocaml) compiler is not installed`. At this point, 

7. Run `./configure -i gmp` -> This is quick. Afterwards, 

7. Run `./configure -i ocaml` -> This takes about 1 hour and then complains that `docopt (docopt) is not installed`.

8. Run `./configure -i docopt` -> This is also fairly quick and then complains that `resultsFile (resultsFile) is not installed`. 

9. Run `./configure -i resultsFile` -> This is also quick and then returns cowsay that `All dependencies installed` but that `/path/to/qp2/build.ninja does not exist`.

10. Run `module load openblas` -> Do not specify amything else (expanse does have different variants of Openblas.)\

11. Run `./configure -c config/gfortran_openblas.cfg` and then 

12. Run `make`-> This takes about 2-hours.