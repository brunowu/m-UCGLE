### Unite and Conquer Multiple GMRES/ERAM LS method implemented with MPI Spawn

Xinzhe WU @ Maison de la Simulation (xinzhe.wu@cea.fr)

Another similar implementation of [UCGLE](https://github.com/brunowu/UCGLE). The previous UCGLE uses a static MPI distributed communication implementation, thus the GMRES and ERAM Component numbers are fixed.

Here, we design a new dynamic distributed communication workflow by MPI_SPWAN. The Spawn functionality allows allocating multiple ERAM and GMRES Components. The multiple ERAM trends to be a [MERAM](https://epubs.siam.org/doi/10.1137/S1064827500366082)(Multiple Explicitly Restarted Arnoldi Method), the exchange of spectral information among the ERAM Components will improve the convergence of resolving no-Hermitian eigenvalue problems.

Meanwhile, multiple GMRES allows to solve simultaneously non-Hermitian linear systems with multiple right-hand sides. It is comparable with deflated [Block GMRES](http://www.sam.math.ethz.ch/~mhg/pub/delhipap.pdf), but with less global communication.

This software is underdevelopment by [Xinzhe Wu](https://brunowu.github.io/)...

#### To run

```bash
cmake .
make
mpirun -np 1 ./test.exe -nb_gmres 1 -proc_gmres 1 -nb_arnoldi 1 -proc_arnoldi 1 -gmres_exec ./gmres.exe -arnoldi_exec ./arnoldi.exe -lsqr_exec ./lsqr.exe --filename=\"utm300.mtx\" --eps-quiet --eps-nodebug --eps-exsitu --eps-sort=\"LM\" \ --eps-nonherm --eps-nev=10 --eps-blockSize=1 --eps-tol=0.1 --eps-no-print-matrix --eps-all-print --eps-numBlocks=20 --eps-maxRestarts=50 --ksp-nodebug --ksp-frequency=1 --ksp-tol=1e-05 --ksp-num-rhs=2 --ksp-block-size=2 --ksp-no-precond --ksp-num-blocks=120 --ksp-fixed-tol --ksp-no-print-matrix --ksp-all-print --ksp-no-dump-data --ksp-lsp-degree=2 --ksp-lsp-latency=1 --ksp-use-lsp > test.txt
```

##### Workflow

![Workflow of UCMGEL](workflow.jpg)

