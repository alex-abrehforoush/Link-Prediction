# PyCCAlg - Correlation-Clustering Algorithms in Python
`Python` implementation of some algorithms for `Correlation Clustering`. Specifically:
* Linear-programming + region-growing O(log n)-approximation algorithms for general weighted graphs
   * `round_demaine` in `src/pyccalg.py`: [Demaine et al.'s](https://www.sciencedirect.com/science/article/pii/S0304397506003227) rounding algorithm
   * `round_charikar` in `src/pyccalg.py`: [Charikar et al.'s](https://www.sciencedirect.com/science/article/pii/S0022000004001424) rounding algorithm
* `kwikcluster` in `src/pyccalg.py`: `KwikCluster` randomized, linear-time algorithm ([Ailon et al., JACM 2008](https://doi.org/10.1145/1411509.1411513)), achieving constant-factor approximation guarantees on complete graphs satisfying certain constraints (e.g., probability constraint and/or triangle-inequality constraint)



### Requirements:

* `Python v3.6+`
* For linear-programming-based algorithms:
   * [`SciPy`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)`v1.6+` and/or [`PuLP`](https://pypi.org/project/PuLP/)
   * `SciPy linprog` comes with various solvers: '*Method `highs-ds` is a wrapper of the C++ high performance dual revised simplex implementation (HSOL). Method `highs-ipm` is a wrapper of a C++ implementation of an interior-point method; it features a crossover routine, so it is as accurate as a simplex solver. Method `highs` chooses between the two automatically. For new code involving linprog, we recommend explicitly choosing one of these three method values instead of `interior-point` (default), `revised simplex`, and `simplex` (legacy)*'. See [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html) for more details.
   * `PuLP` comes with two solvers by default: [`CBC`](https://projects.coin-or.org/Cbc) (linear and integer programming) and [`CHOCO`](https://choco-solver.org/) (constraint programming), but it can connect to many others (e.g., `GUROBI`, `CPLEX`, `SCIP`, `MIPCL`, `XPRESS`, `GLPK9`) if you have them installed
   * Here we use `highs-ipm` with `SciPy linprog` and the default `CBC` with `PuLP` 
   * However, any linear-programming `Python`  (other than `SciPy linprog` or `PuLP`) library can alternatively be used with minimal adaption

### Usage:

``` python src/pyccalg.py -d <DATASET_FILE> [-r <LB,UB>] [-a <PROB>] [-s {'pulp','scipy'}] [-m {'charikar','demaine','kwik'}]```

* Optional arguments: 
   * `-r <LB,UB>`, if you want to generate random edge weights from `[LB,UB]` range
   * `-a <PROB>`, if you want to randomly add edges with probability `PROB`
   * `-m {'charikar','demaine','kwik'}`, to choose the algorithm (default: `'charikar'`). NOTE: `kwikcluster` is always run too
   * `-s {'pulp','scipy'}`, to select the solver to be used (default: `'scipy'` (it seems faster))
* Dataset-file format:
   * First line: `#VERTICES \t #EDGES`
   * One line per edge; every line is a quadruple: `NODE1 \t NODE2 \t POSITIVE_WEIGHT \t NEGATIVE_WEIGHT` (`POSITIVE_WEIGHT` and `NEGATIVE_WEIGHT` are ignored if code is run with `-r` option)
   * Look at `data` folder for some examples

