	 ______ __  __ _____ _      _____ 
	|  ____|  \/  |_   _| |    |_   _|
	| |__  | \  / | | | | |      | |  
	|  __| | |\/| | | | | |      | |  
	| |____| |  | |_| |_| |____ _| |_ 
	|______|_|  |_|_____|______|_____|

emili — Multi-objective search (TPLS / PLS / TPPLS)
===================================================

This document describes the multi-objective components implemented in `multiobjective.cpp` and how to run them from the command line. It follows the same style and rigor as the main README.

Overview
--------

`multiobjective.cpp` implements bi-objective search procedures for the permutation flow shop problem (PFSP) with two objectives, typically Makespan (`PFSP_MS`) and Total Completion Time (`PFSP_TCT`). It provides:

- Time-based local search sweeps (TPLS) with a weight λ in [0, 1]
- Pareto Local Search (PLS) starting from a non-dominated set
- A combined TPPLS flow (TPLS followed by PLS)
- Export of the final Pareto set as a CSV (`graph.csv`) for plotting

New in this variant
-------------------

- Clean file: removed obsolete debug prints and dead code; detailed English comments and function docstrings were added.
- Conditional debug logging: enable logs by compiling with the `ENABLE_MO_DEBUG` macro.
- Safer memory management: cloned solutions are freed at the end of `executeMultiObjective`.
- Deterministic λ generation: `nextLambda()` explores [0, 1] by interval bisection.
- PLS without a time cap: if no time is provided (or `-it 0`), PLS continues until there are no unexplored solutions left in the current Pareto set.

Usage
-----

General syntax (see also the main README):

```
emili INSTANCE_FILE_PATH PROBLEM < ALGORITHM DESCRIPTION > [-it|-ro time] [rnds seed] [ps]
```

For bi-objective PFSP runs, use `PFSP_PLS` followed by the two objectives, then your algorithm template. A typical TPLS+PLS (TPPLS) run with a global time budget:

```
./emili ../instancesDD_for_experiments/100x20_1.txt PFSP_PLS PFSP_MS PFSP_TCT first neh locmin insert -step 10 -it 200
```
```
./emili ../Benchmarks/ta051.txt PFSP_PLS PFSP_MS PFSP_TCT ga 1000 0.8 0.1 random time 60 ga_tournament 2 ga_ox ga_swapmut -step 10 -it 200
```


- `PFSP_PLS PFSP_MS PFSP_TCT` selects the bi-objective PFSP with Makespan and Total Completion Time.
- `first neh locmin insert` is an example algorithm template (initialization and neighborhood).
- `-step 10` sets the number of TPLS steps (the code uses `getTPLSStep()+1`).
- `-it 200` sets a **global** time budget (seconds) shared by TPLS and PLS phases. TPLS consumes time first; PLS uses whatever remains.

### Global time budget semantics (`-it`)

The flag `-it <seconds>` now applies to the **entire multi-objective run**, not only to PLS. Concretely:

1. The TPLS phase iteratively runs local searches while checking elapsed global time.
2. When TPLS finishes (either all steps executed or time exhausted), the initial Pareto front is built.
3. The PLS phase starts only if time remains; otherwise it is skipped with a message.
4. The PLS phase itself also respects the global start time and stops when the total elapsed time reaches the `-it` limit.

Implications:

- When `-it` is small, you may see only TPLS solutions (PLS skipped).
- Increasing `-step` increases TPLS work and can reduce time left for PLS; tune both together.
- If you need to guarantee some PLS exploration, either reduce `-step` or raise `-it`.

No time limit (unbounded PLS exploration)
-----------------------------------------

To run without a global time cap, set `-it 0` (or omit it if the default is zero in your build). In this mode:

```
./emili ../instancesDD_for_experiments/100x20_1.txt PFSP_PLS PFSP_MS PFSP_TCT first neh locmin insert -step 10 -it 0
```

- TPLS runs for the configured number of steps (not time-limited).
- PLS then continues until every Pareto solution discovered has been explored (including new ones added during the process).

Note: Even in unlimited mode, TPLS still stops after its configured number of steps; there is intentionally no infinite TPLS loop.

Output
------

At the end of the run, the set of non-dominated objective pairs is written to `graph.csv` with header `X,Y,Color` and RGB tuples (currently fixed to red):

```
X,Y,Color
123,456,"(1,0,0)"
...
```

Note: When executing from the build directory (`emili/build`), the file is written to `../graph.csv`, i.e., into the project root `emili/graph.csv`.

Enabling debug logs
-------------------

Debug logs are disabled by default. To enable them, define `ENABLE_MO_DEBUG` at compile time. With CMake, you can reconfigure with a custom CXX flag, for example:

```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-DENABLE_MO_DEBUG" ../
make -j
```

This will activate concise progress logs in the multi-objective phases.

Build notes
-----------

Use CMake as in the main README:

```
mkdir -p build
cd build
cmake ../
make -j
```

Tips
----

- The λ schedule in TPLS is either any-time (via `nextLambda()`) or a uniform sweep in [0, 1] across steps.
- PLS adds neighbors that improve at least one objective and maintains a Pareto set sorted by the first objective for efficient filtering.
- If you plot `graph.csv`, you should see a non-increasing trade-off curve (in Y) as X increases.

License
-------

This module is part of `emili` and inherits the same BSD 2-Clause License. See the main README for details and citation guidelines.
