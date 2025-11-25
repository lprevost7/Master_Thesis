/**
 *  multiobjective.h
 *  ---------------------------------------------------------------
 *  Public interface for multi-objective search components (TPLS / PLS / TPPLS)
 *  tailored to the Permutation Flow Shop Problem (PFSP) with two objectives.
 *
 *  Key concepts:
 *  - TPLS (Time-Partitioned Local Search): runs several local-search sweeps
 *    while varying a scalarization weight lambda in [0,1].
 *  - PLS (Pareto Local Search): explores neighborhoods of a non-dominated
 *    set, adding neighbors that improve at least one objective.
 *  - TPPLS: TPLS to generate a seed Pareto set, then PLS expansion.
 *
 *  Notes:
 *  - Most functions clone solutions; the caller (or orchestrator) is
 *    responsible for deallocating returned pointers when appropriate.
 *  - If the PLS time limit is not provided (<= 0), PLS continues until
 *    there are no unexplored solutions remaining in the Pareto set.
 */
#ifndef MULTIOBJECTIVE_H
#define MULTIOBJECTIVE_H

#include <vector>
#include <ctime>

#include "emilibase.h"
#include "pfsp/permutationflowshop.h"

/**
 * ExploredSolution
 * ---------------------------------------------------------------
 * Lightweight wrapper carrying a solution pointer plus a flag indicating
 * whether the solution's neighborhood has already been explored in PLS.
 *
 * Ownership: When inserted into collections managed by PLS utilities,
 * the pointed Solution is typically owned and must be deleted on cleanup.
 */
class ExploredSolution
{
public:
    explicit ExploredSolution(emili::Solution* solution) : Solution(solution) {}

    emili::Solution* Solution; // Raw pointer to a PFSP solution (cloned elsewhere)
    bool explored = false;     // True once its neighborhood has been visited
};

/** Writes the Pareto points (X,Y) to CSV at ../graph.csv with a Color column. */
void writeResultsToFile(const std::vector<int>& ObjValX, const std::vector<int>& ObjValY);

/**
 * Generates the next lambda in [0,1] by interval bisection.
 * First call returns 1.0, then midpoints to spread samples across the range.
 */
double nextLambda();

/** Sorts solutions ascending by the first objective, required by Pareto sweep. */
void sortSolutionsForPLS(std::vector<emili::Solution*>& solVector, emili::pfsp::PFSP_TPLS* problem);

/**
 * Extracts a Pareto front from a list sorted by objective1 (X), sweeping on
 * objective2 (Y). When multiple solutions share X, keeps only the best Y.
 * Returns a vector of cloned solutions (caller must delete them).
 */
std::vector<emili::Solution*> generateParetoFront(std::vector<emili::Solution*>& solVector, emili::pfsp::PFSP_TPLS* problem);

/**
 * Filters a list of ExploredSolution* (sorted by X) to a non-dominated set
 * with respect to Y, counting how many entries remain unexplored.
 * Returns a new vector that takes ownership of the kept entries; the input
 * vector is cleared.
 */
std::vector<ExploredSolution*> generateParetoFrontExplored(std::vector<ExploredSolution*>& Aexpl, emili::pfsp::PFSP_TPLS* problem, int& soltoExplore);

/**
 * Runs n TPLS iterations, collecting cloned solutions into solVector and
 * updating the PFSP_TPLS lambda according to the chosen policy.
 * Returns 0 on success, non-zero on error.
 */
/**
 * Runs up to n TPLS iterations, respecting a global time budget.
 * - startTime: common program start (used to compute elapsed time)
 * - timeLimitSeconds: if > 0, total global time budget; if <= 0, no time cap
 * Returns 0 on success.
 */
int tplsPhase(emili::LocalSearch* ls, clock_t startTime, double timeLimitSeconds, int n, std::vector<emili::Solution*>& solVector);

/**
 * PLS exploration phase. If timeLimitSeconds <= 0, continues until the Pareto
 * set has no unexplored solutions; otherwise runs until time is exhausted or
 * all solutions are explored. Returns 0 on success.
 */
int plsPhase(clock_t startTime, double timeLimitSeconds, std::vector<ExploredSolution*>& Aexpl, emili::pfsp::PFSP_TPLS* problem, emili::LocalSearch* ls);

/**
 * Orchestrates TPLS and optional PLS and writes results to CSV.
 * The flags select which phases to execute.
 */
void executeMultiObjective(emili::LocalSearch* ls, float pls, clock_t startTime, bool is_tpls, bool is_pls, bool is_tppls);

#endif
