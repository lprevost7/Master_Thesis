/**
 *  multiobjective.cpp
 *  ---------------------------------------------------------------
 *  Implements multi-objective search phases (TPLS / PLS / TPPLS) for
 *  Permutation Flow Shop Problems (PFSP) with two objective functions.
 *
 *  Main responsibilities:
 *   - Execute time-based local search (TPLS) iterations while adjusting a
 *     weight (lambda) for the aggregated objective.
 *   - Generate Pareto fronts of solutions for bi-objective evaluation.
 *   - Perform a Pareto Local Search (PLS) exploration over neighborhoods of
 *     non-dominated solutions.
 *   - Persist the final set of objective value pairs to a CSV file for
 *     visualization.
 *
 *  Debug output can be enabled by defining ENABLE_MO_DEBUG at compile time.
 *  (e.g. add -DENABLE_MO_DEBUG to your compiler flags.)
 *
 *  NOTE: Memory ownership – many functions clone solutions. The caller is
 *  responsible for eventual deallocation of returned solution vectors.
 */

#include "multiobjective.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <utility>
#include <ctime>
#include <limits>

#ifndef restrict
#define restrict __restrict__
#endif

#include "moocore/c/hv.h"
#include "moocore/c/config.h"  // for dimension_t

// Conditional logging macro (silent by default unless ENABLE_MO_DEBUG defined)
#ifndef MO_LOG
#ifdef ENABLE_MO_DEBUG
#define MO_LOG(msg) do { std::cout << msg << std::endl; } while (0)
#else
#define MO_LOG(msg) do {} while (0)
#endif
#endif

// Anonymous namespace for internal helpers / state.
namespace
{
int numOfLoop = 0; // Counts TPLS iterations executed in tplsPhase.

// Safely destroys an ExploredSolution entry and its owned solution.
void destroyExploredSolution(ExploredSolution* entry)
{
    if (!entry)
    {
        return;
    }
    delete entry->Solution;
    delete entry;
}
} // namespace

/**
 * writeResultsToFile
 * ---------------------------------------------------------------
 * Writes the final set of objective value pairs (X,Y) to a CSV file.
 * Each line includes a placeholder RGB color tuple for potential plotting.
 */
void writeResultsToFile(const std::vector<int>& ObjValX, const std::vector<int>& ObjValY)
{
    // Overwrite file with header first.
    {
        std::ofstream file("../graph.csv", std::ios::trunc);
        file << "X,Y,Color\n";
    }

    // Append all points.
    {
        std::ofstream file("../graph.csv", std::ios::app);
        for (int i = 0; i < static_cast<int>(ObjValX.size()); i++)
        {
            // Fixed color (red). Could be mapped to dominance depth in future.
            const double colR = 1.0;
            const double colG = 0.0;
            const double colB = 0.0;
            file << ObjValX[i] << "," << ObjValY[i] << ",\"(" << colR << "," << colG << "," << colB << ")\"\n";
        }
    }
}

/**
 * nextLambda
 * ---------------------------------------------------------------
 * Generates the next lambda value by recursively bisecting intervals
 * in [0,1]. The first call returns 1.0; subsequent calls yield midpoints
 * ensuring coverage of the space without repetition until intervals shrink.
 */
double nextLambda()
{
    static std::queue<std::pair<double, double>> intervals; // FIFO of intervals awaiting bisection.
    static int state_lambda = 0;

    if (state_lambda == 0)
    {
        ++state_lambda;
        intervals.push(std::make_pair(0.0, 1.0));
        return 1.0; // Start at extreme to emphasize first objective if weighted.
    }

    auto [left, right] = intervals.front();
    intervals.pop();
    double mid = (left + right) / 2.0;
    intervals.push(std::make_pair(left, mid));
    intervals.push(std::make_pair(mid, right));
    return mid;
}

/**
 * sortSolutionsForPLS
 * ---------------------------------------------------------------
 * Sorts solutions ascending by the first objective. This ordering is
 * required for efficient Pareto front extraction (sweeping along X).
 */
void sortSolutionsForPLS(std::vector<emili::Solution*>& solVector, emili::pfsp::PFSP_TPLS* problem)
{
    if (!problem)
    {
        std::cerr << "Error: Problem is not of type PFSP_TPLS." << std::endl;
        std::exit(1);
    }
    std::sort(solVector.begin(), solVector.end(), [&problem](emili::Solution* a, emili::Solution* b)
              {
                  auto* pfspSolutionA = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(a);
                  auto* pfspSolutionB = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(b);

                  if (!pfspSolutionA || !pfspSolutionB)
                  {
                      std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
                      std::exit(1);
                  }

                  int solA = problem->getValueSolutionProblem1(pfspSolutionA->getJobSchedule());
                  int solB = problem->getValueSolutionProblem1(pfspSolutionB->getJobSchedule());
                  return solA < solB;
              });
}

double computeHypervolume2D(const std::vector<std::pair<double,double>>& points,
                            double refX,
                            double refY)
{
    if (points.empty())
    {
        return 0.0;
    }

    // Flatten the point set as expected by moocore (point-major order).
    const size_t n = points.size();
    std::vector<double> data;
    data.reserve(n * 2);
    for (const auto& p : points)
    {
        data.push_back(p.first);
        data.push_back(p.second);
    }

    // Reference point must be strictly dominated by all points (minimisation).
    double ref[2] = {refX, refY};
    const dimension_t d = 2; // Provided by moocore/config.h via hv.h

    return fpli_hv(data.data(), n, d, ref);
}

/**
 * generateParetoFront
 * ---------------------------------------------------------------
 * Assumes input solutions are sorted ascending by objective1 (X).
 * Performs a linear sweep to build a Pareto set w.r.t objective2 (Y).
 * When multiple solutions share identical X, only the last (best Y)
 * is kept (older clones with same X are discarded).
 */
std::vector<emili::Solution*> generateParetoFront(std::vector<emili::Solution*>& solVector, emili::pfsp::PFSP_TPLS* problem)
{
    std::vector<emili::Solution*> solPLSVector;
    for (int i = 0; i < static_cast<int>(solVector.size()); i++)
    {
        auto* currentSolution = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(solVector[i]);
        if (!currentSolution)
        {
            std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
            std::exit(1);
        }
        if (i == 0)
        {
            solPLSVector.push_back(currentSolution->clone());
        }
        for (int j = 0; j < i; j++)
        {
            auto* prevSolution = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(solVector[j]);
            if (!prevSolution)
            {
                std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
                std::exit(1);
            }

            int currValProb2 = problem->getValueSolutionProblem2(currentSolution->getJobSchedule());
            int prevValProb2 = problem->getValueSolutionProblem2(prevSolution->getJobSchedule());

            if (prevValProb2 < currValProb2)
            {
                break; // Current solution dominated in Y by a previous with same or better X ordering.
            }
            if (j == i - 1)
            {
                int currValProb1 = problem->getValueSolutionProblem1(currentSolution->getJobSchedule());
                while (!solPLSVector.empty())
                {
                    auto* lastSolution = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(solPLSVector.back());
                    if (!lastSolution)
                    {
                        std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
                        std::exit(1);
                    }

                    int lastValProb1 = problem->getValueSolutionProblem1(lastSolution->getJobSchedule());
                    if (currValProb1 != lastValProb1)
                    {
                        break; // Stop removing when X differs.
                    }

                    delete solPLSVector.back();
                    solPLSVector.pop_back();
                }
                solPLSVector.push_back(currentSolution->clone());
            }
        }
    }
    return solPLSVector;
}

/**
 * generateParetoFrontExplored
 * ---------------------------------------------------------------
 * Given a vector of ExploredSolution* sorted by objective1, filters to a
 * non-dominated set (w.r.t objective2). Removes dominated entries and
 * updates the count of still-unexplored solutions.
 */
std::vector<ExploredSolution*> generateParetoFrontExplored(std::vector<ExploredSolution*>& Aexpl, emili::pfsp::PFSP_TPLS* problem, int& soltoExplore)
{
    std::vector<ExploredSolution*> paretoSet;
    int numOfExplored = 0;
    for (auto* entry : Aexpl)
    {
        auto* currentSol = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(entry->Solution);
        if (!currentSol)
        {
            std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
            std::exit(1);
        }

        int currX = problem->getValueSolutionProblem1(currentSol->getJobSchedule());
        int currY = problem->getValueSolutionProblem2(currentSol->getJobSchedule());

        // Dominance check against last element (because of sorting by X).
        if (!paretoSet.empty())
        {
            auto* lastSol = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(paretoSet.back()->Solution);
            if (!lastSol)
            {
                std::cerr << "Error: Problem is not of type PermutationFlowShopSolution." << std::endl;
                std::exit(1);
            }
            int bestY = problem->getValueSolutionProblem2(lastSol->getJobSchedule());
            if (currY > bestY)
            {
                // Current is dominated; release memory and skip.
                destroyExploredSolution(entry);
                continue;
            }
        }

        // Remove same-X solutions with worse (higher) Y.
        while (!paretoSet.empty())
        {
            auto* lastSol = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(paretoSet.back()->Solution);
            if (!lastSol)
            {
                std::cerr << "Error: Problem is not of type PermutationFlowShopSolution." << std::endl;
                std::exit(1);
            }
            int lastX = problem->getValueSolutionProblem1(lastSol->getJobSchedule());
            if (lastX != currX)
            {
                break; // Different X; stop cleanup.
            }
            destroyExploredSolution(paretoSet.back());
            paretoSet.pop_back();
        }
        paretoSet.push_back(entry);
    }
    // Count solutions still requiring neighborhood exploration.
    for (auto* keptEntry : paretoSet)
    {
        if (!keptEntry->explored)
        {
            numOfExplored++;
        }
    }
    Aexpl.clear();
    soltoExplore = numOfExplored;
    return paretoSet;
}

/**
 * tplsPhase
 * ---------------------------------------------------------------
 * Executes n iterations of (timed) local search collecting cloned solutions.
 * Adjusts lambda either by nextLambda() (any-time variant) or linearly spaced
 * values in [0,1] to sweep trade-offs between objectives.
 */
int tplsPhase(emili::LocalSearch* ls, clock_t startTime, double timeLimitSeconds, int n, std::vector<emili::Solution*>& solVector)
{
    const bool noTimeLimit = (timeLimitSeconds <= 0.0);
    for (numOfLoop = 0; numOfLoop < n; numOfLoop++)
    {
        // Respect global time budget.
        double elapsed = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
        if (!noTimeLimit && elapsed >= timeLimitSeconds)
        {
            break;
        }

        // Do not start a per-step timer: run a single untimed search iteration.
        emili::Solution* solution = ls->search();
        solution = ls->getBestSoFar(); // Ensure we take incumbent best.

        emili::pfsp::PFSP_TPLS* problem = dynamic_cast<emili::pfsp::PFSP_TPLS*>(&ls->getInitialSolution().getProblem());
        if (!problem)
        {
            std::cerr << "Error: Problem is not of type PFSP_TPLS 2." << std::endl;
            return 1;
        }
        solVector.push_back(solution->clone());

        // Lambda update strategy.
        if (ls->getAnyTime())
        {
            problem->setLambda(nextLambda());
        }
        else
        {
            if (n <= 1)
            {
                problem->setLambda(1.0);
            }
            else
            {
                double lambdaValue = static_cast<double>(numOfLoop) / static_cast<double>(n - 1);
                problem->setLambda(lambdaValue);
            }
        }
        ls->changeTPLSPreviousBestSoFar(solution);
    }
    return 0;
}

/**
 * plsPhase
 * ---------------------------------------------------------------
 * Performs Pareto Local Search within a given time limit. Randomly selects
 * unexplored Pareto solutions and explores their neighborhood, adding new
 * solutions that improve either objective.
 *
 * If timeLimitSeconds <= 0, the phase runs until there are no more
 * unexplored solutions in the current Pareto set (no time cap).
 */
int plsPhase(clock_t startTime, double timeLimitSeconds, std::vector<ExploredSolution*>& Aexpl, emili::pfsp::PFSP_TPLS* problem, emili::LocalSearch* ls)
{
    const bool noTimeLimit = (timeLimitSeconds <= 0.0);
    const double timeLimit = noTimeLimit ? std::numeric_limits<double>::infinity() : timeLimitSeconds;
    double timeelapsed3 = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
    int step = 0;
    int solToExplore = static_cast<int>(Aexpl.size());
    while (true)
    {
        timeelapsed3 = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
        // Stop when no solutions remain; stop on time only if a limit is set.
        if (solToExplore <= 0 || (!noTimeLimit && timeelapsed3 >= timeLimit))
        {
            break; // Stop when no solutions remain or time exhausted.
        }
        int randomIndex;
        bool randomOk = false;
        while (!randomOk)
        {
            randomIndex = emili::generateRandomNumber() % Aexpl.size();
            if (!Aexpl[randomIndex]->explored)
            {
                randomOk = true;
                Aexpl[randomIndex]->explored = true;
                step++;
            }
        }
        emili::Solution* s = Aexpl[randomIndex]->Solution;
        emili::Neighborhood& sameNeigh = ls->getNeighborhood();
        for (auto it = sameNeigh.begin(s); it != sameNeigh.end(); ++it)
        {
            timeelapsed3 = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
            if (timeelapsed3 >= timeLimit)
            {
                MO_LOG("time : " << timeelapsed3 << " step : " << step);
                break;
            }
            emili::Solution* sPrime = *it;

            emili::pfsp::PermutationFlowShopSolution& sa = dynamic_cast<emili::pfsp::PermutationFlowShopSolution&>(*s);
            emili::pfsp::PermutationFlowShopSolution& sb = dynamic_cast<emili::pfsp::PermutationFlowShopSolution&>(*sPrime);
            int valx1 = problem->getValueSolutionProblem1(sa.getJobSchedule());
            int valy1 = problem->getValueSolutionProblem2(sa.getJobSchedule());
            int valx2 = problem->getValueSolutionProblem1(sb.getJobSchedule());
            int valy2 = problem->getValueSolutionProblem2(sb.getJobSchedule());

            // Accept if it improves at least one objective (weak dominance move).
            if ((valx2 < valx1) || (valy2 < valy1))
            {
                Aexpl.push_back(new ExploredSolution(sPrime->clone()));
            }
        }

        // Maintain sorted order by objective1 and update Pareto set.
        std::sort(Aexpl.begin(), Aexpl.end(), [&problem](ExploredSolution* a, ExploredSolution* b)
                  {
                      auto* pfspSolutionA = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(a->Solution);
                      auto* pfspSolutionB = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(b->Solution);
                      if (!pfspSolutionA || !pfspSolutionB)
                      {
                          std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
                          std::exit(1);
                      }
                      int solA = problem->getValueSolutionProblem1(pfspSolutionA->getJobSchedule());
                      int solB = problem->getValueSolutionProblem1(pfspSolutionB->getJobSchedule());
                      return solA < solB;
                  });
        Aexpl = generateParetoFrontExplored(Aexpl, problem, solToExplore);
    }
    return 0;
}

/**
 * executeMultiObjective
 * ---------------------------------------------------------------
 * Orchestrates multi-objective solving:
 *  1. Run TPLS iterations collecting candidate solutions.
 *  2. Extract initial Pareto front (A).
 *  3. If requested, perform PLS expansion within remaining time budget.
 *  4. Write resulting non-dominated objective pairs to CSV.
 */
void executeMultiObjective(emili::LocalSearch* ls, float pls, clock_t startTime, bool is_tpls, bool is_pls, bool is_tppls)
{
    std::vector<int> ObjValX; // First objective values.
    std::vector<int> ObjValY; // Second objective values.

    std::vector<emili::Solution*> solVector;   // Raw TPLS solutions.
    std::vector<emili::Solution*> solPLSVector; // Initial Pareto front.
    double searchTimeLimit = static_cast<double>(ls->getSearchTime());

    int numOfLoop;
    if (is_tpls)
    {
        MO_LOG("Running TPLS...");
    }
    else if (is_pls)
    {
        MO_LOG("Running PLS...");
    }
    else if (is_tppls)
    {
        MO_LOG("Running TPPLS...");
    }

    numOfLoop = (is_tpls || is_tppls) ? (ls->getTPLSStep() + 1) : 2;
    // Use global time budget for the whole program; do not start per-step timers.
    tplsPhase(ls, startTime, searchTimeLimit, numOfLoop, solVector);
    MO_LOG("Number of solutions before PLS: " << solVector.size());

    emili::pfsp::PFSP_TPLS* problem = dynamic_cast<emili::pfsp::PFSP_TPLS*>(&ls->getInitialSolution().getProblem());
    sortSolutionsForPLS(solVector, problem);
    solPLSVector = generateParetoFront(solVector, problem);
    MO_LOG("Number of solutions after Pareto front: " << solPLSVector.size());

    // Build the explored set using the lightweight wrapper.
    std::vector<ExploredSolution*> Aexpl;
    for (int i = 0; i < static_cast<int>(solPLSVector.size()); i++)
    {
        Aexpl.push_back(new ExploredSolution(solPLSVector[i]->clone()));
    }

    double elapsedTime = static_cast<double>(clock() - startTime) / CLOCKS_PER_SEC;
    double remainingTime = (searchTimeLimit > 0.0) ? std::max(0.0, searchTimeLimit - elapsedTime) : 0.0;

    if (is_pls || is_tppls)
    {
        if (searchTimeLimit > 0.0 && remainingTime <= 0.0)
        {
            MO_LOG("Allocated time (-it) exhausted before PLS phase; stopping.");
        }
        else
        {
            double phaseLimit = (searchTimeLimit > 0.0) ? remainingTime : 0.0;
            plsPhase(startTime, phaseLimit, Aexpl, problem, ls);
        }
    }

    // Collect objective pairs.
    for (int i = 0; i < static_cast<int>(Aexpl.size()); i++)
    {
        auto* currentSolution = dynamic_cast<emili::pfsp::PermutationFlowShopSolution*>(Aexpl[i]->Solution);
        if (!currentSolution)
        {
            std::cerr << "Error: Solution is not of type PermutationFlowShopSolution." << std::endl;
            std::exit(1);
        }
        int val1 = problem->getValueSolutionProblem1(currentSolution->getJobSchedule());
        int val2 = problem->getValueSolutionProblem2(currentSolution->getJobSchedule());
        ObjValX.push_back(val1);
        ObjValY.push_back(val2);
    }

    writeResultsToFile(ObjValX, ObjValY);


    std::vector<std::pair<double,double>> moPoints;
    moPoints.reserve(ObjValX.size());

    int maxX = std::numeric_limits<int>::min();
    int maxY = std::numeric_limits<int>::min();

    for (size_t i = 0; i < ObjValX.size(); ++i)
    {
        int x = ObjValX[i];
        int y = ObjValY[i];
        moPoints.emplace_back(static_cast<double>(x),
                            static_cast<double>(y));

        if (x > maxX) maxX = x;
        if (y > maxY) maxY = y;
    }

    // Choisir le point de référence (un peu pire que le pire point)
    double refX = static_cast<double>(maxX) * 1.01;
    double refY = static_cast<double>(maxY) * 1.01;

    // Calculer l'hypervolume
    double hv = computeHypervolume2D(moPoints, refX, refY);

    // L'afficher pour irace (ou logs)
    std::cout << "Hypervolume (2D, minimisation) : " << hv << std::endl;

    // Print total execution time (like single-objective version) unconditionally.
    double time_elapsed = (double)(clock() - startTime) / CLOCKS_PER_SEC;
    std::cout << "time : " << time_elapsed << std::endl;
    // Also print iteration counter for consistency with single-objective output.
    std::cout << "iteration counter : " << emili::iteration_counter() << std::endl;

    // --- Cleanup cloned solutions to avoid memory leaks ---
    for (auto* s : solVector)
    {
        delete s;
    }
    for (auto* s : solPLSVector)
    {
        delete s;
    }
    for (auto* e : Aexpl)
    {
        destroyExploredSolution(e);
    }
}
