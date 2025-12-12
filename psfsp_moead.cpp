#include "pfsp_moead.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>

#ifndef restrict
#define restrict __restrict__
#endif

#include "moocore/c/hv.h"
#include "moocore/c/config.h"  // for dimension_t

#ifdef _OPENMP
#include <omp.h>
#endif

namespace emili {
namespace pfsp {

/**
 * computeHypervolume2D
 * ---------------------------------------------------------------
 * Computes the 2D hypervolume using moocore's fpli_hv algorithm.
 * Points are expected to be in minimization form.
 * refX, refY: reference point (should be dominated by all points).
 */
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

MOEAD::MOEAD(PFSP_TPLS& prob,
             emili::InitialSolution& initSol,
             GASelection& sel,
             GACrossover& cross,
             GAMutation& mut,
             emili::LocalSearch* improverLS,
             int populationSize,
             int neighborhoodSize,
             float cr,
             float mr,
             float timeLimit)
    : problem(prob),
      init(initSol),
      selection(sel),
      crossover(cross),
      mutation(mut),
      improver(improverLS),
      N(populationSize),
      T(neighborhoodSize),
      crossoverRate(cr),
      mutationRate(mr),
      timeLimitSeconds(timeLimit)
{
    if (N <= 0) {
        std::cerr << "MOEAD: populationSize must be > 0, setting to 10\n";
        N = 10;
    }
    if (T <= 0 || T > N) {
        std::cerr << "MOEAD: neighborhoodSize invalid; setting T = min(10, N)\n";
        T = std::min(10, N);
    }

    lambdas.resize(N);
    neighbors.resize(N);
    population.assign(N, nullptr);
    objectives.resize(N);

    zstar[0] = std::numeric_limits<int>::max();
    zstar[1] = std::numeric_limits<int>::max();
}

MOEAD::~MOEAD()
{
    clearPopulation();
}

void MOEAD::clearPopulation()
{
    for (std::size_t i = 0; i < population.size(); ++i) {
        delete population[i];
        population[i] = nullptr;
    }
}

void MOEAD::initializeWeights()
{
    if (N == 1) {
        lambdas[0] = {1.0, 0.0};
        return;
    }

    for (int i = 0; i < N; ++i) {
        double w1 = static_cast<double>(i) / static_cast<double>(N - 1);
        double w2 = 1.0 - w1;
        lambdas[i] = {w1, w2};
    }
}

void MOEAD::initializeNeighborhoods()
{
    std::vector<std::pair<double,int>> dist(N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double dx = lambdas[i][0] - lambdas[j][0];
            double dy = lambdas[i][1] - lambdas[j][1];
            double d  = std::sqrt(dx*dx + dy*dy);
            dist[j]   = std::make_pair(d, j);
        }
        std::sort(dist.begin(), dist.end(),
                  [](const std::pair<double,int>& a,
                     const std::pair<double,int>& b) {
                      return a.first < b.first;
                  });
        neighbors[i].clear();
        for (int k = 0; k < T; ++k) {
            neighbors[i].push_back(dist[k].second);
        }
    }
}

void MOEAD::evaluateSolution(emili::Solution* s, std::array<int,2>& f) const
{
    auto* pfspSol = dynamic_cast<PermutationFlowShopSolution*>(s);
    if (!pfspSol) {
        std::cerr << "MOEAD::evaluateSolution: wrong solution type\n";
        std::exit(1);
    }
    const std::vector<int>& seq = pfspSol->getJobSchedule();
    // getValueSolutionProblem1/2 expect a non-const reference; create a non-const copy
    std::vector<int> seq_nc(seq);
    f[0] = problem.getValueSolutionProblem1(seq_nc);
    f[1] = problem.getValueSolutionProblem2(seq_nc);
}

void MOEAD::updateIdealPoint(const std::array<int,2>& f)
{
    if (f[0] < zstar[0]) zstar[0] = f[0];
    if (f[1] < zstar[1]) zstar[1] = f[1];
}

double MOEAD::g_tchebycheff(int i, const std::array<int,2>& f) const
{
    double v1 = lambdas[i][0] * std::abs(f[0] - zstar[0]);
    double v2 = lambdas[i][1] * std::abs(f[1] - zstar[1]);
    return (v1 > v2) ? v1 : v2;
}

void MOEAD::initializePopulation()
{
    clearPopulation();
    zstar[0] = std::numeric_limits<int>::max();
    zstar[1] = std::numeric_limits<int>::max();

    // Initialisation simple (NON parallélisée, pour rester safe avec le RNG).
    for (int i = 0; i < N; ++i) {
        emili::Solution* s = init.generateSolution();
        population[i]      = s;
        std::array<int,2> f;
        evaluateSolution(s, f);
        objectives[i] = f;
        updateIdealPoint(f);
    }
}

emili::Solution* MOEAD::makeChild(emili::Solution* p1, emili::Solution* p2)
{
    emili::Solution* child = nullptr;

    float r = emili::generateRealRandomNumber();
    if (r < crossoverRate) {
        child = crossover.crossover(*p1, *p2);
    } else {
        child = p1->clone();
    }

    r = emili::generateRealRandomNumber();
    if (r < mutationRate) {
        mutation.mutate(*child);
    }

    // Partie mémétique : même logique que ton GA.
    if (improver != nullptr) {
        emili::Solution* improved = improver->search(child);
        if (improved && *improved < *child) {
            delete child;
            child = improved;
        } else {
            delete improved;
        }
    }

    return child;
}

void MOEAD::run()
{
    initializeWeights();
    initializeNeighborhoods();
    initializePopulation();

    clock_t start = std::clock();
    const double limit = (timeLimitSeconds > 0.0)
                             ? timeLimitSeconds
                             : std::numeric_limits<double>::infinity();

    while (true) {
        double elapsed = (double)(std::clock() - start) / (double)CLOCKS_PER_SEC;
        if (elapsed >= limit) {
            break;
        }

        // Boucle principale multithreadée : 1 enfant par sous-problème i.
        // NB: on protège la partie "update" avec une section critique.
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            // Sélection de parents dans le voisinage de i
            const std::vector<int>& B = neighbors[i];
            int r1_index = B[ emili::generateRandomNumber() % B.size() ];
            int r2_index = B[ emili::generateRandomNumber() % B.size() ];
            emili::Solution* p1 = population[r1_index];
            emili::Solution* p2 = population[r2_index];

            // Génération de l’enfant
            emili::Solution* child = makeChild(p1, p2);

            // Évaluation de l’enfant
            std::array<int,2> fchild;
            evaluateSolution(child, fchild);

            // Mise à jour globale (point idéal + voisins) protégée
            #pragma omp critical(moead_update)
            {
                // update du point idéal
                updateIdealPoint(fchild);

                // mise à jour des voisins
                for (int jIdx = 0; jIdx < (int)B.size(); ++jIdx) {
                    int j = B[jIdx];
                    double gj_old = g_tchebycheff(j, objectives[j]);
                    double gj_new = g_tchebycheff(j, fchild);

                    if (gj_new <= gj_old) {
                        delete population[j];
                        population[j] = child->clone();
                        objectives[j] = fchild;
                    }
                }
            }

            delete child;
        }
    }
}

std::vector<emili::Solution*> MOEAD::getNonDominatedSet() const
{
    std::vector<emili::Solution*> pareto;

    for (int i = 0; i < N; ++i) {
        const std::array<int,2>& fi = objectives[i];
        bool dominated = false;

        // test si fi est dominé par une solution déjà dans pareto
        for (std::size_t k = 0; k < pareto.size() && !dominated; ++k) {
            auto* pfspSol = dynamic_cast<PermutationFlowShopSolution*>(pareto[k]);
            if (!pfspSol) {
                std::cerr << "MOEAD::getNonDominatedSet: wrong solution type\n";
                std::exit(1);
            }
            std::array<int,2> fp = {
                problem.getValueSolutionProblem1(pfspSol->getJobSchedule()),
                problem.getValueSolutionProblem2(pfspSol->getJobSchedule())
            };

            bool p_better_eq  = (fp[0] <= fi[0]) && (fp[1] <= fi[1]);
            bool p_strict     = (fp[0] <  fi[0]) || (fp[1] <  fi[1]);
            if (p_better_eq && p_strict) {
                dominated = true;
            }
        }

        if (!dominated) {
            // enlever les solutions de pareto dominées par fi
            std::vector<emili::Solution*> newPareto;
            newPareto.reserve(pareto.size() + 1);
            for (std::size_t k = 0; k < pareto.size(); ++k) {
                auto* pfspSol = dynamic_cast<PermutationFlowShopSolution*>(pareto[k]);
                std::array<int,2> fp = {
                    problem.getValueSolutionProblem1(pfspSol->getJobSchedule()),
                    problem.getValueSolutionProblem2(pfspSol->getJobSchedule())
                };

                bool i_better_eq = (fi[0] <= fp[0]) && (fi[1] <= fp[1]);
                bool i_strict    = (fi[0] <  fp[0]) || (fi[1] <  fp[1]);

                if (!(i_better_eq && i_strict)) {
                    newPareto.push_back(pareto[k]);
                }
                // sinon pareto[k] est dominée par fi, on la jette (mais pas delete
                // ici, car elle appartient à population)
            }
            pareto.swap(newPareto);

            // on ajoute un clone de population[i]
            pareto.push_back(population[i]->clone());
        }
    }

    return pareto;
}

void MOEAD::printHypervolume() const
{
    // Récupération du Pareto non dominé (des clones)
    std::vector<emili::Solution*> pareto = getNonDominatedSet();

    // Collecte (f1,f2)
    std::vector<std::pair<double,double>> moPoints;
    moPoints.reserve(pareto.size());

    double maxX = -1e300;
    double maxY = -1e300;

    for (emili::Solution* s : pareto)
    {
        auto* pfspSol = dynamic_cast<PermutationFlowShopSolution*>(s);
        if (!pfspSol)
        {
            delete s;
            continue;
        }

        const std::vector<int>& seq = pfspSol->getJobSchedule();

        // Copie si getValueSolutionProblemX modifie la séquence
        std::vector<int> seq1(seq);
        std::vector<int> seq2(seq);

        double f1 = problem.getValueSolutionProblem1(seq1);
        double f2 = problem.getValueSolutionProblem2(seq2);

        moPoints.emplace_back(f1, f2);

        if (f1 > maxX) maxX = f1;
        if (f2 > maxY) maxY = f2;

        delete s; // libération du clone
    }

    // Si aucun point → hypervolume nul
    if (moPoints.empty())
    {
        std::cout << "Hypervolume (2D,min) : 0\n";
        return;
    }

    // Référence pour hypervolume (un peu plus grand que le max)
    double refX = maxX * 1.01;
    double refY = maxY * 1.01;

    double hv = computeHypervolume2D(moPoints, refX, refY);

    // Impression simple, propre, sans CSV :
    std::cout << "HV " << hv << std::endl;
}

} // namespace pfsp
}
