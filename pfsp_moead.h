#ifndef PFSP_MOEAD_H
#define PFSP_MOEAD_H

#include "emilibase.h"
#include "multiobjective.h"
#include <vector>
#include <array>
#include <limits>
#include <ctime>
#include <fstream>
#include <string>

namespace emili {
namespace pfsp {

class MOEAD
{
public:
    MOEAD(PFSP_TPLS& problem,
          emili::InitialSolution& init,
          GASelection& selection,
          GACrossover& crossover,
          GAMutation& mutation,
          emili::LocalSearch* improver,
          int populationSize,
          int neighborhoodSize,
          float crossoverRate,
          float mutationRate,
          float timeLimitSeconds);

    ~MOEAD();

    /// Lance l’algorithme jusqu’à atteindre la limite de temps.
    void run();

    void printHypervolume() const;

    /// Accès direct à la population interne (à ne pas delete).
    const std::vector<emili::Solution*>& getPopulation() const { return population; }

    /// Renvoie un ensemble de solutions non dominées (clonées).
    /// L’appelant est responsable de delete les pointeurs.
    std::vector<emili::Solution*> getNonDominatedSet() const;

private:
    PFSP_TPLS&              problem;
    emili::InitialSolution& init;
    GASelection&            selection;
    GACrossover&            crossover;
    GAMutation&             mutation;
    emili::LocalSearch*     improver; // peut être nullptr

    int    N;                 // nb de sous-problèmes = taille de la population
    int    T;                 // taille du voisinage dans l’espace des poids
    float  crossoverRate;
    float  mutationRate;
    float  timeLimitSeconds;  // en secondes

    // Données MOEA/D
    std::vector<std::array<double,2>> lambdas;    // vecteurs de poids (w1,w2)
    std::vector<std::vector<int>>     neighbors;  // indices des voisins pour chaque i
    std::vector<emili::Solution*>     population; // x_i
    std::vector<std::array<int,2>>    objectives; // (f1,f2) pour chaque x_i
    std::array<int,2>                 zstar;      // point idéal (min f1, min f2)

    // Méthodes internes
    void clearPopulation();
    void initializeWeights();
    void initializeNeighborhoods();
    void initializePopulation();

    void evaluateSolution(emili::Solution* s, std::array<int,2>& f) const;
    void updateIdealPoint(const std::array<int,2>& f);
    double g_tchebycheff(int subproblemIndex, const std::array<int,2>& f) const;

    /// Génère un enfant à partir de p1,p2 avec crossover/mutation + improver (memetic).
    emili::Solution* makeChild(emili::Solution* p1, emili::Solution* p2);
};

} // namespace pfsp
} // namespace emili

#endif // PFSP_MOEAD_H
