//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef PFSPBUILDER_H
#define  PFSPBUILDER_H 
#include "../generalParser.h"
#include "permutationflowshop.h"
namespace prs
{
void info();

class PfspBuilder: public Builder
{
public:
    PfspBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):
        Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition);
    virtual bool canOpenInstance(char* problem_definition);
    virtual emili::Problem* openInstance();
    virtual emili::Problem* openInstance(emili::Problem* instance1, emili::Problem* instance2);
    virtual emili::LocalSearch* buildAlgo();
    virtual emili::InitialSolution* buildInitialSolution();
    virtual emili::Neighborhood* buildNeighborhood();
    virtual emili::Termination* buildTermination();
    virtual emili::Perturbation* buildPerturbation();
    virtual emili::Acceptance* buildAcceptance();
    virtual emili::TabuMemory* buildTabuTenure();
    virtual emili::GASelection* buildGASelection();
    virtual emili::GACrossover* buildGACrossover();
    virtual emili::GAMutation* buildGAMutation();
    virtual emili::Problem* buildProblem();    
};

}

// PFSP GA operator classes (permutation representation)
namespace emili {
namespace pfsp {

class PfspGATournamentSelection : public emili::GASelection {
protected:
    int k;
public:
    PfspGATournamentSelection(int k_ = 2):k(k_) {}
    virtual int select(const std::vector<emili::Solution*>& population);
};

class PfspGAOXCrossover : public emili::GACrossover {
protected:
    emili::pfsp::PermutationFlowShop& instance;
public:
    PfspGAOXCrossover(emili::pfsp::PermutationFlowShop& inst):instance(inst) {}
    virtual emili::Solution* crossover(const emili::Solution& parent1, const emili::Solution& parent2);
};

class PfspGASwapMutation : public emili::GAMutation {
protected:
    emili::pfsp::PermutationFlowShop& instance;
public:
    PfspGASwapMutation(emili::pfsp::PermutationFlowShop& inst):instance(inst) {}
    virtual void mutate(emili::Solution& individual);
};

}
}

#endif //  PFSPBUILDER_H
