#pragma once

#include "Evaluator.h"

#include <random>
#include <vector>
#include "GeneticAlgorithm.h"
#include "ParameterSearcher.h"
#include "P3.h"

using namespace std;

class COptimizer
{
public:
	COptimizer(CLFLnetEvaluator &cEvaluator);
	~COptimizer();

	void vInitialize();
	void vRunIteration();

	vector<int>* pvGetCurrentBest() { return &v_current_best; }
private:
	CLFLnetEvaluator &c_evaluator;
	AbstractEvolutionaryAlgorithm* evolutionary_algorithm;
	double d_current_best_fitness;
	vector<int> v_current_best;
};//class COptimizer