#pragma once
#include <vector>
using namespace std;
class AbstractEvolutionaryAlgorithm
{

public:
	virtual ~AbstractEvolutionaryAlgorithm();
	virtual void run_iteration() = 0;
	virtual vector<int> get_best_solution() = 0;
	virtual double get_best_fitness() = 0;
};

