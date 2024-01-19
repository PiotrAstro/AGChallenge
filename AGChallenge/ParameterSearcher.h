#pragma once
#include "GeneticAlgorithm.h"
#include "ThreadPool.h"

namespace ParameterSearcher{


	void parameters_search(
		int population_size,
		CLFLnetEvaluator * evaluator,
		vector<float> cross_probs,
		vector<float> mut_probs,
		vector<float> choose_better_in_cross_probs);
	void run_n_iterations_genetic_algorithm(GeneticAlgorithm *genetic_algorithm, int iterations);

	void save_stats_every_n_seconds(vector<AbstractEvolutionaryAlgorithm* > evolutionary_algorithms, int seconds, int max_seconds = -1); // if max_seconds = -1 it means there are no max seconds
}