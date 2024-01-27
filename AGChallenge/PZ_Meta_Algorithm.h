#pragma once
#include "Evaluator.h"

#include <random>
#include <vector>
#include "GeneticAlgorithm.h"
#include "ParameterSearcher.h"
#include "PZ_algorithm.h"
#include "P3.h"


#define RESET_AFTER_N_ITERATIONS_NO_IMPROVEMENT 300

#define RUN_META_ITERATION_EVERY_N_ITERATIONS 3
#define INSERT_INDIVIDUALS_META_EVERY_N_ITERATIONS 5

class PZ_Meta_Algorithm : public AbstractEvolutionaryAlgorithm
{
private:
	PZ_algorithm * pz_algorithm;
	double current_best_fitness_tmp;
	int iteration_no_improvement;

	int self_iteration_counter;

	GeneticAlgorithm* meta_GA;
	CLFLnetEvaluator* evaluator; // this one is passed, so it shouldnt be deleted
	vector<Individual* > best_individuals;
	Individual* best_individual;
	void check_and_actualise_best_individual(Individual & individual);

public:
	PZ_Meta_Algorithm(CLFLnetEvaluator* evaluator);
	~PZ_Meta_Algorithm();
	void run_iteration();
	vector<int> get_best_solution();
	double get_best_fitness();
};

