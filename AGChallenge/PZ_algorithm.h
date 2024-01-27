#pragma once
#include "Individual.h"
#include "Evaluator.h"
#include "RandomValuesHolder.h"
#include "LinkageTree.h"
#include "AbstractEvolutionaryAlgorithm.h"
#include "GeneticAlgorithm.h"
#include "ThreadPool.h"

#define THREADS_NUMBER_ALL 6 // all number of threads used in the algorithm, including the one used for the main algorithm
#define NEW_LOCALY_SEARCHED_INDIVIDUAL_MIN_ITERATIONS 100
#define NEW_LOCALY_SEARCHED_INDIVIDUAL_STOP_AFTER_NO_IMPROVEMENT_ITERATIONS 200
#define NEW_LOCALY_SEARCHED_INDIVIDUAL_MAX_ITERATIONS 10000

#define TAKE_FROM_DIFFERENT_ISLAND_PROB 0.02

#define LOGGING_EVERY_N_ITERATIONS 5
#define VERBOSE false

using namespace std;

class PZ_algorithm : public AbstractEvolutionaryAlgorithm
{
private:
	GeneticAlgorithm* algorithm_best_individuals;
	Individual* best_individual;

	ThreadPool* thread_pool;
	vector<Individual * > best_individuals_from_threads;
	mutex mutex_best_individuals_from_threads;
	vector<bool> threads_finished;
	mutex mutex_threads_finished;
	vector<future<void> > threads_futures;

	int thread_pool_size;
	int current_thread_to_take_from;
	CLFLnetEvaluator* evaluator; // this one is passed, so it shouldnt be deleted

	void check_and_actualise_best_individual(Individual& individual);
public:
	PZ_algorithm(CLFLnetEvaluator* evaluator);
	~PZ_algorithm();
	void run_iteration();
	vector<int> get_best_solution();
	double get_best_fitness();
	Individual* get_best_individual();
	void run_thread_searching_new_optima(int index_in_vectors);
};

