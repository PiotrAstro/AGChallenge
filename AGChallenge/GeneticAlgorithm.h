#pragma once
#include "Individual.h"
#include "Evaluator.h"
#include "PZ_Math.h"
#include "RandomValuesHolder.h"
#include "AbstractEvolutionaryAlgorithm.h"
#include <chrono>

using namespace std;
class GeneticAlgorithm : public AbstractEvolutionaryAlgorithm
{
private:
	int param_population_size;
	float param_cross_prob;
	float param_mut_prob;
	float param_choose_better_in_cross_prob;



	int runned_iterations;
	Individual * best_individual;
	vector<Individual * > population;
	RandomValuesHolder * random_values_holder;
	vector<int> values_range;
	LinkageTree * linkage_tree;
	CLFLnetEvaluator * evaluator; // this one is passed, so it shouldnt be deleted
	chrono::high_resolution_clock::time_point start_time;
	void actualise_best_individual();
	vector<Individual * > generate_crossed_population();
	vector<Individual * > generate_crossed_every_cluster_population();
	vector<Individual* > generate_crossed_random_cluster_population_size();
	void perform_LTGA_routine();
	Individual * get_random_individual_after_fight_of_2();
	int get_random_individual_index_after_fight_of_2();
	vector<vector<int> * > get_population_genotypes();
	void print_iteration_summary();

public:
	GeneticAlgorithm(CLFLnetEvaluator *evaluator, int population_size, float cross_prob, float mut_prob, float choose_better_in_cross_prob);
	~GeneticAlgorithm();
	void run_iteration();
	vector<int> get_best_solution();
	double get_best_fitness();
	double mean_population_entropy();
	vector<double> fitness_quantiles(vector<float> quantiles);
	double GeneticAlgorithm::get_run_time_milis();
	Individual * get_best_individual();

	// getters of parameters
	int get_population_size();
	float get_cross_prob();
	float get_mut_prob();
	float get_choose_better_in_cross_prob();
};


