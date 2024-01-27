#pragma once
#include "Individual.h"
#include "Evaluator.h"
#include "PZ_Math.h"
#include "RandomValuesHolder.h"
#include "AbstractEvolutionaryAlgorithm.h"
#include <chrono>

#define RESET_POPULATION_AFTER_NO_IMPROVEMENT_FOR 100
#define STATS_AND_LINKAGE_TREE_EVERY_N_ITERATIONS 20
#define FHIC_TRY_N_GENE_VALUES 15
#define INTRODUCE_PERCENT_OF_NEW_INDIVIDUALS 0.01
#define SCATTERED_ROUTINE_CROSSES_PER_ITERATION 1000
#define LTGA_ROUTINE_CROSSES_PER_ITERATION 500
#define LTGA_FHIC_EVERY_N_ITERATIONS 40
#define LTGA_CHILD_FHIC_PROB 0.00005
#define LTGA_CHILD_PROB_FHIC_WITH_CROSSED 0.001

#define TREE_CROSS_GREEDY_CHOOSE_RANDOMLY_PROB 0.5
#define TREE_CROSS_RANDOM_CHANGE_GENES_PROB 0.6
#define TREE_CROSS_STOP_DELVING_DEEPER_PROB 0.01

using namespace std;
class GeneticAlgorithm : public AbstractEvolutionaryAlgorithm
{
private:
	int tmp_test_loaded_individual_index = 0;


	int param_population_size;
	float param_cross_prob;
	float param_mut_prob;
	float param_choose_better_in_cross_prob; // in this version it is not used
	int param_linkage_tree_separation_size;
	int param_linkage_tree_min_cluster;
	int param_linkage_tree_max_cluster;
	bool param_use_generic_tree;
	bool param_verbose;


	int last_improvement_iteration;
	int best_individuals_counter;
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
	void perform_scatter_routine();
	Individual * get_random_individual_after_fight_of_2();
	int get_random_individual_index_after_fight_of_2();
	int get_random_individual_index_not_the_only_best();
	vector<vector<int> * > get_population_genotypes();

public:
	GeneticAlgorithm(
		CLFLnetEvaluator* evaluator,
		int population_size,
		float cross_prob,
		float mut_prob,
		float choose_better_in_cross_prob,
		int linkage_tree_separation_size,
		int linkage_tree_min_cluster,
		int linkage_tree_max_cluster,
		bool use_generic_tree,
		bool verbose
	);
	~GeneticAlgorithm();
	void print_iteration_summary();
	void run_iteration();
	void run_iteration_scattered();
	void insert_individuals(vector<Individual* > individuals_to_add);
	void run_iteration_LTGA();
	void introduce_n_new_individuals(int introduce_new_individuals);
	vector<int> get_best_solution();
	double get_best_fitness();
	double mean_population_entropy();
	double static mean_population_entropy(const vector<Individual* > & individuals);
	vector<double> fitness_quantiles(vector<float> quantiles);
	double GeneticAlgorithm::get_run_time_milis();
	Individual * get_best_individual();
	int get_runned_iterations();
	RandomValuesHolder * get_random_values_holder();

	// getters of parameters
	int get_population_size();
	float get_cross_prob();
	float get_mut_prob();
	float get_choose_better_in_cross_prob();
};


