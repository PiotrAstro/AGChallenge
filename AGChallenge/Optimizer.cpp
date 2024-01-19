#include "Optimizer.h"

#include <cfloat>
#include <iostream>
#include <windows.h>
#include "PZ_Math.h"

using namespace std;

COptimizer::COptimizer(CLFLnetEvaluator &cEvaluator)
	: c_evaluator(cEvaluator)
{
}//COptimizer::COptimizer(CEvaluator &cEvaluator)

void COptimizer::vInitialize()
{
	PZ_Math::PZ_Math_Initialize();
	/*vector<float> cross_probs = {0.2, 0.4, 0.6, 0.7, 0.8};
	vector<float> choose_better_in_cross_probs = { 0.3, 0.5, 0.7, 0.85, 0.95};
	vector<float> mut_probs = { 0.0001, 0.0003, 0.001, 0.003, 0.01};

	ParameterSearcher::parameters_search(
		100,
		&c_evaluator,
		cross_probs,
		mut_probs,
		choose_better_in_cross_probs
	);*/

	


	// normal run
	int param_population_size = 10;
	float param_cross_prob = 0.8;
	float param_mut_prob = 0.001;
	float param_choose_better_in_cross_prob = 0.5;
	int param_linkage_tree_separation_size = 100;
	int param_linkage_tree_min_cluster = 3;
	int param_linkage_tree_max_cluster = 99;
	bool param_use_generic_tree = true;


	this->evolutionary_algorithm = new GeneticAlgorithm(
		&c_evaluator,
		param_population_size,
		param_cross_prob,
		param_mut_prob,
		param_choose_better_in_cross_prob,
		param_linkage_tree_separation_size,
		param_linkage_tree_min_cluster,
		param_linkage_tree_max_cluster,
		param_use_generic_tree
	);
	
	//vector<CLFLnetEvaluator> evaluators(4);
	//vector<GeneticAlgorithm * > algorithms(4);
	//for (int i = 0 ; i < evaluators.size(); ++i) {
	//	evaluators[i].bConfigure(c_evaluator.sGetNetName());
	//	algorithms[i] = new GeneticAlgorithm(&evaluators[i], population_size, cross_prob, mut_prob, choose_better_in_cross_prob);
	//}
	//ThreadPool pool(4);
	//
	//std::vector<std::future<void>> results;
	//for (GeneticAlgorithm* algorithm : algorithms) {
	//	// Enqueue tasks and store futures.
	//	results.emplace_back(pool.enqueue(ParameterSearcher::run_n_iterations_genetic_algorithm, algorithm, 100));
	//}

	//// Wait for all tasks to complete and gather results.
	//for (int i = 0; i < algorithms.size(); ++i) {
	//	results[i].get();
	//	PZ_Math::save_vector(algorithms[i]->get_best_solution(), "C:\\Piotr\\2023_studia\\semestr3\\TEP\\AG\\logs\\saved_individuals\\best_solution_" + to_string(i + 4) + ".txt");
	//}
	
	//P3 algorithm:
	
	//int param_FHIC_genes_checked = -1;
	//int param_FHIC_max_number_of_values_per_gene = 3;
	//
	//this->evolutionary_algorithm = new P3(&c_evaluator, param_FHIC_genes_checked, param_FHIC_max_number_of_values_per_gene);
	
	
	d_current_best_fitness = 0;
}//void COptimizer::vInitialize()

void COptimizer::vRunIteration()
{
	this->evolutionary_algorithm->run_iteration();
	if (this->evolutionary_algorithm->get_best_fitness() > d_current_best_fitness)
	{
		d_current_best_fitness = this->evolutionary_algorithm->get_best_fitness();
		cout << d_current_best_fitness << endl;
		v_current_best = this->evolutionary_algorithm->get_best_solution();
	}
}//void COptimizer::vRunIteration()

COptimizer::~COptimizer()
{
	PZ_Math::PZ_Math_Destroy();
	delete evolutionary_algorithm;
}//COptimizer::~COptimizer()
