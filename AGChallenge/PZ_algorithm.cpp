#include "PZ_algorithm.h"

void PZ_algorithm::run_thread_searching_new_optima(int index_in_vectors) {
	// params:
	int param_population_size = 100;
	float param_cross_prob = 1.0;
	float param_mut_prob = 0.001;
	float param_choose_better_in_cross_prob = 0.0; // previously it was 0.5 both values are ok
	int param_linkage_tree_separation_size = 100;
	int param_linkage_tree_min_cluster = 3;
	int param_linkage_tree_max_cluster = 99;
	bool param_use_generic_tree = true;
	bool param_verbose = false;
	// end of params

	int no_improvement_for_n_iterations = 0;
	double best_fitness = 0;

	CLFLnetEvaluator evaluator_copy = CLFLnetEvaluator();
	evaluator_copy.bConfigure(evaluator->sGetNetName());

	GeneticAlgorithm evolutionary_algorithm(
		&evaluator_copy,
		param_population_size,
		param_cross_prob,
		param_mut_prob,
		param_choose_better_in_cross_prob,
		param_linkage_tree_separation_size,
		param_linkage_tree_min_cluster,
		param_linkage_tree_max_cluster,
		param_use_generic_tree,
		param_verbose
	);

	while (
		evolutionary_algorithm.get_runned_iterations() < NEW_LOCALY_SEARCHED_INDIVIDUAL_MIN_ITERATIONS ||
		evolutionary_algorithm.get_runned_iterations() > NEW_LOCALY_SEARCHED_INDIVIDUAL_MAX_ITERATIONS ||
		no_improvement_for_n_iterations > NEW_LOCALY_SEARCHED_INDIVIDUAL_STOP_AFTER_NO_IMPROVEMENT_ITERATIONS)
	{
		evolutionary_algorithm.run_iteration_scattered();

		if (evolutionary_algorithm.get_best_fitness() > best_fitness) {
			best_fitness = evolutionary_algorithm.get_best_fitness();
			no_improvement_for_n_iterations = 0;
			Individual * best_individual_copy = new Individual(*(evolutionary_algorithm.get_best_individual()));
			best_individual_copy->set_evaluator(nullptr);
			best_individual_copy->set_random_values_holder(nullptr);
			// actualising best individual for main thread
			{
				std::lock_guard<std::mutex> lock(this->mutex_best_individuals_from_threads);
				if (this->best_individuals_from_threads[index_in_vectors] != nullptr) {
					delete this->best_individuals_from_threads[index_in_vectors];
				}
				this->best_individuals_from_threads[index_in_vectors] = best_individual_copy;
			}
		}
		else {
			no_improvement_for_n_iterations++;
		}
	}

	{
		std::lock_guard<std::mutex> lock(this->mutex_threads_finished);
		this->threads_finished[index_in_vectors] = true;
	}
}


PZ_algorithm::PZ_algorithm(CLFLnetEvaluator* evaluator)
{
	this->thread_pool_size = THREADS_NUMBER_ALL - 1;
	this->evaluator = evaluator;
	this->thread_pool = new ThreadPool(thread_pool_size);
	this->current_thread_to_take_from = 0;
	
	// params for best individuals algorithm
	int param_population_size = 100;
	float param_cross_prob = 1.0;
	float param_mut_prob = 0.001;
	float param_choose_better_in_cross_prob = 0.0; // previously it was 0.5 both values are ok
	int param_linkage_tree_separation_size = this->evaluator->iGetNumberOfBits();
	int param_linkage_tree_min_cluster = 2;
	int param_linkage_tree_max_cluster = param_linkage_tree_separation_size - 1;
	bool param_use_generic_tree = false; // if set to false it will create real tree
	bool param_verbose = true;

	this->algorithm_best_individuals = new GeneticAlgorithm(
		this->evaluator,
		param_population_size,
		param_cross_prob,
		param_mut_prob,
		param_choose_better_in_cross_prob,
		param_linkage_tree_separation_size,
		param_linkage_tree_min_cluster,
		param_linkage_tree_max_cluster,
		param_use_generic_tree,
		param_verbose
	);
	// end of params


	this->best_individual = new Individual(*this->algorithm_best_individuals->get_best_individual());
	best_individuals_from_threads = vector<Individual*>(thread_pool_size, nullptr);
	threads_finished = vector<bool>(thread_pool_size, false);

	for(int i = 0; i < thread_pool_size; i++) {
		threads_futures.emplace_back(thread_pool->enqueue(&PZ_algorithm::run_thread_searching_new_optima, this, i));
	}
}



void PZ_algorithm::run_iteration() {
	vector<Individual * > individuals_to_add_to_algorithm (0);

	// ckecking for new best individuals and collecting individuals to add to algorithm
	{
		std::lock_guard<std::mutex> lock(this->mutex_best_individuals_from_threads);
		for (Individual * individual : this->best_individuals_from_threads)
		{
			if (individual != nullptr) {
				check_and_actualise_best_individual(*individual);
			}
		}
		if (this->best_individuals_from_threads[current_thread_to_take_from] != nullptr) {
			individuals_to_add_to_algorithm.push_back(new Individual(*this->best_individuals_from_threads[current_thread_to_take_from]));
		}
	}
	current_thread_to_take_from++;
	current_thread_to_take_from %= thread_pool_size;


	// starting new searching if there are any that has finished
	{
		std::lock_guard<std::mutex> lock(this->mutex_threads_finished);
		
		for (int i = 0; i < thread_pool_size; i++) {
			if (threads_finished[i]) {
				threads_finished[i] = false;
				threads_futures[i].get();
				threads_futures[i] = thread_pool->enqueue(&PZ_algorithm::run_thread_searching_new_optima, this, i);
			}
		}
	}
	
	algorithm_best_individuals->run_iteration_LTGA(individuals_to_add_to_algorithm);
	check_and_actualise_best_individual(*(algorithm_best_individuals->get_best_individual()));
}



void PZ_algorithm::check_and_actualise_best_individual(Individual & individual) {
	if (this->best_individual == nullptr) {
		this->best_individual = new Individual(individual);
		this->best_individual->set_evaluator(this->evaluator);
		this->best_individual->set_random_values_holder(this->algorithm_best_individuals->get_random_values_holder());
	}
	else if (individual.get_fitness() > this->best_individual->get_fitness()) {
		delete this->best_individual;
		this->best_individual = new Individual(individual);
		this->best_individual->set_evaluator(this->evaluator);
		this->best_individual->set_random_values_holder(this->algorithm_best_individuals->get_random_values_holder());
	}
}

PZ_algorithm::~PZ_algorithm()
{
	delete thread_pool;
	delete algorithm_best_individuals;
	
	{
		std::lock_guard<std::mutex> lock(this->mutex_best_individuals_from_threads);
		for (Individual * individual : best_individuals_from_threads)
		{
			delete individual;
		}
	}

	delete best_individual;
}

vector<int> PZ_algorithm::get_best_solution() {
	return best_individual->get_genotype();
}

double PZ_algorithm::get_best_fitness() {
	return best_individual->get_fitness();
}

Individual * PZ_algorithm::get_best_individual() {
	return best_individual;
}