#include "PZ_Meta_Algorithm.h"



PZ_Meta_Algorithm::PZ_Meta_Algorithm(CLFLnetEvaluator* evaluator) {
	this->evaluator = evaluator;
	this->pz_algorithm = new PZ_algorithm(this->evaluator);
	this->best_individual = nullptr;


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

	meta_GA = new GeneticAlgorithm(
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


	current_best_fitness_tmp = 0;
	iteration_no_improvement = 0;
	self_iteration_counter = 0;
	check_and_actualise_best_individual(*(this->pz_algorithm->get_best_individual()));
	best_individuals.push_back(new Individual(*(pz_algorithm->get_best_individual())));
}


void PZ_Meta_Algorithm::run_iteration() {
	self_iteration_counter++;
	if (get_best_fitness() < 1.0) {

		this->pz_algorithm->run_iteration();
		if (this->pz_algorithm->get_best_fitness() > current_best_fitness_tmp) {
			iteration_no_improvement = 0;
			current_best_fitness_tmp = this->pz_algorithm->get_best_fitness();
			check_and_actualise_best_individual(*(pz_algorithm->get_best_individual()));
			delete best_individuals[best_individuals.size() - 1];
			best_individuals[best_individuals.size() - 1] = new Individual(*(pz_algorithm->get_best_individual()));
		}
		else {
			iteration_no_improvement++;

			if (iteration_no_improvement > RESET_AFTER_N_ITERATIONS_NO_IMPROVEMENT) {
				current_best_fitness_tmp = 0;
				iteration_no_improvement = 0;
				delete this->pz_algorithm;
				this->pz_algorithm = new PZ_algorithm(evaluator);
				best_individuals.push_back(new Individual(*(pz_algorithm->get_best_individual())));
			}
		}
	}

	if (self_iteration_counter % RUN_META_ITERATION_EVERY_N_ITERATIONS == 0) {
		if ((self_iteration_counter / RUN_META_ITERATION_EVERY_N_ITERATIONS)
			% INSERT_INDIVIDUALS_META_EVERY_N_ITERATIONS == 0)
		{
			vector<Individual* > best_individuals_copy = vector<Individual* >(0);
			for (Individual* individual : best_individuals) {
				best_individuals_copy.push_back(new Individual(*individual));
			}
			meta_GA->insert_individuals(best_individuals_copy);
		}

		meta_GA->run_iteration();
		check_and_actualise_best_individual(*(meta_GA->get_best_individual()));
		//meta_GA->print_iteration_summary();
	}
}


PZ_Meta_Algorithm::~PZ_Meta_Algorithm() {
	if (meta_GA != nullptr) {
		delete meta_GA;
	}
	if (pz_algorithm != nullptr) {
		delete pz_algorithm;
	}

	delete best_individual;
	for (Individual* individual : best_individuals) {
		delete individual;
	}
}

void PZ_Meta_Algorithm::check_and_actualise_best_individual(Individual& individual) {
	if (this->best_individual == nullptr) {
		this->best_individual = new Individual(individual);
		this->best_individual->set_evaluator(this->evaluator);
		this->best_individual->set_random_values_holder(nullptr);
	}
	else if (individual.get_fitness() > this->best_individual->get_fitness()) {
		delete this->best_individual;
		this->best_individual = new Individual(individual);
		this->best_individual->set_evaluator(this->evaluator);
		this->best_individual->set_random_values_holder(nullptr);
	}
}

vector<int> PZ_Meta_Algorithm::get_best_solution() {
	return best_individual->get_genotype();
}

double PZ_Meta_Algorithm::get_best_fitness() {
	return best_individual->get_fitness();
}