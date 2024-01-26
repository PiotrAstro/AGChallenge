#include "P3.h"



Individual * P3::create_localy_optimized_individual(CLFLnetEvaluator * evaluator) {
	this->creating_new_individuals_finished.fetch_add(1);

	// params:
	int param_population_size = 100;
	float param_cross_prob = 1.0;
	float param_mut_prob = 0.001;
	float param_choose_better_in_cross_prob = 0.5;
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
		evolutionary_algorithm.run_iteration();
		if (evolutionary_algorithm.get_best_fitness() > best_fitness) {
			best_fitness = evolutionary_algorithm.get_best_fitness();
			no_improvement_for_n_iterations = 0;
		}
		else {
			no_improvement_for_n_iterations++;
		}
	}

	this->creating_new_individuals_finished.fetch_sub(1);

	while (this->creating_new_individuals_finished.load() > 0) {
		evolutionary_algorithm.run_iteration();
	}

	return new Individual(*evolutionary_algorithm.get_best_individual());
}





P3::P3(CLFLnetEvaluator* evaluator) {
	this->evaluator = evaluator;
	this->threads_number = THREADS_NUMBER;
	this->thread_pool = new ThreadPool(this->threads_number);
	levels = vector<P3Level* >();

	int genotype_size = evaluator->iGetNumberOfBits();
	genes_ranges = vector<int>(genotype_size);
	for (int i = 0; i < genes_ranges.size(); i++) {
		genes_ranges[i] = evaluator->iGetNumberOfValues(i);
	}
	this->random_values_holder = new RandomValuesHolder(genes_ranges, 1);
	best_individual = new Individual(this->evaluator, this->random_values_holder, this->genes_ranges.size());

	// tmp!!!

	vector<Individual* > loaded_individuals = vector<Individual* >(0);
	vector<vector<int>* > original_genotypes = vector<vector<int> * >(0);
	for (int i = 0; i < 26; i++) {
		loaded_individuals.push_back(new Individual(this->evaluator, this->random_values_holder, this->evaluator->iGetNumberOfBits()));
		loaded_individuals[loaded_individuals.size() - 1]->load_from_csv("C:\\Piotr\\2023_studia\\semestr3\\TEP\\AG\\logs\\saved_individuals\\best_solution_" + to_string(i) + ".txt");
		original_genotypes.push_back(loaded_individuals[loaded_individuals.size() - 1]->get_original_genotype());
	}

	levels.push_back(new P3Level(*loaded_individuals[0], &genes_ranges, *random_values_holder));
	delete loaded_individuals[0];
	for (int i = 1; i < loaded_individuals.size(); i++) {
		levels[0]->add_individual(*loaded_individuals[i]);
		delete loaded_individuals[i];
	}

	// end tmp

	creating_new_individuals_finished = 0;

	//for (int i = 0; i < this->threads_number; i++) {
	//	futures_Individuals.emplace_back(
	//		this->thread_pool->enqueue([this]() {
	//				return this->create_localy_optimized_individual(this->evaluator);
	//			}
	//		)
	//	);
	//}
}

void P3::run_iteration() {
	if (levels.size() == 0) {
		vector<Individual * > individuals = vector<Individual * >();
		for (int i = 0; i < futures_Individuals.size(); ++i) {
			individuals.push_back(futures_Individuals[i].get());
		}

		levels.push_back(new P3Level(*individuals[0], &genes_ranges, *random_values_holder));
		delete individuals[0];
		for (int i = 1; i < individuals.size(); i++) {
			run_individual_through_pyramid(individuals[i], true);
			delete individuals[i];
		}
	}
	else {
		//if (creating_new_individuals_finished.load() > 0) {
			Individual * random_individual_from_0_level = new Individual(*levels[0]->get_population()[
				levels[0]->get_random_values_holder()->get_random_individual_index()
			]); // do not delete - it is in the population, so while deleting population it will be deleted

			run_individual_through_pyramid(random_individual_from_0_level, false);
			delete random_individual_from_0_level;
		//}
		//else {
		//	vector<Individual * > individuals = vector<Individual * >();
		//	for (int i = 0; i < futures_Individuals.size(); ++i) {
		//		individuals.push_back(futures_Individuals[i].get());
		//	}

		//	for (int i = 0; i < individuals.size(); i++) {
		//		run_individual_through_pyramid(individuals[i], true);
		//		delete individuals[i];
		//	}

		//	int run_number_of_threads = this->threads_number - 1;
		//	for (int i = 0; i < run_number_of_threads; i++) {
		//		futures_Individuals.emplace_back(
		//			this->thread_pool->enqueue([this]() {
		//					return this->create_localy_optimized_individual(this->evaluator);
		//				}
		//			)
		//		);
		//	}
		//}
	}
}

void P3::run_individual_through_pyramid(Individual* individual, bool add_to_first_level) {
	check_and_actualise_best_individual(individual);

	bool add_to_next_level = add_to_first_level;
	for (int i = 0; i < levels.size(); i++) {
		Individual individual_copy = Individual(*individual);
		individual->set_random_values_holder(levels[i]->get_random_values_holder());
		//individual->mix_self_LinkageTree(levels[i]->get_linkage_tree(), levels[i]->get_population());
		// I changed to take linkage tree from level 0 !!!
		individual->mix_self_LinkageTree(levels[0]->get_linkage_tree(), levels[i]->get_population());

		if (add_to_next_level) {
			levels[i]->add_individual(individual_copy);
		}

		check_and_actualise_best_individual(individual);

		if (individual->get_fitness() > individual_copy.get_fitness()) {
			add_to_next_level = true;
		}
		else {
			add_to_next_level = false;
		}
	}

	if (add_to_next_level) {
		levels.push_back(new P3Level(*individual, &genes_ranges, *random_values_holder));
	}
}

void P3::check_and_actualise_best_individual(Individual* individual) {
	if (best_individual == nullptr) {
		best_individual = new Individual(*individual);
		best_individual->set_random_values_holder(random_values_holder);
	}
	else {
		if (individual->get_fitness() > best_individual->get_fitness()) {
			delete best_individual;
			best_individual = new Individual(*individual);
			best_individual->set_random_values_holder(random_values_holder);


			cout << endl << endl << "levels structure:" << endl;
			for (int i = 0; i < levels.size(); i++) {
				cout << "level " << i << endl;
				cout << "best individual: " << levels[i]->get_best_individual()->get_fitness() << endl;
				cout << "population size: " << levels[i]->get_population().size() << endl;
			}
		}
	}
}

Individual * P3::get_best_individual() {

	return best_individual;
}

vector<int> P3::get_best_solution() {
	return this->best_individual->get_genotype();
}

double P3::get_best_fitness() {
	return this->best_individual->get_fitness();
}

P3::~P3() {
	for (P3Level *level : levels) {
		delete level;
	}
	if (best_individual != nullptr) {
		delete best_individual;
	}
	delete random_values_holder;
	delete thread_pool;
}













// P3 level definitions
P3Level::P3Level(Individual& first_individual, vector<int>* genes_ranges, RandomValuesHolder& random_values_holder) {
	population = vector<Individual*>();
	population.push_back(new Individual(first_individual));
	best_individual = population[0];
	linkage_tree = nullptr;
	this->genes_ranges = genes_ranges;
	this->random_values_holder = new RandomValuesHolder(random_values_holder);
	population[0]->set_random_values_holder(this->random_values_holder);
	is_linkage_tree_actual = false;
}

P3Level::~P3Level() {
	for (Individual* individual : population) {
		delete individual;
	}
	if (linkage_tree != nullptr) {
		delete linkage_tree;
	}
	delete random_values_holder;
}

void P3Level::add_individual(Individual& individual) {
	population.push_back(new Individual(individual));
	Individual * new_individual = population[population.size() - 1];

	random_values_holder->set_population_size(population.size());

	new_individual->set_random_values_holder(random_values_holder);
	if (new_individual->get_fitness() > best_individual->get_fitness()) {
		best_individual = new_individual;
	}

	is_linkage_tree_actual = false;
}

void P3Level::self_actualise_linkage_tree() {
	if (linkage_tree != nullptr) {
		delete linkage_tree;
	}

	vector<vector<int> * > genotypes = vector<vector<int> * >();
	genotypes.reserve(population.size());
	for (int i = 0; i < population.size(); i++) {
		genotypes.push_back(population[i]->get_original_genotype());
	}

	linkage_tree = new LinkageTree(*genes_ranges);
	linkage_tree->build_tree(genotypes, *random_values_holder);
}

LinkageTree& P3Level::get_linkage_tree() {
	if (!is_linkage_tree_actual) {
		self_actualise_linkage_tree();
		is_linkage_tree_actual = true;
	}
	return *linkage_tree;
}

Individual* P3Level::get_best_individual() {
	return best_individual;
}

vector<Individual*>& P3Level::get_population() {
	return population;
}

RandomValuesHolder * P3Level::get_random_values_holder() {
	return random_values_holder;
}