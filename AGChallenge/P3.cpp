#include "P3.h"




Individual P3::create_localy_optimized_individual(CLFLnetEvaluator& evaluator) {
	CLFLnetEvaluator evaluator_copy = CLFLnetEvaluator();
	evaluator_copy.bConfigure(evaluator.sGetNetName());

	int param_population_size = 10;
	float param_cross_prob = 0.8;
	float param_mut_prob = 0.001;
	float param_choose_better_in_cross_prob = 0.5;
	int param_linkage_tree_separation_size = 100;
	int param_linkage_tree_min_cluster = 3;
	int param_linkage_tree_max_cluster = 99;
	bool param_use_generic_tree = true;


	GeneticAlgorithm evolutionary_algorithm(
		&evaluator_copy,
		param_population_size,
		param_cross_prob,
		param_mut_prob,
		param_choose_better_in_cross_prob,
		param_linkage_tree_separation_size,
		param_linkage_tree_min_cluster,
		param_linkage_tree_max_cluster,
		param_use_generic_tree
	);

	while()

}





P3::P3(CLFLnetEvaluator* evaluator) {
	this->evaluator = evaluator;

	best_individual = nullptr;
	levels = vector<P3Level* >();
	int genotype_size = evaluator->iGetNumberOfBits();
	genes_ranges = vector<int> (genotype_size);
	for (int i = 0; i < genes_ranges.size(); i++) {
		genes_ranges[i] = evaluator->iGetNumberOfValues(i);
	}
	random_values_holder = new RandomValuesHolder(genes_ranges, 1);
}

void P3::run_iteration() {
	Individual * new_individual = new Individual(evaluator, random_values_holder, genes_ranges.size());
	//new_individual->FHIC(param_FHIC_max_number_of_values_per_gene, param_FHIC_genes_checked);

	if (levels.size() == 0) {
		levels.push_back(new P3Level(*new_individual, &genes_ranges, *random_values_holder));
		best_individual = new Individual(*new_individual);
	}
	else {
		check_and_actualise_best_individual(new_individual);

		bool add_to_next_level = true;
		for (int i = 0; i < levels.size(); i++) {
			Individual new_individual_copy = Individual(*new_individual);
			new_individual->set_random_values_holder(levels[i]->get_random_values_holder());
			new_individual->mix_self_LinkageTree(levels[i]->get_linkage_tree(), levels[i]->get_population());
			
			if (add_to_next_level) {
				levels[i]->add_individual(new_individual_copy);
			}

			check_and_actualise_best_individual(new_individual);

			if (new_individual->get_fitness() > new_individual_copy.get_fitness()) {
				add_to_next_level = true;
			}
			else {
				add_to_next_level = false;
			}
		}

		if (add_to_next_level) {
			levels.push_back(new P3Level(*new_individual, &genes_ranges, *random_values_holder));
		}
	}
	delete new_individual;
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