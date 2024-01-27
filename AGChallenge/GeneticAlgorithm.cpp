#include "GeneticAlgorithm.h"

GeneticAlgorithm::GeneticAlgorithm(
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
) {
	this->evaluator = evaluator;
	this->param_verbose = verbose;
	this->param_population_size = population_size;
	this->param_cross_prob = cross_prob;
	this->param_mut_prob = mut_prob;
	this->param_choose_better_in_cross_prob = choose_better_in_cross_prob;	
	this->runned_iterations = 0;
	this->start_time = chrono::high_resolution_clock::now();
	this->param_linkage_tree_separation_size = linkage_tree_separation_size;
	this->param_linkage_tree_min_cluster = linkage_tree_min_cluster;
	this->param_linkage_tree_max_cluster = linkage_tree_max_cluster;
	this->param_use_generic_tree = use_generic_tree;
	this->best_individuals_counter = 1;
	this->last_improvement_iteration = 0;

	int genotype_size = evaluator->iGetNumberOfBits();
	this->values_range = vector<int> (genotype_size);
	for (int i = 0; i < evaluator->iGetNumberOfBits(); i++) {
		this->values_range[i] = evaluator->iGetNumberOfValues(i);
	}
	random_values_holder = new RandomValuesHolder(this->values_range, population_size);

	this->population = vector<Individual * > ();
	this->population.reserve(population_size);
	for (int i = 0; i < population_size; i++) {
		this->population.push_back(new Individual(this->evaluator, random_values_holder, genotype_size));
	}
	this->best_individual = new Individual(*this->population[0]);
	
	//actualise_best_individual();
	linkage_tree = new LinkageTree(this->values_range, linkage_tree_separation_size, linkage_tree_min_cluster, linkage_tree_max_cluster);
	linkage_tree->build_generic_tree(*this->random_values_holder);

	if (this->param_verbose) {
		cout << "GeneticAlgorithm created" << endl;
	}
}

double GeneticAlgorithm::get_run_time_milis() {
	chrono::high_resolution_clock::time_point current_time = chrono::high_resolution_clock::now();
	return chrono::duration_cast<chrono::milliseconds>(current_time - this->start_time).count();
}


GeneticAlgorithm::~GeneticAlgorithm() {
	delete this->best_individual;
	for (int i = 0; i < this->param_population_size; i++) {
		delete this->population[i];
	}

	delete this->random_values_holder;
	delete this->linkage_tree;

	if (this->param_verbose) {
		cout << "GeneticAlgorithm deleted" << endl;
	}
}

void GeneticAlgorithm::run_iteration() {
	// for real use:
	run_iteration_scattered();

	// for testing use:
	//vector<Individual * > loaded_individuals = vector<Individual * >();
	//loaded_individuals.push_back(new Individual(this->evaluator, this->random_values_holder, this->evaluator->iGetNumberOfBits()));
	//string file_path = "C:\\Piotr\\2023_studia\\semestr3\\TEP\\AG\\logs\\saved_individuals\\best_solution_" + to_string(tmp_test_loaded_individual_index) + ".txt";
	//loaded_individuals[loaded_individuals.size() - 1]->load_from_csv(file_path);
	//tmp_test_loaded_individual_index++;
	//tmp_test_loaded_individual_index %= 6;
	//run_iteration_LTGA(loaded_individuals);
}

void GeneticAlgorithm::insert_individuals(vector<Individual* > individuals_to_add) {
	for (Individual* individual : individuals_to_add) {
		individual->set_evaluator(this->evaluator);
		individual->set_random_values_holder(this->random_values_holder);
		int individual_index = get_random_individual_index_not_the_only_best();
		delete this->population[individual_index];
		this->population[individual_index] = individual;
	}
}

void GeneticAlgorithm::run_iteration_LTGA() {  // probably scatter inside!!!
	this->runned_iterations++;

	// this one is test
	int introduce_new_individuals = round(INTRODUCE_PERCENT_OF_NEW_INDIVIDUALS * this->param_population_size);
	//int introduce_new_individuals = round((0.1 - mean_population_entropy()) * 20);
	//if (introduce_new_individuals > 20) {
	//	introduce_new_individuals = 20;
	//}
	//else if (introduce_new_individuals < 1) {
	//	introduce_new_individuals = 1;
	//}
	introduce_n_new_individuals(introduce_new_individuals);

	if (this->runned_iterations - this->last_improvement_iteration > RESET_POPULATION_AFTER_NO_IMPROVEMENT_FOR) {
		introduce_n_new_individuals(this->param_population_size);
		this->last_improvement_iteration = this->runned_iterations;
	}

	perform_LTGA_routine();  // probably scatter inside!!!
	//vector<Individual* > new_population = generate_crossed_population();
	//for (int i = 0; i < this->param_population_size; i++) {
	//	delete this->population[i];
	//	new_population[i]->mutate(this->param_mut_prob);
	//	this->population[i] = new_population[i];
	//}

	//if (this->runned_iterations % LTGA_FHIC_EVERY_N_ITERATIONS == 0) {
	//	bool stop = false;
	//	for (int i = 0; i < this->param_population_size && !stop; i++) {
	//		if (this->population[i]->get_fitness() >= this->best_individual->get_fitness()) {
	//			this->population[i]->FHIC(FHIC_TRY_N_GENE_VALUES, -1);
	//			stop = true;
	//		}
	//	}
	//}

	actualise_best_individual();

	if (this->runned_iterations % STATS_AND_LINKAGE_TREE_EVERY_N_ITERATIONS == 0) {
		if (this->param_verbose) {
			print_iteration_summary();
		}

		if (this->param_use_generic_tree) {
			linkage_tree->build_generic_tree(*this->random_values_holder);
		}
		else {
			linkage_tree->build_tree(get_population_genotypes(), *this->random_values_holder);
		}
	}
	else {
		if (this->param_verbose) {
			cout << "Iteration: " << this->runned_iterations << endl;
		}
	}
}

void GeneticAlgorithm::introduce_n_new_individuals(int introduce_new_individuals) {
	for (int i = 0; i < introduce_new_individuals; i++) {
		int individual_index = get_random_individual_index_not_the_only_best();
		delete this->population[individual_index];
		this->population[individual_index] = new Individual(this->evaluator, random_values_holder, this->evaluator->iGetNumberOfBits());
	}
}

void GeneticAlgorithm::run_iteration_scattered() {
	this->runned_iterations++;

	int introduce_new_individuals = round(INTRODUCE_PERCENT_OF_NEW_INDIVIDUALS * this->param_population_size);
	introduce_n_new_individuals(introduce_new_individuals);

	//if (this->runned_iterations - this->last_improvement_iteration > RESET_POPULATION_AFTER_NO_IMPROVEMENT_FOR) {
	//	introduce_n_new_individuals(this->param_population_size);
	//	this->last_improvement_iteration = this->runned_iterations;
	//}

	perform_scatter_routine();
	actualise_best_individual();

	if (this->param_verbose) {
		if (this->runned_iterations % STATS_AND_LINKAGE_TREE_EVERY_N_ITERATIONS == 0) {
			print_iteration_summary();
		}
		else {
			cout << "Iteration: " << this->runned_iterations << endl;
		}
	}
}

void GeneticAlgorithm::perform_LTGA_routine() {  // I probably do here scatter!!!
	for (int i = 0; i < LTGA_ROUTINE_CROSSES_PER_ITERATION; i++) {
		int individual_1 = get_random_individual_index_after_fight_of_2();
		int individual_2 = get_random_individual_index_after_fight_of_2();

		while (individual_2 == individual_1) {
			individual_2 = get_random_individual_index_after_fight_of_2();
		}

		// for fully random use:

		//vector<Individual* > crossed_children = this->population[individual_1]->cross_individual_with_tree_randomly(
		//	this->population[individual_2],
		//	*this->linkage_tree,
		//	TREE_CROSS_RANDOM_CHANGE_GENES_PROB,
		//	TREE_CROSS_STOP_DELVING_DEEPER_PROB
		//);

		//this is test!!! it is actually scatter!!!
		//vector<Individual* > crossed_children = this->population[individual_1]->cross_individual_scattered(this->population[individual_2], this->param_cross_prob);

		//for (Individual* child : crossed_children) {
		//	child->mutate(this->param_mut_prob);
		//}

		//delete this->population[individual_1];
		//delete this->population[individual_2];
		//this->population[individual_1] = crossed_children[0];
		//this->population[individual_2] = crossed_children[1];

		//int worse_individual = this->population[individual_1]->get_fitness() < this->population[individual_2]->get_fitness() ? individual_1 : individual_2;
		//int better_child = crossed_children[0]->get_fitness() > crossed_children[1]->get_fitness() ? 0 : 1;

		//if (this->random_values_holder->get_random_probability() <= LTGA_CHILD_FHIC_PROB) {
		//	crossed_children[better_child]->FHIC(FHIC_TRY_N_GENE_VALUES, -1);
		//}

		//delete this->population[worse_individual];
		//delete crossed_children[1 - better_child];
		//this->population[worse_individual] = crossed_children[better_child];
		
		// for greedy use:

		//Individual* crossed_child = this->population[individual_1]->cross_individual_with_tree_greedy(
		//		this->population[individual_2],
		//		*(this->linkage_tree),
		//		TREE_CROSS_GREEDY_CHOOSE_RANDOMLY_PROB,
		//		TREE_CROSS_STOP_DELVING_DEEPER_PROB
		//	);
		Individual* crossed_child = this->population[individual_1]->cross_individual_scattered_greedy(
			this->population[individual_2],
			this->param_cross_prob,
			(this->random_values_holder->get_random_probability() <= LTGA_CHILD_PROB_FHIC_WITH_CROSSED ? 1.0 : 0.0)
		);
		crossed_child->mutate(this->param_mut_prob);

		int worse_individual = this->population[individual_1]->get_fitness() < this->population[individual_2]->get_fitness() ? individual_1 : individual_2;

		delete this->population[worse_individual];
		this->population[worse_individual] = crossed_child;
	}
}

void GeneticAlgorithm::perform_scatter_routine() {
	for (int i = 0; i < SCATTERED_ROUTINE_CROSSES_PER_ITERATION; i++) {
		int individual_1 = get_random_individual_index_after_fight_of_2();
		int individual_2 = get_random_individual_index_after_fight_of_2();

		while (individual_2 == individual_1) {
			individual_2 = get_random_individual_index_after_fight_of_2();
		}

		vector<Individual* > crossed_children = this->population[individual_1]->cross_individual_scattered(this->population[individual_2], this->param_cross_prob);
		for (Individual* child : crossed_children) {
			child->mutate(this->param_mut_prob);
		}

		//// test!!!
		//delete this->population[individual_1];
		//if (individual_1 != individual_2) {
		//	delete this->population[individual_2];
		//}
		//this->population[individual_1] = crossed_children[0];
		//this->population[individual_2] = crossed_children[1];

		// this was used for best results!!!
		int worse_individual = this->population[individual_1]->get_fitness() < this->population[individual_2]->get_fitness() ? individual_1 : individual_2;
		int better_child = crossed_children[0]->get_fitness() > crossed_children[1]->get_fitness() ? 0 : 1;

		delete this->population[worse_individual];
		delete crossed_children[1 - better_child];
		this->population[worse_individual] = crossed_children[better_child];

		//if(crossed_children[0]->get_fitness() > this->population[individual_1]->get_fitness()) {
		//	delete this->population[individual_1];
		//	this->population[individual_1] = crossed_children[0];
		//}
		//else {
		//	delete crossed_children[0];
		//}

		//if (crossed_children[1]->get_fitness() > this->population[individual_2]->get_fitness()) {
		//	delete this->population[individual_2];
		//	this->population[individual_2] = crossed_children[1];
		//}
		//else {
		//	delete crossed_children[1];
		//}
	}
}

vector< Individual* > GeneticAlgorithm::generate_crossed_population() {
	vector<Individual * > new_population = vector<Individual * > ();
	new_population.reserve(this->param_population_size);

	int current_new_population_index = 0;
	while (current_new_population_index < this->param_population_size) {
		Individual * individual_1 = get_random_individual_after_fight_of_2();
		Individual * individual_2 = get_random_individual_after_fight_of_2();

		vector<Individual* > crossed_children = individual_1->cross_individual_scattered(individual_2, this->param_cross_prob);
		int children_current_index = 0;
		while (current_new_population_index < this->param_population_size && children_current_index < crossed_children.size()) {
			new_population.push_back(crossed_children[children_current_index]);
			current_new_population_index++;
			children_current_index++;
		}
		for (; children_current_index < crossed_children.size(); children_current_index++) {
			delete crossed_children[children_current_index];
		}
	}

	return new_population;
}

vector< Individual* > GeneticAlgorithm::generate_crossed_every_cluster_population() {
	vector<Individual* > new_population = vector<Individual* >();
	new_population.reserve(linkage_tree->get_clusters_ordered().size() * 2);

	for (LinkageCluster* cluster : linkage_tree->get_clusters_ordered()) {
		Individual * individual_1 = get_random_individual_after_fight_of_2();
		Individual * individual_2 = get_random_individual_after_fight_of_2();

		vector<Individual* > crossed_children = individual_1->cross_individual_with_cluster(individual_2, *cluster);
		for (Individual* child : crossed_children) {
			new_population.push_back(child);
		}
	}

	return new_population;
}

vector<Individual* > GeneticAlgorithm::generate_crossed_random_cluster_population_size() {
	vector<Individual* > new_population = vector<Individual* >();
	new_population.reserve(this->param_population_size);

	vector<LinkageCluster *> clusters = vector<LinkageCluster*> (linkage_tree->get_clusters_ordered());
	random_values_holder->shuffle_vector(clusters);

	int current_new_population_index = 0;
	while (current_new_population_index < this->param_population_size) {
		LinkageCluster* cluster = clusters.back();
		clusters.pop_back();
		Individual* individual_1 = get_random_individual_after_fight_of_2();
		Individual* individual_2 = get_random_individual_after_fight_of_2();

		vector<Individual* > crossed_children = individual_1->cross_individual_with_cluster(individual_2, *cluster);
		int children_current_index = 0;
		while (current_new_population_index < this->param_population_size && children_current_index < crossed_children.size()) {
			new_population.push_back(crossed_children[children_current_index]);
			current_new_population_index++;
			children_current_index++;
		}
		for (; children_current_index < crossed_children.size(); children_current_index++) {
			delete crossed_children[children_current_index];
		}
	}

	return new_population;

	return new_population;
}

Individual* GeneticAlgorithm::get_random_individual_after_fight_of_2() {
	return this->population[get_random_individual_index_after_fight_of_2()];
}

int GeneticAlgorithm::get_random_individual_index_after_fight_of_2() {
	int individual_1 = this->random_values_holder->get_random_individual_index();
	int individual_2 = this->random_values_holder->get_random_individual_index();

	bool fight_result = this->param_choose_better_in_cross_prob >= this->random_values_holder->get_random_probability();
	int chosen_individual;
	if (fight_result) {
		if (this->population[individual_1]->get_fitness() > this->population[individual_2]->get_fitness()) {
			chosen_individual = individual_1;
		}
		else {
			chosen_individual = individual_2;
		}
	}
	else {
		bool choose_first = this->random_values_holder->get_random_probability() >= 0.5;
		if (choose_first) {
			chosen_individual = individual_1;
		}
		else {
			chosen_individual = individual_2;
		}
	}

	return chosen_individual;
}

void GeneticAlgorithm::actualise_best_individual() {
	this->best_individuals_counter = 0;

	for (int i = 0; i < this->param_population_size; i++) {
		if (this->population[i]->get_fitness() > this->best_individual->get_fitness()) {
			delete this->best_individual;

			if (this->param_verbose) {
				cout << "New best, FHIC" << endl;
			}

			this->population[i]->FHIC(FHIC_TRY_N_GENE_VALUES, -1);

			if (this->param_verbose) {
				cout << " done FHIC" << endl;
			}

			this->best_individual = new Individual(*this->population[i]);

			this->best_individuals_counter = 1;
			this->last_improvement_iteration = this->runned_iterations;
		}
		else if (this->population[i]->get_fitness() == this->best_individual->get_fitness()) {
			this->best_individuals_counter++;
		}
	}
}

void GeneticAlgorithm::print_iteration_summary() {
	vector<double> quantiles = fitness_quantiles(vector<float>({ 0.25, 0.5, 0.75, 1}));
	cout << "Time: " << get_run_time_milis() / 1000.0 << " s"
		<< "\nIteration: " << this->runned_iterations
		<< "\nBest individual with fitness: " << this->best_individual->get_fitness()
		<< "\nQuantiles:\t(0.25) " << quantiles[0] << "\t(0.5) " << quantiles[1] << "\t(0.75) " << quantiles[2] << "\t(1) " << quantiles[3]
		<< "\nmean entropy: " << mean_population_entropy()
		<< "\n\n";
}

Individual* GeneticAlgorithm::get_best_individual() {
	return this->best_individual;
}

vector<int> GeneticAlgorithm::get_best_solution() {
	return this->best_individual->get_genotype();
}

vector<vector<int> * > GeneticAlgorithm::get_population_genotypes() {
	vector<vector<int> * > population_genotypes = vector<vector<int> * >();
	population_genotypes.reserve(this->param_population_size);
	for (int i = 0; i < this->param_population_size; i++) {
		population_genotypes.push_back(this->population[i]->get_original_genotype());
	}
	return population_genotypes;
}

double GeneticAlgorithm::get_best_fitness() {
	return this->best_individual->get_fitness();
}


double GeneticAlgorithm::mean_population_entropy() {
	return mean_population_entropy(population);
}

double GeneticAlgorithm::mean_population_entropy(const vector<Individual* > & individuals) {
	if (individuals.size() > 0) {
		double entropy_sum = 0;
		int genes_number = individuals[0]->get_original_genotype()->size();

		for (int i = 0; i < genes_number; i++) {
			vector <int> values_count = vector<int>(individuals[0]->get_evaluator()->iGetNumberOfValues(i));
			for (int j = 0; j < individuals.size(); j++) {
				values_count[individuals[j]->get_original_genotype()->at(i)]++;
			}

			for (int j = 0; j < values_count.size(); j++) {
				float probability = (float)values_count[j] / (float)individuals.size();
				entropy_sum -= (double)PZ_Math::prob_mul_log(probability);
			}
		}
		return entropy_sum / (double)genes_number;
	}
	else {
		return 0;
	}
}

vector<double> GeneticAlgorithm::fitness_quantiles(vector<float> quantiles) {
	vector<double> fitnesses = vector<double> (this->param_population_size);
	vector<double> quantiles_values = vector<double>(quantiles.size());

	for (int i = 0; i < this->param_population_size; i++) {
		fitnesses[i] = this->population[i]->get_fitness();
	}

	sort(fitnesses.begin(), fitnesses.end());

	for (int i = 0; i < quantiles.size(); i++) {
		int quantile_index = (int)(quantiles[i] * (double)(this->param_population_size - 1));
		quantiles_values[i] = fitnesses[quantile_index];
	}

	return quantiles_values;
}

int GeneticAlgorithm::get_runned_iterations() {
	return this->runned_iterations;
}

int GeneticAlgorithm::get_population_size() {
	return this->param_population_size;
}

float GeneticAlgorithm::get_cross_prob() {
	return this->param_cross_prob;
}

float GeneticAlgorithm::get_mut_prob() {
	return this->param_mut_prob;
}

float GeneticAlgorithm::get_choose_better_in_cross_prob() {
	return this->param_choose_better_in_cross_prob;
}

int GeneticAlgorithm::get_random_individual_index_not_the_only_best() {
	int individual_index = this->random_values_holder->get_random_individual_index();
	if (this->best_individuals_counter == 1) {
		while (this->population[individual_index]->get_fitness() == this->best_individual->get_fitness()) {
			individual_index = this->random_values_holder->get_random_individual_index();
		}
	}
	else if (this->population[individual_index]->get_fitness() == this->best_individual->get_fitness()) {
		this->best_individuals_counter --;
	}
	return individual_index;
}

RandomValuesHolder* GeneticAlgorithm::get_random_values_holder() {
	return this->random_values_holder;
}