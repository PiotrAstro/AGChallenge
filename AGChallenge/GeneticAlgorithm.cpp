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
	bool use_generic_tree
) {
	this->evaluator = evaluator;
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
	cout << "GeneticAlgorithm deleted" << endl;
}

void GeneticAlgorithm::run_iteration() {
	this->runned_iterations++;
	//vector<Individual * > new_population = generate_crossed_random_cluster_population_size();
	////new_population.reserve(this->population_size + new_population.size());
	//for (Individual* individual : this->population) {
	//	//new_population.push_back(individual);
	//	delete individual;
	//}

	//for (int i = 0; i < new_population.size(); i++) {
	//	new_population[i]->mutate(this->mut_prob);
	//}
	//this->population = new_population;

	////sort(new_population.begin(), new_population.end(), [](Individual* a, Individual* b) {return a->get_fitness() > b->get_fitness(); });
	////for(int i = new_population.size() - 1; i >= population_size; i--) {
	////	delete new_population[i];
	////	new_population.pop_back();
	////}
	// 
	// 
	perform_LTGA_routine();

	//test!!!
	int introduce_new_individuals = round(INTRODUCE_PERCENT_OF_NEW_INDIVIDUALS * this->param_population_size);

	for (int i = 0; i < introduce_new_individuals; i++) {
		int individual_index = this->random_values_holder->get_random_individual_index();
		delete this->population[individual_index];
		this->population[individual_index] = new Individual(this->evaluator, random_values_holder, this->evaluator->iGetNumberOfBits());
	}

	if (this->runned_iterations % STATS_AND_LINKAGE_TREE_EVERY_N_ITERATIONS == 0) {

		//// fo previous great results I didnt do FHIC here!!! test!!!
		//Individual* individual_to_FHIC = this->population[0];
		//for (Individual * individual : this->population) {
		//	if (individual->get_fitness() > individual_to_FHIC->get_fitness()) {
		//		individual_to_FHIC = individual;
		//	}
		//}
		////Individual* individual_to_FHIC = this->population[this->random_values_holder->get_random_individual_index()];
		//individual_to_FHIC->FHIC(30, -1);

		print_iteration_summary();
		//linkage_tree->build_tree(get_population_genotypes(), *this->random_values_holder);

		if (this->param_use_generic_tree) {
			linkage_tree->build_generic_tree(*this->random_values_holder);
		}
		else {
			linkage_tree->build_tree(get_population_genotypes(), *this->random_values_holder);
		}
	}
	//this->population[0]->FHIC(10, -1);

	actualise_best_individual();
}

void GeneticAlgorithm::perform_LTGA_routine() {
	for (LinkageCluster* cluster : linkage_tree->get_clusters_ordered()) {
		int individual_1 = get_random_individual_index_after_fight_of_2();
		int individual_2 = get_random_individual_index_after_fight_of_2();

		vector<Individual* > crossed_children = this->population[individual_1]->cross_individual_with_cluster(this->population[individual_2], *cluster);	
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

		if (crossed_children[better_child]->get_fitness() >= this->population[worse_individual]->get_fitness()) {
			delete this->population[worse_individual];
			delete crossed_children[1 - better_child];
			this->population[worse_individual] = crossed_children[better_child];
		}
		else {
			delete crossed_children[better_child];
			delete crossed_children[1 - better_child];
		}

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
	for (int i = 0; i < this->param_population_size; i++) {
		if (this->population[i]->get_fitness() > this->best_individual->get_fitness()) {
			delete this->best_individual;
			this->population[i]->FHIC(10, -1);
			this->best_individual = new Individual(*this->population[i]);
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
	double entropy_sum = 0;

	for (int i = 0; i < this->values_range.size(); i++) {
		vector <int> values_count = vector<int>(this->values_range[i]);
		for (int j = 0; j < this->param_population_size; j++) {
			values_count[this->population[j]->get_genotype()[i]]++;
		}

		for (int j = 0; j < this->values_range[i]; j++) {
			float probability = (float)values_count[j] / (float)this->param_population_size;
			entropy_sum -= (double)PZ_Math::prob_mul_log(probability);
		}
	}

	return entropy_sum / (double)this->values_range.size();
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


