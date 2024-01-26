#include "Individual.h"
#include <unordered_set>

Individual::Individual(CLFLnetEvaluator* evaluator, RandomValuesHolder* random_values_holder, int genotype_size) {
	this->evaluator = evaluator;
	this->random_values_holder = random_values_holder;
	this->genotype = vector<int> (genotype_size);
	this->fitness = 0;
	this->is_fitness_actual = false;

	for (int i = 0; i < genotype_size; i++) {
		this->genotype[i] = this->random_values_holder->get_random_gene_value(i);
	}
}

Individual::Individual(const Individual & other) {
	this->evaluator = other.evaluator;
	this->random_values_holder = other.random_values_holder;
	this->genotype = vector<int>(other.genotype);
	this->fitness = other.fitness;
	this->is_fitness_actual = other.is_fitness_actual;
}

void Individual::set_random_values_holder(RandomValuesHolder* random_values_holder) {
	this->random_values_holder = random_values_holder;
}

void Individual::set_evaluator(CLFLnetEvaluator* evaluator) {
	this->evaluator = evaluator;
}

void Individual::mutate(float mut_prob) {
	int mutated_genes = 0;
	for (int i = 0; i < this->genotype.size(); i++) {
		if (this->random_values_holder->get_random_probability() <= mut_prob) {
			mutated_genes++;
			this->genotype[i] = this->random_values_holder->get_random_gene_value(i);
		}
	}
	if (mutated_genes > 0) {
		this->is_fitness_actual = false;
	}
}

double Individual::get_fitness() {
	if (!this->is_fitness_actual) {
		this->fitness = this->evaluator->dEvaluate(&this->genotype);
		this->is_fitness_actual = true;
	}
	return this->fitness;
}

vector<int> Individual::get_genotype() {
	return vector<int>(this->genotype);
}

vector<int> * Individual::get_original_genotype() {
	return &this->genotype;
}

Individual::~Individual() {
}

vector<Individual*> Individual::cross_individual(Individual* individual_2, float cross_prob) {
	vector<Individual * > children = vector<Individual * >();
	children.reserve(2);
	children.push_back(new Individual(*this));
	children.push_back(new Individual(*individual_2));

	if (this->random_values_holder->get_random_probability() <= cross_prob) {
		int cross_point = this->random_values_holder->get_random_gene_index();

		for (int i = cross_point; i < this->genotype.size(); i++) {
			children[0]->genotype[i] = individual_2->genotype[i];
			children[1]->genotype[i] = this->genotype[i];
		}

		children[0]->is_fitness_actual = false;
		children[1]->is_fitness_actual = false;
	}

	return children;
}

vector<Individual* > Individual::cross_individual_with_cluster(Individual* individual_2, LinkageCluster& cluster) {
	Individual * child_1 = new Individual(*this);
	Individual * child_2 = new Individual(*individual_2);

	for (int i : cluster.get_indecies()) {
		child_1->genotype[i] = individual_2->genotype[i];
		child_2->genotype[i] = this->genotype[i];
	}

	child_1->is_fitness_actual = false;
	child_2->is_fitness_actual = false;

	vector<Individual* > children = vector<Individual* >();
	children.reserve(2);
	children.push_back(child_1);
	children.push_back(child_2);

	return children;

}

vector<Individual*> Individual::cross_individual_scattered(Individual* individual_2, float cross_prob) {
	vector<Individual* > children = vector<Individual* >();
	children.reserve(2);
	children.push_back(new Individual(*this));
	children.push_back(new Individual(*individual_2));

	if (this->random_values_holder->get_random_probability() <= cross_prob) {
		for (int i = 0; i < this->genotype.size(); i++) {
			if (this->random_values_holder->get_random_probability() > 0.5)
			{
				children[0]->genotype[i] = individual_2->genotype[i];
				children[1]->genotype[i] = this->genotype[i];
			}
		}

		children[0]->is_fitness_actual = false;
		children[1]->is_fitness_actual = false;
	}

	return children;
}

void Individual::load_from_csv(string filename) {

	std::vector<int> new_genotype;
	new_genotype.reserve(this->genotype.size());
	std::ifstream file(filename);
	std::string line;

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		std::string cell;

		while (std::getline(lineStream, cell, ',')) {
			new_genotype.push_back(std::stoi(cell));
		}
	}

	this->genotype = new_genotype;
}

void Individual::FHIC(int max_number_of_values_per_gene, int number_of_genes_to_check) {
	if (number_of_genes_to_check > genotype.size() || number_of_genes_to_check < 0) {
		number_of_genes_to_check = genotype.size();
	}
	
	vector<int> genes_indecies = vector<int>(genotype.size());
	for (int i = 0; i < genotype.size(); i++) {
		genes_indecies[i] = i;
	}

	random_values_holder->shuffle_vector(genes_indecies);

	
	for (int i = 0; i < number_of_genes_to_check; i++) {
		int gene_index = genes_indecies[i];
		int old_gene_value = genotype[gene_index];
		int number_of_values = evaluator->iGetNumberOfValues(gene_index);
		vector<int> possible_values_per_gene = vector<int>(number_of_values);
		for (int i = 0; i < number_of_values; i++) {
			possible_values_per_gene[i] = i;
		}
		PZ_Math::swap_and_remove(possible_values_per_gene, old_gene_value);
		random_values_holder->shuffle_vector(possible_values_per_gene);

		
		double old_fitness = this->get_fitness();
		bool improved = false;

		for (int j = 0; j < max_number_of_values_per_gene && !improved && possible_values_per_gene.size() > 0; j++) {
			int new_gene_value = possible_values_per_gene[possible_values_per_gene.size() - 1];
			possible_values_per_gene.pop_back();

			genotype[gene_index] = new_gene_value;
			is_fitness_actual = false;
			double new_fitness = this->get_fitness();

			if (new_fitness > old_fitness) {
				improved = true;
				fitness = new_fitness;
				is_fitness_actual = true;
			}
			else {
				fitness = old_fitness;
				is_fitness_actual = true;
				genotype[gene_index] = old_gene_value;
			}
		}
	}	
}

vector<Individual* > Individual::cross_individual_with_tree_randomly(
	Individual* individual_2,
	LinkageTree& tree,
	float change_genes_prob,
	float stop_considering_child_clusters_prob)
{
	Individual* child_1 = new Individual(*this);
	Individual* child_2 = new Individual(*individual_2);
	vector<Individual * > children = vector<Individual * >();
	children.push_back(child_1);
	children.push_back(child_2);

	vector<bool> change_genes = vector<bool>(this->genotype.size(), false);
	
	queue< pair<LinkageCluster*, LinkageCluster* > > clusters_to_consider;
	clusters_to_consider.push(
		make_pair(
			tree.get_cluster_of_whole_genotype()->get_child1(),
			tree.get_cluster_of_whole_genotype()->get_child2()
		)
	);

	while (clusters_to_consider.size() > 0) {
		LinkageCluster* clusters_1 = clusters_to_consider.front().first;
		LinkageCluster* clusters_2 = clusters_to_consider.front().second;
		clusters_to_consider.pop();

		// checking if we should randomly choose one of the clusters

		if (this->random_values_holder->get_random_probability() <= change_genes_prob) {
			LinkageCluster* cluster_to_change;
			if (this->random_values_holder->get_random_probability() <= 0.5) {
				cluster_to_change = clusters_1;
			}
			else {
				cluster_to_change = clusters_2;
			}

			for (int gene_index : cluster_to_change->get_indecies()) {
				change_genes[gene_index] = !change_genes[gene_index];
			}
		}

		// adding children to the queue
		if (clusters_1->get_child1() != nullptr && this->random_values_holder->get_random_probability() > stop_considering_child_clusters_prob) {  // it is enough to ckech one, cause either both are null or both are not null
			clusters_to_consider.push(
				make_pair(
					clusters_1->get_child1(),
					clusters_1->get_child2()
				)
			);
		}

		if (clusters_2->get_child1() != nullptr && this->random_values_holder->get_random_probability() > stop_considering_child_clusters_prob) {
			clusters_to_consider.push(
				make_pair(
					clusters_2->get_child1(),
					clusters_2->get_child2()
				)
			);
		}
	}


	//int counter_test_changed = 0;
	for (int i = 0; i < change_genes.size(); i++) {
		if (change_genes[i]) {
			//counter_test_changed++;
			child_1->genotype[i] = individual_2->genotype[i];
			child_2->genotype[i] = this->genotype[i];
		}
	}
	//cout << "changed " << counter_test_changed << endl;

	child_1->is_fitness_actual = false;
	child_2->is_fitness_actual = false;

	return children;
}


bool Individual::mix_self_LinkageTree(LinkageTree& linkage_tree, vector<Individual* >& population) {
	vector<LinkageCluster * > clusters = linkage_tree.get_clusters_ordered();
	vector<int> previous_genotype;
	double initial_fitness = this->fitness;

	for (int i = 0; i < clusters.size() && i < 300; i++) {
		LinkageCluster* cluster = clusters[i];
		previous_genotype = vector<int>(this->genotype);
		double previous_fitness = this->get_fitness();
		vector<int> cluster_indecies = cluster->get_indecies();
		//Individual * chosen_for_cross = population[this->random_values_holder->get_random_individual_index()];
		//
		//for (int gene_index : cluster_indecies) {
		//	this->genotype[gene_index] = chosen_for_cross->genotype[gene_index];
		//}

		//this->is_fitness_actual = false;
		//double new_fitness = this->get_fitness();


		//if (new_fitness <= previous_fitness) {
		//	this->genotype = move(previous_genotype);
		//	this->fitness = previous_fitness;
		//	this->is_fitness_actual = true;
		//}
		//else {
		//	previous_genotype.clear();
		//	this->fitness = new_fitness;
		//	this->is_fitness_actual = true;
		//	cout << "mixing with cluster " << i << "\t\tprevious fitness " << previous_fitness << "\t\tnew fitness " << new_fitness << endl;
		//	cout << "MIXED" << endl << endl << endl;
		//}
		vector<Individual * > population_copy = vector<Individual * >(population);
		this->random_values_holder->shuffle_vector(population_copy);
		bool improved = false;

		while (population_copy.size() > 0 && !improved) {
			Individual* chosen_for_cross = population_copy.back();
			population_copy.pop_back();
			for (int gene_index : cluster_indecies) {
				this->genotype[gene_index] = chosen_for_cross->genotype[gene_index];
			}

			this->is_fitness_actual = false;
			double new_fitness = this->get_fitness();


			if (new_fitness <= previous_fitness) {
				this->genotype = previous_genotype;
				this->fitness = previous_fitness;
				this->is_fitness_actual = true;
			}
			else {
				previous_genotype.clear();
				this->fitness = new_fitness;
				this->is_fitness_actual = true;
				cout << "mixing with cluster " << i << "\t\tprevious fitness " << previous_fitness << "\t\tnew fitness " << new_fitness << endl;
				cout << "MIXED" << endl << endl << endl;
				improved = true;
			}
		}

		if (!improved) {
			this->genotype = move(previous_genotype);
			this->fitness = previous_fitness;
			this->is_fitness_actual = true;
		}
	}

	if (this->fitness != initial_fitness) {
		cout << endl << "changed" << endl;
		return true; // indicates that it has changed
	}
	else {
		cout << endl << "not changed" << endl;
		return false;
	}
}


Individual* Individual::cross_individual_with_tree_greedy(
		Individual* individual_2,
		LinkageTree& tree,
		float choose_randomly_prob,
		float stop_considering_child_clusters_prob)
{
	Individual * child = new Individual(*this);
	queue< pair<LinkageCluster*, LinkageCluster* > > clusters_to_consider;
	clusters_to_consider.push(
		make_pair(
			tree.get_cluster_of_whole_genotype()->get_child1(),
			tree.get_cluster_of_whole_genotype()->get_child2()
		)
	);

	while (clusters_to_consider.size() > 0) {
		LinkageCluster * clusters_1 = clusters_to_consider.front().first;
		LinkageCluster * clusters_2 = clusters_to_consider.front().second;
		clusters_to_consider.pop();

		// checking if we should randomly choose one of the clusters
		if (this->random_values_holder->get_random_probability() <= choose_randomly_prob) {
			LinkageCluster* cluster_from_this;
			LinkageCluster* cluster_from_individual_2;
			if (this->random_values_holder->get_random_probability() <= 0.5) {
				cluster_from_this = clusters_1;
				cluster_from_individual_2 = clusters_2;
			}
			else {
				cluster_from_this = clusters_2;
				cluster_from_individual_2 = clusters_1;
			}

			for (int gene_index : cluster_from_this->get_indecies()) {
				child->genotype[gene_index] = this->genotype[gene_index];
			}

			for (int gene_index : cluster_from_individual_2->get_indecies()) {
				child->genotype[gene_index] = individual_2->genotype[gene_index];
			}

			child->is_fitness_actual = false;
		}
		else {
			Individual* child_tmp_version = new Individual(*child);

			for (int gene_index : clusters_1->get_indecies()) {
				child->genotype[gene_index] = this->genotype[gene_index];
				child_tmp_version->genotype[gene_index] = individual_2->genotype[gene_index];
			}

			for (int gene_index : clusters_2->get_indecies()) {
				child->genotype[gene_index] = individual_2->genotype[gene_index];
				child_tmp_version->genotype[gene_index] = this->genotype[gene_index];
			}
			child->is_fitness_actual = false;
			child_tmp_version->is_fitness_actual = false;
			if (child->get_fitness() < child_tmp_version->get_fitness()) {
				delete child;
				child = child_tmp_version;
			}
			else {
				delete child_tmp_version;
			}
		}



		// adding children to the queue
		if (clusters_1->get_child1() != nullptr && this->random_values_holder->get_random_probability() > stop_considering_child_clusters_prob) {  // it is enough to ckech one, cause either both are null or both are not null
			clusters_to_consider.push(
				make_pair(
					clusters_1->get_child1(),
					clusters_1->get_child2()
				)
			);
		}

		if (clusters_2->get_child1() != nullptr && this->random_values_holder->get_random_probability() > stop_considering_child_clusters_prob) {
			clusters_to_consider.push(
				make_pair(
					clusters_2->get_child1(),
					clusters_2->get_child2()
				)
			);
		}
	}

	return child;
}