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

Individual::Individual(Individual const& other) {
	this->evaluator = other.evaluator;
	this->random_values_holder = other.random_values_holder;
	this->genotype = vector<int>(other.genotype);
	this->fitness = other.fitness;
	this->is_fitness_actual = other.is_fitness_actual;
}

void Individual::set_random_values_holder(RandomValuesHolder* random_values_holder) {
	this->random_values_holder = random_values_holder;
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
		unordered_set<int> checked_values = unordered_set<int>();
		int gene_index = genes_indecies[i];
		int old_gene_value = genotype[gene_index];
		double old_fitness = this->get_fitness();
		checked_values.insert(old_gene_value);
		bool improved = false;

		for (int j = 0; j < max_number_of_values_per_gene && !improved; j++) {
			int new_gene_value = random_values_holder->get_random_gene_value(gene_index);
			auto result = checked_values.insert(new_gene_value);
			if (result.second) {  // if new_gene_value was inserted
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
}


bool Individual::mix_self_LinkageTree(LinkageTree& linkage_tree, vector<Individual* >& population) {
	vector<LinkageCluster * > clusters = linkage_tree.get_clusters_ordered();
	vector<int> previous_genotype;

	for (int i = 0; i < clusters.size() && i < 300; i++) {
		LinkageCluster* cluster = clusters[i];
		previous_genotype = vector<int>(this->genotype);
		double previous_fitness = this->get_fitness();
		vector<int> cluster_indecies = cluster->get_indecies();
		Individual * chosen_for_cross = population[this->random_values_holder->get_random_individual_index()];

		for (int gene_index : cluster_indecies) {
			this->genotype[gene_index] = chosen_for_cross->genotype[gene_index];
		}

		this->is_fitness_actual = false;
		//double new_fitness = this->get_fitness();

		//if (new_fitness < previous_fitness) {
		//	this->genotype = move(previous_genotype);
		//	this->is_fitness_actual = true;
		//}
		//else {
		//	previous_genotype.clear();
		//	this->fitness = new_fitness;
		//	this->is_fitness_actual = true;
		//}
	}
}