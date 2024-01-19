#pragma once
#include "Evaluator.h"
#include "RandomValuesHolder.h"
#include <fstream>
#include "LinkageTree.h"

using namespace std;

class Individual
{
private:
	vector<int> genotype;
	RandomValuesHolder * random_values_holder; // this one is passed, so it shouldnt be deleted
	double fitness;
	bool is_fitness_actual;
	CLFLnetEvaluator * evaluator; // this one is passed, so it shouldnt be deleted
public:
	Individual(CLFLnetEvaluator * evaluator, RandomValuesHolder * random_values_holder, int genotype_size);
	Individual(Individual const &other);
	void set_random_values_holder(RandomValuesHolder * random_values_holder);
	~Individual();

	// P3 operators
	void FHIC(int max_number_of_values_per_gene, int number_of_genes_to_check = -1); // First Improvement Hill Climber:  number_of_genes_to_check=-1 check all genes
	bool mix_self_LinkageTree(LinkageTree& linkage_tree, vector<Individual* >& population); // returns true if individual was changed, false otherwise
	// normal operators
	void mutate(float mut_prob);
	double get_fitness();
	vector<int> get_genotype();
	vector<int> * get_original_genotype();
	vector<Individual * > cross_individual(Individual* individual_2, float cross_prob); //returns array of 2 individuals
	vector<Individual* > cross_individual_scattered(Individual* individual_2, float cross_prob);
	vector<Individual* > cross_individual_with_cluster(Individual* individual_2, LinkageCluster & cluster); //returns array of 2 individuals
};

