#pragma once
#include "Individual.h"
#include "Evaluator.h"
#include "RandomValuesHolder.h"
#include "LinkageTree.h"
#include "AbstractEvolutionaryAlgorithm.h"

class P3Level {
private:
	Individual* best_individual; // do not delete - it is in the population, so while deleting population it will be deleted
	vector<Individual* > population;
	vector<int>* genes_ranges; // do not delete - it is in the population, so while deleting population it will be deleted
	LinkageTree* linkage_tree;
	RandomValuesHolder * random_values_holder;
	bool is_linkage_tree_actual;
public:
	P3Level(Individual & first_individual, vector<int>* genes_ranges, RandomValuesHolder& random_values_holder);
	~P3Level();
	Individual* get_best_individual();
	void add_individual(Individual& individual);
	vector<Individual* >& get_population();
	void self_actualise_linkage_tree();
	LinkageTree& get_linkage_tree();
	RandomValuesHolder * get_random_values_holder();
};


class P3 : public AbstractEvolutionaryAlgorithm 
{
private:
	int param_FHIC_genes_checked;
	int param_FHIC_max_number_of_values_per_gene;

	Individual * best_individual;
	CLFLnetEvaluator* evaluator; // this one is passed, so it shouldnt be deleted
	vector<int> genes_ranges;
	vector<P3Level* > levels;
	RandomValuesHolder * random_values_holder;
	void check_and_actualise_best_individual(Individual* individual);
public:
	P3(CLFLnetEvaluator* evaluator, int param_FHIC_genes_checked, int param_FHIC_max_number_of_values_per_gene);
	~P3();
	void run_iteration();
	vector<int> get_best_solution();
	double get_best_fitness();
	Individual * get_best_individual();
};