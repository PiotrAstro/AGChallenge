#include "RandomValuesHolder.h"

RandomValuesHolder::RandomValuesHolder(vector<int> &values_range, int population_size) {
	this->values_range_distribution = vector<uniform_int_distribution<int> >();
	this->values_range_distribution.reserve(values_range.size());
	for (int i = 0; i < values_range.size(); i++) {
		this->values_range_distribution.push_back(uniform_int_distribution<int>(0, values_range[i] - 1));
	}

	this->gene_range_distribution = uniform_int_distribution<int>(0, values_range.size() - 1);
	this->probability_range_distribution = uniform_real_distribution<float>(0.0, 1.0);
	set_population_size(population_size);

	rand_engine = mt19937(std::chrono::high_resolution_clock::now().time_since_epoch().count() ^
		std::hash<std::thread::id>()(std::this_thread::get_id()));
}

RandomValuesHolder::RandomValuesHolder(const RandomValuesHolder & other) {
	this->values_range_distribution = vector<uniform_int_distribution<int> >(other.values_range_distribution);
	this->gene_range_distribution = uniform_int_distribution<int>(other.gene_range_distribution);
	this->probability_range_distribution = uniform_real_distribution<float>(0.0, 1.0);
	this->individual_range_distribution = uniform_int_distribution<int>(other.individual_range_distribution);
	rand_engine = mt19937(std::chrono::high_resolution_clock::now().time_since_epoch().count() ^
		std::hash<std::thread::id>()(std::this_thread::get_id()));
}

void RandomValuesHolder::set_population_size(int population_size) {
	this->individual_range_distribution = uniform_int_distribution<int>(0, population_size - 1);
}

RandomValuesHolder::~RandomValuesHolder() {
}

double RandomValuesHolder::get_random_probability() {
	return this->probability_range_distribution(this->rand_engine);
}

int RandomValuesHolder::get_random_individual_index() {
	return this->individual_range_distribution(this->rand_engine);
}

int RandomValuesHolder::get_random_gene_index() {
	return this->gene_range_distribution(this->rand_engine);
}

int RandomValuesHolder::get_random_gene_value(int gene_index) {
	return this->values_range_distribution[gene_index](this->rand_engine);
}

int RandomValuesHolder::get_random_int_from_0_to_n(int n) {
	uniform_int_distribution<int> range_to_n = uniform_int_distribution<int> (0, n);
	return range_to_n(this->rand_engine);
}