#pragma once
#include <random>
#include <vector>
#include <iostream>
#include <time.h>

using namespace std;

class RandomValuesHolder
{
	private:
		vector<uniform_int_distribution<int>> values_range_distribution;
		uniform_int_distribution<int> gene_range_distribution;
		uniform_real_distribution<float> probability_range_distribution;
		uniform_int_distribution<int> individual_range_distribution;
		mt19937 rand_engine;
	public:
		RandomValuesHolder(vector<int> &values_range, int population_size);
		RandomValuesHolder(const RandomValuesHolder& other);
		~RandomValuesHolder();
		void set_population_size(int population_size);
		double get_random_probability();
		int get_random_individual_index();
		int get_random_gene_index();
		int get_random_gene_value(int gene_index);

		template <typename T>
		void shuffle_vector(vector<T>& to_shuffle) {
			shuffle(to_shuffle.begin(), to_shuffle.end(), rand_engine);
		}
};

