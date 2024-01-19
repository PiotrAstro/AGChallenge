#pragma once
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "RandomValuesHolder.h"
#include "PZ_Math.h"

#define DEFAULT_SEPARATION_SIZE 100  // considering matrix of 2500x2500 is to big, so I consider (2500 / SEPARATION_SIZE) * (SEPARATION_SIZE * SEPARATION_SIZE)
// therefore I do not catch all distances, but I think that iterative approach will allow different genes to be combined with different ones 

// Originally it was set to 2
#define DEFAULT_MIN_CLUSTER_SIZE 3 // between 1 and inf, it includes clusters of this size, so it is >=
#define DEFAULT_MAX_CLUSTER_SIZE 99 // between 1 and inf, it is included in the range, so if max is set to 100, then 100 will alsobe included



using namespace std;

class LinkageCluster {
private:
	vector<int> linked_indecies;
	LinkageCluster* child1;
	LinkageCluster* child2;
	int id;
public:
	LinkageCluster(vector<int> linked_indecies, int id, LinkageCluster* child1, LinkageCluster* child2);
	~LinkageCluster();
	vector<int>& get_indecies();
	int get_id();
};

class LinkageTree
{
	private:
		int param_separation_size;
		int param_min_cluster_size;
		int param_max_cluster_size;


		vector<LinkageCluster* > usable_clusters_ordered;
		LinkageCluster* cluster_of_whole_genotype;
		vector< int > genes_range;
		int get_distance_id(LinkageCluster* cluster1, LinkageCluster* cluster2);
	public:
		LinkageTree(vector<int> & genes_ranges);
		LinkageTree(vector<int>& genes_ranges, int separation_size, int min_cluster_size, int max_cluster_size);
		~LinkageTree();
		void build_tree(const vector<vector<int> * > & genotypes, RandomValuesHolder& random_values_holder);  // it is effective, it cheks with flag if it has to actualise self
		void build_generic_tree(RandomValuesHolder& random_values_holder); // builds generic tree, simply doesnt take linkage information into account
		vector<LinkageCluster* > & get_clusters_ordered();
};
