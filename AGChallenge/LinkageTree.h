#pragma once
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include "RandomValuesHolder.h"
#include "PZ_Math.h"

using namespace std;

class LinkageCluster {
private:
	vector<int> linked_indecies;

	priority_queue<pair<float, LinkageCluster **>, vector<pair<float, LinkageCluster**> >, greater<> > * min_distance_to_clusters;  // min distance to other clusters, and pointer to them
	int tmp_index_in_vector;
	LinkageCluster* tmp_self_in_priority_queue;
	void heal_tmp_min_distance();

	LinkageCluster* child1;
	LinkageCluster* child2;
	int id;
public:
	LinkageCluster(vector<int> linked_indecies, int id, LinkageCluster* child1, LinkageCluster* child2);
	~LinkageCluster();
	vector<int>& get_indecies();
	LinkageCluster* get_child1();
	LinkageCluster* get_child2();
	int get_id();

	//tmp:
	void clear_distance_data();
	void add_tmp_distance(float tmp_distance, LinkageCluster* tmp_min_distance_to);
	void set_tmp_index_in_vector(int tmp_index_in_vector);
	float get_tmp_min_distance();
	LinkageCluster* get_tmp_min_distance_to();
	int get_tmp_index_in_vector();
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
		LinkageCluster* get_cluster_of_whole_genotype();
		void build_tree(const vector<vector<int> * > & genotypes, RandomValuesHolder& random_values_holder);  // it is effective, it cheks with flag if it has to actualise self
		void build_generic_tree(RandomValuesHolder& random_values_holder); // builds generic tree, simply doesnt take linkage information into account
		vector<LinkageCluster* > & get_clusters_ordered();
};
