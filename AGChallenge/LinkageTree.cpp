#include "LinkageTree.h"


LinkageTree::LinkageTree(vector<int>& genes_ranges) {
	// initial creations
	this->genes_range = vector<int>(genes_ranges);
	
	cluster_of_whole_genotype = nullptr;
	this->param_separation_size = DEFAULT_SEPARATION_SIZE;
	this->param_min_cluster_size = DEFAULT_MIN_CLUSTER_SIZE;
	this->param_max_cluster_size = DEFAULT_MAX_CLUSTER_SIZE;
}

LinkageTree::LinkageTree(vector<int>& genes_ranges, int separation_size, int min_cluster_size, int max_cluster_size) {
	// initial creations
	this->genes_range = vector<int>(genes_ranges);

	cluster_of_whole_genotype = nullptr;
	this->param_separation_size = separation_size;
	this->param_min_cluster_size = min_cluster_size;
	this->param_max_cluster_size = max_cluster_size;
}

void LinkageTree::build_tree(const vector<vector<int> * > & genotypes, RandomValuesHolder& random_values_holder) {
	if (cluster_of_whole_genotype != nullptr) {
		delete cluster_of_whole_genotype;
	}

	// initial actions
	int new_clusters_counter = 0;

	usable_clusters_ordered = vector<LinkageCluster* >();
	vector<LinkageCluster* > separated_clusters_to_finally_merge = vector<LinkageCluster* >();
	vector<int> genes_randomized_indecies = vector<int>(genes_range.size());
	for (int i = 0; i < genes_range.size(); i++) {
		genes_randomized_indecies[i] = i;
	}
	random_values_holder.shuffle_vector(genes_randomized_indecies);

	vector<vector<float> > gene_log_prob = vector<vector<float> >();
	vector<vector<int> > index_gene_value_to_index_here = vector<vector<int> > ();
	gene_log_prob.reserve(genes_range.size());
	for (int i = 0; i < genes_range.size(); i++) {
		gene_log_prob.push_back(vector<float>(0));
		index_gene_value_to_index_here.push_back(vector<int>(0));
		vector<int> counts_here(genes_range[i]);

		for (int individual_id = 0; individual_id < genotypes.size(); individual_id++) {
			counts_here[(*(genotypes[individual_id]))[i]]++;
		}

		for (int gene_value = 0; gene_value < genes_range[i]; gene_value++) {
			if (counts_here[gene_value] > 0) {
				gene_log_prob[i].push_back(PZ_Math::log_0_1(gene_log_prob[i][gene_value] / (float)genotypes.size())); // it handles 0 well, it returns 0 at 0
				index_gene_value_to_index_here[i].push_back(gene_log_prob[i].size() - 1);
			}
			else {
				index_gene_value_to_index_here[i].push_back(0);
			}
		}
	}

	for (int index_start = 0; index_start < genes_range.size(); index_start += this->param_separation_size) {
		int index_end = index_start + this->param_separation_size;
		if (index_end > genes_range.size()) {
			index_end = genes_range.size();
		}

		// creating clusters
		vector<LinkageCluster* > clusters_to_process = vector<LinkageCluster* >();
		clusters_to_process.reserve(this->param_separation_size);
		unordered_map<int, float> distance_to_clusters = unordered_map<int, float>();
		distance_to_clusters.reserve(this->param_separation_size * this->param_separation_size);


		// creating initial clusters with only one gene
		for (int i = index_start; i < index_end; i++) {
			clusters_to_process.push_back(new LinkageCluster(vector<int> {genes_randomized_indecies[i]}, new_clusters_counter, nullptr, nullptr));
			if (this->param_min_cluster_size <= 1) {
				usable_clusters_ordered.push_back(clusters_to_process[clusters_to_process.size() - 1]);
			}
			
			new_clusters_counter++;
		}

		// calulating distance beetwen genes pairs
		for (int i = 0; i < clusters_to_process.size(); i++) {
			for (int j = i + 1; j < clusters_to_process.size(); j++) {
				int distance_id = get_distance_id(clusters_to_process[i], clusters_to_process[j]);
				int i_id = clusters_to_process[i]->get_indecies()[0];
				int j_id = clusters_to_process[j]->get_indecies()[0];
				vector<vector<int> > gene_i_j_counter = vector<vector<int> >(gene_log_prob[i_id].size(), vector<int>(gene_log_prob[j_id].size()));

				for (int individual_id = 0; individual_id < genotypes.size(); individual_id++) {
					int gene_i = index_gene_value_to_index_here[individual_id][(*(genotypes[individual_id]))[i_id]];
					int gene_j = index_gene_value_to_index_here[individual_id][(*(genotypes[individual_id]))[j_id]];

					gene_i_j_counter[gene_i][gene_j]++;
					//gene_i_log_prob[gene_i]++;
					//gene_j_log_prob[gene_j]++;
				}


				//for (int gene_1 = 0; gene_1 < genes_range[i_id]; gene_1++) {
				//	if (gene_i_log_prob[gene_1] != 0) {
				//		gene_i_log_prob[gene_1] = log(gene_i_log_prob[gene_1] / genotypes.size());
				//	}
				//	else {
				//		gene_i_log_prob[gene_1] = 0;  // I can do it, cause later it would be 0 anyway
				//	}
				//}

				//for (int gene_2 = 0; gene_2 < genes_range[j_id]; gene_2++) {
				//	if (gene_j_log_prob[gene_2] != 0) {
				//		gene_j_log_prob[gene_2] = log(gene_j_log_prob[gene_2] / genotypes.size());
				//	}
				//	else {
				//		gene_j_log_prob[gene_2] = 0;  // I can do it, cause later it would be 0 anyway
				//	}
				//}


				float nominator = 0;
				float denominator = 0;

				for (int gene_1 = 0; gene_1 < gene_i_j_counter.size(); gene_1++) {
					for (int gene_2 = 0; gene_2 < gene_i_j_counter[0].size(); gene_2++) {
						float prob = gene_i_j_counter[gene_1][gene_2] / ((float)genotypes.size());
						float log_prob = PZ_Math::log_0_1(prob); // it handles 0 well, it returns 0 at 0

						nominator += prob * (gene_log_prob[i_id][gene_1] + gene_log_prob[j_id][gene_2]);
						denominator += prob * log_prob;
					}
				}

				float distance;

				if (denominator == 0) {
					distance = 1;
				}
				else {
					distance = 2 - (nominator / denominator);
				}

				distance_to_clusters[distance_id] = distance;
			}
		}

		// creating part of linkage tree
		while (clusters_to_process.size() > 1) {
			// searching for the closest clusters
			float min_distance = 1000000000;
			int child1_id = 0;
			int child2_id = 0;
			for (int i = 0; i < clusters_to_process.size(); i++) {
				for (int j = i + 1; j < clusters_to_process.size(); j++) {
					int distance_id = get_distance_id(clusters_to_process[i], clusters_to_process[j]);
					float distance = distance_to_clusters[distance_id];
					if (distance < min_distance) {
						min_distance = distance;
						child1_id = i;
						child2_id = j;
					}
				}
			} // child 2 should be bigger than child 1

			// creating new cluster
			LinkageCluster* child1 = clusters_to_process[child1_id];
			LinkageCluster* child2 = clusters_to_process[child2_id];
			vector<int> new_indecies = vector<int>(child1->get_indecies());
			new_indecies.insert(new_indecies.end(), child2->get_indecies().begin(), child2->get_indecies().end());
			LinkageCluster* new_cluster = new LinkageCluster(new_indecies, new_clusters_counter, child1, child2);
			new_clusters_counter++;

			// removing old clusters
			PZ_Math::swap_and_remove(clusters_to_process, child2_id);
			PZ_Math::swap_and_remove(clusters_to_process, child1_id);

			// updating distances
			for (int i = 0; i < clusters_to_process.size(); i++) {
				int distance_id = get_distance_id(clusters_to_process[i], new_cluster);
				int child1_distance_id = get_distance_id(clusters_to_process[i], child1);
				int child2_distance_id = get_distance_id(clusters_to_process[i], child2);
				float distance_child1 = distance_to_clusters[child1_distance_id];
				float distance_child2 = distance_to_clusters[child2_distance_id];
				//float new_distance = (distance_child1 * child1->get_indecies().size() + distance_child2 * child2->get_indecies().size()) / (new_cluster->get_indecies().size());
				float new_distance;
				if (distance_child1 > distance_child2) {  //currently take max distance
					new_distance = distance_child1;
				}
				else {
					new_distance = distance_child2;
				}
				distance_to_clusters[distance_id] = new_distance;
			}

			// adding new cluster
			if (new_cluster->get_indecies().size() >= this->param_min_cluster_size && new_cluster->get_indecies().size() <= this->param_max_cluster_size) {
				usable_clusters_ordered.push_back(new_cluster);
			}
			clusters_to_process.push_back(new_cluster);
			
		}

		separated_clusters_to_finally_merge.push_back(clusters_to_process[0]);
	}

	//randomly creating final linkage tree
	while (separated_clusters_to_finally_merge.size() > 1) {
		int child1_id = 0;
		int child2_id = 1;

		// creating new cluster
		LinkageCluster* child1 = separated_clusters_to_finally_merge[child1_id];
		LinkageCluster* child2 = separated_clusters_to_finally_merge[child2_id];
		vector<int> new_indecies = vector<int>(child1->get_indecies());
		new_indecies.insert(new_indecies.end(), child2->get_indecies().begin(), child2->get_indecies().end());
		LinkageCluster* new_cluster = new LinkageCluster(new_indecies, new_clusters_counter, child1, child2);

		// removing old clusters
		PZ_Math::swap_and_remove(separated_clusters_to_finally_merge, child2_id);
		PZ_Math::swap_and_remove(separated_clusters_to_finally_merge, child1_id);

		// adding new cluster
		separated_clusters_to_finally_merge.push_back(new_cluster);

		if (new_cluster->get_indecies().size() >= this->param_min_cluster_size && new_cluster->get_indecies().size() <= this->param_max_cluster_size) {
			usable_clusters_ordered.push_back(new_cluster);
		}
		new_clusters_counter++;
	}

	// saving final cluster
	cluster_of_whole_genotype = separated_clusters_to_finally_merge[0];
	
	sort(usable_clusters_ordered.begin(), usable_clusters_ordered.end(),
		[](LinkageCluster* a, LinkageCluster* b) {
			return a->get_indecies().size() < b->get_indecies().size();
		}
	);
}


void LinkageTree::build_generic_tree(RandomValuesHolder& random_values_holder) {
	if (cluster_of_whole_genotype != nullptr) {
		delete cluster_of_whole_genotype;
	}

	// initial actions
	int new_clusters_counter = 0;

	usable_clusters_ordered = vector<LinkageCluster* >();
	vector<int> genes_randomized_indecies = vector<int>(genes_range.size());
	for (int i = 0; i < genes_range.size(); i++) {
		genes_randomized_indecies[i] = i;
	}
	random_values_holder.shuffle_vector(genes_randomized_indecies);

	vector<LinkageCluster* > clusters_to_process = vector<LinkageCluster* >();
	clusters_to_process.reserve(this->genes_range.size());

	for (int i = 0; i < this->genes_range.size(); i++) {
		clusters_to_process.push_back(new LinkageCluster(vector<int> {genes_randomized_indecies[i]}, new_clusters_counter, nullptr, nullptr));
		if (this->param_min_cluster_size <= 1) {
			usable_clusters_ordered.push_back(clusters_to_process[clusters_to_process.size() - 1]);
		}

		new_clusters_counter++;
	}

	//randomly creating final linkage tree
	while (clusters_to_process.size() > 1) {
		int child1_id = 0;
		int child2_id = 1;

		// creating new cluster
		LinkageCluster* child1 = clusters_to_process[child1_id];
		LinkageCluster* child2 = clusters_to_process[child2_id];
		vector<int> new_indecies = vector<int>(child1->get_indecies());
		new_indecies.insert(new_indecies.end(), child2->get_indecies().begin(), child2->get_indecies().end());
		LinkageCluster* new_cluster = new LinkageCluster(new_indecies, new_clusters_counter, child1, child2);

		// removing old clusters
		clusters_to_process.erase(clusters_to_process.begin(), clusters_to_process.begin() + 2);

		// adding new cluster
		clusters_to_process.push_back(new_cluster);

		if (new_cluster->get_indecies().size() >= this->param_min_cluster_size && new_cluster->get_indecies().size() <= this->param_max_cluster_size) {
			usable_clusters_ordered.push_back(new_cluster);
		}
		new_clusters_counter++;
	}

	// saving final cluster
	cluster_of_whole_genotype = clusters_to_process[0];

	sort(usable_clusters_ordered.begin(), usable_clusters_ordered.end(),
		[](LinkageCluster* a, LinkageCluster* b) {
			return a->get_indecies().size() < b->get_indecies().size();
		}
	);
}

LinkageTree::~LinkageTree() {
	if (cluster_of_whole_genotype != nullptr) {
		delete cluster_of_whole_genotype;
	}
}

int LinkageTree::get_distance_id(LinkageCluster* cluster1, LinkageCluster* cluster2) {
	int min_id = cluster1->get_id();
	int max_id = cluster2->get_id();

	if (min_id > max_id) {
		int tmp = min_id;
		min_id = max_id;
		max_id = tmp;
	}

	return min_id * 10000 + max_id;
}

vector<LinkageCluster* > & LinkageTree::get_clusters_ordered() {
	return usable_clusters_ordered;
}











// Linkage cluster definition
LinkageCluster::LinkageCluster(vector<int> indecies, int id, LinkageCluster* child1, LinkageCluster* child2) {
	this->linked_indecies = indecies;
	this->id = id;
	this->child1 = child1;
	this->child2 = child2;
}

LinkageCluster::~LinkageCluster() {
	if (child1 != nullptr) {
		delete child1;
	}
	if (child2 != nullptr) {
		delete child2;
	}
}

vector<int> & LinkageCluster::get_indecies() {
	return linked_indecies;
}

int LinkageCluster::get_id() {
	return this->id;
}