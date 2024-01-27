#include "LinkageTree.h"


LinkageTree::LinkageTree(vector<int>& genes_ranges) {
	// initial creations
	this->genes_range = vector<int>(genes_ranges);
	
	cluster_of_whole_genotype = nullptr;
	this->param_separation_size = genes_ranges.size();
	this->param_min_cluster_size = 2;
	this->param_max_cluster_size = genes_ranges.size() - 1;
}

LinkageTree::LinkageTree(vector<int>& genes_ranges, int separation_size, int min_cluster_size, int max_cluster_size) {
	// initial creations
	this->genes_range = vector<int>(genes_ranges);

	cluster_of_whole_genotype = nullptr;
	this->param_separation_size = separation_size;
	this->param_min_cluster_size = min_cluster_size;
	this->param_max_cluster_size = max_cluster_size;
}

void LinkageTree::clean_self() {
	usable_clusters_ordered.clear();
	if (cluster_of_whole_genotype != nullptr) {
		queue< pair<LinkageCluster*, LinkageCluster* > > clusters_to_consider;
		clusters_to_consider.push(
			make_pair(
				this->cluster_of_whole_genotype->get_child1(),
				this->cluster_of_whole_genotype->get_child2()
			)
		);
		delete cluster_of_whole_genotype;

		while (clusters_to_consider.size() > 0) {
			LinkageCluster* clusters_1 = clusters_to_consider.front().first;
			LinkageCluster* clusters_2 = clusters_to_consider.front().second;
			clusters_to_consider.pop();

			// adding children to the queue
			if (clusters_1 != nullptr) {
				clusters_to_consider.push(
					make_pair(
						clusters_1->get_child1(),
						clusters_1->get_child2()
					)
				);
				delete clusters_1;
			}

			if (clusters_2 != nullptr) {
				clusters_to_consider.push(
					make_pair(
						clusters_2->get_child1(),
						clusters_2->get_child2()
					)
				);
				delete clusters_2;
			}
		}
	}
}

void LinkageTree::build_tree(const vector<vector<int> * > & genotypes, RandomValuesHolder& random_values_holder)  {
	clean_self();

	// initial actions
	int new_clusters_counter = 0;

	usable_clusters_ordered = vector<LinkageCluster* >(0);
	vector<LinkageCluster* > separated_clusters_to_finally_merge = vector<LinkageCluster* >(0);
	vector<int> genes_randomized_indecies = vector<int>(genes_range.size());
	for (int i = 0; i < genes_range.size(); i++) {
		genes_randomized_indecies[i] = i;
	}
	random_values_holder.shuffle_vector(genes_randomized_indecies);

	vector<vector<float> > gene_log_prob = vector<vector<float> >(0);
	vector<float> gene_entropy = vector<float>(0);
	vector<vector<int> > index_gene_value_to_index_here = vector<vector<int> > (0);
	gene_entropy.reserve(genes_range.size());
	gene_log_prob.reserve(genes_range.size());
	index_gene_value_to_index_here.reserve(genes_range.size());
	for (int i = 0; i < genes_range.size(); i++) {
		gene_log_prob.push_back(vector<float>(0));
		gene_entropy.push_back(0);
		index_gene_value_to_index_here.push_back(vector<int>(0));
		vector<int> counts_here(genes_range[i]);
		double entropy_here = 0;

		for (int individual_id = 0; individual_id < genotypes.size(); individual_id++) {
			counts_here[(*(genotypes[individual_id]))[i]]++;
		}

		for (int gene_value = 0; gene_value < genes_range[i]; gene_value++) {
			if (counts_here[gene_value] > 0) {
				float prob = counts_here[gene_value] / ((float)genotypes.size());
				float log_prob = prob * PZ_Math::log_0_1(prob); // it handles 0 well, it returns 0 at 0
				entropy_here -= log_prob;
				gene_log_prob[i].push_back(log_prob); // it handles 0 well, it returns 0 at 0

				index_gene_value_to_index_here[i].push_back(gene_log_prob[i].size() - 1);
			}
			else {
				index_gene_value_to_index_here[i].push_back(-1);
			}
		}

		gene_entropy[i] = entropy_here;
	}

	for (int index_start = 0; index_start < genes_range.size(); index_start += this->param_separation_size) {
		int index_end = index_start + this->param_separation_size;
		if (index_end > genes_range.size()) {
			index_end = genes_range.size();
		}

		// creating clusters
		vector<LinkageCluster* > clusters_to_process = vector<LinkageCluster* > (0);
		clusters_to_process.reserve(this->param_separation_size);
		unordered_map<int, float> distance_to_clusters = unordered_map<int, float>(0);
		distance_to_clusters.reserve(this->param_separation_size * this->param_separation_size);


		// creating initial clusters with only one gene
		for (int i = index_start; i < index_end; i++) {
			clusters_to_process.push_back(new LinkageCluster(vector<int> {genes_randomized_indecies[i]}, new_clusters_counter, nullptr, nullptr));
			clusters_to_process[i]->set_tmp_index_in_vector(i);
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

				//if (i_id == 2188 && j_id == 2342) {
				//	cout << "here" << endl;
				//	int i = 98;
				//	if (i == 98) {
				//		cout << "fefe";
				//	}
				//}

				vector<vector<int> > gene_i_j_counter = vector<vector<int> >(gene_log_prob[i_id].size(), vector<int>(gene_log_prob[j_id].size()));

				for (int individual_id = 0; individual_id < genotypes.size(); individual_id++) {
					int gene_i = index_gene_value_to_index_here[i_id][(*(genotypes[individual_id]))[i_id]];
					int gene_j = index_gene_value_to_index_here[j_id][(*(genotypes[individual_id]))[j_id]];

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


				//float nominator = 0;
				//float denominator = 0;
				double entropy_i_j = 0;

				for (int gene_1 = 0; gene_1 < gene_i_j_counter.size(); gene_1++) {
					for (int gene_2 = 0; gene_2 < gene_i_j_counter[0].size(); gene_2++) {
						float prob = gene_i_j_counter[gene_1][gene_2] / ((float)genotypes.size());
						float log_prob = PZ_Math::log_0_1(prob); // it handles 0 well, it returns 0 at 0

						entropy_i_j -= prob * log_prob;
						//nominator += prob * (gene_log_prob[i_id][gene_1] + gene_log_prob[j_id][gene_2]);
						//denominator += prob * log_prob;
					}
				}

				float distance;

				if (entropy_i_j == 0) {
					distance = 1;
				}
				else {
					distance = 2 - ((gene_entropy[i_id] + gene_entropy[j_id]) / entropy_i_j);
					if (distance < 0) {
						distance = 0;
					}
					else if (distance > 1) {
						distance = 1;
					}
				}

				distance_to_clusters[distance_id] = distance;
				clusters_to_process[i]->add_tmp_distance(distance, clusters_to_process[j]);
				clusters_to_process[j]->add_tmp_distance(distance, clusters_to_process[i]);
			}
		}

		float min_distance = 1000000000;
		int child1_id = 0;
		int child2_id = 0;
		for (int i = 0; i < clusters_to_process.size(); i++) {
			if (clusters_to_process[i]->get_tmp_min_distance() < min_distance) {
				min_distance = clusters_to_process[i]->get_tmp_min_distance();
				child1_id = i;
				child2_id = clusters_to_process[i]->get_tmp_min_distance_to()->get_tmp_index_in_vector();
			}
		}

		// creating part of linkage tree
		while (clusters_to_process.size() > 1) {
			// searching for the closest clusters
			 // child 2 should be bigger than child 1
			if (child1_id > child2_id) {
				int tmp = child1_id;
				child1_id = child2_id;
				child2_id = tmp;
			}

			// creating new cluster
			LinkageCluster* child1 = clusters_to_process[child1_id];
			LinkageCluster* child2 = clusters_to_process[child2_id];
			vector<int> new_indecies = vector<int>(child1->get_indecies());
			new_indecies.insert(new_indecies.end(), child2->get_indecies().begin(), child2->get_indecies().end());
			LinkageCluster* new_cluster = new LinkageCluster(new_indecies, new_clusters_counter, child1, child2);
			new_clusters_counter++;

			// removing old clusters
			PZ_Math::swap_and_remove(clusters_to_process, child2_id);
			clusters_to_process[child2_id]->set_tmp_index_in_vector(child2_id);
			child1->clear_distance_data();

			PZ_Math::swap_and_remove(clusters_to_process, child1_id);
			clusters_to_process[child1_id]->set_tmp_index_in_vector(child1_id);
			child2->clear_distance_data();

			new_cluster->set_tmp_index_in_vector(clusters_to_process.size());

			min_distance = 1000000000;
			child1_id = 0;
			child2_id = 0;

			// updating distances
			for (int i = 0; i < clusters_to_process.size(); i++) {
				int distance_id = get_distance_id(clusters_to_process[i], new_cluster);
				int child1_distance_id = get_distance_id(clusters_to_process[i], child1);
				int child2_distance_id = get_distance_id(clusters_to_process[i], child2);
				float distance_child1 = distance_to_clusters[child1_distance_id];
				float distance_child2 = distance_to_clusters[child2_distance_id];
				float new_distance = (distance_child1 * child1->get_indecies().size() + distance_child2 * child2->get_indecies().size()) / (new_cluster->get_indecies().size());
				//float new_distance;
				//if (distance_child1 > distance_child2) {  //currently take max distance
				//	new_distance = distance_child1;
				//}
				//else {
				//	new_distance = distance_child2;
				//}
				distance_to_clusters[distance_id] = new_distance;

				//if (distance_child1 == distance_child2) {
				//	cout << "error" << endl;
				//}
				//else
				//{
				//	cout << "ok" << endl;
				//	cout << "new distance: " << new_distance << "\t\tdistance_to_child1: " << distance_child1 << "\t\tdistance_to_child2: " << distance_child2
				//		<< endl << "distance id 1: " << child1_distance_id << "\t\tdistance id 2: " << child2_distance_id << endl;
				//}
				
				clusters_to_process[i]->add_tmp_distance(new_distance, new_cluster);
				new_cluster->add_tmp_distance(new_distance, clusters_to_process[i]);


				if (clusters_to_process[i]->get_tmp_min_distance() < min_distance) {
					min_distance = clusters_to_process[i]->get_tmp_min_distance();
					child1_id = i;
					child2_id = clusters_to_process[i]->get_tmp_min_distance_to()->get_tmp_index_in_vector();
				}
			}

			// adding new cluster
			if (new_cluster->get_indecies().size() >= this->param_min_cluster_size && new_cluster->get_indecies().size() <= this->param_max_cluster_size) {
				usable_clusters_ordered.push_back(new_cluster);
			}
			clusters_to_process.push_back(new_cluster);
			//cout << "clusters left: " << clusters_to_process.size() << endl;
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
	clean_self();

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
	clean_self();
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

LinkageCluster* LinkageTree::get_cluster_of_whole_genotype() {
	return cluster_of_whole_genotype;
}









// Linkage cluster definition
LinkageCluster::LinkageCluster(vector<int> indecies, int id, LinkageCluster* child1, LinkageCluster* child2) {
	this->linked_indecies = indecies;
	this->id = id;
	this->child1 = child1;
	this->child2 = child2;

	this->tmp_index_in_vector = -1;
	tmp_self_in_priority_queue = this;
	min_distance_to_clusters = new priority_queue<pair<float, LinkageCluster**>, vector<pair<float, LinkageCluster**> >, greater<> >();
}

LinkageCluster::~LinkageCluster() {
	if (this->min_distance_to_clusters != nullptr) {
		delete this->min_distance_to_clusters;
	}
}

vector<int> & LinkageCluster::get_indecies() {
	return linked_indecies;
}

int LinkageCluster::get_id() {
	return this->id;
}

void LinkageCluster::add_tmp_distance(float tmp_min_distance, LinkageCluster* tmp_min_distance_to) {
	heal_tmp_min_distance();
	this->min_distance_to_clusters->push(make_pair(tmp_min_distance, &(tmp_min_distance_to->tmp_self_in_priority_queue)));
}

void LinkageCluster::set_tmp_index_in_vector(int tmp_index_in_vector) {
	this->tmp_index_in_vector = tmp_index_in_vector;
}

float LinkageCluster::get_tmp_min_distance() {
	return this->min_distance_to_clusters->top().first;
}

void LinkageCluster::heal_tmp_min_distance() {
	while (!this->min_distance_to_clusters->empty() && *(this->min_distance_to_clusters->top().second) == nullptr) {
		min_distance_to_clusters->pop();
	}
}

LinkageCluster* LinkageCluster::get_tmp_min_distance_to() {
	return *(this->min_distance_to_clusters->top().second);
}

int LinkageCluster::get_tmp_index_in_vector() {
	return this->tmp_index_in_vector;
}

void LinkageCluster::clear_distance_data() {
	delete this->min_distance_to_clusters;
	this->min_distance_to_clusters = nullptr;
	this->tmp_index_in_vector = -1;
	this->tmp_self_in_priority_queue = nullptr;
}

LinkageCluster* LinkageCluster::get_child1() {
	return child1;
}

LinkageCluster* LinkageCluster::get_child2() {
	return child2;
}