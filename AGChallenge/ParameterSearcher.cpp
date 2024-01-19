#include "ParameterSearcher.h"

void ParameterSearcher::run_n_iterations_genetic_algorithm(GeneticAlgorithm * genetic_algorithm, int iterations)
{
	for (int i = 0; i < iterations; i++) {
		genetic_algorithm->run_iteration();
	}
}

void ParameterSearcher::parameters_search(
	int population_size,
	CLFLnetEvaluator* evaluator,
	vector<float> cross_probs,
	vector<float> mut_probs,
	vector<float> choose_better_in_cross_probs)
{
	//// this function searches for the best parameters for genetic algorithm and saves to file all results

	//const int max_iterations = 10000;
	//const int show_result_every_n_iterations = 1000;
	//ThreadPool pool(4);


	//vector<GeneticAlgorithm*> genetic_algorithms;
	//string log_file_name = "C:/Piotr/2023_studia/semestr3/TEP/AG/AGChallengePubl_vs2022_konkurs/logs/parameters_search_log.csv";

	//ofstream log_file;
	//log_file.open(log_file_name);

	//if (!log_file.is_open()) {
	//	std::cerr << "Failed to open file: " << log_file_name << std::endl;
	//	return;
	//}

	//log_file << "iteration;best_fitness;cross_prob;mut_prob;choose_better_in_cross_prob" << endl;
	//log_file.flush();
	//try {
	//	for (int i = 0; i < cross_probs.size(); i++) {
	//		for (int j = 0; j < mut_probs.size(); j++) {
	//			for (int k = 0; k < choose_better_in_cross_probs.size(); k++) {
	//				genetic_algorithms.push_back(
	//					new GeneticAlgorithm(
	//						evaluator,
	//						population_size,
	//						cross_probs[i],
	//						mut_probs[j],
	//						choose_better_in_cross_probs[k]
	//					)
	//				);
	//			}
	//		}
	//	}



	//	for (int i = 0; i < max_iterations; i += show_result_every_n_iterations) {
	//		for (int j = 0; j < genetic_algorithms.size(); j++) {
	//			run_n_iterations_genetic_algorithm(genetic_algorithms[j], show_result_every_n_iterations);
	//			//auto future = pool.enqueue(run_n_iterations_genetic_algorithm, genetic_algorithms[j], show_result_every_n_iterations);
	//			//threads.push_back(thread(run_n_iterations_genetic_algorithm, genetic_algorithms[j], show_result_every_n_iterations));
	//		}
	//		//for (int j = 0; j < threads.size(); j++) {
	//		//	threads[j].join();
	//		//}
	//		//threads.clear();

	//		for (int j = 0; j < genetic_algorithms.size(); j++) {
	//			log_file << i << ";"
	//				<< genetic_algorithms[j]->get_best_individual()->get_fitness() << ";"
	//				<< genetic_algorithms[j]->get_cross_prob() << ";"
	//				<< genetic_algorithms[j]->get_mut_prob() << ";"
	//				<< genetic_algorithms[j]->get_choose_better_in_cross_prob() << endl;
	//			log_file.flush();
	//		}
	//	}
	//	for (int j = 0; j < genetic_algorithms.size(); j++) {
	//		delete genetic_algorithms[j];
	//	}
	//}
	//catch (exception e) {
	//	log_file << "exception: " << e.what() << endl;
	//	log_file.flush();
	//}
	//log_file.close();
}




void ParameterSearcher::save_stats_every_n_seconds(vector<AbstractEvolutionaryAlgorithm* > evolutionary_algorithms, int seconds, int max_seconds)
{
	// this function searches for the best parameters for genetic algorithm and saves to file all results

	//ThreadPool pool(4);
	//string log_file_name = "C:/Piotr/2023_studia/semestr3/TEP/AG/AGChallengePubl_vs2022_konkurs/logs/parameters_search_log.csv";

	//ofstream log_file;
	//log_file.open(log_file_name);

	//if (!log_file.is_open()) {
	//	std::cerr << "Failed to open file: " << log_file_name << std::endl;
	//	return;
	//}

	//log_file << "seconds_passed;algorithm_name;best_fitness" << endl;
	//log_file.flush();


	//try {
	//	for (int i = 0; i < max_iterations; i += show_result_every_n_iterations) {
	//		for (int j = 0; j < genetic_algorithms.size(); j++) {
	//			run_n_iterations_genetic_algorithm(genetic_algorithms[j], show_result_every_n_iterations);
	//			auto future = pool.enqueue(run_n_iterations_genetic_algorithm, genetic_algorithms[j], show_result_every_n_iterations);
	//			threads.push_back(thread(run_n_iterations_genetic_algorithm, genetic_algorithms[j], show_result_every_n_iterations));
	//		}
	//		for (int j = 0; j < threads.size(); j++) {
	//			threads[j].join();
	//		}
	//		threads.clear();

	//		for (int j = 0; j < genetic_algorithms.size(); j++) {
	//			log_file << i << ";"
	//				<< genetic_algorithms[j]->get_best_individual()->get_fitness() << ";"
	//				<< genetic_algorithms[j]->get_cross_prob() << ";"
	//				<< genetic_algorithms[j]->get_mut_prob() << ";"
	//				<< genetic_algorithms[j]->get_choose_better_in_cross_prob() << endl;
	//			log_file.flush();
	//		}
	//	}
	//	for (int j = 0; j < genetic_algorithms.size(); j++) {
	//		delete genetic_algorithms[j];
	//	}
	//}
	//catch (exception e) {
	//	log_file << "exception: " << e.what() << endl;
	//	log_file.flush();
	//}
	//log_file.close();
}