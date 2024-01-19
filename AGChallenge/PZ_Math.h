#pragma once
#include <math.h>
#include <vector>
#include <string>
#include <fstream>

#define LOG_0_1_RESOLUTION 10000

using namespace std;
class PZ_Math
{
private:
	static float * log_0_1_table;

public:
	static void PZ_Math_Initialize();
	static void PZ_Math_Destroy();

	static float prob_mul_log(float prob);
	static float log_0_1(float x); // return log between 0 and 1, if > 1 then 0, if <= 0 then 0
	
	template <typename T>
	static void swap_and_remove(vector<T> & vector_given, int index_to_remove) {
		T last_element = vector_given.back();
		vector_given[index_to_remove] = last_element;
		vector_given.pop_back();
	}

	template <typename T>
	static bool save_vector(const vector<T> & vector_given, const string filename) {
		ofstream file(filename);

		if (file.is_open()) {
			// Iterate over the vector and write each element to the file
			for (size_t i = 0; i < vector_given.size(); ++i) {
				file << vector_given[i];
				if (i < vector_given.size() - 1) {
					// Add a comma after all elements except the last one
					file << ",";
				}
			}

			// Close the file
			file.close();
			//std::cout << "Data written to output.csv" << std::endl;
			return true;
		}
		else {
			//std::cerr << "Unable to open file" << std::endl;
			return false;
		}
	}
};

