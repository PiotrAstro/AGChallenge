#include "PZ_Math.h"

float* PZ_Math::log_0_1_table = nullptr;


void PZ_Math::PZ_Math_Initialize()
{
	log_0_1_table = new float[LOG_0_1_RESOLUTION + 1];
	log_0_1_table[0] = log(1 / (float)(10*LOG_0_1_RESOLUTION));
	for (int i = 1; i <= LOG_0_1_RESOLUTION; i++)
	{
		log_0_1_table[i] = log(i / (float)LOG_0_1_RESOLUTION);
	}
}

void PZ_Math::PZ_Math_Destroy()
{
	delete[] log_0_1_table;
}

float PZ_Math::prob_mul_log(float prob) {
	if (prob <= 0 || prob >= 1) {
		return 0;
	}
	else {
		return prob * log(prob);
	}
}


float PZ_Math::log_0_1(float x)
{
	if (x <= 0 || x >= 1)
	{
		return 0;
	}
	else
	{
		int lower = (int) (x * LOG_0_1_RESOLUTION);
		float prob = x * LOG_0_1_RESOLUTION - lower;
		return log_0_1_table[lower] * (1 - prob) + log_0_1_table[lower + 1] * prob;
		//return log_0_1_table[(int) (x * LOG_0_1_RESOLUTION)];
	}
}