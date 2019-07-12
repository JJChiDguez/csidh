#include <math.h>
#include<time.h>

#include "fp.h"
#include "edwards_curve.h"

// Measuring the perfomance
static uint64_t get_cycles()
{
   uint32_t lo, hi;
   asm volatile("rdtsc":"=a"(lo),"=d"(hi));
   return ((uint64_t)hi<<32) | lo;
};

static uint8_t csidh(proj out, const uint8_t sk[], const proj in)
{
	FP_ADD_COMPUTED = 0;
	FP_SQR_COMPUTED = 0;
	FP_MUL_COMPUTED = 0;

	if (!validate(in)) {
		return 0;
	};
	action_evaluation(out, sk, in);
	return 1;
};

unsigned long its = 1024;

int main()
{
	unsigned int i;

	uint64_t add_min = 0xFFFFFFFFFFFFFFFF, add_max = 0, 
	         sqr_min = 0xFFFFFFFFFFFFFFFF, sqr_max = 0, 
	         mul_min = 0xFFFFFFFFFFFFFFFF, mul_max = 0;

	float add_mean = 0, add_variance = 0,
	      sqr_mean = 0, sqr_variance = 0,
	      mul_mean = 0, mul_variance = 0;

	uint64_t add_sample[its],
	         sqr_sample[its],
	         mul_sample[its];

	// ---
	uint8_t key[N];
	proj random_E;
	point_copy(random_E, E);

	for(i = 0; i < its; ++i)
	{

		if (its >= 100 && i % (its / 100) == 0) {
			printf("Doing %lu iterations of action with validation key:\t", its);
            printf("%2lu%%", 100 * i / its);
            fflush(stdout);
            printf("\r\x1b[K");
        }

		random_key(key);
		assert(csidh(random_E, key, random_E));
		
		// ---
		add_sample[i] = FP_ADD_COMPUTED;
		sqr_sample[i] = FP_SQR_COMPUTED;
		mul_sample[i] = FP_MUL_COMPUTED;

		/**************************************/
		if(add_min > add_sample[i])
			add_min = add_sample[i];
		if(sqr_min > sqr_sample[i])
			sqr_min = sqr_sample[i];
		if(mul_min > mul_sample[i])
			mul_min = mul_sample[i];

		/**************************************/
		if(add_max < add_sample[i])
			add_max = add_sample[i];
		if(sqr_max < sqr_sample[i])
			sqr_max = sqr_sample[i];
		if(mul_max < mul_sample[i])
			mul_max = mul_sample[i];
		
		/**************************************/
		add_mean += (float)add_sample[i];
		sqr_mean += (float)sqr_sample[i];
		mul_mean += (float)mul_sample[i];
	};


	add_mean = add_mean / ((float)its * 1.0);
	sqr_mean = sqr_mean / ((float)its * 1.0);
	mul_mean = mul_mean / ((float)its * 1.0);

	for (i = 0; i < its; ++i)
	{
		add_variance += (add_sample[i] - add_mean)*(add_sample[i] - add_mean);
		sqr_variance += (sqr_sample[i] - sqr_mean)*(sqr_sample[i] - sqr_mean);
		mul_variance += (mul_sample[i] - mul_mean)*(mul_sample[i] - mul_mean);
	};

	add_variance = add_variance / ((float)its - 1.0);
	sqr_variance = sqr_variance / ((float)its - 1.0);
	mul_variance = mul_variance / ((float)its - 1.0);
		
	printf("\x1b[01;33mIterations: %lu\x1b[0m\n\n", its);

	printf("\x1b[33mAverage costs:\x1b[0m\n");
	printf("\t %f additions,\n", add_mean);
	printf("\t\x1b[32m %f squarings,\x1b[0m\n", sqr_mean);
	printf("\t\x1b[31m %f multiplications.\x1b[0m\n", mul_mean);

	printf("\n");

	printf("\x1b[33mStandard deviation of the costs:\x1b[0m\n");
	printf("\t %f additions,\n", sqrt(add_variance));
	printf("\t\x1b[32m %f squarings,\x1b[0m\n", sqrt(sqr_variance));
	printf("\t\x1b[31m %f multiplications.\x1b[0m\n", sqrt(mul_variance));

	printf("\n");

	printf("\x1b[33mMinimum costs:\x1b[0m\n");
	printf("\t %lu additions,\n", add_min);
	printf("\t\x1b[32m %lu squarings,\x1b[0m\n", sqr_min);
	printf("\t\x1b[31m %lu multiplications.\x1b[0m\n", mul_min);

	printf("\n");

	printf("\x1b[33mMaximum costs:\x1b[0m\n");
	printf("\t %lu additions,\n", add_max);
	printf("\t\x1b[32m %lu squarings,\x1b[0m\n", sqr_max);
	printf("\t\x1b[31m %lu multiplications.\x1b[0m\n", mul_max);

	return 0;
};
