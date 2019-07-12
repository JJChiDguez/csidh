#include "edwards_curve.h"

void random_key(uint8_t key[])
{
	uint8_t i, tmp;
	for(i = 0; i < N; i++)
	{	
		// Random number from |[0, B]|
		randombytes(&tmp, 1);
		while ( issmaller((int32_t)B[i], (int32_t)tmp) == -1 )	// constant-time comparison
			randombytes(&tmp, 1);
		key[i] = tmp;
	
	};
};

void printf_key(uint8_t key[], char *c)
{
	int i;
	printf("%s := ", c);
	printf("{\t  %2d", key[0]);

	for(i = 1; i < N; i++)
	{
		printf(", %2d", key[i]);
		if( (i % 18) == 17 )
			printf("\n\t\t");
	};

	printf("};\n");

};

/* ----------------------------------------------------------------------------------------------- *
   action_evaluation()
   inputs: a the secret key, the Edwards curve constants A[0]:=a, and A[1]:=(a - d);
   output: the isogenous Edwards curve constants C[0]:=a' and C[1]:=(a' - d') determined by the action 
           evaluated at the secret key and public curve A
   
    NOTE: As far as we've understood how simba works; this next code implements simba approach.
          The action computed by the next code uses only one torsion point T_{+}, which its affine 
          y-coordinate (in the isomorphic Montgomery curve) belongs to Fp.
 * ----------------------------------------------------------------------------------------------- */
void action_evaluation(proj C, const uint8_t key[], const proj A)
{
	// --------------------------------------------------------------------------------------------------------
	// SIMBA parameters
	// Batches
	uint8_t batches[NUMBER_OF_BATCHES][SIZE_OF_EACH_BATCH[0]];
	uint8_t size_of_each_batch[NUMBER_OF_BATCHES];

	for(uint8_t i = 0; i < NUMBER_OF_BATCHES; i++)
		memcpy(batches[i], BATCHES[i], sizeof(uint8_t) * SIZE_OF_EACH_BATCH[i]);

	memcpy(size_of_each_batch, SIZE_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// Complement of each batch
	uint8_t complement_of_each_batch[NUMBER_OF_BATCHES][N];
	uint8_t size_of_each_complement_batch[NUMBER_OF_BATCHES];

	memcpy(complement_of_each_batch, COMPLEMENT_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES * N);
	memcpy(size_of_each_complement_batch, SIZE_OF_EACH_COMPLEMENT_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Copy of public and private data (the private key is modified each iteration)
	uint8_t tmp_e[N];
	memcpy(tmp_e, key, sizeof(uint8_t) * N);	// exponents

	proj current_A[2], current_Tp[2];
	point_copy(current_A[0], A);			// initial Edwards curve constants a and (a -d)
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Vairables required for running SIMBA
	int8_t ec = 0, mask;
	uint16_t count = 0;
	proj G[2], K[(LARGE_L >> 1) + 1], Z;		// Current kernel
	uint8_t finished[N];				// flag that determines if the maximum number of isogeny constructions has been reached
	memset(finished, 0, sizeof(uint8_t) * N);

	int8_t counter[N];				// This variable determines how many isogeny construcctions has been perfomed
	memset(counter, 0, sizeof(int8_t) * N);
	memcpy(counter, B, sizeof(int8_t) * N);		// At the beginning, we must perfomed b_i isogeny constructions for each l_i
	uint64_t isog_counter = 0;			// Total number of isogeny construction perfomed

	uint8_t last_isogeny[NUMBER_OF_BATCHES];
	//index for skipping point evaluations (the last one of each batch)
	memcpy(last_isogeny, LAST_ISOGENY, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	uint32_t bc;
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Main loop
	uint8_t m = 0, i, j;
	uint64_t number_of_batches = NUMBER_OF_BATCHES;
	uint32_t si;

	while (isog_counter < NUMBER_OF_ISOGENIES)
	{
		m = (m + 1) % number_of_batches;
		
		if(count == MY*number_of_batches) {  	//merge the batches after my rounds
			m = 0;
			size_of_each_complement_batch[m] = 0;
			size_of_each_batch[m] = 0;
			number_of_batches = 1;

			for(i = 0; i < N; i++) {
				if( counter[i] == 0 )
				{
					// l_i's reached
					complement_of_each_batch[m][size_of_each_complement_batch[m]] = i;
					size_of_each_complement_batch[m] += 1;
				}
				else
				{
					last_isogeny[0] = i;
					// l_i's not reached
					batches[m][size_of_each_batch[m]] = i;
					size_of_each_batch[m] += 1;
				};
			}
		}

		// Before constructing isogenies, we must to search for suitable point
		elligator(current_Tp[0], current_Tp[1], current_A[0]);

		// Next, it is required to multiply the point by 4 and each l_i that doesn't belong to the current batch
		yDBL(current_Tp[0], current_Tp[0], current_A[0]); // mult. by [2]
		yDBL(current_Tp[0], current_Tp[0], current_A[0]); // mult. by [2]
		// Now, it is required to multiply by the complement of the batch
		for(i = 0; i < size_of_each_complement_batch[m]; i++)
			yMUL(current_Tp[0], current_Tp[0], current_A[0], complement_of_each_batch[m][i]);

		for(i = 0; i < size_of_each_batch[m]; i++)
		{
			if( finished[batches[m][i]] == 1 )
			{ 
				//depends only on randomness
				continue;
			}
			else
			{
				// Now, a degree-(l_{batches[m][i]}) will be constructed
				point_copy(G[0], current_Tp[0]);
				for (j = (i + 1); j < size_of_each_batch[m]; j++)
				{
					if( finished[batches[m][j]] == 0 )
					{
						//depends only on randomness
						yMUL(G[0], G[0], current_A[0], batches[m][j]);
					};
				};

				if ( isinfinity(G[0]) != 1 )
				{
					point_copy(G[1], current_Tp[0]);

					ec = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time
					bc = isequal(ec, 0) & 1;		// Bit that determines if a dummy operation will be perfomed

					fp_cswap(G[0][0], G[1][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
					fp_cswap(G[0][1], G[1][1], bc);		// constant-time swap: dummy or not dummy, that is the question.

					yISOG(K, current_A[1], G[0], current_A[0], batches[m][i]);
					
					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						mask = isequal(L[batches[m][i]], 3);	// Just for catching the case l = 3. This ask is done in constant-time
						si = (L[batches[m][i]] >> 1);		// (l - 1) / 2

						yEVAL(current_Tp[1], current_Tp[0], K, batches[m][i]);

						yADD(Z, K[(si + mask) - 1], G[0], K[(si + mask) - 2]);	// [(l + 1)/2]G[0]
						fp_cswap(Z[0], K[si][0], mask ^ 1);			// constant-time swap: catching degree-3 isogeny case
						fp_cswap(Z[1], K[si][1], mask ^ 1);			// constant-time swap: catching degree-3 isogeny case
						yADD(current_Tp[0], K[si], K[si - 1], G[0]);		// [l]G[0]

						fp_cswap(current_Tp[0][0], current_Tp[1][0], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.
						fp_cswap(current_Tp[0][1], current_Tp[1][1], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.
					};

					fp_cswap(current_A[0][0], current_A[1][0], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.
					fp_cswap(current_A[0][1], current_A[1][1], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.

					tmp_e[batches[m][i]] = ec - (bc ^ 1);
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
				}

				if( counter[batches[m][i]] == 0 )
				{	
					//depends only on randomness
					finished[batches[m][i]] = 1;
					complement_of_each_batch[m][size_of_each_complement_batch[m]] = batches[m][i];
					size_of_each_complement_batch[m] += 1;
				};
			};
		};
		count += 1;
	};
	
	// --------------------------------------------------------------------------------------------------------	
	point_copy(C, current_A[0]);
};
