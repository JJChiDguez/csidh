#include "edwards_curve.h"

void random_key(uint8_t key[])
{
	uint8_t i, tmp, r;
	int8_t exp, sgn;
	for(i = 0; i < N; i++)
	{

		r = B[i] & 0x1;						// B_i mod 2

		// exp is randomly selected from |[ 0, B ]|
		randombytes(&tmp, 1);
		while ( issmaller((int32_t)B[i], (int32_t)tmp) == -1 )	// constant-time comparison
			randombytes(&tmp, 1);

		exp = (int8_t)tmp;

		// Mapping integers from |[ 0, B |] into
		//                                      |[ -B/2, B/2]| if B is even, or
		//                                      |[ -(B+1)/2, (B-1)/2 ]| if B is odd.
		exp = ( (exp << 1) - (B[i] + r) ) >> 1;

		// Mapping into the set |[-B, B]|.
		exp = (exp << 1) + r;
		sgn = exp >> 7;	// sign of exp

		// Next, to write  key[i] = e || ((1 + sgn)/2)
		cmov(&exp, -exp, sgn == -1);
		key[i] = (exp << 1) ^ (1 & (1 + sgn));
	};
};

void printf_key(uint8_t key[], char *c)
{
	int i;
	printf("%s := ", c);
	printf("{\t  %3d", (int)( (2*(key[0] & 0x1) - 1) * (key[0] >> 1) ));

	for(i = 1; i < N; i++)
	{
		printf(", %3d", (int)( (2*(key[i] & 0x1) - 1) * (key[i] >> 1) ) );
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

	proj current_A, current_T[2];
	point_copy(current_A, A);			// initial Edwards curve constants a and (a -d)
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Vairables required for running SIMBA
	int8_t ec = 0;
	uint16_t count = 0;
	proj G[2], K[(LARGE_L >> 1) + 1];		// Current kernel
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
		elligator(current_T[1], current_T[0], current_A);

		// Next, it is required to multiply the point by 4 and each l_i that doesn't belong to the current batch
		// T_{-}
		yDBL(current_T[0], current_T[0], current_A); // mult. by [2]
		yDBL(current_T[0], current_T[0], current_A); // mult. by [2]
		// T_{+}
		yDBL(current_T[1], current_T[1], current_A); // mult. by [2]
		yDBL(current_T[1], current_T[1], current_A); // mult. by [2]
		// Now, it is required to multiply by the complement of the batch
		for(i = 0; i < size_of_each_complement_batch[m]; i++)
		{
			yMUL(current_T[0], current_T[0], current_A, complement_of_each_batch[m][i]);
			yMUL(current_T[1], current_T[1], current_A, complement_of_each_batch[m][i]);
		};

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
				point_copy(G[0], current_T[0]);
				point_copy(G[1], current_T[1]);
                
				ec = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time
				fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                
				for (j = (i + 1); j < size_of_each_batch[m]; j++)
				{
					if( finished[batches[m][j]] == 0 )
					{
						//depends only on randomness
						yMUL(G[0], G[0], current_A, batches[m][j]);
					};
				};

				if ( (isinfinity(G[0]) != 1) && (isinfinity(G[1]) != 1) )	// Depending on randomness
				{
					bc = isequal(ec >> 1, 0) & 1;		// Bit that determine the current isogeny. This ask is done in constant-time

					yISOG(K, current_A, G[0], current_A, batches[m][i]);
					
					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[1]

						yMUL(current_T[1], current_T[1], current_A, batches[m][i]);	// [l]T[1]
					};

					tmp_e[batches[m][i]] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
				}
				else
				{
					// We must perform at most two scalar multiplications by l.
					yMUL(current_T[1], current_T[1], current_A, batches[m][i]);
				};

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                
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
	point_copy(C, current_A);
};
