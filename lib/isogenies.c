#include "edwards_curve.h"

/* ----------------------------------------------------------------------------- *
   yISOG()
   Inputs: the projective Edwards y-coordinate of y(P)=YP/ZP, the Edwards curve 
           constant A[0]:=a and A[1]:=(a - d), and integer number 0 <= i < N;
   Output: degree-l isogenous Edwards curve constants C[0]:=a and C[1]:=(a - d) 
           determined by y(P), and the y-coordinate projective of y(P), y([2]P), 
           ... y([(l-1)/2]P)
 * ----------------------------------------------------------------------------- */
void yISOG(proj Pk[], proj C, const proj P, const proj A, const uint8_t i)
{
	uint8_t mask;
	int64_t bits_l;
	uint64_t j,
                 l = L[i],				// l_i
	         s = l >> 1;				// s <- (l_i - 1) / 2

	bits_l = 0;
	while (l > 0)
	{
		l >>= 1;
		bits_l += 1;
	};	// number of bits of l
	l = L[i];
	// ---
	fp By[2], Bz[2], tmp_0, tmp_1, tmp_d;

	copy(tmp_0, A[0], NUMBER_OF_WORDS);		// a
	fp_sub(tmp_d, A[0], A[1]);			// d
	copy(tmp_1, tmp_d, NUMBER_OF_WORDS);

	copy(By[0], P[0], NUMBER_OF_WORDS); copy(By[1], P[0], NUMBER_OF_WORDS);
	copy(Bz[0], P[1], NUMBER_OF_WORDS); copy(Bz[1], P[1], NUMBER_OF_WORDS);

	point_copy(Pk[0], P);				// P
	yDBL(Pk[1], P, A);				// [2]P
		
	for(j = 2; j < s; j++)
	{
		fp_mul(By[0], By[0], Pk[j - 1][0]);
		fp_mul(Bz[0], Bz[0], Pk[j - 1][1]);
		yADD(Pk[j], Pk[j - 1], P, Pk[j - 2]);	// [j + 1]P

		FP_MUL_COMPUTED += 2;
	};

	mask = isequal(l, 3) ^ 1;		// If l = 3 then we keep with the current values of By[0] and Bz[0]. This ask is done in constant-time
	fp_mul(By[1], By[0], Pk[s - 1][0]);	// This an extra cost for a degree-3 construction
	fp_mul(Bz[1], Bz[0], Pk[s - 1][1]);	// This an extra cost for a degree-3 construction
	fp_cswap(By[0], By[1], mask);		// constant-time swap: dummy or not dummy, that is the question.
	fp_cswap(Bz[0], Bz[1], mask);		// constant-time swap: dummy or not dummy, that is the question.

	// left-to-right method for computing a^l and d^l
	bits_l -= 1;
	for(j = 1; j <= bits_l; j++)
	{
		fp_sqr(tmp_0, tmp_0);
		fp_sqr(tmp_1, tmp_1);
		//if( ( (0x1 << (bits_l - j) ) & l ) != 0x0)
		if( ( (l >> (bits_l - j)) & 1 ) != 0)
		{
			fp_mul(tmp_0, tmp_0, A[0]);
			fp_mul(tmp_1, tmp_1, tmp_d);
			
			FP_MUL_COMPUTED += 2;
		};

		FP_SQR_COMPUTED += 2;
	};

	for(j = 0; j < 3; j++)
	{
		fp_sqr(By[0], By[0]);
		fp_sqr(Bz[0], Bz[0]);

		FP_SQR_COMPUTED += 2;
	};

	fp_mul(C[0], tmp_0, Bz[0]);
	fp_mul(C[1], tmp_1, By[0]);
	fp_sub(C[1], C[0], C[1]);

	FP_ADD_COMPUTED += 2;
	FP_MUL_COMPUTED += 4;
};// Cost ~ (3l + log(l) - 7)M + (l + 2log(l) + 3)S + (3l - 1)a

/* ----------------------------------------------------------------------------- *
   yEVAL()
   Inputs: the projective Edwards y-coordinate of y(Q)=YQ/ZQ, the y-coordinate 
           projective of y(P), y([2]P), ..., y([(l-1)/2]P), and integer number 
           0 <= i < N;
   Output: the image of y(Q) under a degree-L[i] isogeny with kernel generated 
           by y(P).
 * ----------------------------------------------------------------------------- */
void yEVAL(proj R, const proj Q, const proj Pk[], const uint8_t i)
{
	int j;
	fp tmp_0, tmp_1, s_0, s_1;

	proj tmp_Q;
	point_copy(tmp_Q, Q);	// This is for allowing Q <- image of Q

	// Evaluating Q
	fp_mul(s_0, tmp_Q[0], Pk[0][1]);
	fp_mul(s_1, tmp_Q[1], Pk[0][0]);
	// Mapping R into the isomorphic Montgomery curve
	fp_add(R[0], s_0, s_1);
	fp_sub(R[1], s_0, s_1);

	uint64_t s = (L[i] >> 1);
	for(j = 1; j < s; j++)
	{
		// Evaluating Q
		fp_mul(s_0, tmp_Q[0], Pk[j][1]);
		fp_mul(s_1, tmp_Q[1], Pk[j][0]);
		fp_add(tmp_0, s_0, s_1);
		fp_sub(tmp_1, s_0, s_1);
		fp_mul(R[0], R[0], tmp_0);
		fp_mul(R[1], R[1], tmp_1);

		FP_ADD_COMPUTED += 2;
		FP_MUL_COMPUTED += 4;
	};

	fp_sqr(R[0], R[0]);
	fp_sqr(R[1], R[1]);
	// Mapping Q into the isomorphic Montgomery curve
	fp_add(tmp_0, tmp_Q[1], tmp_Q[0]);
	fp_sub(tmp_1, tmp_Q[1], tmp_Q[0]);
	fp_mul(tmp_0, R[0], tmp_0);
	fp_mul(tmp_1, R[1], tmp_1);
	// Mapping R into the Edwards curve
	fp_sub(R[0], tmp_0, tmp_1);
	fp_add(R[1], tmp_0, tmp_1);

	FP_ADD_COMPUTED += 6;
	FP_SQR_COMPUTED += 2;
	FP_MUL_COMPUTED += 4;
};// Cost : 2(l - 1)M + 2S + (3 + l)a

