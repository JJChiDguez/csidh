#include "edwards_curve.h"

/* ------------------------------------------------------------- *
   isinfinity()
   inputs: the projective Edwards y-coordinates of y(P)=YP/ZP;
   output:
            1 if YP == ZP, or
            0 if YP != ZP
 * ------------------------------------------------------------- */
int isinfinity(const proj P)
{
	fp tmp;
	fp_sub(tmp, P[0], P[1]);				// A substraction in order to ask in constant-time: this can be improvement
	return iszero(tmp, NUMBER_OF_WORDS);	// constant-time comparison
};

/* ------------------------------------------------------------- *
   areEqual()
   inputs: the projective Edwards y-coordinates of y(P)=YP/ZP and
           y(Q)=YQ/ZQ;
   output:
            1 if YP*ZQ == ZP*YQ, or
            0 if YP*ZQ != ZP*YQ
 * ------------------------------------------------------------- */
uint8_t areEqual(const proj P, const proj Q)
{
	fp YPZQ, ZPYQ;
	fp_mul(YPZQ, P[0], Q[1]);
	fp_mul(ZPYQ, P[1], Q[0]);
	return (0 == compare(YPZQ, ZPYQ, NUMBER_OF_WORDS));
};
/* ------------------------------------------------------------- *
   point_copy()
   inputs: a projective Edwards y-coordinates of y(P)=YP/ZP;
   output: a copy of the projective Edwards's y-coordinate of y(P)
 * ------------------------------------------------------------- */
void point_copy(proj Q, const proj P)
{
	copy(Q[0], P[0], NUMBER_OF_WORDS);
	copy(Q[1], P[1], NUMBER_OF_WORDS);
};

/* ---------------------------------------------------------------------- *
   This is NOT the stereotypical Edwards y-only doubling, since it make 
   use of the constant ad := a-d

   yDBL()
   inputs: the projective Edwards y-coordinates of y(P)=YP/ZP, and the
           Edwards curve constant A[0]:=a, and A[1]:=(a - d);
   output: the projective Edwards x-coordinates y([2]P)=Y2P/Z2P
 * ---------------------------------------------------------------------- */
void yDBL(proj Q, const proj P, const proj A)
{
	fp tmp_0, tmp_1;
	// Firstly, yDBL is performed like a xDBL but on the isomorphic Montgomery curve
	fp_sqr(tmp_0, P[0]);
	fp_sqr(tmp_1, P[1]);

	fp_mul(Q[1], A[1], tmp_0);
	fp_mul(Q[0], Q[1], tmp_1);
	fp_sub(tmp_1, tmp_1, tmp_0);
	fp_mul(tmp_0, A[0], tmp_1);
	fp_add(Q[1], Q[1], tmp_0);
	fp_mul(tmp_0, Q[1], tmp_1);

	// Lastly, the result is mapping into the Edward's curve
	fp_add(Q[1], Q[0], tmp_0);
	fp_sub(Q[0], Q[0], tmp_0);

	FP_ADD_COMPUTED += 4;
	FP_SQR_COMPUTED += 2;
	FP_MUL_COMPUTED += 4;
};// Cost : 4M + 2S + 4a

/* ---------------------------------------------------------------------- *
   yADD()
   inputs: the projective Edwards y-coordinates of y(P)=YP/ZP, y(Q)=YQ/ZQ, 
           and y(P-Q)=YPQ/ZPQ;
   output: the projective Edwards y-coordinates of y(P+Q)
 * ---------------------------------------------------------------------- */
void yADD(proj R, const proj P, const proj Q, const proj PQ)
{
	fp tmp_0, tmp_1, xD, zD;
	// Firstly, the difference is mapping into the isomorphic Montgomery curve
	fp_add(xD, PQ[1], PQ[0]);
	fp_sub(zD, PQ[1], PQ[0]);

	// Secondly, yADDL is performed like a xADD (similarly to our yDBL)
	fp_mul(tmp_0, P[1], Q[0]);
	fp_mul(tmp_1, P[0], Q[1]);

	fp_sub(R[1], tmp_0, tmp_1);
	fp_add(R[0], tmp_0, tmp_1);

	fp_sqr(R[1], R[1]);
	fp_sqr(R[0], R[0]);

	fp_mul(tmp_0, R[0], zD);
	fp_mul(tmp_1, R[1], xD);

	// Lastly, the result is mapping into the Edward's curve
	fp_sub(R[0], tmp_0, tmp_1);
	fp_add(R[1], tmp_0, tmp_1);

	FP_ADD_COMPUTED += 6;
	FP_SQR_COMPUTED += 2;
	FP_MUL_COMPUTED += 4;
};// Cost : 4M + 2S + 6a

/* ---------------------------------------------------------------------- *
   yMUL()
   inputs: the projective Edwards y-coordinates of y(P)=YP/ZP, the Edwards 
           curve constant A[0]:=a, and A[1]:=(a - d), and an positive 
           integer integer number 0 <= i <= (N - 1);
   output: the projective Edwards y-coordinates y([l_i]P)
 * ---------------------------------------------------------------------- */
void yMUL(proj Q, const proj P, const proj A, uint8_t const i)
{
	proj R[3], T;

	// Initial 3-tuple of points	
	point_copy(R[0], P);		// P
	yDBL(R[1], P, A);		// [2]P
	yADD(R[2], R[1], R[0], P);	// [3]P

	// main loop
	uint8_t j;
	uint32_t tmp = ADDITION_CHAIN[i];
	for(j = 0; j < ADDITION_CHAIN_LENGTH[i]; j++)
	{
		if (isinfinity(R[tmp & 0x1]) == 1)
			yDBL(T, R[2], A);
		else
			yADD(T, R[2], R[(tmp & 0x1) ^ 0x1], R[tmp & 0x1]);
		// updating
		point_copy(R[0], R[(tmp & 0x1) ^ 0x1]);
		point_copy(R[1], R[2]);
		point_copy(R[2], T);
		
		tmp >>= 1;
	};
	point_copy(Q, R[2]);	// At last, R[2] is equal to [l_{i}]P
};// Cost ~ 1.5*Ceil[log_2(l)]*(4M + 2S)

/* ------------------------------------------------------------------------------- *
   elligator()
   Inputs: the Edwards curve constant A[0]:=a, and A[1]:=(a - d), and an integer 
           number i in {0, ..., 15};
   output: the projective Edwards y-coordinates of y(T_{+}) = YT_{+}/ZT_{+} and
           y(T_{-}) = YT_{-}/ZT_{-} such that the y-coordinate of the image of the
           affine points of T_{+} and T_{-} in the Montgomery curve belong to F_p 
           and F_{p^2}\F_p, respectively.
 * ------------------------------------------------------------------------------- */
void elligator(proj T_plus, proj T_minus, const proj A)
{
	set_zero(T_plus[0], NUMBER_OF_WORDS);			// Initial value is zero
	set_zero(T_minus[0], NUMBER_OF_WORDS);			// Initial value is zero

	// u is randomly selecting from {2, ..., (p-1)/2}
	fp u;
	fp_random(u);
	while ( compare(u, (uint64_t *)p_minus_1_halves, NUMBER_OF_WORDS) > 0)
		fp_random(u);

	fp_mul(u, u, R_squared_mod_p);	// mapping u into the Montgomery domain

	// ---
	fp tmp, u2_plus_1, Cu2_minus_1, tmp_0, tmp_1, alpha, beta;
	set_zero(alpha, NUMBER_OF_WORDS);			// 0
	fp_add(beta, alpha, u);						// u
	
	fp_sqr(T_plus[1], u);						// u^2
	fp_add(u2_plus_1, T_plus[1], R_mod_p);		// u^2 + 1
	fp_sub(tmp, T_plus[1], R_mod_p);			// u^2 - 1
	fp_mul(Cu2_minus_1, A[1], tmp);					// C' * (u^2 - 1)

	// The goal is to evaluate in the projective Montgomery curve isomorphic to A
	fp_sub(T_minus[1], A[0], A[1]);					// A' := 2 * (a + d) and C' := (a - d) are
	fp_add(T_minus[1], T_minus[1], A[0]);			// projective constants of the isomorphic
	fp_add(T_minus[1], T_minus[1], T_minus[1]);		// Montgomery curve

	fp_mul(tmp_0, T_minus[1], Cu2_minus_1);			// A' * C' * (u^2 - 1)

	fp_sqr(tmp_1, T_minus[1]);				// (A')^2
	fp_mul(tmp_1, tmp_1, T_plus[1]);		// (A' * u)^2
	fp_sqr(tmp, Cu2_minus_1);				// [C' * (u^2 - 1)]^2
	fp_add(tmp_1, tmp_1, tmp);				// (A' * u)^2 + [C' * (u^2 - 1)]^2
	
	fp_mul(tmp, tmp_0, tmp_1);				// {A' * C' * (u^2 - 1)} * {(A' * u)^2 + [C' * (u^2 - 1)]^2} =? { [C' * (u^2 - 1)]^2 * w}^2
	
	//
	fp_cswap(alpha, beta ,iszero(tmp, NUMBER_OF_WORDS));	// alpha = 0 if A' = 0; alpha = u otherwise
	fp_mul(u2_plus_1, alpha, u2_plus_1);					// u2_plus_1 = 0 if A' != 0; u2_plus_1 = u^3 + u
	fp_mul(alpha, alpha, Cu2_minus_1);						// alpha * C' * (u^2 - 1)

	// the projective y-coordinate of T_{+} or T_{-}
	fp_add(T_plus[0], T_plus[0], T_minus[1]);			// [C' * (u^2 - 1) * v] = A'
	// the projective y-coordinate of T_{-} or T_{+}
	fp_sub(T_minus[0], T_minus[0], T_minus[1]);			// -[C' * (u^2 - 1) * v] = -A'
	fp_mul(T_minus[0], T_minus[0], T_plus[1]);			// -[ (C' * v) + A'] * (u^2) = -A' * (u^2)

	fp_add(T_plus[0], T_plus[0], alpha);				//  A' + alpha*C'*(u^2-1)
	fp_sub(T_minus[0], T_minus[0], alpha);				// -A' * (u^2) - alpha*C'*(u^2-1)

	fp_add(tmp, tmp, u2_plus_1);						// Now, if A'=0 then tmp = u^3 + u
	// Only one legendre symbol computation is required
	uint8_t legendre_symbol = fp_issquare(tmp) &  0x1;
	fp_cswap(T_plus[0], T_minus[0], legendre_symbol ^ 1);		// constant-time swap for determining T_{+} or T_{-}.
	
	// Finally, we mapping the points into the Edward's curves
	// T_{+}
	fp_add(T_plus[1], T_plus[0], Cu2_minus_1);
	fp_sub(T_plus[0], T_plus[0], Cu2_minus_1);
	// T_{-}
	fp_add(T_minus[1], T_minus[0], Cu2_minus_1);
	fp_sub(T_minus[0], T_minus[0], Cu2_minus_1);

	FP_ADD_COMPUTED += 16;
	FP_SQR_COMPUTED += 3;
	FP_MUL_COMPUTED += 8;

	FP_MUL_COMPUTED += 1;	// This multiplication is for mapping the input u<-{2, ..., (p-1)/2} into the Montgomery domain
};// Cost : 1(legendre s.) + 8M + 3S + 16a

/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
void cofactor_multiples(proj P[], const proj A, int8_t lower, int8_t upper)
{
	assert(lower < upper);

	if ( (upper - lower) == 1)
		return;

	int8_t mid = lower + (upper - lower + 1) / 2;

	point_copy(P[mid], P[lower]);
	for (int8_t i = lower; i < mid; ++i)
		yMUL(P[mid], P[mid], A, i);

	for (int8_t i = mid; i < upper; ++i)
		yMUL(P[lower], P[lower], A, i);

	cofactor_multiples(P, A, lower, mid);
	cofactor_multiples(P, A, mid, upper);
};

/* never accepts invalid keys. */
uint8_t validate(const proj A)
{
	do {

		proj P[N];

		fp_random(P[0][0]);
		while ( compare(P[0][0], (uint64_t *)p, NUMBER_OF_WORDS) > 0)
			fp_random(P[0][0]);

		set_zero(P[0][1], NUMBER_OF_WORDS);
		fp_add(P[0][1], P[0][1], R_mod_p);	// Z is set to 1 (in montgomery domain)
		
		yDBL(P[0], P[0], A); // mult. by [2]
		yDBL(P[0], P[0], A); // mult. by [2]

		cofactor_multiples(P, A, 0, N);

		uint16_t bits_of_the_order = 0;
		for (uint8_t i = N - 1; i < N; --i) {

			/* we only gain information if [(p+1)/l] P is non-zero */
			if (isinfinity(P[i]) != 1)
			{
				yMUL(P[i], P[i], A, i);

				if (isinfinity(P[i]) != 1)
				{
					/* P does not have order dividing p+1. */
					return 0;
				}

				bits_of_the_order += BITS_OF_L[i];	// in this case order ~ 2^bits_of_the_order > 4*p^(1/2) ~ 2^BITS_OF_SQRT_OF_P

				if (bits_of_the_order > BITS_OF_4SQRT_OF_P)
				{ 
					/* returns borrow */
					/* order > 4 sqrt(p), hence definitely supersingular */
					return 1;
				}
			}
		}

		/* P didn't have big enough order to prove supersingularity. */
	} while (1);
};

