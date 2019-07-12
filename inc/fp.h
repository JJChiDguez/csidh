#ifndef _Fp_H_
#define _Fp_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "rng.h"

#if defined FP_512
	#define N 74			// Number of l_i's such that l_i | [(p+1)/4]
	#define LOG2_OF_N_PLUS_ONE 8
	#define NUMBER_OF_WORDS 8	// Number of 64-bit words

// If a new field arithmetic will be used, just add it but taking in count that it is required:
//         fp_add(output, input_x, input_y),
//         fp_sub(output, input_x, input_y),
//         fp_mul(output, input_x, input_y),
//         fp_sqr(output, input_x),
//         fp_inv(input), 
//         fp_issquare(input), and
//         fp_random(input);
// This above functions must be implemented allowing that the output variable can be one of the inputs.

#endif

typedef uint64_t fp[NUMBER_OF_WORDS]      __attribute__((aligned(64)));		// 512-bits integer number in Montgomery domain (To be used with the patching)

extern const fp p;
extern const fp R_mod_p;
extern const fp R_squared_mod_p;	// required for mapping a random fp element 2 <= u <= (p-1)/2 into the Montgomery domain
extern const fp p_minus_1_halves;	// (p-1)/2 used for find a random fp element 2 <= u <= (p-1)/2

// All operations are perfomed in the Montgomery domain
void fp_cswap(fp x, fp y, uint8_t c);

void fp_add(fp c, const fp a, const fp b);
void fp_sub(fp c, const fp a, const fp b);
void fp_mul(fp c, const fp a, const fp b);
void fp_sqr(fp b, const fp a);

void fp_inv(fp x);
uint8_t fp_issquare(fp const x);
void fp_random(fp x);	// This function should be modified in order to have a better random function: e.g., shake256.

#define set_zero(x, NUM)\
	memset(x, 0, sizeof(uint64_t) * NUM);

#define set_one(x, NUM) {\
	int i;\
	x[0] = 1;   \
	for (i=1; i < NUM; i++)\
			x[i] = 0;\
}

#define copy(x, y, NUM)\
	memcpy(x, y, sizeof(uint64_t)*NUM);

/* ------------------------------------------------------------- *
   compare()
   inputs: two integer numbers x and y, and the number of 64-bits 
           words of x and y.
   outputs:
            +1 if x > y,
            -1 if x < y, or
             0 if x = y.
 * ------------------------------------------------------------- */
static inline int compare(uint64_t *x, uint64_t *y, int NUM)
{
	int i;
	for (i=NUM-1; i >= 0; i--) 
    {
		if (x[i] != y[i]) 
            return x[i] > y[i] ? 1 : -1; 
	}
	return 0;
}

/* ------------------------------------------------------------- *
   iszero()
   inputs: two integer numbers x and y, and the number of 64-bits 
           words of x and y.
   outputs:
             1 if x = 0,
             0 otherwise
 * ------------------------------------------------------------- */
static inline int iszero(uint64_t *x, int NUM)
{
	int i;
	uint64_t c = 0;
	for (i=NUM-1; i >= 0; i--) 
    	c |= x[i];
	return (c == 0);
}

/* ------------------------------------------------------------- *
   fp_print()
   inputs: an integer number x, the number of 64-bits words of x,
           an integer in {1,0}, and a string.
   Note: This function prints the integer x
 * ------------------------------------------------------------- */
static void fp_print(uint64_t *x, int NUM, int TYPE, char *c)
{
    int i;
    printf("%s := 0x", c);
    for(i=NUM-1; i > -1; i--){
        if(TYPE == 1)
            printf("%.16" PRIX64 " ", x[i]);
        else
            printf("%.16" PRIX64 "", x[i]);
    }
    printf(";\n");
}

/* decision bit b has to be either 0 or 1 */
static void cmov(int8_t *r, const int8_t a, uint32_t b)
{
	uint32_t t;
	b = -b; /* Now b is either 0 or 0xffffffff */
	t = (*r ^ a) & b;
	*r ^= t;
}

/* check if a and b are equal in constant time  */
static uint32_t isequal(uint32_t a, uint32_t b)
{
	//size_t i;
	uint32_t r = 0;
	unsigned char *ta = (unsigned char *)&a;
	unsigned char *tb = (unsigned char *)&b;
	r = (ta[0] ^ tb[0]) | (ta[1] ^ tb[1]) | (ta[2] ^ tb[2]) |  (ta[3] ^ tb[3]);
	r = (-r);
	r = r >> 31;
	return (uint32_t)(1-r);
}

/* get priv[pos] in constant time  */
static uint32_t lookup(size_t pos, int8_t const priv[])
{
	int b;
	int8_t r = priv[0];
	for(size_t i = 1; i < N; i++)
	{
		b = isequal(i, pos);
		cmov(&r, priv[i], b);
	}
	return r;
}

// constant-time comparison: -1 if x < y, 0 otherwise.
static int32_t issmaller(int32_t x, int32_t y)
{
	int32_t xy = x ^ y;
	int32_t c = x - y;
	c ^= xy & (c ^ x);
	return (c >> 31);
}
#endif /* Fp */
