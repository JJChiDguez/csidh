#ifndef _EC_H_
#define _EC_H_

#include <assert.h>
#include<x86intrin.h>
#include<stdint.h>

#include "fp.h"

// projective Edwards y-coordinates (to be used with the patching)
typedef uint64_t proj[2][NUMBER_OF_WORDS]	__attribute__((aligned(64)));
// A curve will be defined as type proj where the first and second entries will be the constants a and (a -d).

uint64_t FP_ADD_COMPUTED,	// Variable used for counting the number of field additions.
	 FP_SQR_COMPUTED,	// Variable used for counting the number of field squarings.
         FP_MUL_COMPUTED;	// Variable used for counting the number of field multiplications.

// Framework to be used: the files required must be in the folder: ./inc/fp$(BITLENGTH_OF_P)/
#include "addc.h"			// Addition chains, Public curve, public points T_{+} and T_{-}, and the list of prime factors l_i's

#if defined WITHDUMMY_1
	#include "simba_withdummy_1.h"	// csidh with dummy operations and using only one torsion point T_{+}
#elif defined WITHDUMMY_2
	#include "simba_withdummy_2.h"	// csidh with dummy operations and using two torsion point T_{+} and T_{-}
#elif defined DUMMYFREE
	#include "simba_dummyfree.h"	// dummy-free csidh using two torsion points T_{+} and T_{-}
#endif

// Functions related with the point arithmetic
int isinfinity(const proj P);			// To determine if a projective y-coordinate point is the infinity
void point_copy(proj Q, const proj P);		// To make a copy of a point
uint8_t areEqual(const proj P, const proj Q);	// To check if two points are equal

void yDBL(proj Q, const proj P, const proj A);
void yADD(proj R, const proj P, const proj Q, const proj PQ);
void yMUL(proj Q, const proj P, const proj A, uint8_t const i);

void elligator(proj T_plus, proj T_minus, const proj A);

void cofactor_multiples(proj P[], const proj A, int8_t lower, int8_t upper);
uint8_t validate(const proj A);

// Functions related with isogenies
void yISOG(proj Pk[], proj C, const proj P, const proj A, const uint8_t i);
void yEVAL(proj R, const proj Q, const proj Pk[], const uint8_t i);

// functions related with the action
void action_evaluation(proj C, const uint8_t key[], const proj A);
void random_key(uint8_t key[]);
void printf_key(uint8_t key[], char *c);

#endif /* Arithmetic and isogenies on Edwards curve */
