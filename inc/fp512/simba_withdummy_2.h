#ifndef _SIMBA_PARAMETERS_H_
#define _SIMBA_PARAMETERS_H_

// SIMBA-(NUMBER_OF_BATCHES)-MY
#define NUMBER_OF_BATCHES 3
#define MY 8

// (each entry corresponds to the number of degree-(l_i) to be required in the action: this the one given in Onuki et al. work)
static int8_t B[] =	{
 2,  3, 3, 3, 3, 3,  3,  3, 
 3,  3, 3, 4, 4, 4,  4,  4, 
 4,  4, 4, 4, 4, 4,  4,  4, 
 4,  4, 5, 5, 5, 5,  5,  5, 
 5,  6, 6, 6, 6, 6,  7,  7, 
 7,  7, 7, 7, 7, 7,  7,  7, 
 7,  7, 8, 9, 9, 9, 10, 10,
10, 10, 9, 8, 8, 8,  7,  7,
 7,  7, 7, 6, 5, 1,  2,  2,
 2,  2
};

// (NUMBER_OF_BATCHES) different subsets (i.e., batches)
static uint8_t BATCH_0[] = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72 };
static uint8_t BATCH_1[] = { 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73 };
static uint8_t BATCH_2[] = { 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71 };

static uint8_t SIZE_OF_EACH_BATCH[NUMBER_OF_BATCHES] = {25, 25, 24};
static uint8_t *BATCHES[NUMBER_OF_BATCHES] = { BATCH_0, BATCH_1, BATCH_2 };

static uint8_t LAST_ISOGENY[NUMBER_OF_BATCHES] = { 72, 73, 71 };
static uint16_t NUMBER_OF_ISOGENIES = 404;

// The complement of each batch
static uint8_t SIZE_OF_EACH_COMPLEMENT_BATCH[NUMBER_OF_BATCHES] = { 49, 49, 50 };
static uint8_t COMPLEMENT_OF_EACH_BATCH[NUMBER_OF_BATCHES][N] = {
{  1,  2,  
   4,  5, 
   7,  8,
  10, 11,
  13, 14,
  16, 17,
  19, 20,
  22, 23,
  25, 26,
  28, 29,
  31, 32,
  34, 35,
  37, 38,
  40, 41,
  43, 44,
  46, 47,
  49, 50,
  52, 53,
  55, 56,
  58, 59,
  61, 62,
  64, 65,
  67, 68,
  70, 71,
  73, 
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74  
},
{  0,  2,
   3,  5,
   6,  8,
   9, 11,
  12, 14,
  15, 17,
  18, 20,
  21, 23,
  24, 26,
  27, 29,
  30, 32,
  33, 35,
  36, 38,
  39, 41,
  42, 44,
  45, 47,
  48, 50,
  51, 53,
  54, 56,
  57, 59,
  60, 62,
  63, 65,
  66, 68,
  69, 71,
  72, 
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74
},
{  0,  1,
   3,  4,
   6,  7,
   9, 10,
  12, 13,
  15, 16,
  18, 19,
  21, 22,
  24, 25,
  27, 28,
  30, 31,
  33, 34,
  36, 37,
  39, 40,
  42, 43,
  45, 46,
  48, 49,
  51, 52,
  54, 55,
  57, 58,
  60, 61,
  63, 64,
  66, 67,
  69, 70,
  72, 73, 
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74
}
};
#endif
