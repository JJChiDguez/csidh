#ifndef _ADDC_H_
#define _ADDC_H_

// The public Edward's curve with coefficients a and b. E is isomorphic to E' / F_p : y^2 = x^3 + x
// a and (a - d) in Montgomery domain represetantion 
static proj E = {
{ 0x767762E5FD1E1599, 0x33C5743A49A0B6F6, 0x68FC0C0364C77443, 0xB9AA1E24F83F56DB, 0x3914101F20520EFB, 0x7B1ED6D95B1542B4, 0x114A8BE928C8828A, 0x3793732BBB24F40},	// a
{ 0xECEEC5CBFA3C2B32, 0x678AE87493416DEC, 0xD1F81806C98EE886, 0x73543C49F07EADB6, 0x7228203E40A41DF7, 0xF63DADB2B62A8568, 0x229517D251910514, 0x6F26E6577649E80}	// a - d
};

// Shortest differential addition chains for each l_i
static uint64_t ADDITION_CHAIN[] = {
0x231, 0x324,  0x10, 0x2D8, 0x140,  0x50,  0x14, 0x108, 
0x101, 0x144, 0x148, 0x122, 0x141,  0x51, 0x12A, 0x109, 
0x118, 0x1A2, 0x181,   0x0, 0x134, 0x194, 0x185, 0x191, 
  0x1, 0x198,  0x82,  0x88,  0x81,  0x8A,  0xC0,  0xA1,  
 0xD0,  0xC2,  0x98,  0xD4,  0xC5,  0x2B,  0x40,  0xE8,
 0xE1,  0x4A,  0x60,   0xC,  0x68,  0x49,   0x0,  0x6C,
  0x4,  0x14,  0x24,   0x9,  0x2C,  0x32,  0x36,   0x1, 
 0x11,  0x18,   0x3,   0x8,   0xA,   0xD,   0x4,   0x5, 
  0x0,   0x1,   0x1,   0x0,   0x0, 0x612, 0x268, 0x312, 
 0xD1, 0x352
};

// Length of the shortest differential addition chain for each l_i
static uint8_t ADDITION_CHAIN_LENGTH[] = {
11, 11, 10, 11, 10, 10, 10, 10,
10, 10, 10, 10, 10, 10, 10, 10,
10, 10, 10,  9, 10, 10, 10, 10, 
 9, 10,  9,  9,  9,  9,  9,  9,
 9,  9,  9,  9,  9,  9,  8,  9, 
 9,  8,  8,  8,  8,  8,  7,  8, 
 7,  7,  7,  7,  7,  7,  7,  6,
 6,  6,  6,  5,  5,  5,  4,  4,
 3,  3,  2,  1,  0, 12, 11, 11,
11, 11
};

// L
static uint32_t L[] = { 
349, 347, 337, 331, 317, 313, 311, 307, 
293, 283, 281, 277, 271, 269, 263, 257, 
251, 241, 239, 233, 229, 227, 223, 211, 
199, 197, 193, 191, 181, 179, 173, 167, 
163, 157, 151, 149, 139, 137, 131, 127, 
113, 109, 107, 103, 101,  97,  89,  83,
 79,  73,  71,  67,  61,  59,  53,  47,
 43,  41,  37,  31,  29,  23,  19,  17,
 13,  11,   7,   5,   3, 587, 373, 367, 
359, 353
};

static uint16_t BITS_OF_L[] = { 
9, 9, 9, 9, 9,  9, 9, 9, 
9, 9, 9, 9, 9,  9, 9, 9, 
8, 8, 8, 8, 8,  8, 8, 8, 
8, 8, 8, 8, 8,  8, 8, 8, 
8, 8, 8, 8, 8,  8, 8, 7, 
7, 7, 7, 7, 7,  7, 7, 7, 
7, 7, 7, 7, 6,  6, 6, 6, 
6, 6, 6, 5, 5,  5, 5, 5, 
4, 4, 3, 3, 2, 10, 9, 9, 
9, 9 
};

#define BITS_OF_4SQRT_OF_P 258
#define LARGE_L 587
// The l_i's are only required for isogeny constructions
#endif /* Addition chains */
