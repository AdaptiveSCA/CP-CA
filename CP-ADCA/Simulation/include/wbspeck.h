#ifndef _HWBSPECK_H_
#define _HWBSPECK_H_
#include "WBMatrix/WBMatrix.h"
#include "math.h"
#include "string.h"

#define u8 uint8_t
#define u16 uint16_t
#define u32 uint32_t
#define u64 uint64_t

#define SPECK_ROUNDS 22
#define SPECK_KEY_LEN 4

#define pt 256
#define tb32 32
#define tb64 64
#define tb128 128

#define R32(x, y, k) (x = ROR32(x, 7), x += y, x ^= k, y = ROL32(y, 2), y ^= x)
#define R64(x, y, k) (x = ROR64(x, 8), x += y, x ^= k, y = ROL64(y, 3), y ^= x)
#define R128(x, y, k) (x = ROR128(x, 8), x += y, x ^= k, y = ROL128(y, 3), y ^= x)

static const double EPS = 1e-6;
void Key_expand32(u16 key[4], u16 roundkey[SPECK_ROUNDS]);
void CP_ADCA32();

u16 ROL32(u16 x, int r);
u16 ROR32(u16 x, int r);
#endif