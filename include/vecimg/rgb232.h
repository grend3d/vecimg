#pragma once
#include <stdint.h>

#define B2TO8 85
#define B3TO8 36

static inline
uint32_t rgb322_to_rgb8(uint8_t c) {
	return 
		  ((((c >> 4) & 7) * B3TO8) << 16)
		| ((((c >> 2) & 3) * B2TO8) << 8)
		| (((c        & 3) * B2TO8));
}

static inline
uint32_t rgb332_to_rgb8(uint8_t c) {
	return 
		  ((((c >> 5) & 7) * B3TO8) << 16)
		| ((((c >> 2) & 7) * B3TO8) << 8)
		| (((c        & 3) * B2TO8));
}

static inline
uint8_t rgb8_to_rgb332(uint32_t px) {
	return 
		  (((px >> 16) & 0xff) / B3TO8) << 5
		| (((px >> 8)  & 0xff) / B3TO8) << 2
		| (((px)       & 0xff) / B2TO8);
}

static inline
uint32_t rgb232_to_rgb8(uint8_t c) {
	return 
		  ((((c >> 5) & 3) * B2TO8) << 16)
		| ((((c >> 2) & 7) * B3TO8) << 8)
		| (((c        & 3) * B2TO8));
}

static inline
uint8_t rgb8_to_rgb232(uint32_t px) {
	return 
		  (((px >> 16) & 0xff) / B2TO8) << 5
		| (((px >> 8)  & 0xff) / B3TO8) << 2
		| (((px)       & 0xff) / B2TO8);
}

static inline
uint32_t rgb223_to_rgb8(uint8_t c) {
	return 
		  ((((c >> 5) & 3) * B2TO8) << 16)
		| ((((c >> 3) & 3) * B2TO8) << 8)
		| (((c        & 7) * B3TO8));
}
