#include <vecimg/rgb232.h>
#include <vecimg/stb_image.h>
#include <vecimg/stb_image_write.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <math.h>

void print_usage(const char *name) {
	fprintf(stderr,
		"Usage: %s [-b block size] [-q quality] [-cd] input output\n"
		"\t-b : block size, mapping input blocks to output pixels\n"
		"\t-e : encode an image [input] to a vectorized image [output]\n"
		"\t-d : decode a vectorized image [input] to a plain image [output]\n"
		"\t-q : 'quality', sets starting granularity for encode, larger is faster",
		name
	);
}

enum mode {
	Invalid,
	Encode,
	Decode,
};

#define MAX(A, B) (((A) > (B))? (A) : (B))
#define MIN(A, B) (((A) < (B))? (A) : (B))

struct pximg {
	int width, height, channels;
	union {
		uint8_t  *data;
		uint32_t *pixels;
	};
};

#define CHANNEL_RED(X)    (X & 0xff)
#define CHANNEL_GREEN(X) ((X & 0xff00)     >> 8)
#define CHANNEL_BLUE(X)  ((X & 0xff0000)   >> 16)
#define CHANNEL_ALPHA(X) ((X & 0xff000000) >> 24)

static inline
float dot(float a[2], float b[2]) {
	return (a[0] * b[0]) + (a[1] * b[1]);
}

uint32_t swap_lower3_bytes(uint32_t px) {
	uint8_t r = px & 0xff;
	uint8_t g = (px >> 8) & 0xff;
	uint8_t b = (px >> 16) & 0xff;
	uint32_t a = px & 0xff000000;
	//uint32_t a = 0xff000000;

	return a
		| (r << 16)
		| (g << 8)
		| (b);
}

static inline
float mix(float a, float b, float amount) {
	return a*amount * b*(1.f - amount);
}


static inline
size_t pximg_stride(struct pximg *img) {
	return img->channels * img->width;
}

static inline
float rotation_error(uint32_t main, uint32_t sec,
                     uint32_t *pixbuf, unsigned blocksize,
                     uint8_t rot, uint8_t off)
{
	float ret = 0.f;

	float rad = M_PI * (rot / 255.f);
	float offset = ((off / 255.f) * 2.f) - 1.f;

	float pnorm[2] = {
		cosf(rad),
		sinf(rad),
	};

	float inc = 2.f / blocksize;
	float by = -1.f;
	for (unsigned y = 0; y < blocksize; y++, by += inc) {
		float bx = -1.f;
		for (unsigned x = 0; x < blocksize; x++, bx += inc) {
			float pos[2] = {bx, by};
			float dist = dot(pnorm, pos) + offset;

			uint32_t px = (dist < 0.f)? main : sec;
			size_t idx = x + y*blocksize;

			int delta_r = (int)CHANNEL_RED(pixbuf[idx])   - (int)CHANNEL_RED(px);
			int delta_g = (int)CHANNEL_GREEN(pixbuf[idx]) - (int)CHANNEL_GREEN(px);
			int delta_b = (int)CHANNEL_BLUE(pixbuf[idx])  - (int)CHANNEL_BLUE(px);

			ret += (delta_r * delta_r) / 65025.f;
			ret += (delta_g * delta_g) / 65025.f;
			ret += (delta_b * delta_b) / 65025.f;
		}
	}

	return ret;
}

typedef struct encoded {
	uint8_t first;
	uint8_t second;
	uint8_t rotation;
	uint8_t offset;
} encoded;

enum swapmode {
	Given,
	Swapped,
	Discard,
};

encoded encode_rotation(uint8_t main, uint8_t sec,
                        uint32_t *pixbuf, unsigned blocksize,
                        unsigned initial_step)
{
	float min = HUGE_VALF;
	uint8_t currot = 0;
	uint8_t curoff = 0;
	uint8_t curmain = 0;
	uint8_t cursec  = 0;

	int off_lower = 0, off_upper = 0x100;
	int rot_lower = 0, rot_upper = 0x100;
	int step = initial_step;

	if (main == sec) {
		// nothing to do if the colors are the same
		return (encoded){main, main, 0, 0};
	}

	uint32_t colora = rgb332_to_rgb8(main); 
	uint32_t colorb = rgb332_to_rgb8(sec); 

	while (step) {
		//unsigned mode = Given;
		for (unsigned mode = Given; mode <= Discard; mode++) {
			for (int off = off_lower; off < off_upper; off += step) {
				for (int rot = rot_lower; rot < rot_upper; rot += step) {
					uint32_t a, b;

					switch (mode) {
						case Given:
							a = colora;
							b = colorb;
							break;

						case Swapped:
							a = colorb;
							b = colora;
							break;

						case Discard:
							//a = colorb;
							//b = colora;
							a = b = colora;
							break;
					}

					float err = rotation_error(a, b, pixbuf, blocksize, rot, off);

					if (err < min) {
						min = err;

						curmain = rgb8_to_rgb332(a);
						cursec  = rgb8_to_rgb332(b);
						currot = rot;
						curoff = off;
					}
				}
			}
		}

		off_lower = MAX(0, (int)curoff - step - step/2);
		rot_lower = MAX(0, (int)currot - step - step/2);

		off_upper = MIN(0x100, (int)curoff + step + step/2);
		rot_upper = MIN(0x100, (int)currot + step + step/2);

		step /= 2;
	}

	/*
	fprintf(stderr, "encoding block with error %g, rotation %u, offset %u\n",
			min, curmin, curoff);
			*/

	return (encoded){
		curmain,
		cursec,
		currot,
		curoff,
	};
}

uint32_t encode_block(struct pximg *in,
                      size_t x,
                      size_t y,
                      unsigned blocksize,
                      unsigned initial_step)
{
	size_t str = pximg_stride(in);

	uint32_t counts[0x100];
	uint32_t pixbuf[blocksize * blocksize];

	memset(counts, 0, sizeof(counts));
	memset(pixbuf, 0, sizeof(pixbuf));

	for (size_t by = 0; by < blocksize; by += 1) {
		for (size_t bx = 0; bx < blocksize; bx += 1) {
			uint32_t px = 0;

			// TODO: alpha
			for (int chan = 0; chan < in->channels && chan < 3; chan++) {
				size_t xidx = x + bx*in->channels + chan;
				size_t yidx = y + by;

				px |= in->data[xidx + str*yidx] << ((2 - chan) * 8);
			}

			uint8_t conv = rgb8_to_rgb332(px);
			counts[conv]++;

			// TODO: converting before correction results in no color splotches,
			//       but less accurate/jagged edges, need to figure out a way
			//       to avoid that
			//       (could try dithering to rgb332 before)
			//pixbuf[bx + blocksize*by] = rgb332_to_rgb8(conv);
			pixbuf[bx + blocksize*by] = px;
		}
	}

	unsigned max = 0;
	unsigned sec = 0;
	unsigned maxcount = 0;
	unsigned seccount = 0;

	for (unsigned i = 0; i < 0x100; i++) {
		// >= to prefer brighter colors
		if (counts[i] >= maxcount) {
			sec = max;
			max = i;

			seccount = maxcount;
			maxcount = counts[i];

		} else if (counts[i] >= seccount) {
			sec = i;
			seccount = counts[i];
		}
	}

	/* ignore small secondary counts TODO: configurable */
	//if (counts[sec] < blocksize/3) {
	if (seccount == 0) {
		sec = max;
	}

	encoded out = encode_rotation(max, sec, pixbuf, blocksize, initial_step);

	return out.first
		| (out.second << 8)
		| (out.rotation << 16)
		| (out.offset << 24);
}

void decode_block(struct pximg *out, size_t x, size_t y,
                  uint32_t px, unsigned blocksize)
{
	size_t str = pximg_stride(out);

	float rad = M_PI * (CHANNEL_BLUE(px) / 255.f);
	float pnorm[2] = {
		cosf(rad),
		sinf(rad),
	};

	uint32_t a = rgb332_to_rgb8(CHANNEL_RED(px));
	uint32_t b = rgb332_to_rgb8(CHANNEL_GREEN(px));

	float inc = 2.0 / blocksize;
	float off = ((CHANNEL_ALPHA(px) / 255.f) * 2.f) - 1.f;

	size_t oy = 0;
	for (float by = -1.f; oy < blocksize; by += inc, oy++) {
		size_t ox = 0;

		for (float bx = -1.f; ox < blocksize; bx += inc, ox++) {
			float pos[2] = {bx, by};
			float dist = dot(pnorm, pos) + off;
			uint32_t px;

			px = (dist < 0.f)? a : b;
			// XXX: no alpha yet
			px |= 0xff000000;

			size_t offx = x * blocksize + ox;
			size_t offy = (y*blocksize + oy) * out->width;

			out->pixels[offx + offy] = swap_lower3_bytes(px);
		}
	}
}

struct pximg *encode(struct pximg *in,
                     unsigned blocksize,
                     unsigned initial_step)
{
	size_t width  = in->width * in->channels;
	size_t height = in->height;
	size_t xinc = blocksize * in->channels;

	// TODO: should check malloc...
	struct pximg *ret = malloc(sizeof(struct pximg));

	ret->width  = (in->width  / blocksize) + !!(in->width  % blocksize);
	ret->height = (in->height / blocksize) + !!(in->height % blocksize);
	ret->channels = 4;
	ret->data = calloc(1, sizeof(uint32_t[ret->width * ret->height]));

	unsigned blocks = ret->width * ret->height;
	unsigned processed = 0;

	for (size_t y = 0, oy = 0; y + blocksize <= height; y += blocksize, oy++) {
		for (size_t x = 0, ox = 0; x + xinc <= width; x += xinc, ox++) {
			uint32_t px = encode_block(in, x, y, blocksize, initial_step);

			ret->pixels[ox + (ret->width * oy)] = px;

			fprintf(stderr, "  [%c] %u/%u          \r", "qpbd"[ox&3], processed++, blocks);
			fflush(stdout);
		}
	}

	putchar('\n');
	return ret;
}

struct pximg *decode(struct pximg *in, unsigned blocksize) {
	// TODO: proper check
	assert(in->channels == 4);

	// TODO: should check malloc...
	struct pximg *ret = malloc(sizeof(struct pximg));

	ret->width  = in->width  * blocksize;
	ret->height = in->height * blocksize;
	ret->channels = 4;
	ret->data = malloc(sizeof(uint32_t[ret->width * ret->height]));

	memset(ret->data, 0xff, sizeof(uint32_t[ret->width * ret->height]));

	size_t str = pximg_stride(in);

	for (size_t y = 0; y < in->height; y++) {
		for (size_t x = 0; x < in->width; x++) {
			decode_block(ret, x, y, in->pixels[x + y*in->width], blocksize);
		}
	}

	return ret;
}

struct pximg load_pximg(const char *filename) {
	struct pximg ret;

	ret.data = stbi_load(filename, &ret.width, &ret.height, &ret.channels, 0);
	return ret;
}

int main(int argc, char *argv[]) {
	int opt;
	int blocksize = 8;
	int initial_step = 16;
	enum mode mode = Invalid;

	const char *fin, *fout;

	while ((opt = getopt(argc, argv, "edpb:q:")) != -1) {
		switch (opt) {
			case 'b':
				blocksize = atoi(optarg);
				break;

			case 'e':
				mode = Encode;
				break;

			case 'd':
				mode = Decode;
				break;

			case 'p':
				// TODO: dump palette
				break;

			case 'q':
				initial_step = atoi(optarg);
				if (initial_step <= 0) {
					fprintf(stderr, "Invalid quality setting %d, must be >= 1\n",
					        initial_step);
					return 1;
				}
				break;

			default:
				fprintf(stderr, "Error: invalid option '-%c'\n", opt);
				print_usage(argv[0]);
				return 1;
		}
	}

	printf("blocksize: %d\n", blocksize);

	if (optind + 2 <= argc) {
		printf("a: %s, b: %s\n", argv[optind], argv[optind + 1]);
		fin  = argv[optind];
		fout = argv[optind + 1];

	} else {
		fprintf(stderr, "Error: no input/output\n");
		return 1;
	}

	if (mode == Invalid) {
		fprintf(stderr, "Error: encode/decode mode not set.\n");
		print_usage(argv[0]);
		return 1;
	}

	struct pximg input = load_pximg(fin);
	FILE *fp = fopen(fout, "w");

	if (!fp || !input.data) {
		fprintf(stderr, "couldn't do the thing\n");
		return 1;
	}

	if (mode == Encode) {
		struct pximg *foo = encode(&input, blocksize, initial_step);

		if (foo) {
			stbi_write_png(fout, foo->width, foo->height, foo->channels,
			               foo->data, foo->width * foo->channels);
		}

	} else if (mode == Decode) {
		struct pximg *foo = decode(&input, blocksize);

		if (foo) {
			stbi_write_png(fout, foo->width, foo->height, foo->channels,
			               foo->data, foo->width * foo->channels);

		} else {
			fprintf(stderr, "borked\n");
			// TODO: error
		}
	}

	return 0;
}
