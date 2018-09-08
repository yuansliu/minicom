#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "kvec.h"
#include "minicom.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/*static inline uint64_t hash64_(uint64_t key)
{
	key = (~key + (key << 21)); // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)); // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)); // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}
*/
/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param p      minimizers; p->a[i].x is the 2k-bit hash value;
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
/*void mm_sketch_lh(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			kmer[0] = (kmer[0] << 2 | c);           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			if (++l >= k)
				info.x = hash64_(kmer[z]), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, *p, min);
}
*/
void mm_sketch_lh_ori(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			if (++l >= k)
				info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		} else l = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k) kv_push(mm128_t, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1) kv_push(mm128_t, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, *p, min);
}

/*void mm_sketch_one(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), kmer[2] = {0,0};
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);
	// if (debug) fprintf(stderr, "mask: %lu\n", mask);
	assert(len > 0 && k > 0);
	// w += f;
	// for (i = 0, l = 0; i < len; i++) {
	for (i = 0, l = 0; i < len; i++) {
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);

		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c);           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		// if (debug) fprintf(stderr, "%lu %lu\n", kmer[0], kmer[1]);
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		// if (debug) fprintf(stderr, "%lu %lu %d\n", kmer[0], kmer[1], z);

		if (++l >= k) {
			info.x = hash64_(kmer[z]), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		if (info.x < min.x) {
			min = info;
		}
		// if (debug) fprintf(stderr, "%lu, %lu\n", min.x, min.y);
	}
	// if (debug) fprintf(stderr, "%lu\n", UINT64_MAX);
	// if (min.x != UINT64_MAX) {
		// kv_push(mm128_t, *p, min);
	// }
	*p = min;
	// if (debug) fprintf(stderr, "\n");

}
*/
/*void mm_sketch_range(const char *str, int begin, int end, int k, uint32_t rid, mm128_t *p) // [begin, end]
{
	uint64_t shift1 = 2 * (k - 1), kmer[2] = {0,0};
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	assert(end - begin + 1 >= k && k > 0);
	for (i = begin, l = 0; i <= end; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c);           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand

		if (++l >= k) {
			info.x = hash64_(kmer[z]), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	*p = min;
}*/

void mm_sketch_two(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);
	// if (debug) fprintf(stderr, "mask: %llu\n", mask);
	// len -= 3;
	assert(len > 0 && k > 0);
	// w += f;
	// debug = 0;
	// for (i = 0, l = 0; i < len; i++) {
	for (i = 0, l = 0; i < len; i++) {
		// if (i >= k - 1) debug = 1;
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);

		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer

		// if (debug) fprintf(stderr, "%llu %llu\n", kmer[0], kmer[1]);

		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand

		// if (debug) fprintf(stderr, "%llu %llu %llu\n", kmer[0], kmer[1], kmer[z]);

		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		// if (debug) fprintf(stderr, "%llu %llu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		if (info.x < min.x) {
			min = info;
		}
		// if (debug) fprintf(stderr, "%lu, %lu\n", min.x, min.y);
	}
	// if (debug) fprintf(stderr, "%lu\n", UINT64_MAX);
	// if (min.x != UINT64_MAX) {
		// kv_push(mm128_t, *p, min);
	// }
	*p = min;
	// fprintf(stderr, "%llu, %llu\n", min.x, min.y);
	// exit(0);
	// if (debug) fprintf(stderr, "\n");

}

/*void mm_sketch_nth_min(const char *str, int len, int k, uint32_t rid, mm128_v *p, int num)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, s;
	// mm128_t min = { UINT64_MAX, UINT64_MAX };

	mm128_t *buf;
	buf = (mm128_t*)alloca(num * 16);
	// memset(buf, 0xff, num * 16);
	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);
	// if (debug) fprintf(stderr, "mask: %lu\n", mask);
	assert(len > 0 && k > 0);
	int maxbufidx;
	// w += f;
	// for (i = 0, l = 0; i < len; i++) {
	for (i = 0, l = 0, s = 0; i < len; ++i) {
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);

		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c)  & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		// if (debug) fprintf(stderr, "%lu %lu\n", kmer[0], kmer[1]);
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		// if (debug) fprintf(stderr, "%lu %lu %d\n", kmer[0], kmer[1], z);

		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;

			if (s < num) {
				buf[s++] = info;
			} else {
				maxbufidx = 0;
				for (j = 1; j < num; ++j) {
					if (buf[j].x > buf[maxbufidx].x) {
						maxbufidx = j;
					}
				} //find the maxinum buf
				if (info.x < buf[maxbufidx].x) {
					buf[maxbufidx] = info;
				}
			}
		}
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu, %lu\n", min.x, min.y);
	}
	for (int i = 0; i < num; ++i) {
		kv_push(mm128_t, *p, buf[i]);
	}
	// if (debug) fprintf(stderr, "\n");

}

void mm_sketch_substring_nth(const char *str, int b, int len, int k, uint32_t rid, mm128_v *p, int num)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l, s, j;

	mm128_t *buf;
	buf = (mm128_t*)alloca(num * 16);
	int maxbufidx;

	assert(len > 0 && k > 0);
	
	for (i = b, l = 0, s = 0; i < len; i++) {
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);
		// fprintf(stderr, "%c", str[i]);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		// if (debug) fprintf(stderr, "%lu %lu\n", kmer[0], kmer[1]);
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		// if (debug) fprintf(stderr, "%lu %lu %d\n", kmer[0], kmer[1], z);

		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			if (s < num) {
				buf[s++] = info;
			} else {
				maxbufidx = 0;
				for (j = 1; j < num; ++j) {
					if (buf[j].x > buf[maxbufidx].x) {
						maxbufidx = j;
					}
				} //find the maxinum buf
				if (info.x < buf[maxbufidx].x) {
					buf[maxbufidx] = info;
				}
			}
		}
	}
	// if (debug) fprintf(stderr, "%lu\n", UINT64_MAX);
	// if (min.x != UINT64_MAX) {
		// kv_push(mm128_t, *p, min);
	// }
	// fprintf(stderr, "\n");
	// *p = min;
	for (int i = 0; i < num; ++i) {
		kv_push(mm128_t, *p, buf[i]);
	}
	// if (debug) fprintf(stderr, "\n");

}

void mm_sketch_split_nth(const char *str, int len, int k, uint32_t rid, mm128_v *p, int avg) { // str is divided to n parts, the length of each part is avg
	// fprintf(stderr, "len = %d, k = %d\n", len, k);
	int begin, end = k - 2;
	avg = avg - k + 1;
	kv_init(*p);

	for (; ; ) {
		begin = end - k + 2;
		end += avg;
		if (end >= len) end = len - 1;
		if (end <= begin + k) break;
		
		// fprintf(stderr, "[%d, %d]\n", begin, end);
		if (end > begin + k) {
			mm_sketch_substring_nth(str, begin, end + 1, k, rid, p, 3);
		}

		// fprintf(stderr, "%lld\n", temp_p.x);
		// fprintf(stderr, "%lu %lu %lu\n--\n", temp_p.x, temp_p.y >> 1, temp_p.y & 1);
	}
	// *p = res;
	// kv_push(mm128_t, *p, res);
}

void mm_sketch_substring(const char *str, int b, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);
	// if (debug) fprintf(stderr, "mask: %lu\n", mask);
	// len -= 3;
	assert(len > 0 && k > 0);
	// w += f;
	// for (i = 0, l = 0; i < len; i++) {
	for (i = b, l = 0; i < len; i++) {
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);
		// fprintf(stderr, "%c", str[i]);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		// if (debug) fprintf(stderr, "%lu %lu\n", kmer[0], kmer[1]);
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		// if (debug) fprintf(stderr, "%lu %lu %d\n", kmer[0], kmer[1], z);

		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		if (info.x < min.x) {
			min = info;
		}
		// if (debug) fprintf(stderr, "%lu, %lu\n", min.x, min.y);
	}
	// if (debug) fprintf(stderr, "%lu\n", UINT64_MAX);
	// if (min.x != UINT64_MAX) {
		// kv_push(mm128_t, *p, min);
	// }
	// fprintf(stderr, "\n");
	*p = min;
	// if (debug) fprintf(stderr, "\n");

}

void mm_sketch_nth_one(const char *str, int len, int k, uint32_t rid, mm128_v *p, int n) { // str is divided to n parts
	// fprintf(stderr, "len = %d, k = %d\n", len, k);
	int avg = (len - k + 1)/n;
	// fprintf(stderr, "avg = %d\n", avg);
	mm128_t res = {0, 0}, temp_p;
	int begin, end = k - 2;
	for (int times = 0; times < n; ++times) {
		begin = end - k + 2;
		end += avg;
		if (times + 1 == n) end = len - 1;
		// fprintf(stderr, "[%d, %d]\n", begin, end + 1);

		mm_sketch_substring(str, begin, end + 1, k, rid, &temp_p);

		// fprintf(stderr, "%lld\n", temp_p.x);
		// fprintf(stderr, "%lu %lu %lu\n--\n", temp_p.x, temp_p.y >> 1, temp_p.y & 1);

		if (times == 0 || (times > 0 && temp_p.x > res.x)) {
			res = temp_p;
		}
	}
	// *p = res;
	kv_push(mm128_t, *p, res);
}

void mm_sketch_nth(const char *str, int len, int k, uint32_t rid, mm128_v *p, int n) { // str is divided to n parts
	// fprintf(stderr, "len = %d, k = %d\n", len, k);
	int avg = (len - k + 1)/n;
	// fprintf(stderr, "avg = %d\n", avg);
	mm128_t temp_p;
	int begin, end = k - 2;
	for (int times = 0; times < n; ++times) {
		begin = end - k + 2;
		end += avg;
		if (times + 1 == n) end = len - 1;
		// fprintf(stderr, "[%d, %d]\n", begin, end + 1);

		mm_sketch_substring(str, begin, end + 1, k, rid, &temp_p);

		// fprintf(stderr, "%lld\n", temp_p.x);
		// fprintf(stderr, "%lu %lu %lu\n--\n", temp_p.x, temp_p.y >> 1, temp_p.y & 1);

		// if (times == 0 || (times > 0 && temp_p.x > res.x)) {
		// 	res = temp_p;
		// }
		kv_push(mm128_t, *p, temp_p);
	}
	// *p = res;
}

void mm_sketch_one_ori(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), kmer[2] = {0,0};
	int i, l;
	mm128_t min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);
	// if (debug) fprintf(stderr, "mask: %lu\n", mask);
	assert(len > 0 && k > 0);
	// w += f;
	for (i = 0, l = 0; i < len; i++) {
	// for (i = 2, l = 0; i < len; i++) {
		// if (rid == 13148457)  fprintf(stderr, "%2d, \n", i);
		// if (debug) fprintf(stderr, "%2d:\n", i);

		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c);           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		// if (debug) fprintf(stderr, "%lu %lu\n", kmer[0], kmer[1]);
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		// if (debug) fprintf(stderr, "%lu %lu %d\n", kmer[0], kmer[1], z);

		if (++l >= k) {
			info.x = hash64_(kmer[z]), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		// if (debug) fprintf(stderr, "%lu %lu\n", info.x, info.y);
		if (info.x < min.x) {
			min = info;
		}
		// if (debug) fprintf(stderr, "%lu, %lu\n", min.x, min.y);
	}
	// if (debug) fprintf(stderr, "%lu\n", UINT64_MAX);
	// if (min.x != UINT64_MAX) {
		// kv_push(mm128_t, *p, min);
	// }
	*p = min;
	// if (debug) fprintf(stderr, "\n");

}
*/

/*
void mm_sketch_one_max(const char *str, int len, int k, uint32_t rid, mm128_t *p)
{
	uint64_t shift1 = 2 * (k - 1), kmer[2] = {0,0};
	int i, l;
	mm128_t min = { 0, 0 };

	assert(len > 0 && k > 0);
	for (i = 0, l = 0; i < len; i++) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { 0, 0 };
		int z;
		kmer[0] = (kmer[0] << 2 | c);           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand

		if (++l >= k) {
			info.x = hash64_(kmer[z]), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x > min.x) {
			min = info;
		}
	}
	*p = min;
}

void mm_sketch0(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "in mm_sketch(): %s\n", str);

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = 0; i < w - 1; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = w - k, l = 0; i < 2*w - k - 1; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = len - w + 1, l = 0; i < len; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = len - 2*w + k + 1, l = 0; i < len - w + k; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

}

void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "%s\n", str);

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = 0, l = 0; i < w - 1; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = len - w + 1, l = 0; i < len; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

}


void mm_sketch4(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, l;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	// int debug = 1;
	// if (debug) fprintf(stderr, "%s\n", str);

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (i = l = 0; i < w - 1; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = w - k, l = 0; i < 2*w - k - 1; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = len - w + 1, l = 0; i < len; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

	min.x = min.y = UINT64_MAX;
	for (i = len - 2*w + k + 1, l = 0; i < len - w + k; i++) {
		// if (debug) fprintf(stderr, "%2d, ", i);
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		int z;
		kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
		if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
		z = kmer[0] < kmer[1]? 0 : 1; // strand
		if (++l >= k) {
			info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
		}
		if (info.x < min.x) {
			min = info;
		}
	}
	if (min.x != UINT64_MAX) {
		kv_push(mm128_t, *p, min);
	}
	// if (debug) fprintf(stderr, "\n");

}

void mm_sketchAll(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int bg, i, l;
	mm128_t *buf, min = { UINT64_MAX, UINT64_MAX };

	assert(len > 0 && w > 0 && k > 0);
	buf = (mm128_t*)alloca(w * 16);
	memset(buf, 0xff, w * 16);

	for (bg = 0; bg <= len - w; ++bg) {
		min.x = min.y = UINT64_MAX;
		kmer[0] = 0, kmer[1] = 0;
		for (i = bg, l = 0; i < bg + w; ++i) {
			int c = seq_nt4_table[(uint8_t)str[i]];
			mm128_t info = { UINT64_MAX, UINT64_MAX };
			int z;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			if (++l >= k) {
				info.x = hash64(kmer[z], mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
			if (info.x < min.x) {
				min = info;
			}
		}
		if (min.x != UINT64_MAX && !(p->n >= 1 && min.x == p->a[p->n - 1].x && min.y == p->a[p->n - 1].y) ) {
			kv_push(mm128_t, *p, min);
		}
	}
}*/
