#ifndef MM_BREADS_H
#define MM_BREADS_H

#include <stdbool.h>
#include <bitset>
#include <assert.h>
#include <unistd.h>
#include "bseq.h"
// #include "config.h"
#include <iostream>
#include <fstream>
#include "khash.h"
#include <algorithm>

typedef struct {
 	uint64_t x, y;
} mm128_t;

typedef struct {
 	int32_t x, y, z; // x for reads id; y for begin position; z for dir
} mm64_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;
typedef struct { size_t n, m; uint64_t *a; } uint64_v;
typedef struct { size_t n, m; uint32_t *a; } uint32_v;
typedef struct { size_t n, m; int32_t *a; } int32_v;
typedef struct { size_t n, m; char *a; } char_v;
typedef struct { size_t n, m; bool *a; } bool_v;
typedef struct { size_t n, m; mm64_t *a; } mm64_v;
typedef struct { size_t n; uint8_t a; } uint8bit_v;
   
// typedef uint64_v cluster_t;
 
typedef struct {
	int allA, allT, allN;
	uint32_v allA_id, allT_id, allAN_id, allTN_id, allN_id;
} sp_reads_t;

struct mm_idx_t;

typedef struct {
	int bpos;//begin postion
	uint32_t id;//tnum|tid
} next_t;

// typedef struct { size_t n, m; mapping_t *a; } link_t;
// typedef struct { size_t n, m; link_t *a; } link_v;

typedef struct {
	size_t n, m;
	uint64_t *a; //(uint32_t)id | (uint32_t)pos
	char *ref;
	// link_t next;
	// next_t next;
	// int bpos;	//for next cluster
	// uint32_t id;// tnum|tid //for next cluster
	int ennum; // number of encoded characters 
} cluster_t;

#define cluster_destroy(v) do { \
		(v).n = (v).m = 0;\
		free((v).ref);\
		free((v).a);\
	} while (0)

typedef struct {
	size_t n, m;
	cluster_t *a;
} cluster_v; 

typedef mm128_v cluster_bucket_t;

typedef struct {
    int n_seq, seq_len;
	int b, f, w, rw, k; //sk;//f is the position of first minimizer
	// uint32_v *a;//an array size of 256
	// uint16_t *equals;
	// uint8_t *flag;
	bool *sg_flag, **c_flag;
	bseq1_t *seq;
	sp_reads_t *sp;
	cluster_bucket_t **B;
	// cluster_bucket_t *B0;
	struct mm_idx_t **mi;
	struct mm_idx_t *sgmi;// single mi_idx_t
	struct mm_idx_t *cmi;// single mi_idx_t
	cluster_v **clusters;
	// cluster_v *clusters2;
	// link_v link;
	// uint64_v *next;
	// bool_v *pre;
	// uint64_v *pre;
	// cluster_v *clusters_left;
	// cluster_v *clusters_right;
	int single;
	uint32_v sg;
	FILE *fpa, *fpt, *fpn;
	uint32_v fpA_id, fpT_id, fpN_id;
	FILE *singleFile;
	uint32_v singleFile_id;
	FILE *Nfile; //single reads contain N
	uint32_v Nfile_id;
} reads_t;

typedef struct {
	mm128_v a;   // (minimizer, position) array
	int32_t n;
	uint64_t *p;
	void *h;
} mm_idx_bucket_t;

typedef struct mm_idx_t{
	uint32_t n;  // number of reference sequences
	mm_idx_bucket_t *B;
	// uint32_t max_occ;
	// float freq_thres;
	// int32_t *len;    // length of each reference sequence
	// char **name; // TODO: if this uses too much RAM, switch one concatenated string
} mm_idx_t; 

// extern mm128_v *mm_threads_idx;
// extern uint32_v *sketch_idx;
extern unsigned char seq_nt4_table[256];
extern const char invert_code_rule[4];
extern int clusize[20];
extern int diff_threshold;
extern int maxthr;
extern int cbthreshold;
extern int first_mininum;
extern int thr_step;
extern int max_rounds;
extern char *folder;
extern int maxmatch;
extern int rw;
extern int maxsearch;
extern int idxv;
extern int half_val;

#define getmax(a,b) ((a)>(b)?(a):(b))
#define getmin(a,b) ((a)<(b)?(a):(b))
#define getmax3(a,b,c) getmax(getmax(a,b),c)
#define getmax4(a,b,c,d) getmax(getmax3(a,b,c),d)

extern reads_t *reads;//global
extern int mm_verbose;
extern double mm_realtime0;
extern uint32_t num_locks;
extern int n_threads;
extern FILE *atfiles;
extern pthread_mutex_t atfiles_mutex;
// #ifdef __cplusplus
// extern "C" {
// #endif
	
double cputime();
double realtime();

void reverse_complement(char *seq, int seq_len);
// #ifdef _PE
void pre_process(const char *fn, const char *fn0, int f, int w, int k, int b, int n_threads);
// #else
void pre_process(const char *fn, int f, int w, int k, int b, int n_threads);
// #endif
int cmp(const void *a, const void *b);
// compute minimizers
void mm_sketch(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);
void mm_sketchAll(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);
void mm_sketch_mid(const char *str, int len, int f, int w, int k, uint32_t rid, mm128_v *p);
void mm_sketch_one(const char *str, int len, int k, uint32_t rid, mm128_t *p);
void mm_sketch_two(const char *str, int len, int k, uint32_t rid, mm128_t *p);
void mm_sketch_substring(const char *str, int b, int len, int k, uint32_t rid, mm128_t *p);
void mm_sketch_one_ori(const char *str, int len, int k, uint32_t rid, mm128_t *p);
void mm_sketch_lh(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);
void mm_sketch_lh_ori(const char *str, int len, int w, int k, uint32_t rid, mm128_v *p);
void mm_sketch_nth(const char *str, int len, int k, uint32_t rid, mm128_v *p, int n);
void mm_sketch_nth_one(const char *str, int len, int k, uint32_t rid, mm128_v *p, int n);
void mm_sketch_nth_min(const char *str, int len, int k, uint32_t rid, mm128_v *p, int num);
void mm_sketch_split_nth(const char *str, int len, int k, uint32_t rid, mm128_v *p, int avg);
void mm_sketch_substring_nth(const char *str, int b, int len, int k, uint32_t rid, mm128_v *p, int num);
void mm_sketch_range(const char *str, int begin, int end, int k, uint32_t rid, mm128_t *p);

//minimizer indexing
mm_idx_t *mm_idx_init(int b);
// void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a);
void mm_idx_destroy(mm_idx_t *mi);

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n);
void mm_idx_generation(int n_threads, mm_idx_t *mi);
void combine_cluster(int n_threads, reads_t *reads, int *index);
void search_sg_reads(int n_threads, reads_t *reads, int index);
void sg_idx(int n_threads, reads_t *reads, int threshold, int k, int num);
void sg_extend(int n_threads, reads_t *reads, int index);
void cluster_dump(int n_threads, reads_t *reads, int index);
void cluster_dump_idx(int n_threads, reads_t *reads, int index);
void cluster_dump_pe(int n_threads, reads_t *reads, int index);
void filter_cluster(int n_threads, reads_t *reads, int *index_);
void cluster_idx(int n_threads, reads_t *reads, int k, int index);
void extend_hash(int n_thread, reads_t *reads, int index, int threshold);
void extend1_hash(int n_threads, reads_t *reads, int index);
// void sg_realign(int n_threads, reads_t *reads, int kk, int index);
void sg_realign(int n_threads, reads_t *reads, int k, int index, int beg_threshold, int max_threshold);
void second_reads(int n_threads, reads_t *reads);
void kt_for_second_bucket(int n_threads, reads_t *reads, long n, int cindex);
void updateSingle();
void realign_hash(int n_thread, reads_t *reads, int index, int max_threshold);
// void realign_hash(int n_thread, reads_t *reads, int index, int beg_threshold, int max_threshold);

int cmpcluster(const void *a_, const void *b_);
int cmpcluster2(const void *a_, const void *b_);
// #ifdef ORDER
// int cmpcluster3(const void *a_, const void *b_);
// #endif

void update_reference(reads_t *reads, cluster_t *p, int n);
// void print_encode(cluster_t *p, FILE *fp);
// void print_encode(cluster_t *p, FILE *fpref, FILE *fppos, FILE *fpdir, FILE *fpdif);
// void print_encode(cluster_t *p, FILE *fpref, FILE *fppos, FILE *fpdir, FILE *fpdif_pos, FILE *fpdif_char);
// void print_encode(cluster_t *p, FILE *fpref, FILE *fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_pos, FILE *fpdif_char);
// double match_pro(char *str0, char *str1, int i_, int j_);
int match_pro(char *str0, char *str1, int i_, int j_);
bool encode_byte(char *seq, char *ref, int pos, int dir);


void radix_sort_128x(mm128_t *beg, mm128_t *end);
void radix_sort_64(uint64_t *beg, uint64_t *end);

// #ifdef __cplusplus
// }
// #endif

//v - ; fp is ofstream; x is the value; n is the number, if v.n == n output v.a, and a.n = 0;
#define DNA_push(v, fp, x) do {									\
		(v).a += ((x) << (2*(v).n));										\
		++ (v).n;										\
		if ((v).n == 4) {										\
			fp.write((char*)&(v).a, sizeof(uint8_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

#define bit_push(v, fp, x) do {									\
		(v).a += ((x) << ((v).n));										\
		++ (v).n;										\
		if ((v).n == 8) {										\
			fp.write((char*)&(v).a, sizeof(uint8_t));							\
			(v).a = (v).n = 0;							\
		}															\
	} while (0)

#endif

