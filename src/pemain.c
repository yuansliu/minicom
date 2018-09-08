#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
// #include "minimap.h"
#include "bseq.h"
#include "kvec.h"
#include "breads.h"
#include "config.h"

int mm_verbose = 3;
double mm_realtime0;
char *folder;

#define MM_VERSION "0.0-Thu-12-Dec"

// const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

void liftrlimit()
{
#ifdef __linux__
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
#endif
}

reads_t *reads;//global
uint32_t num_locks = 0x1000000;
// int diff_threshold = 2;
int diff_threshold = 4;
// int first_mininum = 100; // 4 6
int first_mininum = 4; // 4 6
int n_threads;
int maxmatch;
int maxsearch = 500;


// FILE *atfiles;
// pthread_mutex_t atfiles_mutex;

int main(int argc, char *argv[])
{
	int i;
 
	liftrlimit();
	mm_realtime0 = realtime();

	// atfiles = fopen("atfiles.seq", "w");
	// pthread_mutex_init(&atfiles_mutex, 0);
	folder = argv[2];

	// char str[100] = "SRR445718.fastq";
	// char str[100] = "SRR490961.fastq";
	// char str[100] = "SRR490961_subset.fastq";
	// char str[100] = "subset.fastq";
	// fprintf(stderr, "%s\n", argv[1]);
	// fprintf(stderr, "%s\n", folder);
	// char str[100] = "head100.fastq";
	// char str[100] = "SRR034940_1.fastq";
	// char str[100] = "SRR1294122.fastq";

	maxmatch = readlen/2;
	n_threads = num_thr;
	// int k = 32, w = 80, b = 14, f = 9;
	// int k = 32, w = 70, b = 14, f = 15;
	// int k = 32, w = 60, b = 14, f = 19;
	// int k = 32, w = 100, b = 14, f = 0;
	// if (readlen)
	// int k = 17;

	int k = 31;
	if (readlen < 80) {
		// k = 17;
		k = 19;
	}
	int w = readlen - k/2;
	// int k = 32;
	// int w = 100;
	fprintf(stderr, "k: %d; w: %d\n", k, w);
	// for 100 length
	// int k = 31;
	// int w = 80, 
	int b = 14, f = 0;
 
	pre_process(argv[1], f, w, k, b, n_threads, readlen); //int w, int k, int b,
	
	// // fprintf(stderr, "xxx\n");
	// fprintf(stderr, "%d, %d\n", reads->n_seq, reads->seq_len);
	// fprintf(stderr, "%d, %d\n", reads->sp->allA, reads->sp->allT);

	// if (mm_verbose >= 3)
	// 	fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
	// 			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), reads->n_seq);
	
	// if (fp)  bseq_close(fp);

	// fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
	// fprintf(stderr, "[M::%s] CMD:", __func__);
	// for (i = 0; i < argc; ++i)
	// 	fprintf(stderr, " %s", argv[i]);
	// fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());

	// fclose(atfiles);
	
	return 0;
}


