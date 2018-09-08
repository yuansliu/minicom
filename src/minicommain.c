#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "bseq.h"
#include "kvec.h"
#include "breads.h"
#include "config.h"

int mm_verbose = 3;
double mm_realtime0;
char *folder;
int idxv;

#define MM_VERSION "0.0-18-Aug"

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

// int diff_threshold = 4; //default
// int first_mininum = 4; // 4 6

// int diff_threshold = 8; //default
// int diff_threshold = 4; //default
// int first_mininum = 10; // 4 6

#ifdef _PE
int half_val;
#endif
//*****// for SRR870667_1******//
// int diff_threshold = 18;
// int first_mininum = 20; // 4 6
//***********//
//*****// for SRR870667_2******//
// int diff_threshold = 4;
// int first_mininum = 20; // 4 6
//***********//

//*****// for SRR327342_1******//
//*****// for SRR635193_1******//
// int diff_threshold = 4;
// int first_mininum = 10; // 4 6
//***********//

int diff_threshold = 4; //
int maxthr; //
int cbthreshold;
int first_mininum = 6; // 4 6
int max_rounds = 35;
int thr_step;
// int diff_threshold = 15; // for ERR174310_1
// int first_mininum = 10; // 4 6

// int diff_threshold = 10; // for ERR174310
// int first_mininum = 20; // 4 6


// int first_mininum = 100; // 4 6
// int first_mininum = 10; // 4 6
int n_threads;
int maxmatch;
int maxsearch = 500;
int rw;
// int maxsearch = 1000;

int main(int argc, char *argv[])
{
	int i;

	liftrlimit();
	mm_realtime0 = realtime();

	maxmatch = readlen/2;
	n_threads = num_thr;
	// maxmatch = readlen*2/3;
	
	int k = 31;
	if (readlen < 80) {
		k = 17;
		// if (readlen > 70) {
		// 	k = 17;
		// 	// k = 19;
		// 	// k = 21;
		// } else 
		// if (readlen > 60) {
		// 	k = 16;
		// } 
		// else {
		// 	k = 13;
		// }
	}
	// echo "#define inik $inik" >> src/config.h
	// echo "#define inithr $threshold" >> src/config.h
	// echo "#define ininumdict $num_dict" >> src/config.h
	// echo "#define iniw $w" >> src/config.h
	
	if (inik > 0) {
		k = inik;
	}
	if (inithr > 0) {
		diff_threshold = inithr;
	}
	rw = iniw;
	if (inim > 0) {
		first_mininum = inim;
	}
	if (inicbthr > 0) {
		cbthreshold = inicbthr;
	} else {
		cbthreshold = 2*diff_threshold;
	}
	if (inimaxrounds > 0 && inimaxrounds < max_rounds) {
		max_rounds = inimaxrounds;
	}
	thr_step = diff_threshold;

	if (inistep > 0) {
		thr_step = inistep;
	} else 
	if (thr_step > 10) {
		thr_step = 5;
	}

	// maxthr = readlen*2/3;
	maxthr = readlen/2;
	if (inimaxthr > 0) {
		maxthr = inimaxthr;
	}
	// k = 25;// for SRR870667_1
	// k = 23;// for SRR870667_2 //final
	// k = 17; // for SRR327342_1
	// k = 17; // for SRR327342_2
	// k = 17; // for SRR635193_1

	// --------
	// k = 25;// for ERR174310_1
	// k = 29;// for ERR174310
	// k = 17;// SRR554369_1 //final
	// k = 17;// SRR554369_2 //final
	// k = 17; //MH0001.081026 // final
	// k = 17;
	// k = 17; // SRR689233
	// k = 13;
	// k = 16;
	
	// fprintf(stderr, "k: %d\n", k);
	// fprintf(stderr, "diff_threshold: %d\n", diff_threshold);
	// fprintf(stderr, "cbthreshold: %d\n", cbthreshold);
	// fprintf(stderr, "first_mininum: %d\n", first_mininum);
	// fprintf(stderr, "thr_step: %d\n", thr_step);
	// fprintf(stderr, "maxthr: %d\n", maxthr);

	int w = readlen - k/2;
	// int k = 32;
	// int w = 100;
	// fprintf(stderr, "k: %d; w: %d\n", k, w);
	// for 100 length
	// int k = 31;
	// int w = 80, 
	int b = 14, f = 0;
 	
 	// if (k < b) {
 	// 	b = k;
 	// }
	
#ifdef _PE
	// pre_process(argv[1], argv[2], f, w, k, b, n_threads, readlen); //int w, int k, int b,
	// fprintf(stderr, "Two reads files: %s and %s\n", argv[1], argv[2]);
	
	folder = argv[3];
	pre_process(argv[1], argv[2], f, w, k, b, n_threads); //int w, int k, int b,
	cluster_dump_pe(n_threads, reads, idxv);
	// fprintf(stderr, "Compressed file: %s.minicom\n", argv[3]);
#else
	// fprintf(stderr, "%s %s\n", argv[1], argv[2]);
	folder = argv[2];
	// pre_process(argv[1], f, w, k, b, n_threads, readlen); //int w, int k, int b,
	pre_process(argv[1], f, w, k, b, n_threads); //int w, int k, int b,
	cluster_dump(n_threads, reads, idxv);
#endif

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


