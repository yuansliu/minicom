#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "kvec.h"
#include "breads.h"
#include "khash.h"

/******************
 * Pre-process *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"
#include "config.h"

void kt_for_bucket(int n_threads, reads_t *reads, long n);
void kt_for_reads(int n_threads, reads_t *reads, long n);

// int clusize[20];

void reverse_complement(char *seq, int seq_len) {
	// char rc_tablle[256];
	char *rc_tablle = (char*)alloca(256 * sizeof(char));
	rc_tablle['A'] = 'T';
	rc_tablle['T'] = 'A';
	rc_tablle['G'] = 'C';
	rc_tablle['C'] = 'G';
	rc_tablle['N'] = 'N';
	// uint32_v *n_pos = (uint32_v*)calloc(1, sizeof(uint32_v));
	char *res = (char*)alloca((seq_len + 1) * sizeof(char));
	for (int i = seq_len - 1, j = 0; i >= 0; --i, ++j) {
		res[j] = rc_tablle[(uint8_t)seq[i]];
	}
	res[seq_len] = '\0';
	strcpy(seq, res);
}

#ifdef _PE
void pre_process(const char *fn, const char *fn0, int f, int w, int k, int b, int n_threads) {
#else
void pre_process(const char *fn, int f, int w, int k, int b, int n_threads) {
#endif
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return;
	reads = (reads_t*)calloc(1, sizeof(reads_t));
	memset(reads, 0, sizeof(reads_t)); 

	reads->seq_len = readlen;
	
	#ifdef _PE
	fprintf(stdout, "Loading the first file %s ...\n", fn);
	#else
	fprintf(stdout, "Loading the file %s ...\n", fn);
	#endif
	reads->seq = bseq_read(fp, &reads->n_seq, readlen);
	bseq_close(fp);

	#ifdef _PE
	int pre_n_seq = reads->n_seq;
	fp = bseq_open(fn0);
	if (fp == 0) return;
	fprintf(stdout, "Loading the second file %s ...\n", fn0);

	bseq_read_second(reads->seq, fp, &reads->n_seq, readlen);
	bseq_close(fp);
	half_val = reads->n_seq - pre_n_seq;

	if ((half_val<<1) != reads->n_seq) {
		fprintf(stdout, "%s and %s contain different number of reads.\n\n", fn, fn0);
		// exit(1);
	}
	#endif

	// fprintf(stderr, "reads->n_seq: %d\n", reads->n_seq);
	
	// reads->flag = (uint8_t*)calloc(reads->n_seq, sizeof(uint8_t));
	// memset(reads->flag, 0, reads->n_seq * sizeof(uint8_t));

	// reads->flag = (uint8_t*)calloc(reads->n_seq, sizeof(uint8_t));
	// memset(reads->flag, 0, reads->n_seq * sizeof(uint8_t));

	// fprintf(stderr, "finish bseq_read() ...\n");
 
	// n_threads = 1;
	reads->f = f; reads->k = k; reads->b = b;
	// reads->w = w; 
	if (readlen >= 70) {
		reads->rw = readlen/2 - reads->k;
		// if (reads->rw >= 5) reads->rw -= 5;

		/*if (readlen == 108) { // for SRR870667_1
			reads->rw = 45 - reads->k;
			// reads->rw = 15;
		}*/
		// reads->rw = 5;
		// if (reads->rw <= 0) reads->rw = 1;
	} else { 
		reads->rw = 3;
	// reads->rw = 3; //small rw can align more reads to the contigs as in contig merging stage having a more strict rule.
	}
	// if (readlen == 44) reads->rw = 3;

	if (rw > 0) {
		reads->rw = rw;
	}

	// fprintf(stderr, "reads->rw (reads window): %d\n", reads->rw);
	// reads->sk = 5;

	fprintf(stdout, "Contigs generation ...\n");

	reads->B = (cluster_bucket_t**)calloc(2, sizeof(cluster_bucket_t*));
	reads->B[0] = (cluster_bucket_t*)calloc(1<<b, sizeof(cluster_bucket_t));
	reads->B[1] = (cluster_bucket_t*)calloc(1<<b, sizeof(cluster_bucket_t));
	for (int i = 0; i < (1<<b); i++) {
		kv_init(reads->B[0][i]);
	}
	// kv_init(reads->B[0]);

	/*reads->B = (cluster_bucket_t*)calloc(1<<b, sizeof(cluster_bucket_t));
	kv_init(*reads->B);
	reads->B0 = (cluster_bucket_t*)calloc(1<<b, sizeof(cluster_bucket_t));
	kv_init(*reads->B0);*/

	reads->sp = (sp_reads_t*)calloc(1, sizeof(sp_reads_t));
	memset(reads->sp, 0, sizeof(sp_reads_t));

	kv_init(reads->sg);
	kv_init(reads->fpA_id);
	kv_init(reads->fpT_id);
	kv_init(reads->fpN_id);
	kv_init(reads->Nfile_id);
	kv_init(reads->singleFile_id);

	double mm_realtime1 = realtime();

	// n_threads = 1;
	// fprintf(stderr, "before kt_for_reads()...\n");
	kt_for_reads(n_threads, reads, reads->n_seq);
	// fprintf(stderr, "after kt_for_reads()...\n");

	reads->mi = (mm_idx_t**)calloc(2, sizeof(mm_idx_t*));
	reads->mi[0] = mm_idx_init(reads->b);
	// fprintf(stderr, "after init(f, w, k, b)...\n");

	// kt_sort(n_threads, reads_sort_pipeline, reads);
	// kt_sort(n_threads, reads);
	// fprintf(stderr, "after reads_sort_pipeline()...\n");
	//-------------

	// n_threads = 1;
	reads->clusters = (cluster_v**)calloc(2, sizeof(cluster_v*));
	for (int i = 0; i < 2; ++i) {
		reads->clusters[i] = (cluster_v*)calloc(n_threads, sizeof(cluster_v));
	}
	for (int j = 0; j < n_threads; ++j) {
		kv_init(reads->clusters[0][j]);
		kv_resize(cluster_t, reads->clusters[0][j], 1 << 10);
	}
	

	reads->single = 0;
	// n_threads = 1;
	kt_for_bucket(n_threads, reads, 1<<reads->b);

	// fprintf(stderr, "sg size: %lu\n", reads->sg.n);
	if (reads->sg.n <= 5000000) {
		maxmatch = readlen*2/3;
		maxsearch = 2000;
	}

	idxv = 0;
	
	// fprintf(stderr, "n_threads: %d\n", n_threads);

	combine_cluster(n_threads, reads, &idxv);
	// 
	// fprintf(stderr, "after combine_cluster()... idxv: %d\n", idxv);

	reads->sg_flag = (bool*)calloc(reads->sg.n, sizeof(bool));
	memset(reads->sg_flag, 0, reads->sg.n * sizeof(bool));

	double mm_realtime2 = realtime();
	fprintf(stdout, "[Stage 1] Real time: %.3f sec;\n", mm_realtime2 - mm_realtime1);
	
	fprintf(stdout, "Realignment ... (please wait, it requires more time than stage 1)\n");

	// fprintf(stderr, "------------------------------------\n");
	int pre_reads_number_in_cluster = 0;
	bool threshold_flag = true;
	// threshold_flag = false;
	// int step = diff_threshold;
	// if (step >= 8) step = 5;
 	// for (int threshold = thresh; threshold_flag && threshold <= 84; threshold += 5) {
 	for (int threshold = diff_threshold; threshold_flag ; threshold += thr_step) {
 	// for (int threshold = thresh; threshold_flag && threshold <= readlen*2/3; threshold += 20) {
 	// for (int threshold = thresh; threshold_flag && threshold <= thresh; threshold += 20) {
 	// for (int threshold = thresh; threshold <= 34; threshold += 5) {
 		// if (threshold > 12) break;
		if (threshold > maxthr) break;
	 	updateSingle();
	 	realign_hash(n_threads, reads, idxv, threshold);
		do {	
			int tot = 0;
			for (int k = 0; k < n_threads; ++k) tot += reads->clusters[idxv][k].n;
			// fprintf(stderr, "cluster number: %d\n", tot);
			int cluster_reads = 0;
		
			for (int k = 0; k < n_threads; ++k) {
				for (int i = 0; i < reads->clusters[idxv][k].n; ++i) {
					cluster_reads += reads->clusters[idxv][k].a[i].n;
					
				}
			}
			// fprintf(stderr, "cluster contains reads number: %d\n", cluster_reads);
			// fprintf(stderr, "cluster_reads - pre_reads_number_in_cluster: %d\n", cluster_reads - pre_reads_number_in_cluster);
			int max_thr = 1000;
			if (reads->sg.n > 1000000 && readlen >= 68) {
				max_thr = 10000;
			}

			if (cluster_reads - pre_reads_number_in_cluster < max_thr) {
				threshold_flag = false;
				break;
			}
			pre_reads_number_in_cluster = cluster_reads;
		} while (0);
		// fprintf(stderr, "threshold = %d\n", threshold);
		// if (threshold >= 12) break;
	}
	updateSingle();
	double mm_realtime3 = realtime();
	fprintf(stdout, "[Stage 2] Real time: %.3f sec;\n", mm_realtime3 - mm_realtime2);
	
	// fprintf(stdout, "remaining reads: %d\n", reads->sg.n);

	// cluster_dump(n_threads, reads, idxv);
	// cluster_dump(1, reads, idxv);
}

void updateSingle() {
	int nn = 0;
	for (int i = 0; i < reads->sg.n; ++i) {
		if (!reads->sg_flag[i]) {
			reads->sg.a[nn++] = reads->sg.a[i];
		}
	}
	reads->sg.n = nn;
	// update sg_flag
	free(reads->sg_flag);
	reads->sg_flag = (bool*)calloc(reads->sg.n, sizeof(bool));
	memset(reads->sg_flag, 0, reads->sg.n * sizeof(bool));
}
