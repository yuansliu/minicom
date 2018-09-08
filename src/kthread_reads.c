#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "breads.h"
#include "kvec.h"

/************
 * kt_reads_for() *
 ************/

struct kt_reads_for_t;

typedef struct {
	struct kt_reads_for_t *t;
	long i;
} ktf_reads_worker_t;

typedef struct kt_reads_for_t {
	int n_threads;
	long n;
	ktf_reads_worker_t *w;
	reads_t *reads;
	pthread_mutex_t *mutex, sg_mutex;
	pthread_mutex_t mutexAllA, mutexAllT, mutexAllN, AAfile, TTfile, NNfile;
	pthread_mutex_t mutexSingleFile;
} kt_reads_for_t;

static inline long steal_reads_work(kt_reads_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void process_reads(kt_reads_for_t *t, uint32_t rid) {
	mm128_t minimizer;
	int mask = (1<<reads->b) - 1, lock;

	bseq1_t *seq = &t->reads->seq[rid];
	char *cur_seq = seq->seq;
	sp_reads_t *s = t->reads->sp;
	// fprintf(stderr, "rid: %ld\n", rid);
	// if (rid == 1662) fprintf(stderr, "%s\n", cur_seq);
	uint32_v *n_pos = (uint32_v*)calloc(1, sizeof(uint32_v));
	int seq_len = reads->seq_len;

	char *en_str = (char*)alloca((seq_len + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));

	int cntA = 0, cntT = 0, cntG = 0, cntC = 0, cntN = 0;
	for (int i = 0; i < seq_len; ++i) {
		if (cur_seq[i] == 'A') {
			cntA++;
		} else 
		if (cur_seq[i] == 'T') {
			cntT++;
		} else
		if (cur_seq[i] == 'G') {
			cntG++;
		} else
		if (cur_seq[i] == 'C') {
			cntC++;
		} else
		if (cur_seq[i] == 'N') {
			cntN++;
			kv_push(uint32_t, *n_pos, i);
		}
	}
	// if (rid == 1662) fprintf(stderr, "%d, %d, %d, %d\n", cntA, cntT, cntG, cntC);
	if (cntN > 0) {
		seq->n_pos = n_pos;
	} else {
		seq->n_pos = NULL;
		free(n_pos);
	}
	// int threshold = seq_len + diff_threshold;
	// int threshold = seq_len - diff_threshold;
	// int threshold = seq_len + 100;
	if (cntA == seq_len) {
		// fprintf(stderr, "allA\n");
		pthread_mutex_lock(&t->mutexAllA);
		s->allA++;
		kv_push(uint32_t, s->allA_id, rid);
		pthread_mutex_unlock(&t->mutexAllA);

		// pthread_mutex_lock(&atfiles_mutex);
		// fprintf(atfiles, "%s\n", cur_seq);
		// pthread_mutex_unlock(&atfiles_mutex);
	} else 
	if (cntT == seq_len) {
		// fprintf(stderr, "allT\n");
		pthread_mutex_lock(&t->mutexAllT);
		s->allT++;
		kv_push(uint32_t, s->allT_id, rid);
		pthread_mutex_unlock(&t->mutexAllT);

		// pthread_mutex_lock(&atfiles_mutex);
		// fprintf(atfiles, "%s\n", cur_seq);
		// pthread_mutex_unlock(&atfiles_mutex);
	} else
	if (cntN == seq_len) {
		pthread_mutex_lock(&t->mutexAllN);
		s->allN++;
		kv_push(uint32_t, s->allN_id, rid);
		pthread_mutex_unlock(&t->mutexAllN);

	} else
	if (cntT + cntG + cntC + cntN <= diff_threshold) {
		pthread_mutex_lock(&t->AAfile);
		kv_push(uint32_t, reads->fpA_id, rid); // most of bases are A
		pthread_mutex_unlock(&t->AAfile);
	} else 
	if (cntA + cntG + cntC + cntN <= diff_threshold) {
		pthread_mutex_lock(&t->TTfile);
		kv_push(uint32_t, reads->fpT_id, rid); // T
		pthread_mutex_unlock(&t->TTfile);
	} else 
	if (cntA + cntT + cntG + cntC <= diff_threshold) { //most of them are N
		pthread_mutex_lock(&t->NNfile);
		kv_push(uint32_t, reads->fpN_id, rid);	// N	
		pthread_mutex_unlock(&t->NNfile);
	} else { //not all A T G C
			// fprintf(stderr, "not all A T\n");
			
			// fprintf(stderr, "begin mm_sketch_mid()...\n");
			// mm_sketch_one(cur_seq, readlen, reads->k, rid, &minimizer);
		/*if (cntN > 0) {		
			if (cntN <= 0.4 * readlen) {
				if (cntN > 0) {
					// return;
					int max_value = getmax4(cntA, cntT, cntG, cntC);
					char replace = 'A';
					if (max_value == cntA) {
						cntA += cntN;
					} else 
					if (max_value == cntT) {
						cntT += cntN;
						replace = 'T';
					} else 
					if (max_value == cntG) {
						cntG += cntN;
						replace = 'G';
					} else 
					if (max_value == cntC) {
						cntC += cntN;
						replace = 'C';
					}
					for (int i = 0; i < cntN; i++) {
						cur_seq[n_pos->a[i]] = replace;
					} 
				} 
				pthread_mutex_lock(&t->sg_mutex);
			 	kv_push(uint32_t, t->reads->sg, rid);
			 	// sg_idx = t->reads->sg.n - 1;
				pthread_mutex_unlock(&t->sg_mutex);
			} else {
				pthread_mutex_lock(&t->mutexSingleFile);
				// fprintf(reads->Nfile, "%s\n", cur_seq);
				// fprintf(reads->Nfile_id, "%s\n", cur_seq);
				kv_push(uint32_t, reads->Nfile_id, rid);
				pthread_mutex_unlock(&t->mutexSingleFile);
			}
		} else {
			if (reads->k < 32) { 
				mm_sketch_two(cur_seq, readlen, reads->k, rid, &minimizer);
			} else {
				mm_sketch_one_ori(cur_seq, readlen, reads->k, rid, &minimizer);
			} 
			// fprintf(stderr, "after mm_sketch_mid()...\n");
			lock = minimizer.x&mask;
			mm128_v *p = &reads->B[0][lock];//mm_idx_bucket_t

			pthread_mutex_lock(&t->mutex[lock]);
			kv_push(mm128_t, *p, minimizer);
			pthread_mutex_unlock(&t->mutex[lock]);
		}*/ 
		if (cntN <= 0.4 * seq_len) {
			if (cntN > 0) {
				// return;
				int max_value = getmax4(cntA, cntT, cntG, cntC);
				char replace = 'A';
				if (max_value == cntA) {
					cntA += cntN;
				} else 
				if (max_value == cntT) {
					cntT += cntN;
					replace = 'T';
				} else 
				if (max_value == cntG) {
					cntG += cntN;
					replace = 'G';
				} else 
				if (max_value == cntC) {
					cntC += cntN;
					replace = 'C';
				}
				for (int i = 0; i < cntN; i++) {
					cur_seq[n_pos->a[i]] = replace;
				} 
			} 
			//set reads->k always < 32; the biggest value is 31.
			// if (reads->k < 32) { //set reads->k always < 32; the biggest value is 31.
			mm_sketch_two(cur_seq, seq_len, reads->k, rid, &minimizer);
			// } else {
				// mm_sketch_one_ori(cur_seq, readlen, reads->k, rid, &minimizer);
			// } 
			// fprintf(stderr, "after mm_sketch_mid()...\n");
			lock = minimizer.x&mask;
			mm128_v *p = &reads->B[0][lock];//mm_idx_bucket_t

			pthread_mutex_lock(&t->mutex[lock]);
			kv_push(mm128_t, *p, minimizer);
			pthread_mutex_unlock(&t->mutex[lock]);
		} else {
			pthread_mutex_lock(&t->mutexSingleFile);
			// fprintf(reads->Nfile, "%s\n", cur_seq);
			// fprintf(reads->Nfile_id, "%s\n", cur_seq);
			kv_push(uint32_t, reads->Nfile_id, rid);
			pthread_mutex_unlock(&t->mutexSingleFile);
		}
			
	}
	// fprintf(stderr, "rid: %ld\n", rid);
	// fprintf(stderr, "over process_reads()...\n");
}

static void *ktf_reads_worker(void *data) {
	ktf_reads_worker_t *w = (ktf_reads_worker_t*)data;
	kt_reads_for_t *t = w->t;
	long i;

	for (;;) {
		i = __sync_fetch_and_add(&w->i, t->n_threads);
		if (i >= t->n) break;
		process_reads(t, i);
	}
	while ((i = steal_reads_work(t)) >= 0)
		process_reads(t, i);
	pthread_exit(0);
}

void kt_for_reads(int n_threads, reads_t *reads, long n) {

	int i;
	kt_reads_for_t t;
	pthread_t *tid;
	t.reads = reads, t.n_threads = n_threads, t.n = n;
	t.w = (ktf_reads_worker_t*)alloca(n_threads * sizeof(ktf_reads_worker_t));
	
	int num_mutex = (1<<reads->b);
	t.mutex = (pthread_mutex_t*)alloca(num_mutex * sizeof(pthread_mutex_t));
	for (int i = 0; i < num_mutex; ++i) {
		pthread_mutex_init(&t.mutex[i], 0);
	}
	pthread_mutex_init(&t.mutexAllA, 0);
	pthread_mutex_init(&t.mutexAllT, 0);
	pthread_mutex_init(&t.mutexAllN, 0);
	pthread_mutex_init(&t.mutexSingleFile, 0);
	pthread_mutex_init(&t.AAfile, 0);
	pthread_mutex_init(&t.TTfile, 0);
	pthread_mutex_init(&t.NNfile, 0);
	pthread_mutex_init(&t.sg_mutex, 0);

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = i; //
	// for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_reads_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);

}
