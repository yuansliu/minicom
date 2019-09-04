#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "breads.h"
#include "kvec.h"

/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct {
	struct kt_for_t *t;
	long i;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	long n, index;
	ktf_worker_t *w;
	reads_t *reads;
	pthread_mutex_t *mutex;
	pthread_mutex_t *lock;
	pthread_mutex_t fp_mutex;
	pthread_mutex_t sg_mutex;
	FILE *single_fp;
	bool last_rounds;
	int kmer;
} kt_for_t;

static inline long steal_bucket_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

int cmpcluster(const void *a_, const void *b_) {
	// idpos_t _a = *(idpos_t*)a_, _b = *(idpos_t*)b_;
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;

	int rid_a = a>>32, rid_b = b>>32;
	int pos_a = (uint32_t)a>>1, pos_b = (uint32_t)b>>1;
	// int dir_a = a&1, dir_b = b&1;
	if (a&1) {
		pos_a = reads->seq_len - pos_a + reads->k - 2;
	}
	if (b&1) {
		pos_b = reads->seq_len - pos_b + reads->k - 2;
	}
	if (pos_a == pos_b) {
		return rid_a - rid_b;
	} 
	return pos_b - pos_a;
	// return strcmp(reads->seq[].seq, reads->seq[*(int*)b].seq);
}

const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

// int outputnum;
// for dgb debug "run SRR870667_1.fastq SRR870667_1_comp"
// void construct_ref(kt_for_t *t, reads_t *reads, cluster_t *p, int n, int index, bool last_rounds, FILE *fp) {
void construct_ref(kt_for_t *t, reads_t *reads, cluster_t *p, int n, int index, bool last_rounds) {
	// uint32_t **count_table;
	int readlen = reads->seq_len;
	int count_table_len = readlen<<1 + 1;
	// uint32_t **count_table = (uint32_t**)calloc(4, sizeof(uint32_t*));
	uint32_t *count_table[4];
	for (int i = 0; i < 4; ++i) {
		count_table[i] = (uint32_t*)calloc(count_table_len, sizeof(uint32_t));
		// memset(count_table[i], 0, count_table_len * sizeof(uint32_t));
	} 
	char *temp_str = (char*)calloc((readlen + 1), sizeof(char));
	int rid, pos, dir, pos_0 = readlen;

	// return; //for debug
// if (t->kmer != 12) {
	for (int k = 0; k < n; ++k) { // p->a[k]
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		dir = y&1;

		strcpy(temp_str, reads->seq[rid].seq);

		if (dir) {
			pos = readlen - pos + reads->k - 2;
			reverse_complement(temp_str, readlen);
		}
		if (k == 0) pos_0 = pos;

		for (int s = 0; s < readlen; ++s) {
			++count_table[seq_nt4_table[(uint8_t)temp_str[s]]][pos_0 - pos + s];
		}
		p->a[k] = (uint64_t)y >>32 << 32 | ((uint64_t)(pos_0 - pos) << 1) | dir;
	}
// }

	char *ref = (char*)calloc(count_table_len, sizeof(char));
	int max_count;
	for (int s = 0; s < count_table_len; ++s) {
		ref[s] = 'A';
		max_count = count_table[0][s];
		for (int k = 1; k < 4; ++k) {
			if (count_table[k][s] > max_count) {
				max_count = count_table[k][s];
				ref[s] = invert_code_rule[k];
			}
		}
		// fprintf(stderr, "%d,\n", max_count);
		if (max_count == 0) {
			ref[s] = '\0';
			break;
		}
	}
	// bool debug = false;
	// if (t->kmer == 12) debug = true;
	// if (debug) {
	// 	fprintf(stderr, "p->n: %d\n", p->n);
	// 	fprintf(stderr, "ref: %s\n", ref);
	// 	// for (int k = 0; k < p->n; ++k) { // p->a[k]
	// 	for (int k = 0; k < p->n; ++k) { // p->a[k]
	// 		uint64_t y = p->a[k];
	// 		rid = y>>32;
	// 		pos = (uint32_t)y>>1;
	// 		dir = y&1;

	// 		fprintf(stderr, "%d, %d : %s\n", pos, dir, reads->seq[rid].seq);
	// 	}
	// 	fprintf(stderr, "---\n");
	// 	// exit(0);
	// }
		
// if (t->kmer != 12) {
	mm128_t minimizer;
	// cluster_t *pp = (cluster_t*)calloc(1, sizeof(cluster_t));
	// kv_init(*pp);
	// kv_resize(uint64_t, *pp, n);
	
	// fprintf(stderr, "%s\n", ref);
	// kv_cluster_push_ref(p, ref);
	int ref_len = strlen(ref);
	p->ref = (char*)calloc(ref_len + 1, sizeof(char));
	strcpy(p->ref, ref);

	int sg_idx;

	//encode reads
	// pp->ennum = 0;
	// char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	// char *int_str = (char*)alloca(10 * sizeof(char));

	int mask = (1<<reads->b) - 1, lock;
	// int temp_str_len = readlen, sg_idx;
	// if (debug) {
	// 	fprintf(stderr, "---- %d %s\n", n, ref);
	// }
	p->n = 0;
	for (int k = 0; k < n; ++k) {
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		dir = y&1;

		strcpy(temp_str, reads->seq[rid].seq);

		if (dir) {
			reverse_complement(temp_str, readlen);
		}
		// if (debug && k < 5) fprintf(stderr, "%s\n", temp_str);

		// int en_str_len = 0;
		// int eq_char_num = 0;
		int dif_num = 0;
		for (int rj = pos, tj = 0; tj < readlen; ++rj, ++tj) {
			if (ref[rj] != temp_str[tj]) {
				++dif_num;
			} // else ++eq_char_num;
		}

		// if (debug) fprintf(stderr, "dif_num: %d; diff_threshold: %d\n", dif_num, diff_threshold);

		if (dif_num <= diff_threshold) {
			// kv_push(uint64_t, *pp, p->a[k]);
			kv_push(uint64_t, *p, p->a[k]);
			// if (debug) fprintf(stderr, "p->n: %d\n", p->n);
			// p->ennum += en_str_len;
		} else { // treat it as single reads;
			// p->a[k];
			rid = p->a[k] >> 32;

			if (last_rounds) {
				pthread_mutex_lock(&t->sg_mutex);
			 	kv_push(uint32_t, t->reads->sg, rid);
			 	// sg_idx = t->reads->sg.n - 1;
				pthread_mutex_unlock(&t->sg_mutex);
				//
			} else {
				mm_sketch_two(reads->seq[rid].seq, readlen, t->kmer, rid, &minimizer);
				// mm_sketch_one_ori(reads->seq[rid].seq, readlen, reads->k, rid, &minimizer);
				lock = minimizer.x&mask;
				mm128_v *p = &reads->B[index^1][lock];//mm_idx_bucket_t

				pthread_mutex_lock(&t->mutex[lock]);
				kv_push(mm128_t, *p, minimizer);
				pthread_mutex_unlock(&t->mutex[lock]);
			}
		}
		
	}
	// if (t->kmer == 12) {
	// 	fprintf(stderr, "555xxxx\n");
	// 	fprintf(stderr, "%d\n", p->n);
	// }
	// for (int i = 0; i < p->n; ++i) {
	// 	fprintf(stderr, "%lld, ", p->a[i]);
	// }
	// fprintf(stderr, "\n");
	// free(p->a);

	// if (t->kmer == 12){
	// 	for (int k = 0; k < 4; k++) {
	// 		for (int s = 0; s < 5; s++) {
	// 			fprintf(stderr, "%d,", count_table[k][s]);
	// 		}
	// 		fprintf(stderr, "--\n");
	// 	}
	// }

	// p->n = pp->n;
	// p->a = pp->a;
	// free(pp);
	// if (t->kmer == 12) fprintf(stderr, "666xxxx\n");
	// update reference of this cluster
	// if (debug) {
	// 	fprintf(stderr, "p->n after update: %d\n", p->n);
	// }
	if (p->n > 0) {
		// if (t->kmer == 12) fprintf(stderr, "111xxxx p->n > 0\n");

		for (int k = 0; k < n; ++k) { // p->a[k]
			uint64_t y = p->a[k];
			rid = y>>32;
			pos = (uint32_t)y>>1;
			dir = y&1;
			// if (rid == 6366077) debug = true;
		}

		for (int i = 0; i < 4; ++i) {
			memset(count_table[i], 0, count_table_len * sizeof(uint32_t));
		}
		// if (debug&&false) {
			// for (int i = 0; i < 4; i++) {
			// 	for (int j = 0; j < count_table_len; ++j) {
			// 		fprintf(stderr, "%d, ", count_table[i][j]);
			// 	} 
			// 	fprintf(stderr, "\n");
			// }
		// }

		int rend = 0;

		for (int i = 0; i < p->n; ++i) {
			uint64_t y = p->a[i];
			rid = y>>32;
			pos = (uint32_t)y>>1;
			dir = y&1;
			strcpy(temp_str, reads->seq[rid].seq);
			if (dir) {
				reverse_complement(temp_str, readlen);
			}
			
			// if (t->kmer == 29 && n <= 4) { fprintf(stderr, "%s %d %d\n", temp_str, pos, dir);} //for debug

			for (int s = 0; s < readlen; ++s) {
				++count_table[seq_nt4_table[(uint8_t)temp_str[s]]][pos + s];
			}

			if (pos + readlen > rend) {
				rend = pos + readlen;
			}
		}

		// if (t->kmer == 29 && n <= 4) {
		// 	outputnum++;
		// 	fprintf(stderr, "---\n");
		// } // for debug
		// if (outputnum >= 50) exit(0); // for debug

		// if (debug&&false) {	
			// for (int i = 0; i < 4; i++) {
			// 	for (int j = 0; j < 7; ++j) {
			// 		fprintf(stderr, "%d,", count_table[i][j]);
			// 	} 
			// 	fprintf(stderr, "\n");
			// }
		// }
		int s;
		for (s = 0; s < ref_len; ++s) {
			max_count = count_table[0][s];
			for (int k = 1; k < 4; ++k) {
				if (count_table[k][s] > max_count) {
					max_count = count_table[k][s];
				}
			}
			// fprintf(stderr, "%d,\n", max_count);
			if (max_count != 0) {
				break;
			}
		}
		int sv = s, r;
		// if (debug) fprintf(stderr, "ref_len: %d\n", ref_len);
		for (r = 0; s < rend; ++s, ++r) {
			p->ref[r] = 'A';
			max_count = count_table[0][s];
			for (int k = 1; k < 4; ++k) {
				if (count_table[k][s] > max_count) {
					max_count = count_table[k][s];
					p->ref[r] = invert_code_rule[k];
				}
			}
			/*if (debug&&false) {
				for (int k = 0; k < 4; ++k) {
					fprintf(stderr, "%d ", count_table[k][s]);
				}
				// fprintf(stderr, "");
				fprintf(stderr, "max_count: %d; s: %d\n", max_count, s);
			}*/
			// if (0 == max_count) {
				/*if (debug) fprintf(stderr, "s: %d\n", s);
				if (debug) fprintf(stderr, "r: %d\n", r);*/
				// break;
			// }
		}
		p->ref[r] = '\0';
		
		for (int k = 0; k < p->n; ++k) { // p->a[k]
			uint64_t y = p->a[k];
			rid = y>>32;
			pos = (uint32_t)y>>1;
			dir = y&1;

			p->a[k] = (uint64_t)y >>32 << 32 | ((uint64_t)(pos - sv) << 1) | dir;
		}
		// if (t->kmer == 12) fprintf(stderr, "222xxxx\n");
	}
	// if (t->kmer == 12) fprintf(stderr, "111YYYYY\n");
	// fprintf(fp, "%s %d\n", ref, n);
	free(temp_str);
// }
	// if (t->kmer == 12) fprintf(stderr, "222YYYYY\n");

	// if (t->kmer == 12){
	// 	for (int k = 0; k < 4; k++) {
	// 		for (int s = 0; s < 5; s++) {
	// 			fprintf(stderr, "%d,", count_table[k][s]);
	// 		}
	// 		fprintf(stderr, "\n");
	// 	}
	// 	// fprintf(stderr, "\n");
	// }

	for (int i = 0; i < 4; ++i) {
		// if (t->kmer == 12) fprintf(stderr, "i = %d\n", i);
		free(count_table[i]);
	}
	// free(count_table);
	// exit(0)
	free(ref);
	// if (t->kmer == 12) fprintf(stderr, "111ZZZZ\n");
}

// static void process_bucket(reads_t *reads, long i, int tid) {
// static void process_bucket(kt_for_t *t, long i, int tid, int index, bool last_rounds, FILE *fp) {
static void process_bucket(kt_for_t *t, long i, int tid, int index, bool last_rounds) {
	int j, start_a, n;
	// idxhash_t *h;
	// mm_idx_t *mi = (mm_idx_t*)g;
	cluster_bucket_t *b = &t->reads->B[index][i];
	if (b->n == 0) return;

	// if (t->kmer == 12) fprintf(stderr, "1111\n");
	mm128_t minimizer;
	// sort by minimizer
	radix_sort_128x(b->a, b->a + b->n);
	// create the hash table
	int rid, lock;
	// mm128_t minimizer;
	int mask = (1<<reads->b) - 1;
	// fprintf(stderr, "func process_bucket()\n");
 	// int sg_idx;
	for (j = 1, n = 1, start_a = 0; j <= b->n; ++j) {
		if (j == b->n || b->a[j].x != b->a[j-1].x) {
			assert(j - start_a == n);
			// if ((!last_rounds && n <= 5) || n < 2 ) {
			if (n < 2) {
				// __sync_fetch_and_add(&reads->single, n);
				// __sync_fetch_and_add(&reads->single, 1);
				// __sync_fetch_and_add(&clusize[n], 1);
				for (int s = 0; s < n; ++s) {
					// rid = b->a[j-1].y>>32;
					rid = b->a[start_a + s].y >> 32;
										
					pthread_mutex_lock(&t->sg_mutex);
				 	kv_push(uint32_t, t->reads->sg, rid);
				 	// sg_idx = t->reads->sg.n - 1;
					pthread_mutex_unlock(&t->sg_mutex);
					/*if (last_rounds) {
						pthread_mutex_lock(&t->sg_mutex);
					 	kv_push(uint32_t, t->reads->sg, rid);
					 	// sg_idx = t->reads->sg.n - 1;
						pthread_mutex_unlock(&t->sg_mutex);
						//
					} else {
						mm_sketch_two(reads->seq[rid].seq, readlen, t->kmer, rid, &minimizer);
						// mm_sketch_one_ori(reads->seq[rid].seq, readlen, reads->k, rid, &minimizer);
						lock = minimizer.x&mask;
						mm128_v *p = &reads->B[index^1][lock];//mm_idx_bucket_t

						pthread_mutex_lock(&t->mutex[lock]);
						kv_push(mm128_t, *p, minimizer);
						pthread_mutex_unlock(&t->mutex[lock]);
					}*/
				}
			} else {
				cluster_t *p;
				// pthread_mutex_lock(&t->mutex);
				kv_pushp(cluster_t, t->reads->clusters[0][tid], &p);
				kv_init(*p);

				kv_resize(uint64_t, *p, n);
				for (int k = 0; k < n; ++k) {
					// tmp_idpos.a = b->a[start_a + k].y;
					kv_push(uint64_t, *p, b->a[start_a + k].y);
				}
				qsort(p->a, n, sizeof(uint64_t), cmpcluster);
				//save to fp
				// construct_ref(t, t->reads, p, n, index, last_rounds, fp);
				// if (t->kmer == 12)
				construct_ref(t, t->reads, p, n, index, last_rounds);

				// if (t->kmer == 12) 
					// fprintf(stderr, "222222\n");

				if (p->n > 1) { // p->ref
					mm128_v mini;
					kv_init(mini);
					int len = strlen(p->ref);
					// fprintf(stderr, "%s---\n", p->ref);
					//set reads->k always < 32; the biggest value is 31.
					// if (reads->k < 32) {
					mm_sketch_lh_ori(p->ref, len, reads->rw, reads->k, ((t->reads->clusters[0][tid].n - 1)<<8) + tid, &mini);
					// } else {
						// mm_sketch_lh(p->ref, len, reads->rw, reads->k, ((t->reads->clusters[0][tid].n - 1)<<8) + tid, &mini);
					// }

					for (int k = 0; k < mini.n && k < first_mininum; ++k) {
					// for (int k = 0; k < mini.n && k < first_mininum; ++k) {
					// for (int k = 0; k < mini.n; ++k) { // ???? number of minimizers
						// fprintf(stderr, "(%lu, %lu)\n", mini.a[i].x, mini.a[i].y>>1);
						// fprintf(stderr, "(%lu, %lu)\n", mini.a[k].x&mask, mini.a[k].y>>1);
						lock = mini.a[k].x&mask;
						mm128_v *pp = &reads->mi[0]->B[lock].a;

						pthread_mutex_lock(&t->lock[lock]);
						kv_push(mm128_t, *pp, mini.a[k]);
						pthread_mutex_unlock(&t->lock[lock]);
					}
					kv_destroy(mini);
				} else {
					if (p->n == 1) {
						// only one reads in a cluster  // this before version Sunday Aug. 5 11:00 AM
						// pthread_mutex_lock(&t->sg_mutex);
					 // 	kv_push(uint32_t, t->reads->sg, p->a[0] >> 32);
						// pthread_mutex_unlock(&t->sg_mutex);
						if (last_rounds) {
							pthread_mutex_lock(&t->sg_mutex);
						 	kv_push(uint32_t, t->reads->sg, p->a[0] >> 32);
							pthread_mutex_unlock(&t->sg_mutex);
							//
						} else {
							rid = p->a[0] >> 32;
							mm_sketch_two(reads->seq[rid].seq, reads->seq_len, t->kmer, rid, &minimizer);
							// mm_sketch_one_ori(reads->seq[rid].seq, readlen, reads->k, rid, &minimizer);
							lock = minimizer.x&mask;
							mm128_v *p = &reads->B[index^1][lock];//mm_idx_bucket_t

							pthread_mutex_lock(&t->mutex[lock]);
							kv_push(mm128_t, *p, minimizer);
							pthread_mutex_unlock(&t->mutex[lock]);
						}
					}
					p->n = 0;
					--t->reads->clusters[0][tid].n;
				} 
			}
			start_a = j, n = 1;
		} else {++n;}
	}
	// deallocate and clear b->a
	free(b->a);
	b->n = b->m = 0, b->a = 0;
}

static void *ktf_bucket_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	kt_for_t *t = w->t;
	// char name[100]; sprintf(name, "encode_%ld.txt", w - w->t->w);
	// FILE *fp = fopen(name, "w");
	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, t->n_threads);
		if (i >= t->n) break;
		// w->t->func(w->t->data, i, w - w->t->w); 
		process_bucket(t, i, w - w->t->w, t->index, t->last_rounds);
	}
	while ((i = steal_bucket_work(t)) >= 0)
		// w->t->func(w->t->data, i, w - w->t->w);
		process_bucket(t, i, w - w->t->w, t->index, t->last_rounds);
	// fclose(fp);
	pthread_exit(0);
}

void kt_for_bucket_index(int n_threads, reads_t *reads, long n, int index, int kmer, bool last_rounds)
{
	int i;
	kt_for_t t;
	pthread_t *tid;
	t.reads = reads, t.n_threads = n_threads, t.n = n, t.index = index, t.last_rounds = last_rounds, t.kmer = kmer;
	
	// fprintf(stderr, "kmer: %d\n", t.kmer);

	t.w = (ktf_worker_t*)alloca(n_threads * sizeof(ktf_worker_t));
	// pthread_mutex_init(&t.mutex, 0);
	int num_mutex = (1<<reads->b);
	t.mutex = (pthread_mutex_t*)alloca(num_mutex * sizeof(pthread_mutex_t));
	for (int i = 0; i < num_mutex; ++i) {
		pthread_mutex_init(&t.mutex[i], 0);
	}
	//for single reads
	t.lock = (pthread_mutex_t*)alloca(num_mutex * sizeof(pthread_mutex_t));
	for (int i = 0; i < num_mutex; ++i) {
		pthread_mutex_init(&t.lock[i], 0);
	}
	pthread_mutex_init(&t.fp_mutex, 0);
	pthread_mutex_init(&t.sg_mutex, 0);

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = i; //
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_bucket_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}

void kt_for_bucket(int n_threads, reads_t *reads, long n) {
	// fprintf(stderr, "begin kt_for_bucket(int n_threads, reads_t *reads, long n)\n");

	// fprintf(stderr, "reads->k: %d; readlen: %d\n", reads->k, reads->seq_len);

	int index = 0;
	// , max_rounds = 10;
	// int index = 0, max_rounds = 6;
	// int index = 0, max_rounds = 1;
	// fprintf(stderr, "max_rounds: %d\n", max_rounds);
	int pre_cluster_reads = 0;
	int last_rounds = 0;

	int max_bucket = 1 << reads->b;
	// for (int k = 1; k <= max_rounds; ++k) {
	// max_rounds = 5;
	for (int k = 1; ; k++) {
		reads->single = 0;
		for (int i = 0; i < max_bucket; i++) {
			kv_init(reads->B[index^1][i]);
		}
		
		if (reads->k - k <= 9) ++last_rounds;
		if (k == max_rounds - 1) ++last_rounds;
		// kt_for_bucket_index(n_threads, reads, n, index, reads->k - (2*k - 1), !(k < max_rounds));
		// fprintf(stderr, "k = %d\n", k);
		// fprintf(stderr, "reads->k - k: %d\n", reads->k - k);
		// if (reads->k - k == 12) 
			// n_threads = 1;
		// if (k == 2) {n_threads = 1; outputnum = 0;}
		kt_for_bucket_index(n_threads, reads, n, index, reads->k - k, last_rounds);

		if (last_rounds) ++last_rounds;

		for (int i = 0; i < max_bucket; i++) {
			kv_destroy(reads->B[index][i]);
		}
		index ^= 1;

		// do {	
			int tot = 0;
			for (int k = 0; k < n_threads; ++k) tot += reads->clusters[0][k].n;
			// fprintf(stderr, "cluster number: %d\n", tot);
			int cluster_reads = 0;
		
			for (int k = 0; k < n_threads; ++k) {
				for (int i = 0; i < reads->clusters[0][k].n; ++i) {
					cluster_reads += reads->clusters[0][k].a[i].n;
					
				}
			}
			// fprintf(stderr, "cluster contains reads number: %d\n", cluster_reads);
			if (cluster_reads - pre_cluster_reads < 100) {
				// break;
				++last_rounds;
				/* code */
			}
			pre_cluster_reads = cluster_reads;
		// } while (0);
		// exit(0);
		if (last_rounds > 1) break;
	}

	for (int i = 0; i < max_bucket; i++) {
		kv_destroy(reads->B[index][i]);
	}
	// fprintf(stderr, "over kt_for_bucket(int n_threads, reads_t *reads, long n)\n");
}