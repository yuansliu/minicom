#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "breads.h"
#include "kvec.h"
#include "config.h"

/************
 * kt_for() *
 ************/

/*typedef struct {
	int tid, n_threads;
	reads_t *reads;
} ktf_find_next_worker_t;*/

struct kt_cb_for_t;

typedef struct {
	struct kt_cb_for_t *t;
	long i, n; //i < n;
	int tid;
} ktf_cb_worker_t;

typedef struct kt_cb_for_t {
	int n_threads, index, win, win_step;//win is the length of window
	ktf_cb_worker_t *w;
	reads_t *reads;
	bool **flag;
	pthread_mutex_t **mutex;
	pthread_mutex_t *lock;
} kt_cb_for_t;

int match_pro(char *str0, char *str1, int i_, int j_) {
	int tot = 0, match = 0;
	int i = i_, j = j_;
	while (str0[i] != '\0' && str1[j] != '\0') {
		++tot;
		if (str0[i] == str1[j]) ++match;
		++i, ++j;
	}
	i = i_ - 1, j = j_ - 1;
	while (i >= 0 && j >= 0) {
		++tot;
		if (str0[i] == str1[j]) ++match;
		--i, --j;
	}
	// return 1.0*match/tot;
	return tot - match;
}

int cmpcluster2(const void *a_, const void *b_) {
	// idpos_t _a = *(idpos_t*)a_, _b = *(idpos_t*)b_;
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;

	// int rid_a = a>>32, rid_b = b>>32;
	int pos_a = (uint32_t)a>>1, pos_b = (uint32_t)b>>1;
	int dir_a = a&1, dir_b = b&1;

	if (pos_a == pos_b) {
		// return rid_a - rid_b; // for order-preserving mode
		return dir_a - dir_b;
	}
	// return pos_b - pos_a;
	return pos_a - pos_b;
	// return strcmp(reads->seq[].seq, reads->seq[*(int*)b].seq);
}

#ifdef ORDER 
int cmpcluster3(const void *a_, const void *b_) {
	// idpos_t _a = *(idpos_t*)a_, _b = *(idpos_t*)b_;
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;

	int rid_a = a>>32, rid_b = b>>32;
	int pos_a = (uint32_t)a>>1, pos_b = (uint32_t)b>>1;
	// int dir_a = a&1, dir_b = b&1;

	if (pos_a == pos_b) {
		return rid_a - rid_b; // for order-preserving mode
	}
	return pos_a - pos_b;
	// return strcmp(reads->seq[].seq, reads->seq[*(int*)b].seq);
}
#endif

#ifdef _PE
int cmpcluster3(const void *a_, const void *b_) {
	// idpos_t _a = *(idpos_t*)a_, _b = *(idpos_t*)b_;
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;

	int rid_a = a>>32, rid_b = b>>32;
	int pos_a = (uint32_t)a>>1, pos_b = (uint32_t)b>>1;
	// int dir_a = a&1, dir_b = b&1;

	if (pos_a == pos_b) {
		return rid_a - rid_b; // for order-preserving mode
	}
	return pos_a - pos_b;
	// return strcmp(reads->seq[].seq, reads->seq[*(int*)b].seq);
}
#endif

void construct_ref2(reads_t *reads, cluster_t *p) {

	qsort(p->a, p->n, sizeof(uint64_t), cmpcluster2);
	
	int tot_len = (((uint32_t)p->a[p->n - 1] >> 1)) + (readlen << 1);

	// fprintf(stderr, "tot_len: %d\n", tot_len);
	uint32_t *count_table[4];
	int count_table_len = tot_len + 1;
	for (int i = 0; i < 4; ++i) {
		count_table[i] = (uint32_t*)calloc(count_table_len, sizeof(uint32_t));
		// memset(count_table[i], 0, count_table_len * sizeof(uint32_t));
	}
	char *temp_str = (char*)calloc((readlen + 1), sizeof(char));
	int rid, pos, dir;

	int rend = 0;
	for (int k = 0; k < p->n; ++k) { // p->a[k]
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		dir = y&1;
		strcpy(temp_str, reads->seq[rid].seq);
		if (dir) {
			reverse_complement(temp_str, readlen);
		}
		for (int s = 0; s < readlen; ++s) {
			++count_table[seq_nt4_table[(uint8_t)temp_str[s]]][pos + s];
		}
		if (pos + readlen > rend) {
			rend = pos + readlen;
		}
	}
	char *ref = (char*)calloc(count_table_len, sizeof(char));
	int max_count;
	for (int s = 0; s < rend; ++s) {
		ref[s] = 'A';
		max_count = count_table[0][s];
		for (int k = 1; k < 4; ++k) {
			if (count_table[k][s] > max_count) {
				max_count = count_table[k][s];
				ref[s] = invert_code_rule[k];
			}
		}
		// fprintf(stderr, "%d,\n", max_count);
		// if (max_count == 0) {
		// 	ref[s] = '\0';
		// 	break;
		// }
	}
	ref[rend] = '\0';
	// fprintf(stderr, "after count_table()...\n");
	p->ref = (char*)calloc(strlen(ref) + 1, sizeof(char));
	strcpy(p->ref, ref);

	p->ennum = 0;
	//encode reads
	char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));

	// int pre_pos = 0;
	// fprintf(fp, "%d %s\n", n, ref);
	for (int k = 0; k < p->n; ++k) {
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		dir = y&1;
		strcpy(temp_str, reads->seq[rid].seq);
		if (dir) {
			reverse_complement(temp_str, readlen);
		}
		// if (rid == 6366077) fprintf(stderr, "rid == 6366077 in kthread_cb.c\n");
		// fprintf(fp, "%s\n", temp_str);
		int en_str_len = 0;
		int eq_char_num = 0;
		int dif_num = 0;
		for (int rj = pos, tj = 0; tj < readlen; ++rj, ++tj) {
			if (ref[rj] != temp_str[tj]) {
				if (eq_char_num > 0) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				}
				en_str[en_str_len++] = temp_str[tj];
				++dif_num;
			} else ++eq_char_num;
		}
		if (en_str_len == 0) {
			en_str[en_str_len++] = '0';
		}
		en_str[en_str_len] = '\0';
		p->ennum +=  en_str_len;

		// if (dif_num > 5) {

		// }
		// if (en_str_len > 10) {
		// 	fprintf(fp, "%s\n", temp_str[k]);
		// }
		// fprintf(fp, "%d %d %s\n", pos[getmax(k-1, 0)] - pos[k], dir[k], en_str);
		// fprintf(fp, "%d %d %s\n", (int)((uint32_t)p->a[k]>>1) - pre_pos, (int)p->a[k]&1, en_str);
		// pre_pos = (uint32_t)p->a[k]>>1;
		//pos[0] - pos[k];
		// fprintf(fp, "%s\n", en_str);
	}

	free(temp_str);
	for (int i = 0; i < 4; ++i) {
		free(count_table[i]);
	}
	free(ref);
}

static void find_next(kt_cb_for_t *t, int tid_, int i_, int tid_ori) {
	// fprintf(stderr, "%d, %d\n", tid_, i_);
	//t->index
	mm128_v mini;
	kv_init(mini);
	cluster_t *p = &reads->clusters[t->index][tid_].a[i_];
	
	int len = strlen(p->ref);

	// if (tid_ == 0 && i_ == 34972)
	// 	fprintf(stderr, "ref: %s\n", p->ref);

	uint32_t rid_ori = (i_<<8) + tid_;
	// if (reads->k < 32) {
	mm_sketch_lh_ori(p->ref, len, t->win, reads->k, rid_ori, &mini);
	// } else {
	// 	mm_sketch_lh(p->ref, len, t->win, reads->k, rid_ori, &mini);
	// }

	bool flagEncode = false;
	uint64_t y;
	int n, lock;
	// mm128_t minimizer;
	int mask = (1<<reads->b) - 1;

	// if (tid_ == 0 && i_ == 34972)
	// 	fprintf(stderr, "n: %d\n", n);

		
	bool debug = false;
	for (int ii = 0; ii < p->n; ++ii) {
		uint64_t y = p->a[ii];
		int rid = y>>32;
		int pos = (uint32_t)y>>1;
		int dir = y&1;
		// strcpy(temp_str, reads->seq[rid].seq);
		// if (dir) {
		// 	reverse_complement(temp_str, readlen);
		// }
		// fprintf(stderr, "%d, %d, %d\n", rid, pos, dir);
		// fprintf(stderr, "%s\n", temp_str);
		// if (rid == 6366077) debug = true;
		// if (strcmp(reads->seq[rid].seq, "TTGGGCAGAAAATGGACCAAGGAAGAATTGTCAAAAATACAGCTATGATGAGAGAAGCTGCAAGAAAAATAGAAGAAAGGCATTTTTCCAACCATGCAAC") == 0) {
		// 	debug = true;
		// }
	}

	for (int j = 0; j < mini.n && !flagEncode; ++j) {
		const uint64_t *r;
		r = mm_idx_get(reads->mi[t->index], mini.a[j].x, &n);

		uint32_t pos_ori = (uint32_t) mini.a[j].y >> 1;
		uint32_t dir_ori = mini.a[j].y&1;

		for (int k = 0; k < n && !flagEncode; ++k) {
			uint32_t rid = (uint32_t) (r[k] >> 32);
			// if (r[k] != mini.a[j].y && rid != rid_ori) { // 
			// if (rid != (uint32_t)(mini.a[j].y >> 32)) { // 
			if (rid != rid_ori) { // not self
				// fprintf(stderr, "%lu\n", r[k]);
				int cnum = rid&((1<<8)-1);
				int cid = rid>>8;
 
				uint32_t pos = (uint32_t) r[k] >> 1;
				uint32_t dir = r[k]&1;
				// && !(tid_ == cnum && i_ == cid) 
				if (dir == dir_ori && !t->flag[tid_][i_] && !t->flag[cnum][cid] && 
					// match_pro(p->ref, reads->clusters[cnum].a[cid].ref, pos_ori, pos) <= 4) {
					// match_pro(p->ref, reads->clusters[t->index][cnum].a[cid].ref, pos_ori, pos) <= diff_threshold) {
					// match_pro(p->ref, reads->clusters[t->index][cnum].a[cid].ref, pos_ori, pos) <= diff_threshold*2) {
					match_pro(p->ref, reads->clusters[t->index][cnum].a[cid].ref, pos_ori, pos) <= cbthreshold) {

// fprintf(stderr, "%s\n%s\n", p->ref, reads->clusters[t->index][cnum].a[cid].ref);
// exit(0);
					//combine to one cluster
					cluster_t *cori = &reads->clusters[t->index][cnum].a[cid];

					cluster_t temp_p;//result
					// kv_pushp(cluster_t, reads->clusters[t->index^1][tid_ori], &tp);
					kv_init(temp_p);
					kv_resize(uint64_t, temp_p, p->n + cori->n);

					if (pos_ori >= pos) {//am.bpos = pos_ori - pos;
						for (int i = 0; i < p->n; ++i) {
							kv_push(uint64_t, temp_p, p->a[i]);
						}
						for (int i = 0; i < cori->n; ++i) {
							// ct->a[j] (rid-32 bits | pos | dir-1 bit)
							y = cori->a[i];
							// if (((uint64_t)((uint32_t)y >> 1) + (uint64_t)(pos_ori - pos)) == 283 && strcmp("TTGGGCAGAAAATGGACCAAGGAAGAATTGTCAAAAATACAGCTATGATGAGAGAAGCTGCAAGAAAAATAGAAGAAAGGCATTTTTCCAACCATGCAAC", reads->seq[(uint32_t) y >> 32].seq)) {
							// 	fprintf(stderr, "11111 y>>1: %u, %d\n", (uint32_t)y >> 1, (pos_ori - pos));
							// 	fprintf(stderr, "11111 pos_ori: %u; pos: %u\n", pos_ori, pos);
							// }
							y = (uint64_t)y >> 32 << 32 | (((uint64_t)((uint32_t)y >> 1) + (uint64_t)(pos_ori - pos)) << 1) | (y&1);
							kv_push(uint64_t, temp_p, y);
						}
					} else {
						for (int i = 0; i < cori->n; ++i) {
							kv_push(uint64_t, temp_p, cori->a[i]);
						}
						for (int i = 0; i < p->n; ++i) {
							y = p->a[i];
							y = (uint64_t)y >> 32 << 32 | (((uint64_t)((uint32_t)y >> 1) + (uint64_t)(pos - pos_ori)) << 1) | (y&1);
							kv_push(uint64_t, temp_p, y);
						}
					}

					construct_ref2(reads, &temp_p);

					// if (temp_p.ennum + strlen(temp_p.ref) <= p->ennum + cori->ennum + strlen(p->ref) + strlen(cori->ref)) {
						pthread_mutex_lock(&t->mutex[tid_][i_]);
						// pthread_mutex_lock(&t->mutex[cnum][cid]);

						while ( pthread_mutex_trylock(&t->mutex[cnum][cid]) ) { /* Test if already locked   */
							pthread_mutex_unlock(&t->mutex[tid_][i_]);  /* Free resource to avoid deadlock */
							sleep(1);
							pthread_mutex_lock(&t->mutex[tid_][i_]);
						}

						if (!t->flag[tid_][i_] && !t->flag[cnum][cid]) {
							t->flag[tid_][i_] = true;
							t->flag[cnum][cid] = true;
							flagEncode = true;
						}
						pthread_mutex_unlock(&t->mutex[cnum][cid]);
						pthread_mutex_unlock(&t->mutex[tid_][i_]);

						if (flagEncode) {
							cluster_t *tp;//result
							kv_pushp(cluster_t, reads->clusters[t->index^1][tid_ori], &tp);
							kv_init(*tp);
							kv_resize(uint64_t, *tp, temp_p.n);

							kv_copy(uint64_t, *tp, temp_p);
							tp->ref = (char*)calloc(strlen(temp_p.ref) + 1, sizeof(char));
							strcpy(tp->ref, temp_p.ref);
							tp->ennum = temp_p.ennum;

							// fprintf(stderr, "new ref: %s\n", tp->ref);
							mm128_v miniv;
							kv_init(miniv);
							
							int len = strlen(tp->ref);
							// fprintf(stderr, "%s\n", p->ref);
							// if (reads->k < 32) {
							mm_sketch_lh_ori(tp->ref, len, t->win + t->win_step, reads->k, ((t->reads->clusters[t->index^1][tid_ori].n - 1)<<8) + tid_ori, &miniv);
							// } else {
							// 	mm_sketch_lh(tp->ref, len, t->win + t->win_step, reads->k, ((t->reads->clusters[t->index^1][tid_ori].n - 1)<<8) + tid_ori, &miniv);
							// }
							// fprintf(stderr, "mini.n: %lu\n", mini.n);
							for (int j1 = 0; j1 < miniv.n && j1 < first_mininum; ++j1) {
							// for (int j1 = 0; j1 < miniv.n; ++j1) { // ???? number of minimizers
								// fprintf(stderr, "(%lu, %lu)\n", mini.a[i].x, mini.a[i].y>>1);
								// fprintf(stderr, "(%lu, %lu)\n", mini.a[k].x&mask, mini.a[k].y>>1);
								lock = miniv.a[j1].x&mask;
								mm128_v *pp = &reads->mi[t->index^1]->B[lock].a;

								pthread_mutex_lock(&t->lock[lock]);
								kv_push(mm128_t, *pp, miniv.a[j1]);
								pthread_mutex_unlock(&t->lock[lock]);
							}
							kv_destroy(miniv);
							// print_encode(tp, fp);
						}
					// }

					free(temp_p.a);
					free(temp_p.ref);					
				}
				// if (i == 2) fprintf(stderr, "%s\n", reads->clusters[cnum].a[cid].ref);
			} //else printf("the same one\n");
		}
	}

	kv_destroy(mini);
}

static void cp_cluster(kt_cb_for_t *t, int tid_, int i_, int tid_ori) {
	// fprintf(stderr, "%d, %d\n", tid_, i_);
	cluster_t *p = &reads->clusters[t->index][tid_].a[i_];
	//t->index
	int lock;
	int mask = (1<<reads->b) - 1;

	cluster_t *tp;
	kv_pushp(cluster_t, reads->clusters[t->index^1][tid_ori], &tp);
	kv_init(*tp);
	kv_resize(uint64_t, *tp, p->n);
	tp->ref = (char*)calloc(strlen(p->ref) + 1, sizeof(char));
	strcpy(tp->ref, p->ref);	
	tp->ennum = p->ennum;
	kv_copy(uint64_t, *tp, *p);

	mm128_v mini;
	kv_init(mini);
	int len = strlen(tp->ref);
	// fprintf(stderr, "%s\n", p->ref);
	// if (reads->k < 32) {
	mm_sketch_lh_ori(tp->ref, len, t->win + t->win_step, reads->k, ((t->reads->clusters[t->index^1][tid_ori].n - 1)<<8) + tid_ori, &mini);
	// } else {
	// 	mm_sketch_lh(tp->ref, len, t->win + t->win_step, reads->k, ((t->reads->clusters[t->index^1][tid_ori].n - 1)<<8) + tid_ori, &mini);
	// }
	// fprintf(stderr, "mini.n: %lu\n", mini.n);
	for (int j = 0; j < mini.n && j < first_mininum; ++j) {
		// fprintf(stderr, "(%lu, %lu)\n", mini.a[i].x, mini.a[i].y>>1);
		// fprintf(stderr, "(%lu, %lu)\n", mini.a[k].x&mask, mini.a[k].y>>1);
		lock = mini.a[j].x&mask;
		mm128_v *pp = &reads->mi[t->index^1]->B[lock].a;

		pthread_mutex_lock(&t->lock[lock]);
		kv_push(mm128_t, *pp, mini.a[j]);
		pthread_mutex_unlock(&t->lock[lock]);
	}
	kv_destroy(mini);
}

static inline long steal_cb_work(kt_cb_for_t *t, int *tid_, int *i_)
{
	int i, max_i = -1;
	long max = -1;
	for (i = 0; i < t->n_threads; ++i) {
		if (t->w[i].n - t->w[i].i > max) {
			max = t->w[i].n - t->w[i].i, max_i = i;
		}
		// if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	}
	if (max <= 0) {
		return -1;
	}
	*i_ = __sync_fetch_and_add(&t->w[max_i].i, 1);
	*tid_ = t->w[max_i].tid;
	// k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	if (*i_ >= t->w[max_i].n) return -1;
	return 1;
}

static void *ktf_find_next_worker(void *data)
{
	ktf_cb_worker_t *w = (ktf_cb_worker_t*)data;
	kt_cb_for_t *t = w->t;

	long i, cur_tid = w - w->t->w;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, 1);
		if (i >= w->n) break;
		// cb_link(t, i, w - w->t->w, fp);
		// find and combine next
		if (!t->flag[w->tid][i]) find_next(t, w->tid, i, cur_tid);
	}
	int tid_, i_;
	while (steal_cb_work(t, &tid_, &i_) > 0) {
		find_next(t, tid_, i_, cur_tid);
	}

	// fprintf(stderr, "tid: %ld over\n", w - w->t->w);
	pthread_exit(0);
}


static void *ktf_find_next_worker2(void *data)
{
	ktf_cb_worker_t *w = (ktf_cb_worker_t*)data;
	kt_cb_for_t *t = w->t;

	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, 1);
		if (i >= w->n) break;
		// cb_link(t, i, w - w->t->w, fp);
		// find and combine next
		if (!t->flag[w->tid][i]) cp_cluster(t, w->tid, i, w - w->t->w);
	}
	// while (steal_cb_work(t, w - w->t->w) > 0) {
	// 	// cb_link(t, i, w - w->t->w, fp);
	// 	// find and combine next
	// 	if (!t->flag[w->tid][w->i]) find_next(t, w->tid, w->i, w - w->t->w, fp);
	// }

	// fprintf(stderr, "tid: %ld over\n", w - w->t->w);
	pthread_exit(0);
}

void kt_find_next_for(int n_threads, reads_t *reads, int index, int win, int win_step)
{
	int i;	
	kt_cb_for_t t;
	pthread_t *tid;
	t.reads = reads, t.n_threads = n_threads, t.index = index, t.win = win, t.win_step = win_step;
	t.flag = (bool**)calloc(n_threads, sizeof(bool*));
	for (int k = 0; k < n_threads; ++k) {
		t.flag[k] = (bool*)calloc(reads->clusters[index][k].n, sizeof(bool));
		// fprintf(stderr, "%lu\n", reads->clusters[index][k].n);
		memset(t.flag[k], 0, reads->clusters[index][k].n * sizeof(bool));
	}

	// t.w = (ktf_cb_worker_t*)alloca(n_threads * sizeof(ktf_cb_worker_t));
	t.w = (ktf_cb_worker_t*)calloc(n_threads, sizeof(ktf_cb_worker_t));
	// t.mutex = (pthread_mutex_t**)alloca(n_threads * sizeof(pthread_mutex_t*));
	t.mutex = (pthread_mutex_t**)calloc(n_threads, sizeof(pthread_mutex_t*));
	for (int i = 0; i < n_threads; ++i) {
		// t.mutex[i] = (pthread_mutex_t*)alloca(reads->clusters[index][i].n * sizeof(pthread_mutex_t));
		t.mutex[i] = (pthread_mutex_t*)calloc(reads->clusters[index][i].n, sizeof(pthread_mutex_t));
		for (int j = 0; j < reads->clusters[index][i].n; ++j) {
			pthread_mutex_init(&t.mutex[i][j], 0);
		}
	} 

	int num_mutex = (1<<reads->b);
	// t.lock = (pthread_mutex_t*)alloca(num_mutex * sizeof(pthread_mutex_t));
	t.lock = (pthread_mutex_t*)calloc(num_mutex, sizeof(pthread_mutex_t));
	for (int i = 0; i < num_mutex; ++i) {
		pthread_mutex_init(&t.lock[i], 0);
	}

	// tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = 0, t.w[i].tid = i, t.w[i].n = reads->clusters[index][i].n;

	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_find_next_worker, &t.w[i]);
	// for (i = 0; i < 1; ++i) pthread_create(&tid[i], 0, ktf_find_next_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	// FILE *fp = fopen("2.txt", "w");
	// for (int tid = 0; tid < n_threads; ++tid) {
	// 	for (int i = 0; i < reads->clusters[index][tid].n; ++i) {
	// 		if (!t.flag[tid][i]) find_next(&t, tid, i, tid, fp);
	// 	}
	// }
	// fclose(fp);

	// fprintf(stderr, "after pthread_join()..\n");

	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = 0, t.w[i].tid = i, t.w[i].n = reads->clusters[index][i].n;

	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_find_next_worker2, &t.w[i]);
	// for (i = 0; i < 1; ++i) pthread_create(&tid[i], 0, ktf_find_next_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);

	for (int k = 0; k < n_threads; ++k) {
		free(t.flag[k]);
		free(t.mutex[k]);
	}
	free(t.mutex);
	free(t.lock);
	free(t.w);
	free(tid);
	free(t.flag);
}

void combine_cluster(int n_threads, reads_t *reads, int *index_) {
	// fprintf(stderr, "begin combine_cluster()...\n");
	int index = *index_;
	int win = reads->rw, win_step = 0; //70 75 80 85 90 95 100
	// for (int i = 0; i < 10; ++i) {
	int pre_tot_clusters = 0;
	// for (int i = 0; i < 15; ++i) {
	for (int i = 0; ; ++i) {
		// fprintf(stderr, "index: %d\n", index);
		// n_threads = 1;
		mm_idx_generation(n_threads, reads->mi[index]);
		for (int j = 0; j < n_threads; ++j) {
			kv_init(reads->clusters[index^1][j]);
			kv_resize(cluster_t, reads->clusters[index^1][j], 1 << 10);
		}

		// fprintf(stderr, "after mm_idx_generation()...\n");
		reads->mi[index^1] = mm_idx_init(reads->b);

		// fprintf(stderr, "before kt_find_next_for()\n");
		// fprintf(stderr, "win: %d; win_step: %d\n", win, win_step);

		kt_find_next_for(n_threads, reads, index, win, win_step);
		// fprintf(stderr, "after kt_find_next_for()\n");

		mm_idx_destroy(reads->mi[index]);

		for (int j = 0; j < n_threads; ++j) {
			for (int k = 0; k < reads->clusters[index][j].n; ++k) {
				cluster_destroy(reads->clusters[index][j].a[k]);
			}
		}

		index ^= 1;

		/*do {	
			int tot = 0;
			for (int k = 0; k < n_threads; ++k) tot += reads->clusters[index][k].n;
			fprintf(stderr, "cluster number: %d\n", tot);
			int cluster_reads = 0;
		
			for (int k = 0; k < n_threads; ++k) {
				for (int i = 0; i < reads->clusters[index][k].n; ++i) {
					cluster_reads += reads->clusters[index][k].a[i].n;
					
				}
			}
			fprintf(stderr, "cluster contains reads number: %d\n", cluster_reads);
			// fprintf(stderr, "-----------------++++++++++++++++++\n");
		} while (0);*/

		int tot_clusters = 0;
		for (int k = 0; k < n_threads; ++k) tot_clusters += reads->clusters[index][k].n;

		// fprintf(stderr, "pre_tot_clusters: %d; tot_clusters: %d;\n", pre_tot_clusters, tot_clusters);
		if (abs(pre_tot_clusters - tot_clusters) < 100) break;
		pre_tot_clusters = tot_clusters;
	}
	mm_idx_destroy(reads->mi[index]);
	
	*index_ = index;

	//output contig length
	/*char fn[1000];
	FILE *ftmp = fopen("filename.txt", "r");
	fscanf(ftmp, "%s", fn);
	fclose(ftmp);
	FILE *ffp = fopen(fn, "w");
	do {	
		int tot = 0;
		for (int k = 0; k < n_threads; ++k) tot += reads->clusters[index][k].n;
		fprintf(ffp, "%d\n", tot);
		// int cluster_reads = 0;
	
		// for (int k = 0; k < n_threads; ++k) {
		// 	for (int i = 0; i < reads->clusters[index][k].n; ++i) {
		// 		cluster_reads += reads->clusters[index][k].a[i].n;
				
		// 	}
		// }
		// fprintf(stderr, "cluster contains reads number: %d\n", cluster_reads);
		// fprintf(stderr, "-----------------++++++++++++++++++\n");
	} while (0);

	for (int k = 0; k < n_threads; ++k) {
		for (int i = 0; i < reads->clusters[index][k].n; ++i) {
			// cluster_reads += reads->clusters[index][k].a[i].n;
			fprintf(ffp, "%d\n", strlen(reads->clusters[index][k].a[i].ref));
		}
	}
	fclose(ffp);*/
}
