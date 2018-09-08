#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "breads.h"
#include "kvec.h"
#include "khash.h"

/************
 * kt_for() *
 ************/

struct kt_idx_for_t;

typedef struct {
	struct kt_idx_for_t *t;
	long i;
} ktf_idx_worker_t;

typedef struct kt_idx_for_t {
	int n_threads;
	long n;
	ktf_idx_worker_t *w;
	void (*func)(void*,long,int);
	void *data;
} kt_idx_for_t;

static inline long idx_steal_work(kt_idx_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void *ktf_idx_worker(void *data)
{
	ktf_idx_worker_t *w = (ktf_idx_worker_t*)data;
	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->w);
	}
	while ((i = idx_steal_work(w->t)) >= 0)
		w->t->func(w->t->data, i, w - w->t->w);
	pthread_exit(0);
}

void kt_idx_for(int n_threads, void (*func)(void*,long,int), void *data, long n)
{
	// fprintf(stderr, "n_threads (in kt_idx_for): %d\n", n_threads);
	int i;
	kt_idx_for_t t;
	pthread_t *tid;
	t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
	t.w = (ktf_idx_worker_t*)alloca(n_threads * sizeof(ktf_idx_worker_t));
	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i)
		t.w[i].t = &t, t.w[i].i = i;
	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_idx_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

mm_idx_t *mm_idx_init(int b) {
	mm_idx_t *mi;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<reads->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>reads->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<reads->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	free(mi);
}

static void worker_post(void *g, long i, int tid)
{
	// fprintf(stderr, "i: %ld\n", i);
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	// b->n is the number of minimizers appearing >1 times
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);//uint64_t *p; // position array for minimizers appearing >1 times

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>reads->b<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	free(b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
void mm_idx_generation(int n_threads, mm_idx_t *mi)
{
	kt_idx_for(n_threads, worker_post, mi, 1<<reads->b);
}