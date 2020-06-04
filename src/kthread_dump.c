#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "breads.h"
#include "kvec.h"
#include "config.h"

/************
 * kt_dump_for() *
 ************/

struct kt_dump_for_t;

typedef struct {
	struct kt_dump_for_t *t;
	long i, n; //i < n;
} ktf_dump_worker_t;

typedef struct kt_dump_for_t {
	int n_threads, index;//win is the length of window
	ktf_dump_worker_t *w;
	reads_t *reads;
	FILE *fffppp;
	int num0, num1;
	pthread_mutex_t mutex;
} kt_dump_for_t;

#ifdef ORDER
// void print_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_pos, FILE *fpdif_char, std::ofstream& fpids, FILE *fffppp) {
// void print_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_pos, FILE *fpdif_char, std::ofstream& fpids) {
void print_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_char, std::ofstream& fpids) {
	qsort(p->a, p->n, sizeof(uint64_t), cmpcluster3);

	// bool debug = false;

	char *temp_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));

	int pre_pos = 0;
	uint32_t rid, pre_rid, temp_rid;
	int pos, dir;
	// fprintf(fpref, "%s\n", p->ref);
	int ref_len = strlen(p->ref);
	for (int i = 0; i < ref_len; ++i) {
		DNA_push(refbin, fpref, seq_nt4_table[(uint8_t)p->ref[i]]);
	}
	// fprintf(stderr, "%s\n", p->ref);
	// fprintf(fp, "%lu %s\n", p->n, p->ref);
	// fprintf(stderr, "%lu ", p->n);
	// if (flag) fprintf(stderr, "p->n: %lu\n", p->n);

	// fprintf(fppos, "%lu ", p->n);
	uint32_t num = p->n;
	uint16_t posbin;
    fppos.write((char*)&num, sizeof(uint32_t));

	for (int k = 0; k < p->n; ++k) {
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		// if (debug) fprintf(stderr, "pos: %d\n", pos);
		dir = y&1;

		bseq1_t *seq = &reads->seq[rid];
		strcpy(temp_str, seq->seq);
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				temp_str[n_pos->a[i]] = 'N';
			}
		}
		/*if (pos >= ref_len) {
			fprintf(stderr, "ref: %s\n", p->ref);
			fprintf(stderr, "seq: %s\n", seq->seq);
			exit(0);
		}*/
		
		// fprintf(fffppp, "%s\n", temp_str);

		if (dir) {
			reverse_complement(temp_str, readlen);
		}
		int en_str_len = 0;
		int eq_char_num = 0;
		for (int rj = pos, tj = 0; tj < readlen; ++rj, ++tj) {
			if (p->ref[rj] != temp_str[tj]) {
				if (eq_char_num > 1) {
	                sprintf(int_str, "%d", eq_char_num);
	                for (char *tk = int_str; *tk != '\0'; ++tk) {
	                    en_str[en_str_len++] = *tk;
	                }
	                eq_char_num = 0;
	            } else {
	                for (int i = tj - eq_char_num; i < tj; ++i) {
	                    en_str[en_str_len++] = temp_str[i];
	                }
	                eq_char_num = 0;
	            }
				en_str[en_str_len++] = temp_str[tj];
			} else ++eq_char_num;
		}
		if (en_str_len == 0) {
			en_str[en_str_len++] = '0';
		}
		en_str[en_str_len] = '\0';
		fprintf(fpdif_char, "%s\n", en_str);
		// fprintf(fppos, "%d ", (int)((uint32_t)p->a[k]>>1) - pre_pos);
		// fprintf(fppos, "%d ", (int)(pos - pre_pos));
		posbin = (uint16_t)(pos - pre_pos);
		fppos.write((char*)&posbin, sizeof(uint16_t));

		//---------------------------------------------// for order
		if (k == 0 || posbin > 0) { //the first reads or different begin position
			fpids.write((char*)&rid, sizeof(uint32_t));
			// fprintf(fffppp, "%u\n", rid);
		} else
		if (posbin == 0) { //with the same begin position
			temp_rid = rid - pre_rid;
			fpids.write((char*)&temp_rid, sizeof(uint32_t));
			// fprintf(fffppp, "%u\n", temp_rid);
			// fprintf(fffppp, "%u\n", rid);
		}
		pre_rid = rid;
		//---------------------------------------------// 
		// fprintf(fpdir, "%d", dir);
		bit_push(dirbin, fpdir, dir);
		// fprintf(fpdir, "%d", (int)p->a[k]&1);

		pre_pos = pos;
		// pre_pos = (uint32_t)p->a[k]>>1;
		// pos[0] - pos[k];
		// fprintf(fp, "%s\n", en_str);
	}
}
#else
// void print_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_pos, FILE *fpdif_char, FILE *fffppp) {
// void print_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_pos, FILE *fpdif_char) {
void print_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_char) {
	qsort(p->a, p->n, sizeof(uint64_t), cmpcluster2);
	
	bool debug = false;

	// int readlen = reads->readlen;

	char *temp_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));

	int pre_pos = 0;
	uint32_t rid;
	int pos, dir;
	// fprintf(fpref, "%s\n", p->ref);
	int ref_len = strlen(p->ref);
	for (int i = 0; i < ref_len; ++i) {
		DNA_push(refbin, fpref, seq_nt4_table[(uint8_t)p->ref[i]]);
	}
	// fprintf(stderr, "%s\n", p->ref);
	// fprintf(fp, "%lu %s\n", p->n, p->ref);
	// fprintf(stderr, "%lu ", p->n);
	// if (flag) fprintf(stderr, "p->n: %lu\n", p->n);

	// fprintf(fppos, "%lu ", p->n);
	uint32_t num = p->n;
	uint16_t posbin;
	fppos.write((char*)&num, sizeof(uint32_t));

	for (int k = 0; k < p->n; ++k) {
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		// if (debug) fprintf(stderr, "pos: %d\n", pos);
		dir = y&1;

		bseq1_t *seq = &reads->seq[rid];
		strcpy(temp_str, seq->seq);
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				temp_str[n_pos->a[i]] = 'N';
			}
		}
		// if (pos >= ref_len) {
		// 	fprintf(stderr, "ref: %s\n", p->ref);
		// 	fprintf(stderr, "seq: %s\n", seq->seq);
		// 	exit(0);
		// }
		
		// fprintf(fffppp, "%s\n", temp_str);

		if (dir) {
			reverse_complement(temp_str, readlen);
		}
		int en_str_len = 0;
		int eq_char_num = 0;
		for (int rj = pos, tj = 0; tj < readlen; ++rj, ++tj) {
			if (p->ref[rj] != temp_str[tj]) {
				if (eq_char_num > 1) {
	                sprintf(int_str, "%d", eq_char_num);
	                for (char *tk = int_str; *tk != '\0'; ++tk) {
	                    en_str[en_str_len++] = *tk;
	                }
	                eq_char_num = 0;
	            } else {
	                for (int i = tj - eq_char_num; i < tj; ++i) {
	                    en_str[en_str_len++] = temp_str[i];
	                }
	                eq_char_num = 0;
	            }
				en_str[en_str_len++] = temp_str[tj];
			} else ++eq_char_num;
		}
		if (en_str_len == 0) {
			en_str[en_str_len++] = '0';
		}
		en_str[en_str_len] = '\0';
		fprintf(fpdif_char, "%s\n", en_str);
		// fprintf(fppos, "%d ", (int)((uint32_t)p->a[k]>>1) - pre_pos);
		// fprintf(fppos, "%d ", (int)(pos - pre_pos));
		posbin = (uint16_t)(pos - pre_pos);
		fppos.write((char*)&posbin, sizeof(uint16_t));

		// fprintf(fpdir, "%d", dir);
		bit_push(dirbin, fpdir, dir);
		// fprintf(fpdir, "%d", (int)p->a[k]&1);

		pre_pos = pos;
		// pre_pos = (uint32_t)p->a[k]>>1;
		// pos[0] - pos[k];
		// fprintf(fp, "%s\n", en_str);
	}
}
#endif

static void *ktf_dump_worker(void *data)
{
	ktf_dump_worker_t *w = (ktf_dump_worker_t*)data;
	kt_dump_for_t *t = w->t;

	char name[100]; 
	sprintf(name, "%s/ref.bin.%ld", folder, w - w->t->w);
	// FILE *fpref = fopen(name, "w");
	std::ofstream fpref(name, std::ios::binary);
	
	sprintf(name, "%s/beg_pos.bin.%ld", folder, w - w->t->w);
	// FILE *fppos = fopen(name, "w");
	std::ofstream fppos(name, std::ios::binary);
	
	sprintf(name, "%s/dir.bin.%ld", folder, w - w->t->w);
	std::ofstream fpdir(name, std::ios::binary);
	
	// sprintf(name, "%s/dif_pos.txt.%ld", folder, w - w->t->w);
	// FILE *fpdif_pos = fopen(name, "w");

	sprintf(name, "%s/dif_char.txt.%ld", folder, w - w->t->w);
	FILE *fpdif_char = fopen(name, "w");
	// fprintf(stderr, "%s\n", name);

	// sprintf(name, "%s/ids.txt.%ld", folder, w - w->t->w);
	// FILE *fffppp = fopen(name, "w");

	#ifdef ORDER
	sprintf(name, "%s/ids.bin.%ld", folder, w - w->t->w);
	std::ofstream fpids(name, std::ios::binary);
	#endif

	uint8bit_v refbin, dirbin;
	refbin.n = refbin.a = 0;
	dirbin.n = dirbin.a = 0;

	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, 1);
		if (i >= w->n) break;
		cluster_t *p = &reads->clusters[t->index][w - w->t->w].a[i];	
		// if (t->last_rounds) 
		// update_reference(reads, p, p->n); 
		// if (strcmp(p->ref, "TGAAGAGTGACAGGTTGGCTCAACTTTTTTCTTTTTTTTTTTTAATTGCCTTGTATTGTAAGTATTCTTCCCTGCAGTCCAAGTGACTTTTCATTTTTTGTTTTA") ==0) 
			// fprintf(stderr, "name: %s\n", name);
		// print_encode(t, p, fpref, refbin, fppos, fpdir, dirbin, fpdif_pos, fpdif_char);
		#ifdef ORDER
		// print_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_pos, fpdif_char, fpids, fffppp);
		// print_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_pos, fpdif_char, fpids);
		print_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_char, fpids);
		#else
		// print_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_pos, fpdif_char, fffppp);
		// print_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_pos, fpdif_char);
		print_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_char);
		#endif
		// if (i > 100) break;
	}
	// while (steal_cb_work(t, w - w->t->w) > 0) {
	// 	// cb_link(t, i, w - w->t->w, fp);
	// 	// find and combine next
	// 	if (!t->flag[w->tid][w->i]) find_next(t, w->tid, w->i, w - w->t->w, fp);
	// }
	// fclose(fpref);
	// fclose(fppos);
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	
	if (refbin.n > 0) {
		fpref.write((char*)&refbin.a, sizeof(uint8_t));
	}
	fpref.close();
	fppos.close();
	fpdir.close();
	#ifdef ORDER
	fpids.close();
	#endif

	// fclose(fpdif_pos);
	fclose(fpdif_char);

	// fclose(fffppp);

	// fprintf(stderr, "tid: %ld over\n", w - w->t->w);
	pthread_exit(0);
}

void kt_dump_for(int n_threads, reads_t *reads, int index)
{
	int i;	
	kt_dump_for_t t;
	pthread_t *tid;
	t.reads = reads, t.n_threads = n_threads, t.index = index;

	t.num0 = t.num1 = 0;
	pthread_mutex_init(&t.mutex, 0);
	// t.fffppp = fopen("compressed.seq", "w");

	t.w = (ktf_dump_worker_t*)alloca(n_threads * sizeof(ktf_dump_worker_t));

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = 0, t.w[i].n = reads->clusters[index][i].n;

	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_dump_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	// FILE *fp = fopen("2.txt", "w");
	// for (int tid = 0; tid < n_threads; ++tid) {
	// 	for (int i = 0; i < reads->clusters[index][tid].n; ++i) {
	// 		if (!t.flag[tid][i]) find_next(&t, tid, i, tid, fp);
	// 	}
	// }
	// fclose(fp);
	// fprintf(stderr, "t.num0: %d; t.num1: %d\n", t.num0, t.num1);
	// fclose(t.fffppp);
	// fprintf(stderr, "after pthread_join()..\n");
	// 
}

int cmp(const void *a, const void *b) {
	int rid0 = *(int*)a;
	int rid1 = *(int*)b;
	return strcmp(reads->seq[rid0].seq, reads->seq[rid1].seq);
}

void cluster_dump(int n_threads, reads_t *reads, int index) {
	// fprintf(stderr, "begin kt_dump_for()\n");
	// n_threads = 1;
	kt_dump_for(n_threads, reads, index);	

	char name[100]; 
	sprintf(name, "%s/info.txt", folder);
	FILE *fp_info = fopen(name, "w");
	// sp_reads_t *s = reads->sp;
	// int readlen = reads->readlen;
	
	fprintf(fp_info, "%d %d\n", readlen, n_threads);
	fprintf(fp_info, "%d %d %d\n", reads->sp->allA, reads->sp->allT, reads->sp->allN);
	#ifdef ORDER
	fprintf(fp_info, "%u\n", reads->n_seq); //number of
	#endif
	fclose(fp_info);

	// reads->singleFile = fopen("compfiles/single_enc.fasta", "w");
	sprintf(name, "%s/single.seq", folder);
	std::ofstream singleFile(name, std::ios::binary);
	uint8bit_v singlebin;
	singlebin.a = singlebin.n = 0;

	char *temp_str = (char*)alloca((readlen + 1) * sizeof(char));
	int rid;
	for (int i = 0; i < reads->sg.n; ++i) {
		if (!reads->sg_flag[i]) {
			rid = reads->sg.a[i];
			bseq1_t *seq = &reads->seq[rid];
			strcpy(temp_str, seq->seq);
			// put N back to reads;
			uint32_v *n_pos = (uint32_v*) seq->n_pos;
			if (n_pos != NULL) {
				reads->sg_flag[i] = true;
				for (int i = 0; i < n_pos->n; ++i) {
					temp_str[n_pos->a[i]] = 'N';
				}
				// fprintf(reads->Nfile, "%s\n", temp_str);
				kv_push(uint32_t, reads->Nfile_id, rid);
			} else {
				// fprintf(fp, "@SRR490961\n");
				// fprintf(reads->singleFile, "%s\n", temp_str);
				// fprintf(fp, "+\n");
				// fprintf(fp, "%s\n", temp_str);
				#ifdef ORDER
				kv_push(uint32_t, reads->singleFile_id, rid);
				#else
				for (int j = 0; j < readlen; ++j) {
					DNA_push(singlebin, singleFile, seq_nt4_table[(uint8_t)temp_str[j]]);
				}
				#endif
			}
		}
	}

#ifdef ORDER
	std::sort(reads->fpA_id.a, reads->fpA_id.a + reads->fpA_id.n);
	std::sort(reads->fpT_id.a, reads->fpT_id.a + reads->fpT_id.n);
	std::sort(reads->fpN_id.a, reads->fpN_id.a + reads->fpN_id.n);
	std::sort(reads->Nfile_id.a, reads->Nfile_id.a + reads->Nfile_id.n);
	std::sort(reads->singleFile_id.a, reads->singleFile_id.a + reads->singleFile_id.n);

	sp_reads_t *s = reads->sp;
	std::sort(s->allA_id.a, s->allA_id.a + s->allA_id.n);
	std::sort(s->allT_id.a, s->allT_id.a + s->allT_id.n);
	std::sort(s->allN_id.a, s->allN_id.a + s->allN_id.n);

	uint32_t temp;

	sprintf(name, "%s/allA.ids.bin", folder);
	std::ofstream allAids(name, std::ios::binary);
	for (int i = 0; i < s->allA_id.n; ++i) {
		if (i == 0) {
			allAids.write((char*)&s->allA_id.a[i], sizeof(uint32_t));
		} else {
			temp = s->allA_id.a[i] - s->allA_id.a[i-1];
			allAids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	allAids.close();

	sprintf(name, "%s/allT.ids.bin", folder);
	std::ofstream allTids(name, std::ios::binary);
	for (int i = 0; i < s->allT_id.n; ++i) {
		if (i == 0) {
			allTids.write((char*)&s->allT_id.a[i], sizeof(uint32_t));
		} else {
			temp = s->allT_id.a[i] - s->allT_id.a[i-1];
			allTids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	allTids.close();

	sprintf(name, "%s/allN.ids.bin", folder);
	std::ofstream allNids(name, std::ios::binary);
	for (int i = 0; i < s->allN_id.n; ++i) {
		if (i == 0) {
			allNids.write((char*)&s->allN_id.a[i], sizeof(uint32_t));
		} else {
			temp = s->allN_id.a[i] - s->allN_id.a[i-1];
			allNids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	allNids.close();

	sprintf(name, "%s/AA.ids.bin", folder);
	std::ofstream fpAids(name, std::ios::binary);
	for (int i = 0; i < reads->fpA_id.n; ++i) {
		if (i == 0) {
			fpAids.write((char*)&reads->fpA_id.a[i], sizeof(uint32_t));
		} else {
			temp = reads->fpA_id.a[i] - reads->fpA_id.a[i-1];
			fpAids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	fpAids.close();

	sprintf(name, "%s/TT.ids.bin", folder);
	std::ofstream fpTids(name, std::ios::binary);
	for (int i = 0; i < reads->fpT_id.n; ++i) {
		if (i == 0) {
			fpTids.write((char*)&reads->fpT_id.a[i], sizeof(uint32_t));
		} else {
			temp = reads->fpT_id.a[i] - reads->fpT_id.a[i-1];
			fpTids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	fpTids.close();

	sprintf(name, "%s/NN.ids.bin", folder);
	std::ofstream fpNids(name, std::ios::binary);
	for (int i = 0; i < reads->fpN_id.n; ++i) {
		if (i == 0) {
			fpNids.write((char*)&reads->fpN_id.a[i], sizeof(uint32_t));
		} else {
			temp = reads->fpN_id.a[i] - reads->fpN_id.a[i-1];
			fpNids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	fpNids.close();

	sprintf(name, "%s/Nfile.ids.bin", folder);
	std::ofstream fpNfileids(name, std::ios::binary);
	for (int i = 0; i < reads->Nfile_id.n; ++i) {
		if (i == 0) {
			fpNfileids.write((char*)&reads->Nfile_id.a[i], sizeof(uint32_t));
		} else {
			temp = reads->Nfile_id.a[i] - reads->Nfile_id.a[i-1];
			fpNfileids.write((char*)&temp, sizeof(uint32_t));
		}
	}
	fpNfileids.close();

	for (int i = 0; i < reads->singleFile_id.n; ++i) {
		rid = reads->singleFile_id.a[i];
		bseq1_t *seq = &reads->seq[rid];
		strcpy(temp_str, seq->seq);

		for (int j = 0; j < readlen; ++j) {
			DNA_push(singlebin, singleFile, seq_nt4_table[(uint8_t)temp_str[j]]);
		}
	}

	sprintf(name, "%s/singleFile.ids.bin", folder);
	std::ofstream fpsingleFileids(name, std::ios::binary);
	// FILE *fsfspp = fopen("special.ids.bin", "w");
	for (int i = 0; i < reads->singleFile_id.n; ++i) {
		if (i == 0) {
			fpsingleFileids.write((char*)&reads->singleFile_id.a[i], sizeof(uint32_t));
			// fprintf(fsfspp, "%u %u\n", reads->singleFile_id.a[i], reads->singleFile_id.a[i]);
		} else {
			temp = reads->singleFile_id.a[i] - reads->singleFile_id.a[i-1];
			fpsingleFileids.write((char*)&temp, sizeof(uint32_t));
			// fprintf(fsfspp, "%u %u\n", reads->singleFile_id.a[i], temp);
		}
	}
	// fclose(fsfspp);
	fpsingleFileids.close();
#endif

	if (singlebin.n > 0) {
		singleFile.write((char*)&singlebin.a, sizeof(uint8_t));
	}
	singleFile.close();
	// fclose(reads->singleFile);

	//
	sprintf(name, "%s/AA.txt", folder);
	reads->fpa = fopen(name, "w");
	sprintf(name, "%s/TT.txt", folder);
	reads->fpt = fopen(name, "w");
	sprintf(name, "%s/NN.txt", folder);
	reads->fpn = fopen(name, "w");
	sprintf(name, "%s/single_N.seq", folder);
	reads->Nfile = fopen(name, "w");

	char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));
	char *cur_seq;
	int en_str_len = 0, eq_char_num = 0;

	for (int i = 0; i < reads->fpA_id.n; ++i) {
		bseq1_t *seq = &reads->seq[reads->fpA_id.a[i]];
		cur_seq = seq->seq;
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				cur_seq[n_pos->a[i]] = 'N';
			}
		}
		// store as special file "A.thread_id"
		// 1. a bit stream; 2. list of position with delta encoding; 
		// 3. number of N;  4. list of position of N
		en_str_len = eq_char_num = 0;
		for (int tj = 0; tj < readlen; ++tj) {
			if ('A' != cur_seq[tj]) {
				if (eq_char_num > 0) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				}
				en_str[en_str_len++] = cur_seq[tj];
			} else ++eq_char_num;
		}
		if (en_str_len == 0) {
			en_str[en_str_len++] = '0';
		}
		en_str[en_str_len] = '\0';
		fprintf(reads->fpa, "%s\n", en_str); 
	}

	for (int i = 0; i < reads->fpT_id.n; ++i) {
		bseq1_t *seq = &reads->seq[reads->fpT_id.a[i]];
		cur_seq = seq->seq;
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				cur_seq[n_pos->a[i]] = 'N';
			}
		}
		en_str_len = eq_char_num = 0;
		for (int tj = 0; tj < readlen; ++tj) {
			if ('T' != cur_seq[tj]) {
				if (eq_char_num > 0) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				}
				en_str[en_str_len++] = cur_seq[tj];
			} else ++eq_char_num;
		}
		if (en_str_len == 0) {
			en_str[en_str_len++] = '0';
		}
		en_str[en_str_len] = '\0';
		fprintf(reads->fpt, "%s\n", en_str);
	}

	for (int i = 0; i < reads->fpN_id.n; ++i) {
		bseq1_t *seq = &reads->seq[reads->fpN_id.a[i]];
		cur_seq = seq->seq;
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				cur_seq[n_pos->a[i]] = 'N';
			}
		}
		en_str_len = eq_char_num = 0;
		for (int tj = 0; tj < readlen; ++tj) {
			if ('N' != cur_seq[tj]) {
				if (eq_char_num > 0) {
					sprintf(int_str, "%d", eq_char_num);
					for (char *tk = int_str; *tk != '\0'; ++tk) {
						en_str[en_str_len++] = *tk;
					}
					eq_char_num = 0;
				}
				en_str[en_str_len++] = cur_seq[tj];
			} else ++eq_char_num;
		}
		if (en_str_len == 0) {
			en_str[en_str_len++] = '0';
		}
		en_str[en_str_len] = '\0';
		fprintf(reads->fpn, "%s\n", en_str);
	}

	for (int i = 0; i < reads->Nfile_id.n; ++i) {
		rid = reads->Nfile_id.a[i];
		bseq1_t *seq = &reads->seq[rid];
		strcpy(temp_str, seq->seq);
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				temp_str[n_pos->a[i]] = 'N';
			}
		}
		fprintf(reads->Nfile, "%s\n", temp_str);
	}

	fclose(reads->fpa);
	fclose(reads->fpt);
	fclose(reads->fpn);
	fclose(reads->Nfile);
	free(reads->sg_flag);
}
