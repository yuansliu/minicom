#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "breads.h"
#include "kvec.h"
#include "config.h"
#include <omp.h>

/************
 * kt_dump_pe_for() *
 ************/

struct kt_dump_pe_for_t;

typedef struct {
	struct kt_dump_pe_for_t *t;
	long i, n; //i < n;
} ktf_dump_pe_worker_t;

typedef struct kt_dump_pe_for_t {
	int n_threads, index;//win is the length of window
	ktf_dump_pe_worker_t *w;
	reads_t *reads;
	FILE *fffppp;
	int num0, num1;
	pthread_mutex_t mutex;
} kt_dump_pe_for_t;

// int half_val, mpvid;
int mpvid;
int *mpv;

void print_pe_encode(cluster_t *p, std::ofstream& fpref, uint8bit_v& refbin, std::ofstream& fppos, std::ofstream& fpdir, uint8bit_v& dirbin, FILE *fpdif_char, FILE *fpids) {
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
	
	uint32_t num = p->n;
	uint16_t posbin;
    fppos.write((char*)&num, sizeof(uint32_t));

	for (int k = 0; k < p->n; ++k) {
		uint64_t y = p->a[k];
		rid = y>>32;
		pos = (uint32_t)y>>1;
		// if (debug) fprintf(stderr, "pos: %d\n", pos);
		dir = y&1;

		if (rid < half_val) {
			fprintf(fpids, "0 %u\n", rid);
		} else {
			fprintf(fpids, "1 %u\n", rid);
		}

		bseq1_t *seq = &reads->seq[rid];
		strcpy(temp_str, seq->seq);
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				temp_str[n_pos->a[i]] = 'N';
			}
		}
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

		bit_push(dirbin, fpdir, dir);

		pre_pos = pos;
		// pre_pos = (uint32_t)p->a[k]>>1;
		// pos[0] - pos[k];
		// fprintf(fp, "%s\n", en_str);
	}
}

static void *ktf_dump_pe_worker(void *data)
{
	ktf_dump_pe_worker_t *w = (ktf_dump_pe_worker_t*)data;
	kt_dump_pe_for_t *t = w->t;

	char name[100]; 
	sprintf(name, "%s/ref.bin.%ld", folder, w - w->t->w);
	// FILE *fpref = fopen(name, "w");
	std::ofstream fpref(name, std::ios::binary);
	
	sprintf(name, "%s/beg_pos.bin.%ld", folder, w - w->t->w);
	// FILE *fppos = fopen(name, "w");
	std::ofstream fppos(name, std::ios::binary);
	
	sprintf(name, "%s/dir.bin.%ld", folder, w - w->t->w);
	std::ofstream fpdir(name, std::ios::binary);
	
	sprintf(name, "%s/dif_char.txt.%ld", folder, w - w->t->w);
	FILE *fpdif_char = fopen(name, "w");

	sprintf(name, "%s/ids.txt.%ld", folder, w - w->t->w);
	FILE *fpids = fopen(name, "w"); //store reads from which file and its postion in reads

	// #ifdef ORDER
	// sprintf(name, "%s/ids.bin.%ld", folder, w - w->t->w);
	// std::ofstream fpids(name, std::ios::binary);
	// #endif

	uint8bit_v refbin, dirbin; //filebin is 0 1, record reads from which file
	refbin.n = refbin.a = 0;
	dirbin.n = dirbin.a = 0;
	// filebin.n = filebin.a = 0;

	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, 1);
		if (i >= w->n) break;
		cluster_t *p = &reads->clusters[t->index][w - w->t->w].a[i];	
		print_pe_encode(p, fpref, refbin, fppos, fpdir, dirbin, fpdif_char, fpids);
	}
	
	if (dirbin.n > 0) { 
		fpdir.write((char*)&dirbin.a, sizeof(uint8_t));
	}
	
	if (refbin.n > 0) {
		fpref.write((char*)&refbin.a, sizeof(uint8_t));
	}

	fpref.close();
	fppos.close();
	fpdir.close();
	fclose(fpids);
	fclose(fpdif_char);

	pthread_exit(0);
}

void kt_dump_pe_for(int n_threads, reads_t *reads, int index)
{
	int i;	
	kt_dump_pe_for_t t;
	pthread_t *tid;
	t.reads = reads, t.n_threads = n_threads, t.index = index;

	t.num0 = t.num1 = 0;
	pthread_mutex_init(&t.mutex, 0);
	// t.fffppp = fopen("compressed.seq", "w");

	t.w = (ktf_dump_pe_worker_t*)alloca(n_threads * sizeof(ktf_dump_pe_worker_t));

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (i = 0; i < n_threads; ++i) 
		t.w[i].t = &t, t.w[i].i = 0, t.w[i].n = reads->clusters[index][i].n;

	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_dump_pe_worker, &t.w[i]);
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

void cluster_dump_pe(int n_threads, reads_t *reads, int index) {
	char name[100]; 
	sprintf(name, "%s/info.txt", folder);
	FILE *fp_info = fopen(name, "w");
	// sp_reads_t *s = reads->sp;
	// int readlen = reads->readlen;
	
	fprintf(fp_info, "%d %d\n", readlen, n_threads);

	// half_val = reads->n_seq >> 1;
	fprintf(fp_info, "%d\n", half_val); //>= half_val is reads from the second file

	mpv = (int*)calloc(half_val, sizeof(int));
	memset(mpv, 0, half_val * sizeof(int));
	mpvid = 0;

	fprintf(fp_info, "%d %d %d\n", reads->sp->allA, reads->sp->allT, reads->sp->allN);
	fclose(fp_info);
	// fprintf(stderr, "begin kt_dump_pe_for()\n");
	kt_dump_pe_for(n_threads, reads, index);	

	int a, b;
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
				kv_push(uint32_t, reads->singleFile_id, rid);
			}
		}
	}

	// #ifdef ORDER
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
	{
		//get the order
		for (int i = 0; i < s->allA_id.n; ++i) {// sprintf(name, "%s/allA.ids.bin", folder);
			if (s->allA_id.a[i] < half_val) {
				mpv[s->allA_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < s->allT_id.n; ++i) {// sprintf(name, "%s/allT.ids.bin", folder);
			if (s->allT_id.a[i] < half_val) {
				mpv[s->allT_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < s->allN_id.n; ++i) {// sprintf(name, "%s/allN.ids.bin", folder);
			if (s->allN_id.a[i] < half_val) {
				mpv[s->allN_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < reads->fpA_id.n; ++i) {// sprintf(name, "%s/AA.ids.bin", folder);
			if (reads->fpA_id.a[i] < half_val) {
				mpv[reads->fpA_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < reads->fpT_id.n; ++i) { // sprintf(name, "%s/TT.ids.bin", folder);
			if (reads->fpT_id.a[i] < half_val) {
				mpv[reads->fpT_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < reads->fpN_id.n; ++i) {// sprintf(name, "%s/NN.ids.bin", folder);
			if (reads->fpN_id.a[i] < half_val) {
				mpv[reads->fpN_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < reads->Nfile_id.n; ++i) {// sprintf(name, "%s/Nfile.ids.bin", folder); singleton contain N
			if (reads->Nfile_id.a[i] < half_val) {
				mpv[reads->Nfile_id.a[i]] = mpvid ++;
			}
		}
		for (int i = 0; i < reads->singleFile_id.n; ++i) {// sprintf(name, "%s/singleFile.ids.bin", folder);
			if (reads->singleFile_id.a[i] < half_val) {
				mpv[reads->singleFile_id.a[i]] = mpvid ++;
			}
		}

		//
		for (int i = 0; i < n_threads; ++i) {
			sprintf(name, "%s/ids.txt.%ld", folder, i);
			FILE *fpin = fopen(name, "r");
			while (fscanf(fpin, "%d%d", &a, &b) != EOF) {
				if (a == 0) {
					mpv[b] = mpvid ++;
				}
			}
			fclose(fpin);
		}
	}

	/// write data  //peids.bin.
	sprintf(name, "%s/peids.bin.sp", folder);
	std::ofstream fpids(name, std::ios::binary);

	sprintf(name, "%s/file.bin.sp", folder);
	std::ofstream fpfile(name, std::ios::binary);
	uint8bit_v filebin; //filebin is 0 1, record reads from which file
	filebin.n = filebin.a = 0;
	
	{
		for (int i = 0; i < s->allA_id.n; ++i) {
			temp = s->allA_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < s->allT_id.n; ++i) {
			temp = s->allT_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < s->allN_id.n; ++i) {
			temp = s->allN_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < reads->fpA_id.n; ++i) {
			temp = reads->fpA_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < reads->fpT_id.n; ++i) {
			temp = reads->fpT_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < reads->fpN_id.n; ++i) {
			temp = reads->fpN_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < reads->Nfile_id.n; ++i) {
			temp = reads->Nfile_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}

		for (int i = 0; i < reads->singleFile_id.n; ++i) {
			rid = reads->singleFile_id.a[i];
			bseq1_t *seq = &reads->seq[rid];
			strcpy(temp_str, seq->seq);

			for (int j = 0; j < readlen; ++j) {
				DNA_push(singlebin, singleFile, seq_nt4_table[(uint8_t)temp_str[j]]);
			}
		}
		for (int i = 0; i < reads->singleFile_id.n; ++i) {
			temp = reads->singleFile_id.a[i];
			if (temp < half_val) {
				bit_push(filebin, fpfile, 0);
			} else {
				bit_push(filebin, fpfile, 1);
				temp = mpv[temp - half_val];
				fpids.write((char*)&temp, sizeof(uint32_t));
			}
		}
		
		if (singlebin.n > 0) {
			singleFile.write((char*)&singlebin.a, sizeof(uint8_t));
		}
		singleFile.close();
	}

	if (filebin.n > 0) {
		fpfile.write((char*)&filebin.a, sizeof(uint8_t));
	}
	fpfile.close();
	fpids.close();

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

	//
	omp_set_num_threads(n_threads);
	#pragma omp parallel for
	for (int i = 0; i < n_threads; ++i) {
		char name[100]; 
		sprintf(name, "%s/peids.bin.%d", folder, i);
		std::ofstream fpidsi(name, std::ios::binary);

		sprintf(name, "%s/file.bin.%d", folder, i);
		std::ofstream fpfilei(name, std::ios::binary);
		uint8bit_v filebini; //filebin is 0 1, record reads from which file
		filebini.n = filebini.a = 0;

		sprintf(name, "%s/ids.txt.%ld", folder, i);
		FILE *fpin = fopen(name, "r");
		int a, b;
		uint32_t temp;
		while (fscanf(fpin, "%d%d", &a, &b) != EOF) {
			if (a == 0) {
				bit_push(filebini, fpfilei, 0);
			} else {
				bit_push(filebini, fpfilei, 1);
				temp = mpv[b - half_val];
				fpidsi.write((char*)&temp, sizeof(uint32_t));
			}
		}
		fclose(fpin);
		
		if (filebini.n > 0) {
			fpfilei.write((char*)&filebini.a, sizeof(uint8_t));
		}
		fpfilei.close();
		fpidsi.close();
	}

}
