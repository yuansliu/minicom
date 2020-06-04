#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
// #include "config.h"

void reverse_complement(char *seq, int seq_len) {
	char rc_tablle[256];
	rc_tablle['A'] = 'T';
	rc_tablle['T'] = 'A';
	rc_tablle['G'] = 'C';
	rc_tablle['C'] = 'G';
	rc_tablle['N'] = 'N';
	// uint32_v *n_pos = (uint32_v*)calloc(1, sizeof(uint32_v));
	char *res = (char*)alloca((seq_len + 1) * sizeof(char));
	for (int i = seq_len - 1, j = 0; i >= 0; --i, ++j) {
		res[j] = rc_tablle[seq[i]];
	}
	res[seq_len] = '\0';
	strcpy(seq, res);
}

const char invert_code_rule[4] = {'A', 'C', 'G', 'T'}; //decoding rule //A - 0; C - 1; G - 2; T - 3;

struct DNAPOOL {
	int n, a[4];
};

struct BITPOOL {
	int n, a[8];
};

int readlen, n_threads, half_val, *per;
char *folder;

int getDir(std::ifstream& fpdir, BITPOOL& curdirpool) {
	if (curdirpool.n >= 8) {
		uint8_t dirbin;
		fpdir.read((char*)&dirbin, sizeof(uint8_t));
		for (int i = 0; i < 8; ++i) {
		// for (int i = 7; i >= 0; --i) {
			curdirpool.a[i] = dirbin&1;
			dirbin >>= 1;
		}
		curdirpool.n = 0;
	}
	int res = curdirpool.a[curdirpool.n];
	++curdirpool.n;
	return res;
}

int getDNAcode(std::ifstream& fpif, DNAPOOL& curdnapool) {
	// fprintf(stderr, "yyyyyy\n");
	if (curdnapool.n >= 4) {
		uint8_t bin;
		// fprintf(stderr, "zzzz\n");
		fpif.read((char*)&bin, sizeof(uint8_t));
		// fprintf(stderr, "%d\n", bin);
		if (fpif.eof()) return -1;

		for (int i = 0; i < 4; ++i) {
			curdnapool.a[i] = bin&3;
			bin >>= 2;
		}
		curdnapool.n = 0;
	}
	int res = curdnapool.a[curdnapool.n];
	++curdnapool.n;
	return res;
}

int getFileBit(std::ifstream& fpfile, BITPOOL& curfilepool) {
	if (curfilepool.n >= 8) {
		uint8_t filebin;
		fpfile.read((char*)&filebin, sizeof(uint8_t));
		for (int i = 0; i < 8; ++i) {
		// for (int i = 7; i >= 0; --i) {
			curfilepool.a[i] = filebin&1;
			filebin >>= 1;
		}
		curfilepool.n = 0;
	}
	int res = curfilepool.a[curfilepool.n];
	++curfilepool.n;
	return res;
}

void getRef(std::ifstream& fpref, DNAPOOL& curdnapool, char *ref, int len) {
	int code, ref_len;
	for (ref_len = strlen(ref); ref_len < len; ++ref_len) {
		code = getDNAcode(fpref, curdnapool);
		if (code < 0) break;
		ref[ref_len] = invert_code_rule[code];
	}
	ref[ref_len] = '\0';
}

uint32_t n_seq;
typedef struct { size_t n, m; char *a; } char_v;
char **reads;
FILE *fp_order;

// int leftcnt;

void decompress_order(int nth) {
	FILE *fpdif_pos = NULL, *fpdif_char = NULL;
	char name[100]; 
	sprintf(name, "%s/ref.bin.%d", folder, nth);
	std::ifstream fpref(name, std::ios::binary);

	sprintf(name, "%s/beg_pos.bin.%d", folder, nth);
	std::ifstream fppos(name, std::ios::binary);
	
	sprintf(name, "%s/dir.bin.%d", folder, nth);
	std::ifstream fpdir(name, std::ios::binary);
	
	sprintf(name, "%s/dif_char.txt.%d", folder, nth);
	fpdif_char = fopen(name, "r");
	if (fpdif_char == NULL) {
		fprintf(stderr, "Open fpdif_char file is error!\n");
	}

	sprintf(name, "%s/ids.bin.%d", folder, nth);
	std::ifstream fpids(name, std::ios::binary);

	// int readlen = 100;
	char *ref = (char*)calloc(1<<25, sizeof(char));
	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	char *temp_char = (char*)alloca(readlen * sizeof(char)), *temp;
	char *sub_str = (char*)alloca((readlen + 1) * sizeof(char));
	int pos, pre_pos, len, eq_char_num, temp_char_len;
	int dir;
	// fprintf(stderr, "xxxx\n");
	uint32_t seq_num;

	int times = 0;
	// bool flag = false;
	uint8_t seq_num_bin[4];
	uint16_t posbin;

	DNAPOOL curdnapool;
	curdnapool.n = 4;

	BITPOOL curdirpool;
	curdirpool.n = 8;

	while (true) {
		// while (fscanf(fppos, "%d", &seq_num) != EOF) {
		// fprintf(stderr, "yyy\n");
		fppos.read((char*)&seq_num, sizeof(uint32_t));
		if (fppos.eof()) break;
		// if (flag) fprintf(stderr, "seq_num: %d\n", seq_num);

		pre_pos = 0;
		// fscanf(fpref, "%s", ref);
		// getRef(ref, curdnapool, len);
		ref[0] = '\0';

		// if (flag) fprintf(stderr, "ref: %s\n", ref);
		// if (strcmp(ref, "TCAACTTTCGATGGCAGTCGCCGTGCCTACCATGGAGACCACGGGTAACGGAGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTAACACCACCAAGGAA") ==0) {
		// 	flag = true;
		// }
		uint32_t idx, pre_idx = 0;
		for (uint32_t tt = 0; tt < seq_num; ++tt) {
			++times;
			// fscanf(fppos, "%d", &pos);
			fppos.read((char*)&posbin, sizeof(uint16_t));

			// if (flag) fprintf(stderr, "pos: %d\n", pos);
			// if (flag) fprintf(stderr, "pos --: %d\n", pos);
			pos = posbin;

			fpids.read((char*)&idx, sizeof(uint32_t));
			if (pos == 0) {
				idx += pre_idx;
			} 
			pre_idx = idx;

			pos += pre_pos;
			pre_pos = pos;

			getRef(fpref, curdnapool, ref, pos + readlen);

			// fscanf(fpdir, "%c", &dir);
			dir = getDir(fpdir, curdirpool);

			fscanf(fpdif_char, "%s", temp_char); 
			// temp_char_len = strlen(temp_char);
			// if (flag) fprintf(stderr, "%s\n", temp_char);
			temp_char_len = strlen(temp_char);
			len = 0;
			eq_char_num = 0;
			for (int i = 0; i < temp_char_len; ++i) {
				if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
					if (eq_char_num > 0) {
						for (int j = 0; j < eq_char_num; ++j) {
							seq[len] = ref[pos + len];
							++ len;
						}
						eq_char_num = 0;
					}
					seq[len++] = temp_char[i];
				} else { // '0' <= temp_char[i] <= '9'
					eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
				}
			}
			for (; len < readlen; ++len) {
				seq[len] = ref[pos + len];
			}
			seq[len] = '\0';
			if (dir) {
				reverse_complement(seq, readlen);
			}
			// fprintf(result, "%s\n", seq);
			// fprintf(fp_order, "%u\n", idx);
			strcpy(reads[idx], seq);
			// if (flag) fprintf(stderr, "%s\n", seq);
			// if (times > 60 ) exit(0);
			// if (flag && seq_num == 0) {
			// if (flag) {
			// 	exit(0);
			// }
		}
	}
	free(ref);

	fpref.close();
	fppos.close();
	fpdir.close();

	fclose(fpdif_char);
}

void decomp_single_order() {
	// std::ifstream fpsingle("%s/single.seq");
	char name[100]; 
	sprintf(name, "%s/single.seq", folder);
	std::ifstream fpsingle(name, std::ios::binary);
	sprintf(name, "%s/singleFile.ids.bin", folder);
	std::ifstream fpsingleFileids(name, std::ios::binary);
	// FILE *fpSingleOut = fopen("decompfiles/single_dec.fasta", "w");
	// FILE *fpSingleOut = fopen("%s/single.fasta", "w");
	DNAPOOL curdnapool;
	curdnapool.n = 4;
	int code;
	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	bool flag = true;
	// fprintf(stderr, "xxxx\n");
	// fprintf(stderr, "xxxxxx\n");
	uint32_t idx, pre_idx = 0;
	// int xxxx = 0;
	while (flag) {
		// fprintf(stderr, "xxxx\n");
		for (int i = 0; flag && i < readlen; ++i) {
			code = getDNAcode(fpsingle, curdnapool);
			// fprintf(stderr, "%d\n", code);
			if (code < 0) {
				flag = false;
				break;
			}
			seq[i] = invert_code_rule[code];
		}
		if (!flag) break;
		seq[readlen] = '\0';
		fpsingleFileids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		// fprintf(stderr, "idx: %d\n", idx);
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
		// if (xxxx ++ > 10) exit(0);
		// fprintf(fpSingleOut, "%s\n", seq);
	}
	// fclose(fpSingleOut);
	fpsingleFileids.close();

	// fprintf(stderr, "xxxxxx\n");
	sprintf(name, "%s/single_N.seq", folder);
	FILE *fp = fopen(name, "r");
	sprintf(name, "%s/Nfile.ids.bin", folder);
	std::ifstream Nfileids(name, std::ios::binary);
	pre_idx = 0;
	while (fscanf(fp, "%s", seq) != EOF) {
		Nfileids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
	}
	fclose(fp);
	Nfileids.close();
	// []
}

void decomp_AATTNN_order() {
	char name[100]; 
	sprintf(name, "%s/info.txt", folder);
	FILE *fp_in = fopen(name, "r");
	// FILE *fp_out = fopen("decompfiles/aatt.fasta", "w");
	fscanf(fp_in, "%d%d", &readlen, &n_threads);

	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	char *temp_char = (char*)alloca(readlen * sizeof(char)), *temp;
	int len, eq_char_num, temp_char_len;

	int a_num, t_num, n_num;
	fscanf(fp_in, "%d%d%d", &a_num, &t_num, &n_num);
	fscanf(fp_in, "%u", &n_seq);
	fclose(fp_in);

    reads = (char**)calloc(n_seq + 10, sizeof(char*));
    for (int i = 0; i < n_seq; ++i) {
	    reads[i] = (char*) calloc(readlen + 1, sizeof(char));
	    reads[i][0] = '\0';
	}
	uint32_t idx, pre_idx = 0;

	sprintf(name, "%s/allA.ids.bin", folder);
	std::ifstream allAids(name, std::ios::binary);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'A';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < a_num; ++i) {
		// fprintf(fp_out, "%s\n", seq);
		allAids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
	}
	allAids.close();
	
	sprintf(name, "%s/allT.ids.bin", folder);
	std::ifstream allTids(name, std::ios::binary);
	pre_idx = 0;
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'T';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < t_num; ++i) {
		// fprintf(fp_out, "%s\n", seq);
		allTids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
	}
	allTids.close();

	sprintf(name, "%s/allN.ids.bin", folder);
	std::ifstream allNids(name, std::ios::binary);
	pre_idx = 0;
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'N';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < n_num; ++i) {
		// fprintf(fp_out, "%s\n", seq);
		allNids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
	}
	allNids.close();

	sprintf(name, "%s/AA.ids.bin", folder);
	std::ifstream fpAids(name, std::ios::binary);
	pre_idx = 0;

	sprintf(name, "%s/AA.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'A';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'A';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';

		fpAids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);

		// fprintf(fp_out, "%s\n", seq);
	}
	fclose(fp_in);
	fpAids.close();

	sprintf(name, "%s/TT.ids.bin", folder);
	std::ifstream fpTids(name, std::ios::binary);
	pre_idx = 0;
	sprintf(name, "%s/TT.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'T';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'T';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';

		fpTids.read((char*)&idx, sizeof(uint32_t));
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
		// fprintf(fp_out, "%s\n", seq);
	}
	fclose(fp_in);
	fpTids.close();

	sprintf(name, "%s/NN.ids.bin", folder);
	std::ifstream fpNids(name, std::ios::binary);
	pre_idx = 0;
	sprintf(name, "%s/NN.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'N';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			// if(temp_char_len == 3) fprintf(stderr, "len: %d\n", len);
			seq[len] = 'N';
		}
		// if (flag) 
		seq[len] = '\0';

		fpNids.read((char*)&idx, sizeof(uint32_t));
		// fprintf(stderr, "decomp_AATTNN\n");
		idx += pre_idx;
		pre_idx = idx;
		// fprintf(fp_order, "%u\n", idx);
		strcpy(reads[idx], seq);
	}
	fclose(fp_in);
	fpNids.close();
	// fclose(fp_out);
}

void decompress_nonorder(int nth, FILE *result) {
	FILE *fpdif_pos = NULL, *fpdif_char = NULL;
	char name[100]; 
	sprintf(name, "%s/ref.bin.%d", folder, nth);
	std::ifstream fpref(name, std::ios::binary);

	sprintf(name, "%s/beg_pos.bin.%d", folder, nth);
	std::ifstream fppos(name, std::ios::binary);
	
	sprintf(name, "%s/dir.bin.%d", folder, nth);
	std::ifstream fpdir(name, std::ios::binary);
	
	sprintf(name, "%s/dif_char.txt.%d", folder, nth);

	fpdif_char = fopen(name, "r");
	if (fpdif_char == NULL) {
		fprintf(stderr, "Open fpdif_char file is error!\n");
	}

	// int readlen = 100;
	char *ref = (char*)calloc(1<<25, sizeof(char));
	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	char *temp_char = (char*)alloca(readlen * sizeof(char)), *temp;
	char *sub_str = (char*)alloca((readlen + 1) * sizeof(char));
	int pos, pre_pos, len, eq_char_num, temp_char_len;
	uint32_t seq_num;
	int dir;
	// fprintf(stderr, "xxxx\n");
	
	int times = 0;
	bool flag = false;
	uint16_t posbin;

	DNAPOOL curdnapool;
	curdnapool.n = 4;

	BITPOOL curdirpool;
	curdirpool.n = 8;

	while (true) {
		// while (fscanf(fppos, "%d", &seq_num) != EOF) {
		// fprintf(stderr, "yyy\n");
		fppos.read((char*)&seq_num, sizeof(uint32_t));
		if (fppos.eof()) break;

		// if (flag) fprintf(stderr, "seq_num: %d\n", seq_num);

		pre_pos = 0;
		// fscanf(fpref, "%s", ref);
		// getRef(ref, curdnapool, len);
		ref[0] = '\0';

		// if (flag) fprintf(stderr, "ref: %s\n", ref);
		// if (strcmp(ref, "TCAACTTTCGATGGCAGTCGCCGTGCCTACCATGGAGACCACGGGTAACGGAGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTAACACCACCAAGGAA") ==0) {
		// 	flag = true;
		// }
		for (uint32_t tt = 0; tt < seq_num; ++tt) {
			++times;
			// fscanf(fppos, "%d", &pos);
			fppos.read((char*)&posbin, sizeof(uint16_t));

			// if (flag) fprintf(stderr, "pos: %d\n", pos);
			// if (flag) fprintf(stderr, "pos --: %d\n", pos);
			pos = posbin;
			pos += pre_pos;
			pre_pos = pos;

			getRef(fpref, curdnapool, ref, pos + readlen);

			// fscanf(fpdir, "%c", &dir);
			dir = getDir(fpdir, curdirpool);

			fscanf(fpdif_char, "%s", temp_char); 
			// temp_char_len = strlen(temp_char);
			// if (flag) fprintf(stderr, "%s\n", temp_char);
			temp_char_len = strlen(temp_char);
			len = 0;
			eq_char_num = 0;
			for (int i = 0; i < temp_char_len; ++i) {
				if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
					if (eq_char_num > 0) {
						for (int j = 0; j < eq_char_num; ++j) {
							seq[len] = ref[pos + len];
							++ len;
						}
						eq_char_num = 0;
					}
					seq[len++] = temp_char[i];
				} else { // '0' <= temp_char[i] <= '9'
					eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
				}
			}
			for (; len < readlen; ++len) {
				seq[len] = ref[pos + len];
			}
			seq[len] = '\0';
			if (dir) {
				reverse_complement(seq, readlen);
			}
			fprintf(result, "%s\n", seq);
			// if (flag) fprintf(stderr, "%s\n", seq);
			// if (times > 60 ) exit(0);
			// if (flag && seq_num == 0) {
			// if (flag) {
			// 	exit(0);
			// }
		}
	}
	free(ref);

	fpref.close();
	fppos.close();
	fpdir.close();

	fclose(fpdif_char);
}

void decomp_single_nonorder() {
	// std::ifstream fpsingle("%s/single.seq");
	char name[100]; 
	sprintf(name, "%s/single.seq", folder);
	std::ifstream fpsingle(name, std::ios::binary);
	
	sprintf(name, "%s/single_dec.fasta", folder);
	FILE *fpSingleOut = fopen(name, "w");
	// FILE *fpSingleOut = fopen("%s/single.fasta", "w");

	DNAPOOL curdnapool;
	curdnapool.n = 4;
	int code;
	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	bool flag = true;
	// fprintf(stderr, "xxxx\n");
	while (flag) {
		// fprintf(stderr, "xxxx\n");
		for (int i = 0; flag && i < readlen; ++i) {
			code = getDNAcode(fpsingle, curdnapool);
			// fprintf(stderr, "%d\n", code);
			if (code < 0) {
				flag = false;
				break;
			}
			seq[i] = invert_code_rule[code];
		}
		if (!flag) break;
		seq[readlen] = '\0';
		fprintf(fpSingleOut, "%s\n", seq);
	}
	fclose(fpSingleOut);
	// []
}

void decomp_AATTNN_nonorder() {
	char name[100]; 
	sprintf(name, "%s/info.txt", folder);
	FILE *fp_in = fopen(name, "r");
	
	sprintf(name, "%s/aatt.fasta", folder);
	FILE *fp_out = fopen(name, "w");

	fscanf(fp_in, "%d%d", &readlen, &n_threads);

	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	char *temp_char = (char*)alloca(readlen * sizeof(char)), *temp;
	int len, eq_char_num, temp_char_len;

	int a_num, t_num, n_num;
	fscanf(fp_in, "%d", &a_num);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'A';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < a_num; ++i) {
		fprintf(fp_out, "%s\n", seq);
	}

	fscanf(fp_in, "%d", &t_num);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'T';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < t_num; ++i) {
		fprintf(fp_out, "%s\n", seq);
	}

	fscanf(fp_in, "%d", &n_num);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'N';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < n_num; ++i) {
		fprintf(fp_out, "%s\n", seq);
	}
	fclose(fp_in);

	sprintf(name, "%s/AA.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'A';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'A';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';
		fprintf(fp_out, "%s\n", seq);
	}
	fclose(fp_in);

	sprintf(name, "%s/TT.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'T';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'T';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';
		fprintf(fp_out, "%s\n", seq);
	}
	fclose(fp_in);

	sprintf(name, "%s/NN.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'N';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'N';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';
		fprintf(fp_out, "%s\n", seq);
	}
	fclose(fp_in);

	fclose(fp_out);
}

void decompress_pe(int nth, FILE *result) {
	FILE *fpdif_pos = NULL, *fpdif_char = NULL;
	char name[100]; 
	sprintf(name, "%s/ref.bin.%d", folder, nth);
	std::ifstream fpref(name, std::ios::binary);

	sprintf(name, "%s/beg_pos.bin.%d", folder, nth);
	std::ifstream fppos(name, std::ios::binary);
	
	sprintf(name, "%s/dir.bin.%d", folder, nth);
	std::ifstream fpdir(name, std::ios::binary);
	
	sprintf(name, "%s/dif_char.txt.%d", folder, nth);

	fpdif_char = fopen(name, "r");
	if (fpdif_char == NULL) {
		fprintf(stderr, "Open fpdif_char file is error!\n");
	}
	
	sprintf(name, "%s/file.bin.%d", folder, nth);
	std::ifstream fpfile(name, std::ios::binary);
	BITPOOL curfilepool;
	curfilepool.n = 8;

	sprintf(name, "%s/peids.bin.%d", folder, nth);
	std::ifstream fpids(name, std::ios::binary);

	// int readlen = 100;
	char *ref = (char*)calloc(1<<25, sizeof(char));
	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	char *temp_char = (char*)alloca(readlen * sizeof(char)), *temp;
	char *sub_str = (char*)alloca((readlen + 1) * sizeof(char));
	int pos, pre_pos, len, eq_char_num, temp_char_len;
	uint32_t seq_num;
	int dir;
	// fprintf(stderr, "xxxx\n");
	
	int times = 0;
	// bool flag = false;
	uint8_t seq_num_bin[4];
	uint16_t posbin;

	DNAPOOL curdnapool;
	curdnapool.n = 4;

	BITPOOL curdirpool;
	curdirpool.n = 8;

	int filevalue;
	uint32_t idx;

	while (true) {
		// while (fscanf(fppos, "%d", &seq_num) != EOF) {
		// fprintf(stderr, "yyy\n");
		fppos.read((char*)&seq_num, sizeof(uint32_t));
		if (fppos.eof()) break;

		// if (flag) fprintf(stderr, "seq_num: %d\n", seq_num);

		pre_pos = 0;
		// fscanf(fpref, "%s", ref);
		// getRef(ref, curdnapool, len);
		ref[0] = '\0';

		// if (flag) fprintf(stderr, "ref: %s\n", ref);
		// if (strcmp(ref, "TCAACTTTCGATGGCAGTCGCCGTGCCTACCATGGAGACCACGGGTAACGGAGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTAACACCACCAAGGAA") ==0) {
		// 	flag = true;
		// }
		for (uint32_t tt = 0; tt < seq_num; ++tt) {
			++times;
			// fscanf(fppos, "%d", &pos);
			fppos.read((char*)&posbin, sizeof(uint16_t));

			// if (flag) fprintf(stderr, "pos: %d\n", pos);
			// if (flag) fprintf(stderr, "pos --: %d\n", pos);
			pos = posbin;
			pos += pre_pos;
			pre_pos = pos;

			getRef(fpref, curdnapool, ref, pos + readlen);

			// fscanf(fpdir, "%c", &dir);
			dir = getDir(fpdir, curdirpool);

			fscanf(fpdif_char, "%s", temp_char); 
			// temp_char_len = strlen(temp_char);
			// if (flag) fprintf(stderr, "%s\n", temp_char);
			temp_char_len = strlen(temp_char);
			len = 0;
			eq_char_num = 0;
			for (int i = 0; i < temp_char_len; ++i) {
				if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
					if (eq_char_num > 0) {
						for (int j = 0; j < eq_char_num; ++j) {
							seq[len] = ref[pos + len];
							++ len;
						}
						eq_char_num = 0;
					}
					seq[len++] = temp_char[i];
				} else { // '0' <= temp_char[i] <= '9'
					eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
				}
			}
			for (; len < readlen; ++len) {
				seq[len] = ref[pos + len];
			}
			seq[len] = '\0';
			if (dir) {
				reverse_complement(seq, readlen);
			}

			filevalue = getFileBit(fpfile, curfilepool);
			if (filevalue) { //1 right ends
				fpids.read((char*)&idx, sizeof(uint32_t));
				strcpy(reads[idx], seq);
			} else { //0 left ends
				fprintf(result, "%s\n", seq);
				// ++leftcnt;
			}

			// if (flag) fprintf(stderr, "%s\n", seq);
			// if (times > 60 ) exit(0);
			// if (flag && seq_num == 0) {
			// if (flag) {
			// 	exit(0);
			// }
		}
	}
	free(ref);

	fpref.close();
	fppos.close();
	fpdir.close();
	fpids.close();
	fpfile.close();
	fclose(fpdif_char);
}

/*void decomp_single_pe() {
	char name[100]; 
	sprintf(name, "%s/single.seq", folder);
	std::ifstream fpsingle(name, std::ios::binary);
	
	sprintf(name, "%s/single_dec.fasta", folder);
	FILE *fpSingleOut = fopen(name, "w");
	// FILE *fpSingleOut = fopen("%s/single.fasta", "w");

	DNAPOOL curdnapool;
	curdnapool.n = 4;
	int code;
	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	bool flag = true;
	int filevalue;
	uint32_t idx;
	// fprintf(stderr, "xxxx\n");
	while (flag) {
		// fprintf(stderr, "xxxx\n");
		for (int i = 0; flag && i < readlen; ++i) {
			code = getDNAcode(fpsingle, curdnapool);
			// fprintf(stderr, "%d\n", code);
			if (code < 0) {
				flag = false;
				break;
			}
			seq[i] = invert_code_rule[code];
		}
		if (!flag) break;
		seq[readlen] = '\0';
		
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			// fprintf(fp_out, "%s\n", seq);
			fprintf(fpSingleOut, "%s\n", seq);
			// ++leftcnt;
		}
	}
	fclose(fpSingleOut);
}*/

void decomp_AATTNN_pe() {
	// leftcnt = 0;
	char name[100]; 
	sprintf(name, "%s/info.txt", folder);
	FILE *fp_in = fopen(name, "r");	

	fscanf(fp_in, "%d%d", &readlen, &n_threads);
	fscanf(fp_in, "%d", &half_val);

	// fprintf(stderr, "readlen: %d; n_threads: %d\n", readlen, n_threads);
	// fprintf(stderr, "half_val: %d\n", half_val);

	per = (int*)calloc(half_val, sizeof(int));
	memset(per, 0, half_val * sizeof(int));

	reads = (char**)calloc(half_val + 10, sizeof(char*));
    for (int i = 0; i < half_val; ++i) {
	    reads[i] = (char*) calloc(readlen + 1, sizeof(char));
	    reads[i][0] = '\0';
	}
	sprintf(name, "%s/file.bin.sp", folder);
	std::ifstream fpfile(name, std::ios::binary);
	BITPOOL curfilepool;
	curfilepool.n = 8;

	sprintf(name, "%s/peids.bin.sp", folder);
	std::ifstream fpids(name, std::ios::binary);

	//---------
	sprintf(name, "%s/aatt.fasta", folder);
	FILE *fp_out = fopen(name, "w");

	char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	char *temp_char = (char*)alloca(readlen * sizeof(char)), *temp;
	int len, eq_char_num, temp_char_len;

	int a_num, t_num, n_num;
	int filevalue;
	uint32_t idx;

	fscanf(fp_in, "%d", &a_num);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'A';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < a_num; ++i) {
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}

	fscanf(fp_in, "%d", &t_num);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'T';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < t_num; ++i) {
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}

	fscanf(fp_in, "%d", &n_num);
	for (int i = 0; i < readlen; ++i) {
		seq[i] = 'N';
	}
	seq[readlen] = '\0';
	for(int i = 0; i < n_num; ++i) {
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}
	fclose(fp_in);

	//----------
	sprintf(name, "%s/AA.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'A';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'A';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';
		// fprintf(fp_out, "%s\n", seq);
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}
	fclose(fp_in);

	sprintf(name, "%s/TT.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'T';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'T';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';
		// fprintf(fp_out, "%s\n", seq);
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}
	fclose(fp_in);

	sprintf(name, "%s/NN.txt", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", temp_char) != EOF) {
		temp_char_len = strlen(temp_char);
		len = 0;
		eq_char_num = 0;
		for (int i = 0; i < temp_char_len; ++i) {
			if (temp_char[i] >= 'A' && temp_char[i] <= 'Z') {
				if (eq_char_num > 0) {
					for (int j = 0; j < eq_char_num; ++j) {
						seq[len] = 'N';
						++ len;
					}
					eq_char_num = 0;
				}
				seq[len++] = temp_char[i];
			} else {
				eq_char_num = eq_char_num * 10 + temp_char[i] - '0';
			}
		}
		for (; len < readlen; ++len) {
			seq[len] = 'N';
		}
		// if (flag) fprintf(stderr, "len: %d\n", len);
		seq[len] = '\0';
		// fprintf(fp_out, "%s\n", seq);
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}
	fclose(fp_in);

	sprintf(name, "%s/single_N.seq", folder);
	fp_in = fopen(name, "r");
	while (fscanf(fp_in, "%s", seq) != EOF) {
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			fprintf(fp_out, "%s\n", seq);
			// ++leftcnt;
		}
	}

	fclose(fp_out);

	// char name[100]; 
	sprintf(name, "%s/single.seq", folder);
	std::ifstream fpsingle(name, std::ios::binary);
	
	sprintf(name, "%s/single_dec.fasta", folder);
	FILE *fpSingleOut = fopen(name, "w");
	// FILE *fpSingleOut = fopen("%s/single.fasta", "w");

	DNAPOOL curdnapool;
	curdnapool.n = 4;
	int code;
	// char *seq = (char*)alloca((readlen + 1) * sizeof(char));
	bool flag = true;
	// int filevalue;
	// uint32_t idx;
	// fprintf(stderr, "xxxx\n");
	while (flag) {
		// fprintf(stderr, "xxxx\n");
		for (int i = 0; flag && i < readlen; ++i) {
			code = getDNAcode(fpsingle, curdnapool);
			// fprintf(stderr, "%d\n", code);
			if (code < 0) {
				flag = false;
				break;
			}
			seq[i] = invert_code_rule[code];
		}
		if (!flag) break;
		seq[readlen] = '\0';
		
		filevalue = getFileBit(fpfile, curfilepool);
		if (filevalue) { //1 right ends
			fpids.read((char*)&idx, sizeof(uint32_t));
			strcpy(reads[idx], seq);
		} else { //0 left ends
			// fprintf(fp_out, "%s\n", seq);
			fprintf(fpSingleOut, "%s\n", seq);
			// ++leftcnt;
		}
	}
	fclose(fpSingleOut);

	fpids.close();
	fpfile.close();
}

int main(int argc, char *argv[]) {
	// fp_order = fopen("order.idx", "w");
	// for (int i = 0; i < argc; ++i) {
	// 	fprintf(stderr, "%s\n", argv[i]);
	// }
	folder = argv[1];
	int n_thc = 0;

	// argv[5]
	for (int i = 0; i < strlen(argv[5]); ++i) {
		n_thc = n_thc*10 + argv[5][i] - '0';
	}
	// fprintf(stderr, "number of threads: %d\n", n_thc);
	omp_set_num_threads(n_thc);

	// fprintf(stderr, "%s\n", folder);
	if (argv[3][0] == 'f') { //no pe
		if (argv[4][0] == 't') { //order
			decomp_AATTNN_order();
			// fprintf(stderr, "decomp_AATTNN_order()...over\n");
			decomp_single_order();
			// fprintf(stderr, "decomp_single_order()...over\n");

			#pragma omp parallel for
			for (int n = 0; n < n_threads; ++n) {
				decompress_order(n);
			}
			FILE *fpresult;
			// char name[100]; 
			// sprintf(name, "decom_final.reads");
			fpresult = fopen(argv[2], "w");
			for (int i = 0; i < n_seq; ++i) {
				/*if (strlen(reads[i]) == 0) {
					fprintf(stderr, "---000000\n");
				}*/
				fprintf(fpresult, "%s\n", reads[i]);
			}
			fclose(fpresult);
		} else {
			decomp_AATTNN_nonorder();
			// fprintf(stderr, "decomp_AATTNN_nonorder()...over\n");
			decomp_single_nonorder();
			// fprintf(stderr, "decomp_single_nonorder()...over\n");
			
			// fprintf(stderr, "n_threads: %d\n", n_threads);

			#pragma omp parallel for
			for (int n = 0; n < n_threads; ++n) {
				FILE *fpresult;
				char name[100]; 
				// sprintf(name, "%s/result_%d.seq", n);
				sprintf(name, "%s/result_%d.seq", folder, n);
				fpresult = fopen(name, "w");
				decompress_nonorder(n, fpresult);
				fclose(fpresult);
			}
		}
	} else { //pe
		decomp_AATTNN_pe();
		// fprintf(stderr, "decomp_AATTNN_pe()...over\n");
		// decomp_single_pe();
		// fprintf(stderr, "decomp_single_pe()...over\n");

		#pragma omp parallel for
		for (int n = 0; n < n_threads; ++n) {
			FILE *fpresult;
			char name[100]; 
			sprintf(name, "%s/result_%d.seq", folder, n);
			fpresult = fopen(name, "w");
			decompress_pe(n, fpresult);
			fclose(fpresult);
		}

		char name[100]; 
		sprintf(name, "%s/catsh.sh", folder);
		FILE *fpcatsh = fopen(name, "w");
		fprintf(fpcatsh, "cat %s/aatt.fasta %s/single_dec.fasta", folder, folder);

		for (int n = 0; n < n_threads; ++n) {
			sprintf(name, "%s/result_%d.seq", folder, n);
			fprintf(fpcatsh, " %s", name);
		}
		fprintf(fpcatsh, " > %s", argv[2]);
		fclose(fpcatsh);
		// cat $decomp/aatt.fasta $decomp/single_dec.fasta $decomp/result_*.seq > $result
		// fprintf(stderr, "leftcnt: %d\n", leftcnt);
		FILE *fpresult = fopen(argv[6], "w");
		for (int i = 0; i < half_val; ++i) {
			fprintf(fpresult, "%s\n", reads[i]);
		}
		fclose(fpresult);
	}
}
