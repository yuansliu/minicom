#include "bbhashdict.h"

char revinttochar[4] = {'A','G','C','T'};//used in bitsettostring
char inttochar[] = {'A','C','G','T'};
char chartorevchar[128];//A-T etc for reverse complement
int chartoint[128];//A-0,C-1 etc. used in updaterefcount
int *dict_start;
int *dict_end; 
int numdict_s = 3;

// std::string outdir = "output/"; ///????????? need mkdir
std::string outdir = std::string(output);
std::string uuid = std::string(uniqid);


uint32_t numreads;
std::bitset<2*readlen> **basemask;//bitset for A,G,C,T at each position 
std::bitset<2*readlen> *positionmask;//bitset for each position (1 at two bits and 0 elsewhere)
std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)
std::bitset<2*readlen> allAbitset, allTbitset;
omp_lock_t allA_lock;
omp_lock_t allT_lock;

u_int64_t boomphf::printPt( pthread_t pt) {
	  unsigned char *ptc = (unsigned char*)(void*)(&pt);
		u_int64_t res =0;
	  for (size_t i=0; i<sizeof(pt); i++) {
		  res+= (unsigned)(ptc[i]);
	  }
		return res;
	}

void bbhashdict::findpos(int64_t *dictidx, uint32_t &startposidx) {
	dictidx[0] = startpos[startposidx];
	auto endidx = startpos[startposidx+1];
	if(read_id[endidx-1] == numreads)//means exactly one read has been removed
		dictidx[1] = endidx-1;
	else if(read_id[endidx-1] == numreads+1)//means two or more reads have been removed (in this case second last entry stores the number of reads left)
		dictidx[1] = dictidx[0] + read_id[endidx-2];
	else
		dictidx[1] = endidx;//no read deleted
	return;
}

void bbhashdict::remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current) {
	auto size = dictidx[1] - dictidx[0];
	if(size == 1)//just one read left in bin
	{
		empty_bin[startposidx] = 1;
		return; //need to keep one read to check during matching
	}
	uint32_t pos = std::lower_bound(read_id+dictidx[0],read_id+dictidx[1],current)-(read_id+dictidx[0]);
	
	std::move(read_id+dictidx[0]+pos+1,read_id+dictidx[1],read_id+dictidx[0]+pos);
	auto endidx = startpos[startposidx+1];
	if(dictidx[1] == endidx)//this is first read to be deleted
		read_id[endidx-1] = numreads;
	else if(read_id[endidx-1] == numreads)//exactly one read has been deleted till now
	{
		read_id[endidx-1] = numreads + 1;
		read_id[endidx-2] = size - 1;//number of reads left in bin
	}
	else//more than two reads have been deleted
		read_id[endidx-2]--;

	return;	
}

std::bitset<2*readlen> stringtobitset(char *s) {
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void bitsettostring(std::bitset<2*readlen> b, char *s) { //only for extend
	unsigned long long ull, rem;
	for(int i = 0; i < 2*readlen/64+1; i++)
	{	
		ull = (b&mask64).to_ullong();
		b>>=64;
		for(int j = 32*i  ; j < 32*i+32 && j < readlen ; j++)
		{
			s[j] = revinttochar[ull%4];	
			ull/=4;
		}
	}
	s[readlen] = '\0';
	return;
}

std::bitset<2*readlen> chartobitset(char *s) {
	std::bitset<2*readlen> b;
	for(int i = 0; i < readlen; i++)
		b |= basemask[i][s[i]];
	return b;
}

void reverse_complement_(char* s, char* s1) {
	for(int j = 0; j < readlen; j++)
		s1[j] = chartorevchar[s[readlen-j-1]];
	return;
}

void generatemasks(std::bitset<2*readlen> *mask, std::bitset<2*readlen> *revmask) {
	for(int i = 0; i < maxmatch; i++)
	{	
		mask[i].reset();
		revmask[i].reset();
		for(int j = 0; j < 2*readlen - 2*i; j++)
			mask[i][j] = 1;
		for(int j = 2*i; j < 2*readlen; j++)
			revmask[i][j] = 1; 	
	}
	return;
}

void generateindexmasks(std::bitset<2*readlen> *mask1, int num_dict) {//masks for dictionary positions
	for(int j = 0; j < num_dict; j++)
		mask1[j].reset();
	for(int j = 0; j < num_dict; j++)
		for(int i = 2*dict_start[j]; i < 2*(dict_end[j]+1); i++)
			mask1[j][i] = 1;
	return;
}

void singleRead2bitset(reads_t *reads, std::bitset<2*readlen> *read, int threshold) {
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	uint32_t i, stop;	
	//doing initial setup and first read
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
	if(tid == omp_get_num_threads()-1)
		stop = numreads;

	char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));

	char *temp_str = (char*)alloca((readlen + 1) * sizeof(char));

	while(i < stop) {
		bseq1_t *seq = &reads->seq[reads->sg.a[i]];
		// char *cur_seq = reads->seq[reads->sg.a[i]].seq;
		read[i] = stringtobitset(seq->seq);

		strcpy(temp_str, seq->seq);
		// put N back to reads;
		uint32_v *n_pos = (uint32_v*) seq->n_pos;
		if (n_pos != NULL) {
			for (int i = 0; i < n_pos->n; ++i) {
				temp_str[n_pos->a[i]] = 'N';
			}
		}

		if (basediff(read[i]^allAbitset) <= threshold) { // put it in all A
			int en_str_len = 0;
			int eq_char_num = 0;
			for (int tj = 0; tj < readlen; ++tj) {
				if ('A' != temp_str[tj]) {
					if (eq_char_num > 0) {
						sprintf(int_str, "%d", eq_char_num);
						for (char *tk = int_str; *tk != '\0'; ++tk) {
							en_str[en_str_len++] = *tk;
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

			if (en_str_len <= readlen*0.4) {
				reads->sg_flag[i] = true;
				omp_set_lock(&allA_lock);
				// fprintf(reads->fpa, "%s\n", en_str);
				kv_push(uint32_t, reads->fpA_id, reads->sg.a[i]);
				omp_unset_lock(&allA_lock);

				/*pthread_mutex_lock(&atfiles_mutex);
				fprintf(atfiles, "%s\n", cur_seq);
				pthread_mutex_unlock(&atfiles_mutex);*/
			}
		} else 
		// if ((read[i]^allTbitset).count() <= threshold ) { // put it in all T
		if (basediff(read[i]^allTbitset) <= threshold ) { // put it in all T
			// sotre as file
			int en_str_len = 0;
			int eq_char_num = 0;
			for (int tj = 0; tj < readlen; ++tj) {
				if ('T' != temp_str[tj]) {
					if (eq_char_num > 0) {
						sprintf(int_str, "%d", eq_char_num);
						for (char *tk = int_str; *tk != '\0'; ++tk) {
							en_str[en_str_len++] = *tk;
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

			if (en_str_len <= readlen*0.4) {
				reads->sg_flag[i] = true;
				omp_set_lock(&allT_lock);
				// fprintf(reads->fpt, "%s\n", en_str);
				kv_push(uint32_t, reads->fpT_id, reads->sg.a[i]);
				omp_unset_lock(&allT_lock);

				/*pthread_mutex_lock(&atfiles_mutex);
				fprintf(atfiles, "%s\n", cur_seq);
				pthread_mutex_unlock(&atfiles_mutex);*/
			}
		}
		i++;
	}
	}
	return;
}

void generateAllATbitset() {
	// fprintf(stderr, "begin generateAllATbitset\n");
	char *str = (char*) alloca((readlen + 1) * sizeof(char));
	// fprintf(stderr, "begin generateAllATbitset 1111\n");
	str[readlen] = '\0';
	for (int i = 0; i < readlen; ++i) str[i] = 'A';
	// fprintf(stderr, "begin generateAllATbitset 2222\n");
	allAbitset = stringtobitset(str);
	// fprintf(stderr, "begin generateAllATbitset 333\n");
	for (int i = 0; i < readlen; ++i) str[i] = 'T';
	// fprintf(stderr, "begin generateAllATbitset 4444\n");
	allTbitset = stringtobitset(str);
	// fprintf(stderr, "after generateAllATbitset\n");

	omp_init_lock(&allA_lock);
	omp_init_lock(&allT_lock);
}

int basediff(std::bitset<2*readlen> bitstr) {
	/*int res = 0;
	for (int i = 0; i < readlen; ++i) {
		if (bitstr[2*i] == 1 || bitstr[2*i + 1] == 1 ) ++ res;
	}
	return res;*/
	return bitstr.count();
}

void freeglobalarrays() {
	for (int i = 0; i < readlen; ++i) {
		free(basemask[i]);
	}
	free(positionmask);
	free(basemask);

	delete[] dict_start;
	delete[] dict_end;
}

