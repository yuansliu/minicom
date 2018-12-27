#include "bbhashdict.h"

void constructdictionary_realign(std::bitset<2*readlen> *read, bbhashdict *dict) {
	// std::bitset<2*readlen> mask[numdict_s];
	std::bitset<2*readlen> *mask = (std::bitset<2*readlen>*)alloca(numdict_s * sizeof(std::bitset<2*readlen>));
	generateindexmasks(mask, numdict_s);

	double mm_realtime0;

	// fprintf(stderr, "begin constructdictionary_realign()...\n");
	for(int j = 0; j < numdict_s; j++)
	{
		uint64_t *ull = new uint64_t[numreads];
		// fprintf(stderr, "*** parallel begin ***\n");

		mm_realtime0 = realtime();

		#pragma omp parallel
		{
		std::bitset<2*readlen> b;
		int tid = omp_get_thread_num();
		std::ofstream foutkey(outdir+uuid+std::string("keys.bin.")+std::to_string(tid),std::ios::binary);
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		//compute keys and write to file and store in ull
		for(; i < stop; i++)
		{
			b = read[i]&mask[j];
			ull[i] = (b>>2*dict_start[j]).to_ullong();
			foutkey.write((char*)&ull[i], sizeof(uint64_t));
		}
		foutkey.close();
		}//parallel end
		// fprintf(stderr, "*** parallel end ***\n");
		// if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
		// mm_realtime0 = realtime();

		//deduplicating ull
		std::sort(ull, ull+numreads);
		uint32_t k = 0;
		for (uint32_t i = 1; i < numreads; i++) 
		        if (ull[i] != ull[k])         
				ull[++k] = ull[i];
		dict[j].numkeys = k+1;
		//construct mphf
		auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(ull), static_cast<const u_int64_t*>(ull+dict[j].numkeys));
		double gammaFactor = 5.0;//balance between speed and memory
		dict[j].bphf = new boomphf::mphf<u_int64_t,hasher_t>(dict[j].numkeys,data_iterator,n_threads,gammaFactor,true,false);
	
		delete[] ull;

		// if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
		// mm_realtime0 = realtime();

		//compute hashes for all reads
		#pragma omp parallel
		{ 
		int tid = omp_get_thread_num();	
		std::ifstream finkey(outdir+uuid+std::string("keys.bin.")+std::to_string(tid),std::ios::binary);
		std::ofstream fouthash(outdir+uuid+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
		uint64_t currentkey, currenthash;
		uint32_t i, stop;
		i = uint64_t(tid)*numreads/omp_get_num_threads();
		stop = uint64_t(tid+1)*numreads/omp_get_num_threads();
		if(tid == omp_get_num_threads()-1)
			stop = numreads;
		for(; i < stop; i++)
		{
			finkey.read((char*)&currentkey, sizeof(uint64_t));
			currenthash = (dict[j].bphf)->lookup(currentkey);
			fouthash.write((char*)&currenthash, sizeof(uint64_t));
		}
		finkey.close();
		remove((outdir+uuid+std::string("keys.bin.")+std::to_string(tid)).c_str());
		fouthash.close();
		}//parallel end
	}
	// fprintf(stderr, "middle constructdictionary_realign()...\n");
		// if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
		// mm_realtime0 = realtime();

	// omp_set_num_threads(std::min(numdict_s,n_threads));
	#pragma omp parallel
	{
	#pragma omp for	
	for (int j = 0; j < numdict_s; ++j) {
		//fill startpos by first storing numbers and then doing cumulative sum
		dict[j].startpos = new uint32_t[dict[j].numkeys+1]();//1 extra to store end pos of last key
		uint64_t currenthash;
		for(int tid = 0; tid < n_threads; tid++)
		{
			std::ifstream finhash(outdir+uuid+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].startpos[currenthash+1]++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
		}
	
		dict[j].empty_bin = new bool[dict[j].numkeys]();
		for(uint32_t i = 1; i < dict[j].numkeys; i++)
			dict[j].startpos[i] =  dict[j].startpos[i] +  dict[j].startpos[i-1];
	
		// if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
		// mm_realtime0 = realtime();

		//insert elements in the dict array
		dict[j].read_id = new uint32_t[numreads];
		uint32_t i = 0;
		for(int tid = 0; tid < n_threads; tid++)
		{
			std::ifstream finhash(outdir+uuid+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j),std::ios::binary);
			finhash.read((char*)&currenthash,sizeof(uint64_t));
			while(!finhash.eof())
			{
				dict[j].read_id[dict[j].startpos[currenthash]++] = i;
				i++;
				finhash.read((char*)&currenthash,sizeof(uint64_t));
			}
			finhash.close();
			remove((outdir+uuid+std::string("hash.bin.")+std::to_string(tid)+'.'+std::to_string(j)).c_str());
		}
		
		// if (mm_verbose >= 3) fprintf(stderr, "[M::%s::%.3f*%.2f] \n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
		// mm_realtime0 = realtime();

		//correcting startpos array modified during insertion
		for(int64_t i = dict[j].numkeys; i >= 1 ; i--)	
			dict[j].startpos[i] = dict[j].startpos[i-1];
		dict[j].startpos[0] = 0;
	}//for end
	}//parallel end
	// fprintf(stderr, "end constructdictionary_realign()...\n");
	return;
}

void setglobalarrays_realign() {
	chartorevchar['A'] = 'T';
	chartorevchar['C'] = 'G';
	chartorevchar['G'] = 'C';
	chartorevchar['T'] = 'A';
	chartoint['A'] = 0;
	chartoint['C'] = 1;
	chartoint['G'] = 2;
	chartoint['T'] = 3;

	// if (readlen > 50) {
		int len_t = 17;
		if (readlen <= 80) len_t = 11;
		// int len_t = 31;
		// int len_t = reads->k;
		numdict_s = readlen / len_t;

		/*if (numdict_s <= 3) {
			maxsearch = 2000;

			len_t = readlen / 3;
			numdict_s = readlen / len_t;
		}*/

		// if (numdict_s > 4) {
			// numdict_s = 4;
		// }
		if (ininumdict > 1 && ininumdict < numdict_s) {
			numdict_s = ininumdict;
		}
		// numdict_s = 3;

		dict_start = new int[numdict_s];
		dict_end = new int[numdict_s];
		/*dict_start[0] = dict1_start;
		dict_end[0] = dict1_end;
		dict_start[1] = dict2_start;
		dict_end[1] = dict2_end;*/

		/*dict_start[0] = 0;
		dict_end[0] = 20;
		dict_start[1] = 21;
		dict_end[1] = 41;*/

		/*dict_start[0] = 0;
		dict_end[0] = 15;
		dict_start[1] = 16;
		dict_end[1] = 31*/;
		
		/*		
		dict_start[1] = 21;
		dict_end[1] = 41;
		dict_start[2] = 42;
		dict_end[2] = 62;		*/
		
		if (ininumdict > 0 && ininumdict < numdict_s) {
			dict_start[0] = readlen/2 - (len_t * numdict_s)/2;
		} else {
			dict_start[0] = 0;
		}
		dict_end[0] = dict_start[0] + len_t - 1;
		for (int i = 1; i < numdict_s; ++i) {
			dict_start[i] = dict_end[i-1] + 1;
			dict_end[i] = dict_start[i] + len_t - 1;
		}
		// for (int i = 0; i < numdict_s; ++i) {
		// 	fprintf(stderr, "%d, %d\n", dict_start[i], dict_end[i]);
		// }
		// fprintf(stderr, "----------------\n");

		/*dict_start[0] = 0;
		dict_end[0] = 19;
		for (int i = 1; i < numdict_s; ++i) {
			dict_start[i] = dict_end[i-1] + 1;
			dict_end[i] = dict_start[i] + 20;
		}*/
		/*dict_start[0] = 0;
		dict_end[0] = 31;
		dict_start[1] = 32;
		dict_end[1] = 63;*/
				// #define dict1_start 18
				// #define dict1_end 49
				// #define dict2_start 50
				// #define dict2_end 81
		// dict_start[0] = 18;
		// dict_end[0] = 47;
		// dict_start[1] = 50;
		// dict_end[1] = 79;
	// } else {
	// 	numdict_s = 2;
	// 	dict_start = new int[numdict_s];
	// 	dict_end = new int[numdict_s];
		
	// 	dict_start[0] = 0;
	// 	dict_end[0] = 20*readlen/50;
	// 	dict_start[1] = 20*readlen/50 + 1;
	// 	dict_end[1] = 41*readlen/50;
	// }
	
	for(int i = 0; i < 64; i++)
		mask64[i] = 1;

	// std::bitset<2*readlen> basemask[readlen][128]
	// std::bitset<2*readlen> positionmask[readlen]
	basemask = (std::bitset<2*readlen>**)calloc(readlen, sizeof(std::bitset<2*readlen>*));
	positionmask = (std::bitset<2*readlen>*)calloc(readlen, sizeof(std::bitset<2*readlen>));

	for(int i = 0; i < readlen; i++) {
		basemask[i] = (std::bitset<2*readlen>*)calloc(128, sizeof(std::bitset<2*readlen>));
		basemask[i]['A'][2*i] = 0;
		basemask[i]['A'][2*i+1] = 0;
		basemask[i]['C'][2*i] = 0;
		basemask[i]['C'][2*i+1] = 1;
		basemask[i]['G'][2*i] = 1;
		basemask[i]['G'][2*i+1] = 0;
		basemask[i]['T'][2*i] = 1;
		basemask[i]['T'][2*i+1] = 1;
		positionmask[i][2*i] = 1;
		positionmask[i][2*i+1] = 1;
	}		
	return;
}

struct kt_realign_hash_for_t;

typedef struct {
	struct kt_realign_hash_for_t *t;
	long i, n; //i < n;
	// int tid;
} ktf_realign_hash_worker_t;

typedef struct kt_realign_hash_for_t {
	int n_threads, index, threshold;//win is the length of window
	ktf_realign_hash_worker_t *w;
	reads_t *reads;
	std::bitset<2*readlen> *read, *mask, *revmask, *mask1;
	bbhashdict *dict;
	pthread_mutex_t *dict_lock;
	pthread_mutex_t *read_lock;
} kt_realign_hash_for_t;

bool encode_byte(char *seq, char *ref, int pos, int dir) {
	char *temp_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *en_str = (char*)alloca((readlen + 1) * sizeof(char));
	char *int_str = (char*)alloca(10 * sizeof(char));
	strcpy(temp_str, seq);
	if (dir) {
		reverse_complement(temp_str, readlen);
	}
	int en_str_len = 0;
	int eq_char_num = 0;
	for (int rj = pos, tj = 0; tj < readlen; ++rj, ++tj) {
		if (ref[rj] != temp_str[tj]) {
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
            }
			en_str[en_str_len++] = temp_str[tj];
		} else ++eq_char_num;
	}
	if (en_str_len == 0) {
		en_str[en_str_len++] = '0';
	}
	en_str[en_str_len] = '\0';
	return en_str_len <= readlen*0.4;
}

static void realign_hash_search(kt_realign_hash_for_t *t, int i_, int tid_) {
	cluster_t *p = &reads->clusters[t->index][tid_].a[i_];
	qsort(p->a, p->n, sizeof(uint64_t), cmpcluster2);
	std::bitset<2*readlen> ref, revref, b;
	int64_t *dictidx = (int64_t*)alloca(2 * sizeof(int64_t));//to store the start and end index (end not inclusive) in the dict read_id array
	uint32_t startposidx;//index in startpos
	bool flag = false;
	uint32_t current, k, rid;
	uint64_t ull, y;

	int pre_pn = p->n;
	// uint64_t y = p->a[p->n - 1]; // the last one
	// uint32_t rid = (uint32_t)(y >> 32);
	// int pos = (uint32_t)y >> 1;
	std::list<uint32_t> *deleted_rids = new std::list<uint32_t> [numdict_s];
	// b = stringtobitset(reads->seq[rid].seq);
	char *s1 = (char*)alloca((readlen + 1) * sizeof(char));
	int ref_len = strlen(p->ref) - readlen + 1;

	// bool debug = false;
	// if (strcmp(p->ref, "CCGTCACCCGGGGTCCCCAGGGTAGGCACGGCGAATACCATCGAAAGTTGATAGGGCAGCCGTTCGAATGGGTCGTCGCCGCCACGGGGGGCGTGCGATCGG") == 0) {
		// debug = true;
	// }

	/*if (debug) {
		fprintf(stderr, "p->n: %d\n", p->n);
		fprintf(stderr, "p->ref: %s\n", p->ref);
		for (int k = 0; k < p->n; ++k) { // p->a[k]
			uint64_t y = p->a[k];
			int rid = y>>32;
			int pos = (uint32_t)y>>1;
			int dir = y&1;

			fprintf(stderr, "%s\n", reads->seq[rid].seq);
		}
		fprintf(stderr, "---\n");
		exit(0);
	}*/

	for (int jj = 0; jj < ref_len; ++jj) {
		ref = stringtobitset(p->ref + jj);
		reverse_complement_(p->ref + jj, s1);
		revref = stringtobitset(s1);
		flag = false;
		int j = 0; // equivalent to  for (int j = 0; j < maxmatch; ++j) {
		//find forward match
		for (int l = 0; l < numdict_s; ++l) {
			if (dict_end[l] + j >= readlen) {
				continue;
			}
			// fprintf(stderr, "l: %d\n", l);

			b = ref & t->mask1[l];
			ull = (b >> 2*dict_start[l]).to_ullong();
			// fprintf(stderr, "%lu\n", );
			startposidx = t->dict[l].bphf->lookup(ull);
			if (startposidx >= t->dict[l].numkeys)//not found
				continue;
			//check if any other thread is modifying same dictpos
			if (pthread_mutex_trylock(&t->dict_lock[startposidx & 0xFFFFFF])) {
				continue;
			}
			// pthread_mutex_lock(&t->dict_lock[startposidx & 0xFFFFFF]);
			t->dict[l].findpos(dictidx, startposidx);
			// fprintf(stderr, "dictidx: %u\n", dictidx);
			if (t->dict[l].empty_bin[startposidx]) { //bin is empty
				pthread_mutex_unlock(&t->dict_lock[startposidx & 0xFFFFFF]);
				continue;
			}
			uint64_t ull1 = ((t->read[t->dict[l].read_id[dictidx[0]]] & t->mask1[l]) >> 2*dict_start[l]).to_ullong();
			if (ull == ull1) { //checking if ull is actually the key for this bin
				// fprintf(stderr, "begin enumurate...\n");
				for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--) {
					auto sg_id = t->dict[l].read_id[i];
					rid = reads->sg.a[sg_id];
					// if ((ref^(t->read[sg_id]&t->mask[j])).count() <= t->threshold  && (t->threshold <= 24 || encode_byte(reads->seq[rid].seq, p->ref, jj, 0))) { 
					// if ((ref^(t->read[sg_id]&t->mask[j])).count() <= t->threshold  && encode_byte(reads->seq[rid].seq, p->ref, jj, 0)) { 
					if (basediff(ref^(t->read[sg_id]&t->mask[j])) <= t->threshold  && encode_byte(reads->seq[rid].seq, p->ref, jj, 0)) { 
					// if ((ref^(t->read[sg_id]&t->mask[j])).count() <= t->threshold) { 
						pthread_mutex_lock(&t->read_lock[sg_id & 0xFFFFFF]);
						if (!reads->sg_flag[sg_id]) {
							reads->sg_flag[sg_id] = true;
							// if (reads->sg.a[sg_id] == 6366077) fprintf(stderr, "rid == 6366077 in kthread_hash_realign.c\n");
							flag = true;
						}
						pthread_mutex_unlock(&t->read_lock[sg_id & 0xFFFFFF]);
						if (flag) {
							flag = false;
							// rid = reads->sg.a[sg_id];
							// fprintf(stderr, "find!!! %s\n", reads->seq[rid].seq);
							// if (rid == 6366077) debug = true;

							y = (uint64_t)rid << 32 | ((uint64_t)(jj) << 1) | 0; 
							kv_push(uint64_t, *p, y);
							for(int l1 = 0; l1 < numdict_s; l1++) {
								deleted_rids[l1].push_back(sg_id);
							}
						}
					}
				}
			}
			pthread_mutex_unlock(&t->dict_lock[startposidx & 0xFFFFFF]);
			
			//delete from dictionaries
			for (int l1 = 0; l1 < numdict_s; ++l1) {
				for(auto it = deleted_rids[l1].begin(); it != deleted_rids[l1].end();) {
					b = t->read[*it] & t->mask1[l1];
					ull = (b >> 2*dict_start[l1]).to_ullong();
					startposidx = t->dict[l1].bphf->lookup(ull);
					if (pthread_mutex_trylock(&t->dict_lock[startposidx & 0xFFFFFF])) {
						++it;
						continue;
					}
					// pthread_mutex_lock(&t->dict_lock[startposidx & 0xFFFFFF]);
					t->dict[l1].findpos(dictidx, startposidx);
					t->dict[l1].remove(dictidx, startposidx, *it);
					it = deleted_rids[l1].erase(it);
					pthread_mutex_unlock(&t->dict_lock[startposidx & 0xFFFFFF]);
				}
			}
		}
		if (flag) continue;	
		//find reverse match
		for (int l = 0; l < numdict_s; l++) {
			if (dict_start[l] <= j) continue;
			b = revref&t->mask1[l];
			ull = (b>>2*dict_start[l]).to_ullong();
			startposidx = t->dict[l].bphf->lookup(ull);
			if (startposidx >= t->dict[l].numkeys)//not found
				continue;
			//check if any other thread is modifying same dictpos
			// pthread_mutex_lock(&t->dict_lock[startposidx & 0xFFFFFF]);
			if (pthread_mutex_trylock(&t->dict_lock[startposidx & 0xFFFFFF])) {
				continue;
			}
			t->dict[l].findpos(dictidx,startposidx);
			if (t->dict[l].empty_bin[startposidx]) {//bin is empty
				pthread_mutex_unlock(&t->dict_lock[startposidx & 0xFFFFFF]);
				continue;
			}
			uint64_t ull1 = ((t->read[t->dict[l].read_id[dictidx[0]]] & t->mask1[l])>>2*dict_start[l]).to_ullong();
			if (ull == ull1) { //checking if ull is actually the key for this bin
				for (int64_t i = dictidx[1] - 1 ; i >= dictidx[0] && i >= dictidx[1] - maxsearch; i--) {
					auto sg_id = t->dict[l].read_id[i];
					rid = reads->sg.a[sg_id];
					if ((revref^(t->read[sg_id]&t->revmask[j])).count() <= t->threshold && (t->threshold <= 24 || encode_byte(reads->seq[rid].seq, p->ref, jj, 1))) {	
					// if ((revref^(t->read[sg_id]&t->revmask[j])).count() <= t->threshold && encode_byte(reads->seq[rid].seq, p->ref, jj, 1)) {	
					// if ((revref^(t->read[sg_id]&t->revmask[j])).count() <= t->threshold) {	
						pthread_mutex_lock(&t->read_lock[sg_id & 0xFFFFFF]);
						if (!reads->sg_flag[sg_id]) {
							reads->sg_flag[sg_id] = true;
							flag = true;
						}
						pthread_mutex_unlock(&t->read_lock[sg_id & 0xFFFFFF]);
						if (flag) {
							flag = false;
							// if (rid == 6366077) debug = true;
							// rid = reads->sg.a[sg_id];
							y = (uint64_t)rid << 32 | ((uint64_t)(jj) << 1) | 1; 
							kv_push(uint64_t, *p, y);
							for(int l1 = 0; l1 < numdict_s; l1++) {
								deleted_rids[l1].push_back(sg_id);
							}
						}
					}
				}
			}
			pthread_mutex_unlock(&t->dict_lock[startposidx & 0xFFFFFF]);
			
			//delete from dictionaries
			for (int l1 = 0; l1 < numdict_s; ++l1) {
				for(auto it = deleted_rids[l1].begin(); it != deleted_rids[l1].end();) {
					b = t->read[*it] & t->mask1[l1];
					ull = (b >> 2*dict_start[l1]).to_ullong();
					startposidx = t->dict[l1].bphf->lookup(ull);
					if (pthread_mutex_trylock(&t->dict_lock[startposidx & 0xFFFFFF])) {
						++it;
						continue;
					}
					// pthread_mutex_lock(&t->dict_lock[startposidx & 0xFFFFFF]);
					t->dict[l1].findpos(dictidx, startposidx);
					t->dict[l1].remove(dictidx, startposidx, *it);
					it = deleted_rids[l1].erase(it);
					pthread_mutex_unlock(&t->dict_lock[startposidx & 0xFFFFFF]);
				}
			}
		}

	}

	delete[] deleted_rids;
	// update_reference(reads, p, p->n);
}


static void *ktf_realign_hash_worker(void *data)
{
	ktf_realign_hash_worker_t *w = (ktf_realign_hash_worker_t*)data;
	kt_realign_hash_for_t *t = w->t;
	// fprintf(stderr, "tid: %ld\n", w - w->t->w);

	long i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, 1);
		if (i >= w->n) break;
		realign_hash_search(t, i, w - w->t->w);
	}
	pthread_exit(0);
}

void kt_realign_hash_for(int n_threads, reads_t *reads, int index, int max_threshold, std::bitset<2*readlen> *read, bbhashdict *dict)
{
	int i;	
	kt_realign_hash_for_t t;
	pthread_t *tid;
	t.reads = reads, t.n_threads = n_threads, t.index = index, t.threshold = max_threshold, t.dict = dict, t.read = read;

	t.w = (ktf_realign_hash_worker_t*)calloc(n_threads, sizeof(ktf_realign_hash_worker_t));
	t.dict_lock = (pthread_mutex_t*)calloc(num_locks, sizeof(pthread_mutex_t));
	t.read_lock = (pthread_mutex_t*)calloc(num_locks, sizeof(pthread_mutex_t));
	for (int j = 0; j < num_locks; ++j) {
		pthread_mutex_init(&t.dict_lock[j], 0);
		pthread_mutex_init(&t.read_lock[j], 0);
	}

	t.mask = (std::bitset<2*readlen>*)calloc(maxmatch, sizeof(std::bitset<2*readlen>));
	t.revmask = (std::bitset<2*readlen>*)calloc(maxmatch, sizeof(std::bitset<2*readlen>));
	generatemasks(t.mask, t.revmask);
	t.mask1 = (std::bitset<2*readlen>*)calloc(numdict_s, sizeof(std::bitset<2*readlen>));
	generateindexmasks(t.mask1, numdict_s);

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));

	// fprintf(stderr, "n_threads: %d\n", n_threads);

	for (i = 0; i < n_threads; ++i) {
		t.w[i].t = &t, t.w[i].i = 0, t.w[i].n = reads->clusters[index][i].n;
	}
	// fprintf(stderr, "before pthread_create()...\n");

	for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktf_realign_hash_worker, &t.w[i]);
	for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);

	// fprintf(stderr, "after pthread_join()...\n");
	free(t.w);
	free(t.dict_lock);
	free(t.read_lock);

	free(t.mask);  
	free(t.mask1);
	free(t.revmask);
}

void realign_hash(int n_threads, reads_t *reads, int index, int max_threshold) { // index is clusters[index]
	
	numreads = reads->sg.n;
	// fprintf(stderr, "numreads: %d\n", numreads);

	// std::cerr << outdir+uuid+std::string("keys.bin.") << "\n";

	omp_set_num_threads(n_threads);
	setglobalarrays_realign();

	generateAllATbitset();

	std::bitset<2*readlen> *read = new std::bitset<2*readlen> [numreads];
	singleRead2bitset(reads, read, max_threshold);
	bbhashdict dict[numdict_s];
	// fprintf(stderr, "Constructing dictionaries\n");
	constructdictionary_realign(read, dict);
	// fprintf(stderr, "begin realign reads\n");
	// n_threads = 1;
	kt_realign_hash_for(n_threads, reads, index, max_threshold, read, dict);
	
	freeglobalarrays();
	delete[] read;

	// fprintf(stderr, "end realign_hash(reads_t *reads, int index)\n************-------\n");
}



