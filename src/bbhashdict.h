#ifndef MM_BBHASHDICT_H
#define MM_BBHASHDICT_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <list>
#include <algorithm>
#include <omp.h>
#include <atomic>
#include "breads.h"
#include "kvec.h"
#include "BooPHF.h"
#include "config.h"

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

class bbhashdict
{
	public:
	boophf_t * bphf;
	uint32_t numkeys;
	uint32_t *startpos;
	uint32_t *read_id;
	bool *empty_bin;
	void findpos(int64_t *dictidx, uint32_t &startposidx);
	void remove(int64_t *dictidx, uint32_t &startposidx, uint32_t current);
	bbhashdict()
	{
		bphf = NULL;
		startpos = NULL;
		read_id = NULL;
	}
	~bbhashdict()
	{
		delete[] startpos;
		delete[] read_id;
		delete bphf;
	}	
};

extern char revinttochar[4];//used in bitsettostring
extern char inttochar[4];
extern char chartorevchar[128];//A-T etc for reverse complement
extern int chartoint[128];//A-0,C-1 etc. used in updaterefcount
extern int *dict_start;
extern int *dict_end;
extern uint32_t numreads;
extern std::string outdir;
extern std::string uuid;
extern int numdict_s;

// std::bitset<2*readlen> basemask[readlen][128];//bitset for A,G,C,T at each position 
extern std::bitset<2*readlen> **basemask;//bitset for A,G,C,T at each position 
//used in stringtobitset, chartobitset and bitsettostring
// std::bitset<2*readlen> positionmask[readlen];//bitset for each position (1 at two bits and 0 elsewhere)
extern std::bitset<2*readlen> *positionmask;//bitset for each position (1 at two bits and 0 elsewhere)
//used in bitsettostring
extern std::bitset<2*readlen> mask64;//bitset with 64 bits set to 1 (used in bitsettostring for conversion to ullong)
extern std::bitset<2*readlen> allAbitset, allTbitset;
extern omp_lock_t allA_lock;
extern omp_lock_t allT_lock;

std::bitset<2*readlen> stringtobitset(char *s);
void bitsettostring(std::bitset<2*readlen> b, char *s);
std::bitset<2*readlen> chartobitset(char *s);
void generateindexmasks(std::bitset<2*readlen> *mask1, int num_dict); //masks for dictionary positions
void generatemasks(std::bitset<2*readlen> *mask, std::bitset<2*readlen> *revmask);
void reverse_complement_(char *s, char *s1);
void singleRead2bitset(reads_t *reads, std::bitset<2*readlen> *read, int threshold);
void generateAllATbitset();
void freeglobalarrays();
int basediff(std::bitset<2*readlen> bitstr);

#endif
