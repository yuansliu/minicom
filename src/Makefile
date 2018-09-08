CC=			g++
CFLAGS=		-O3 -Wno-unused-function -std=c++11 -w -march=native -fopenmp
# CFLAGS=		-g -O3 -Wno-unused-function -std=c++11 -w -march=native -fopenmp
CPPFLAGS=
INCLUDES=	-I.
OBJS=		bseq.o misc.o preprocess.o sketch.o bbhashdict.o kthread_reads.o kthread_bucket.o kthread_idx.o kthread_cb.o kthread_dump.o kthread_hash_realign.o
PEOBJS=		bseq.o misc.o preprocess.o sketch.o bbhashdict.o kthread_reads.o kthread_bucket.o kthread_idx.o kthread_cb.o kthread_dump.o kthread_dump_pe.o kthread_hash_realign.o
# PROG=		example
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

# all:$(PROG)

# extra:all $(PROG_EXTRA)

minicomsg:minicommain.o libminicomsg.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminicomsg $(LIBS)

minicompe:minicommain.o libminicompe.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminicompe $(LIBS)

decompress: decompress.o
		$(CC) $(CFLAGS) $< -o $@ -L. $(LIBS)

libminicomsg.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

libminicompe.a:$(PEOBJS)
		$(AR) -csru $@ $(PEOBJS)

# sdust:sdust.c kdq.h kvec.h kseq.h sdust.h
# 		$(CC) -D_SDUST_MAIN $(CFLAGS) $< -o $@ -lz

clean:
		rm -fr *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

# depend:
# 		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bseq.o: bseq.h kseq.h
# example.o: minicom.h bseq.h kseq.h
minicommain.o: bseq.h kseq.h kvec.h breads.h kdq.h khash.h config.h
# main.o: minicom.h bseq.h
misc.o: minicom.h bseq.h ksort.h
sketch.o: kvec.h minicom.h bseq.h
bbhashdict.o: bbhashdict.h breads.h bseq.h kvec.h config.h
preprocess.o: breads.h bseq.h kvec.h config.h
reordersecondfile.o: breads.h bseq.h
# readfiles.o: breads.h bseq.h
breads.o: breads.h bseq.h config.h
kthread_reads.o: breads.h bseq.h
kthread_bucket.o: breads.h bseq.h
kthread_idx.o: breads.h khash.h kvec.h
kthread_cb.o: breads.h khash.h kvec.h config.h
kthread_cluster.o: breads.h khash.h kvec.h
kthread_dump.o: breads.h khash.h kvec.h config.h
kthread_dump_pe.o: breads.h khash.h kvec.h config.h
kthread_hash_realign.o: bbhashdict.h breads.h khash.h kvec.h BooPHF.h config.h
decompress.o: config.h
