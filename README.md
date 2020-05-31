# ***Notice***
1. As "alloca" is used in the source code. It needs more stack memory. Plase check the stack size before running `minicom` using the command `ulimit -s`. If it is very small, please change it by `ulimit -s unlimited`. We are grateful to Prof. Sebastian Deorowicz for reporting the bug.

2. Tomasz Kowalski and Szymon Grabowski reported that minicom fails to correctly decompress several datasets in they work [PgRC](https://www.biorxiv.org/content/10.1101/710822v1). We fixed it in the commit `28a736a`. The bug is caused when construting contigs. We re-ran some datasets and found that the compression ratio is not affected much (see following table).

| Dataset | Compressed file size reported with bug | Latest compressed file size|
| ------------- | :-------------: | :-------------: |
| SRR445718 | 126064640 | 126085120 |
| SRR445724 | 245585920 | 245585920 |
| SRR445726 | 215572480 | 215582720 |
| SRR490961 | 122552320 | 122552320 |
| SRR490976 | 146636800 | 146575360 |

# minicom

Minicom is a tool for compressing short reads in FASTQ. The minicom program is written in C++11 and works on Linux. It is availble under an open-source license.

Note: Minicom only compresses DNA sequences in the FASTQ file. It does not support to compress the whole FASTQ file.

## Download & Install

	git clone https://github.com/yuansliu/minicom.git
	cd minicom
	sh install.sh

In the script `install.sh`, it downloads the tools *bsc* and *p7zip*. Please make sure the two tools can be ran on your machine.

## Usage
To compress:

    ./minicom -r IN.fastq 					
    ./minicom -r IN.fastq -p 				#order-preserving mode
	./minicom -1 IN_1.fastq -2 IN_2.fastq 			#preserving paired-end information		

Minicom creates compressed file `IN_comp.minicom`, `IN_comp_order.minicom` and `IN_comp_pe.minicom` respectively.

Options for compression

	-h 		print help message
	-t 		number of threads, default: 24 
	-r 		reads file 
	-1 		reads file of paired-end reads [used for paired-end]
	-2 		reads file of paired-end reads [used for paired-end]
	-k 		length of k-mer, k <= 31, default: 31
	-e 		difference threshold, default: 4
	-m 		number of minimizer, default 6
	-w 		window size, default L/2-k, where L is the length of reads
	-s 		number of indexed substring, default L/17 for L > 80; otherwise L/11
	-S 		step size of threshold, default S = e
	-E 		maximum threshold, default L/2
	-p 		order-preserving mode [only used for compressing single FASTQ file]

To decompress:

	./minicom -d IN.minicom [-t number of threads, default 24]

<!-- ## Status -->
<!-- Under review -->

## Citation
Yuansheng Liu, Zuguo Yu, Marcel E. Dinger, Jinyan Li; Index suffix-prefix overlaps by (w, k)-minimizer to generate long contigs for reads compression. *Bioinformatics*, 35(12):2066-2074, 2019.

### Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>

---
