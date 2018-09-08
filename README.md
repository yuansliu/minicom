# minicom

Minicom is a tool for compressing short reads in FASTQ. The minicom program is written in C++11 and works on Linux. It is availble under an open-source license.

## Download & Install

	git clone https://github.com/yuansliu/minicom.git
	cd minicom
	sh install.sh

In the script `install.sh`, it downloads the tools *bsc* and *p7zip*. Please make sure the two tools can be ran on your machine.
    
## Usage
To compress:

    ./minicom -r IN.fastq 					
    ./minicom -r IN.fastq -p 				#order-preserving mode
	./minicom -1 IN_1.fastq -2 IN_2.fastq 	#preserving paired-end information		

Minicom creates compressed file `IN_comp.minicom`, `IN_comp_order.minicom` and `IN_comp_pe.minicom` respectively.

Options for compression

	-h 			print help message
	-t 			number of threads, default: 24 
	-r 			reads file 
	-1 			reads file of paired-end reads [used for paired-end]
	-2 			reads file of paired-end reads [used for paired-end]
	-k 			length of k-mer, k <= 31, default: 31
	-e 			difference threshold, default: 4
	-m 			number of minimizer, default 6
	-w 			window size, default L/2-k, where L is the length of reads
	-s 			number of indexed substring, default L/17 for L > 80; otherwise L/11
	-S 			step size of threshold, default S = e
	-E 			maximum threshold, default L/2
	-p 			order-preserving mode [only used for compressing single FASTQ file]
	-d 			a compressed file .minicom [only for decompression]

To decompress:

	./minicom -d IN.minicom [-t number of threads, default 24]

## Status
Under review

## Citation
Yuansheng Liu, Zuguo Yu, Marcel E. Dinger, Jinyan Li; Index suffix-prefix overlaps by (w, k)-minimizer to generate long contigs for reads compression. 2018.

### Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>

---

I thank my families and apologize to my babies.

<img src="mybabies.jpeg">
