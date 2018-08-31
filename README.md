# minicom
---
Minicom is a tool for compressing short reads in FASTQ. The minicom program is written in C++11 and works on Linux. It is availble under an open-source license.

## Download & Install

		git clone https://github.com/yuansliu/minicom.git
		cd minicom
		sh install.sh

In the script \`install.sh', it downloads the tools *bsc* and *p7zip*. Please make sure the two tools can be ran on your machine.
    
## Usage
To compress:

    	./minicom -r IN.fastq 

or order-preserving mode

    	./minicom -r IN.fastq -p

or with a paired-end reads		
		
		./minicom -1 IN_1.fastq -2 IN_2.fastq

Minicom create compressed file *IN_comp.minicom*, *IN_comp_order.minicom* and *IN_comp_pe.minicom* respectively.

To decompress:

		./minicom -d IN.minicom

Options for compression

		-h 			print help message
		-t 			number of threads, default: 24 
		-r 			reads file 
		-1 			reads file of paired-end reads
		-2 			reads file of paired-end reads
		-k 			length of k-mer, k <= 31, default: 31
		-e 			difference threshold, default: 4
		-s 			number of indexed substring
		-w 			window size, default L/2-k, where L is the length of reads
		-m 			number of minimizer, default 6
		-S 			step size of threshold, default S = e
		-E 			maximum threshold, default L/2
		-p 			order-preserving mode [only used for compressing single FASTQ file]
		-d 			a compressed file .minicom [only for decompression]

When decompressing, only two options `-d` and `-t` can be used. Other options used for compressing.

## Status
Under review

## Citation
Yuansheng Liu, Zuguo Yu, Jinyan Li; Index suffix-prefix overlaps by (w, k)-minimizer to generate long contigs for reads compression. 2018.

### Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>
