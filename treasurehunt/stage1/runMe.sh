#!/bin/bash

set -e

../trinity_msg_decoder --encode msg1.txt  >  target.stage1.seq.fa

~/GITHUB/trinityrnaseq/util/misc/simulate_illuminaPE_from_transcripts.pl --transcripts target.stage1.seq.fa  --require_proper_pairs --frag_length 50 --read_length 30

cp reads.simPE_R30_F50_FR.left.fa reads_stage1.left.fa
cp reads.simPE_R30_F50_FR.right.fa reads_stage1.right.fa

~/GITHUB/trinityrnaseq/Trinity --seqType fa --left reads_stage1.left.fa --right reads_stage1.right.fa --CPU 1 --max_memory 1G --bypass_java_version_check

../trinity_msg_decoder --decode trinity_out_dir/Trinity.fasta | sort


