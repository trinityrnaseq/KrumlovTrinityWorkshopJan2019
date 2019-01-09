#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use FindBin;
use Cwd;
use Carp;

use lib ("$FindBin::Bin/PerlLib");
use IniReader;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $trinotate_conf_file = "$FindBin::Bin/conf.txt";

my $usage = <<__EOUSAGE__;

#################################################################################
#
#  Optional:
#
#  --autopilot         automatically run the pipeline end-to-end
#
#  --conf <string>     configuration file to use.  Default: $trinotate_conf_file
#
#################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $AUTO_MODE = 0;


&GetOptions( 'help|h' => \$help_flag,
             'autopilot' => \$AUTO_MODE,
             
             'conf=s' => \$trinotate_conf_file,
             
    );

if ($help_flag) {
    die $usage;
}


unless ($trinotate_conf_file =~ m|^/|) {
    $trinotate_conf_file = cwd() . "/$trinotate_conf_file";
}

my $trinity_dir = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME set to Trinity installation directory";
$ENV{PATH} .= ":$trinity_dir";  ## adding it to our PATH setting.

my $trinotate_dir = $ENV{TRINOTATE_HOME} or die "Error, need env var TRINOTATE_HOME set to Trinotate installation directory (note this is different than Trinity) ";

unless ($ENV{TRANSDECODER_HOME}) {
    $ENV{TRANSDECODER_HOME} = "$ENV{HOME}/GITHUB/TransDecoder";
}

if (! -d "$ENV{TRANSDECODER_HOME}") {
    die "Error, cannot locate TransDecoder dir at  env var setting:  $ENV{TRANSDECODER_HOME} ";
}
$ENV{PATH} = "$trinity_dir/trinity-plugins/BIN/:$ENV{TRANSDECODER_HOME}:$ENV{PATH}";


my $OS_type = `uname`;

my $workdir = cwd();

## first check for tools needed.

my @tools = qw (Trinity
    bowtie2
    samtools
    igv.sh
    TransDecoder.LongOrfs
    TransDecoder.Predict
);
 
{
    my $missing_tool_flag = 0;
    foreach my $tool (@tools) {
        my $path = `which $tool`;
        unless ($path =~ /\w/) {
            print STDERR "Error, cannot find path to: $tool\n";
            $missing_tool_flag = 1;
        }
    }
    
    if ($missing_tool_flag) {
        die "\n\nTools must be in PATH setting before proceeding.\n\n";
    }

}


my $checkpoints_dir = $workdir . "/__TrinDemo_checkpoints_dir";
unless (-d $checkpoints_dir) {
    mkdir $checkpoints_dir or die "Error, cannot mkdir $checkpoints_dir";
}


my $STEP_COUNT = 0; # incremented at each process_cmd


##############
# run Trinity.

my $run_Trinity_cmd = "$trinity_dir/Trinity --seqType fq "
    . " --samples_file data/samples.txt "
    . " --CPU 2 --max_memory 2G --min_contig_length 150 "; 


&process_cmd($run_Trinity_cmd, "$checkpoints_dir/trinity.ok");

# Examine top of Trinity.fasta file
&process_cmd("head trinity_out_dir/Trinity.fasta", "$checkpoints_dir/head_trinity.ok");

# count the number of transcripts assembled.
&process_cmd("grep '>' trinity_out_dir/Trinity.fasta | wc -l ", "$checkpoints_dir/count_trans.ok");


# Get Trinity stats:
&process_cmd("$trinity_dir/util/TrinityStats.pl trinity_out_dir/Trinity.fasta", "$checkpoints_dir/trin_stats.ok");


## representation of reads by the assembly
&process_cmd("bowtie2-build trinity_out_dir/Trinity.fasta trinity_out_dir/Trinity.fasta", "$checkpoints_dir/bowtie2_build_read_assess.ok");

&process_cmd("bowtie2 --local --no-unal -x trinity_out_dir/Trinity.fasta -q -1 data/wt_SRR1582651_1.fastq -2 data/wt_SRR1582651_2.fastq | samtools view -Sb - | samtools sort -o - - > bowtie2.bam", "$checkpoints_dir/bowtie2_align_reads_assessss.ok");


## Examine read alignments in IGV

&process_cmd("samtools index bowtie2.bam", "$checkpoints_dir/idx_bowtie_alignments.ok");

my $igv_cmd = "igv.sh -g trinity_out_dir/Trinity.fasta bowtie2.bam";
if ($AUTO_MODE) {
    $igv_cmd .= " & ";
}
&process_cmd($igv_cmd, "$checkpoints_dir/igv_trinity_reads.ok")
;

###########################################
## assess number of full-length transcripts

&process_cmd("blastx -query trinity_out_dir/Trinity.fasta -db data/mini_sprot.pep -out blastx.outfmt6 -evalue 1e-20 -num_threads 2 -max_target_seqs 1 -outfmt 6", "$checkpoints_dir/blastx_for_full_length.ok");

&process_cmd("head blastx.outfmt6 | column -t", "$checkpoints_dir/check_blast_output.ok");

&process_cmd("$trinity_dir/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 trinity_out_dir/Trinity.fasta data/mini_sprot.pep", "$checkpoints_dir/tophat_blast_cov_stats.ok");


###################################
## Abundance estimation using salmon
###################################


my $align_estimate_command = "$trinity_dir/util/align_and_estimate_abundance.pl --seqType fq "
    . " --samples_file data/samples.txt "
    . " --transcripts trinity_out_dir/Trinity.fasta "
    . " --est_method salmon "
    . " --trinity_mode  "
    . " --prep_reference "
    . " --output_dir salmon";

&process_cmd($align_estimate_command, "$checkpoints_dir/salmon_expr_estimate.ok");

&process_cmd("ls -ltr", "$checkpoints_dir/see_expr_output_dir_results.ok");

&process_cmd("ls -ltr wt_SRR1582651", "$checkpoints_dir/show_expr_sample_out_dir.ok");

&process_cmd("head wt_SRR1582651/quant.sf | column -t", "$checkpoints_dir/show_quant_sf.ok");


my $get_quant_file_listing_cmd = 'find wt_* GSNO_* -name "quant.sf" | tee quant_files.list';
&process_cmd($get_quant_file_listing_cmd, "$checkpoints_dir/salmon_quant_files.list.ok");


my $make_matrix_cmd = "$trinity_dir/util/abundance_estimates_to_matrix.pl --est_method salmon --out_prefix Trinity --name_sample_by_basedir --quant_files quant_files.list --gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map";

&process_cmd($make_matrix_cmd, "$checkpoints_dir/get_matrix.ok");


# look at the output
&process_cmd("head  Trinity.isoform.counts.matrix | column -t ", "$checkpoints_dir/head.trinity_isoform_counts.ok");

&process_cmd("head Trinity.isoform.TMM.EXPR.matrix | column -t", "$checkpoints_dir/head.trinity_expr.ok");


## Examine the E90N50 statistic
&process_cmd("$trinity_dir//util/misc/contig_ExN50_statistic.pl Trinity.isoform.TMM.EXPR.matrix trinity_out_dir/Trinity.fasta > ExN50.stats", "$checkpoints_dir/ExNstats.ok");

&process_cmd("cat ExN50.stats", "$checkpoints_dir/cat_ExNstats.ok");

## Plot the values
&process_cmd("$trinity_dir/util/misc/plot_ExN50_statistic.Rscript ExN50.stats", "$checkpoints_dir/plot_ExN50.ok");

&show("ExN50.stats.plot.pdf");



##############
## DE analysis
##############


## run edgeR
&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Trinity.isoform.counts.matrix --samples_file data/samples.txt --method DESeq2 --output DESeq2_trans", "$checkpoints_dir/run.DESeq2_trans.ok");

# take a look at what edgeR generated:
&process_cmd("ls -ltr DESeq2_trans/", "$checkpoints_dir/ls.DESeq2_trans.dir.ok");

&process_cmd("head DESeq2_trans/Trinity.isoform.counts.matrix.GSNO_vs_wt.DESeq2.DE_results | column -t ", "$checkpoints_dir/head.DESeq2_trans.DE_results.ok");

&show("DESeq2_trans/Trinity.isoform.counts.matrix.GSNO_vs_wt.DESeq2.DE_results.MA_n_Volcano.pdf");

&change_dir("DESeq2_trans", "$checkpoints_dir/cd.DESeq2.ok");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../Trinity.isoform.TMM.EXPR.matrix --samples ../data/samples.txt -P 1e-3 -C 2",
             "$checkpoints_dir/analyze_diff_expr.ok");


&process_cmd("wc -l diffExpr.P1e-3_C2.matrix", "$checkpoints_dir/wc_diff_expr_matrix.ok"); # number of DE transcripts + 1

&show("diffExpr.P1e-3_C2.matrix.log2.centered.genes_vs_samples_heatmap.pdf");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData",
             "$checkpoints_dir/cut_clusters_tree.ok");

&show("diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60/my_cluster_plots.pdf");


## Now run the DE analysis for the genes

&change_dir("../", "$checkpoints_dir/cd_back_to_wd_after_DESeq2.ok");

&process_cmd("$trinity_dir/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Trinity.gene.counts.matrix --samples_file data/samples.txt  --method DESeq2 --output DESeq2_gene", "$checkpoints_dir/gene_DE_analysis.ok");

&process_cmd("ls -ltr DESeq2_gene/", "$checkpoints_dir/ls_DESeq2_gene_dir.ok");

#########################
## Time now for Trinotate
#########################

&process_cmd("pwd", "$checkpoints_dir/ensure_pwd_pre_Trinotate.ok");


# run the Trinotate bioinf processing and generate report
&run_Trinotate_demo(); # note, this cd's into Trinotate

#######################################
## Prep TrinotateWeb w/ Expression Data
#######################################

&change_dir("../Trinotate", "$checkpoints_dir/cd_back_to_Trinotate_from_DESeq2.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl --sqlite Trinotate.sqlite --transcript_mode --samples_file ../data/samples.txt --count_matrix ../Trinity.isoform.counts.matrix --fpkm_matrix ../Trinity.isoform.TMM.EXPR.matrix",
             "$checkpoints_dir/Trinotate.load_trans_expr_data.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl --sqlite Trinotate.sqlite --transcript_mode --samples_file ../data/samples.txt --DE_dir ../DESeq2_trans",
    "$checkpoints_dir/Trinotate.load_trans_DE_data.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_transcript_clusters.pl --sqlite Trinotate.sqlite --group_name DE_all_vs_all --analysis_name diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60 ../DESeq2_trans/diffExpr.P1e-3_C2.matrix.RData.clusters_fixed_P_60/*matrix",
             "$checkpoints_dir/Trinotate.load_expr_clusters.ok");

# load in the gene results

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl  --sqlite Trinotate.sqlite --gene_mode  --samples_file ../data/samples.txt --count_matrix ../Trinity.gene.counts.matrix --fpkm_matrix ../Trinity.gene.TMM.EXPR.matrix", "$checkpoints_dir/Trinotate.load_gene_expr_data.ok");

&process_cmd("$trinotate_dir/util/transcript_expression/import_expression_and_DE_results.pl --sqlite Trinotate.sqlite --gene_mode  --samples_file ../data/samples.txt --DE_dir ../DESeq2_gene ", "$checkpoints_dir/Trinotate.load_gene_DE_data.ok");


print STDERR "\n\n\tCommand-line Demo complete.  Congratulations! :)  Now explore your data via TrinotateWeb\n\n\n\n";



exit(0);

####
sub process_cmd {
    my ($cmd, $checkpoint) = @_;

    unless ($checkpoint) {
        die "Error, need checkpoint file defined";
    }
    
    $STEP_COUNT++;
    
    if (-e $checkpoint) { return; }

    
    unless ($AUTO_MODE) {
        
        my $response = "";
        while ($response !~ /^[YN]/i) {
            print STDERR "\n\n"
                . "###############################################\n"
                . "CMD: $cmd\n"
                . "###############################################\n\n"
                . "Execute (Y/N)? ";

            $response = <STDIN>;
        }

        if ($response =~ /^N/i) {
            print STDERR "\t *** Exiting on demand. ****\n\n"
                . "Goodbye. \n\n\n";
            exit(0);
        }
    }
    
    print STDERR "\tNow running:\n\t\t$cmd\n\n\n";
    
    my $ret = system($cmd);
    if ($ret) {
        confess  "Error, cmd: $cmd died with ret $ret";
    }
    
    system("touch $checkpoint");
    
    return;
}


sub show {
    my ($image) = @_;

    my $cmd;

    if ($OS_type =~ /linux/i) {
        ## use xpdf
        $cmd = "xpdf $image";
    }
    else {
        ## rely on ImageMagick:
        $cmd = "open $image";
    }
    
    if ($AUTO_MODE) {
        $cmd .= " & ";
    }
    
    &process_cmd($cmd, "$checkpoints_dir/view." . basename($image) . ".ok");

    return;
}



####
sub get_fq_files_listings {
    my (%samples) = @_;
    
    my @left_fq_files;
    my @right_fq_files;

    foreach my $condition_fq_lists_href (values %samples) {
        
        foreach my $replicates_fq_lists_aref (values %$condition_fq_lists_href) {
            
            my ($left_fq, $right_fq) = @$replicates_fq_lists_aref;
            push (@left_fq_files, "$left_fq");
            push (@right_fq_files, "$right_fq");
        }
    }

    return(\@left_fq_files, \@right_fq_files);
}


####
sub change_dir {
    my ($dest_dir, $checkpoint) = @_;

    
    eval {
        &process_cmd("cd $dest_dir", $checkpoint);
    };
    if ($@) {
        print STDERR "\n ** Note: if you see an error message above about not being able to cd, just ignore it... it's a weird behavior of this demo script. Rest assured we've \'cd $dest_dir\' just fine.   :)\n\n";
        system("touch $checkpoint");
    }

    # now do it in the script. :)
    chdir("$dest_dir") or die "Error, could not cd to $dest_dir"; 
 
    return;
}

####
sub run_Trinotate_demo {
   
    &process_cmd("mkdir Trinotate", "$checkpoints_dir/mkdir_Trinotate.ok");

    
    &change_dir("Trinotate", "$checkpoints_dir/cd_Trinotate.ok");
    
    my $ini_reader = new IniReader($trinotate_conf_file);
    
    my @sections = $ini_reader->get_section_headings();
    @sections = grep { $_ ne 'GLOBALS' } @sections;

    my %globals = $ini_reader->get_section_hash('GLOBALS');
    $globals{TRANSCRIPTS_FASTA} = "../trinity_out_dir/Trinity.fasta";
    $globals{GENE_TO_TRANS_MAP} = "../trinity_out_dir/Trinity.fasta.gene_trans_map";
    $globals{CPU} = 2;
    $globals{TRINOTATE_HOME} = $trinotate_dir;
    
    ## get command structs
    my @cmd_structs;
    foreach my $section (@sections) {
        my %keyvals = $ini_reader->get_section_hash($section);
        $keyvals{__SECTION__} = $section;
        
        if ($keyvals{RUN} =~ /^T/i) {
            push (@cmd_structs, \%keyvals);
        }
    }

    @cmd_structs = sort {$a->{RANK}<=>$b->{RANK}} @cmd_structs;

    foreach my $cmd_struct (@cmd_structs) {
        my $CMD = $cmd_struct->{CMD};
        $CMD = &substitute_tokens($CMD, \%globals);
        
        my $section_name = $cmd_struct->{__SECTION__};
        my $checkpoint_file = "$checkpoints_dir/Trinotate.$section_name.ok";
        
        &process_cmd($CMD, $checkpoint_file);
        
    }
    
    return;
}

####
sub substitute_tokens {
    my ($cmd, $globals_href) = @_;

    my %token_templates;
    while ($cmd =~ /(\{__\S+__\})/g) {
        my $token_template = $1;
        
        $token_templates{$token_template}++;
    }

    if (%token_templates) {
        foreach my $token_template (keys %token_templates) {
            $token_template =~ /\{__(\S+)__\}/ or die "Error, not able to parse token template: $token_template";
            my $token_name = $1;

            my $replacement_val = $globals_href->{$token_name};
            unless (defined $replacement_val) {
                die "Error, unable to identify global value for token name: $token_name of cmd: $cmd";
            }
            $cmd =~ s/$token_template/$replacement_val/g;
        }
    }

    return($cmd);
}
