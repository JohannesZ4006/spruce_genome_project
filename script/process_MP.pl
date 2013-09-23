#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;


##########################################################################################################
##########                                                                                      ##########
##########                 MATE-PAIRS-APPLICATION PIPELINE                                      ##########
##########                                                                                      ##########
##########################################################################################################

## 2013-09-17: changed by Johannes, based on http://sprucewiki.scilifelab.se/spruceassembly/node/1720:
## -a TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG -a CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA
### Original:
### adapter
## -a CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -a AGATGTGTATAAGAGACAG
### references /bubo/nobackup/uppnex/reference/biodata/biodata.stage/genomes/Hsapiens/hg19/bwa/



my @reads_1             = ();
my @reads_2             = (); 
my @IDs                 = (); ## headers name of all outputs files
my $histWidth           = 10000;
my $noAdaptorRemoval    = 0;
my $noDuplicatesRemoval = 0;
my $noClean             = 0;
my $reference           = "";
my $help                = 0;

my $scriptDir           = "/bubo/home/h18/johann/scripts_spruce_qc/";


GetOptions(
    'reads=s'           => \@reads_1,
    'reference=s'       => \$reference,
    'hist-width=s'      => \$histWidth,
    'no-adapt-removal'  => \$noAdaptorRemoval,
    'no-remove-dupl'    => \$noDuplicatesRemoval,
    'no-clean'          => \$noClean,
    'help|?'            => \$help
) or pod2usage(2);



pod2usage(1) if $help;
pod2usage("--reads cannot be left empty: specify at least one  dataset")  if (@reads_1 == 0);


#my check if files exists and they are provided with the correct mate pair files
my $referenceID = "";
checkDataConsistency();
## NOW @IDs contains the base name

mkdir "toBeDeleted";
my @ToBeDeleted = ();

if($reference eq "") {
    MPpipeline_noRef();
} else {
    MPpipeline_Ref();
}

foreach my $file (@ToBeDeleted) {
    system("mv $file toBeDeleted");
}

exit 0;


sub MPpipeline_Ref {
    foreach my $i (0..@IDs-1) {
        #### compute number of reads prior to duplicate removal and quality trimming
       
        print "$reads_1[$i]\n";
        print "$reads_2[$i]\n";
        print "$IDs[$i]\n";
 
        my $read_1   = $reads_1[$i];
        my $read_2   = $reads_2[$i];
        my $ID       = $IDs[$i];

        my ($numReadsFile1, $numReadsFile2) = countReads($read_1, $read_2);
        print "number of original reads pairs is $numReadsFile1\n";
        
        ## align original reads
        alignReads($read_1, $read_2, $ID);
        
        if($noAdaptorRemoval == 0) { #remove adoptor and re-align
            ### NOW TRIM ADAPTER
            ($read_1, $read_2, $ID) = removeLinker($read_1, $read_2, $ID);
            ##now re-align
            alignReads($read_1, $read_2, $ID);
        }
    }
    
}




sub MPpipeline_noRef {
    
    foreach my $i (0..@IDs-1) {
        #### compute number of reads prior to duplicate removal and quality trimming
        my @ToBeDeleted = ();
        print "$reads_1[$i]\n";
        print "$reads_2[$i]\n";
        print "$IDs[$i]\n";
        
        my $read_1   = $reads_1[$i];
        my $read_2   = $reads_2[$i];
        my $ID       = $IDs[$i];
        my $outputID = $IDs[$i];
        my $tailID   = "";
        
        my ($numReadsFile1, $numReadsFile2) = countReads($read_1, $read_2);
        print "number of original reads pairs is $numReadsFile1\n";
        
        
        #create BAM files
        print "running FastqToSam.jar ...";
        system("java -jar -Xmx16g \$PICARD_HOME/FastqToSam.jar FASTQ=$read_1  FASTQ2=$read_2 OUTPUT=$outputID\.bam QUALITY_FORMAT=Standard SAMPLE_NAME=$ID LIBRARY_NAME=$ID > output_FastqToSam.std 2> output_FastqToSam.err") == 0 or die("Failed while running FastqToSam.jar: $!");
        print "done\n";
        
        print "running EstimateLibraryComplexity.jar ...";
        system("java -jar -Xmx16g \$PICARD_HOME/EstimateLibraryComplexity.jar INPUT=$outputID\.bam O=$outputID\_libComplexity\.txt MAX_RECORDS_IN_RAM=5000000 > output_EstimateLibraryComplexity.std 2> output_FastqToDam.err") == 0 or die("Failed while running EstimateLibraryComplexity.jar: $!");
        print "done\n";
        
        ###Print Table and Histogram
        my $line;
        open(PICARD_RES, "$outputID\_libComplexity\.txt") or die("problems while opening $outputID\_libComplexity\.txt: $!");
        $line = <PICARD_RES>; #->  ## net.sf.picard.metrics.StringHeader
        $line = <PICARD_RES>; #->  # net.sf.picard.sam.EstimateLibraryComplexity INPUT=[Small_1_121113_BC1C10ACXX_P205_101B_index2.bam] OUTPUT=Small_1_121113_BC1C10ACXX_P205_101B_index2_libComplexity.txt ....
        $line = <PICARD_RES>; #->  ## net.sf.picard.metrics.StringHeader
        $line = <PICARD_RES>; #->  # Started on: Fri Feb 08 14:12:10 CET 2013
        $line = <PICARD_RES>; #->
        $line = <PICARD_RES>; #->  ## METRICS CLASS        net.sf.picard.sam.DuplicationMetrics
        $line = <PICARD_RES>; #->  LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED ....
        my @results = (split /\t/,$line);
        print "$line";
        $line = <PICARD_RES>; #->  Small_1_121113_BC1C10ACXX_P205_101B_index2      0       40999   0       0       642     642     0.015659
        print "$line";
        $line = <PICARD_RES>; #->
        $line = <PICARD_RES>; #->  ## HISTOGRAM    java.lang.Integer
        $line = <PICARD_RES>; #->  duplication_group_count Small_1_121113_BC1C10ACXX_P205_101B_index2
        
        open(HIST, ">$outputID\_libComplexity\.hist");
        print HIST "count\tNumReads\n";
        while($line = <PICARD_RES>) {
            print HIST $line;
        }
        close HIST;
        close PICARD_RES;
        
        
        push @ToBeDeleted, "$outputID\.bam";
        
        # if($noDuplicatesRemoval == 0) {
        #    $tailID .= ".noDup";
        #    system("/proj/a2012043/software/clc-assembly-cell-beta-4.0.6-linux_64/remove_duplicates -p -r -i $read_1 $read_2 -o $outputID\.fq -s $outputID\_CLCstats.txt > output_CLCremoveDuplicates.std 2> output_CLCremoveDuplicates.err")  == 0 or die("not able to generate hitogram : $!");
            ## now split the file
        #
        #    $read_1 = "$outputID\_1$tailID\.fq";
        #    $read_2 = "$outputID\_2$tailID\.fq";
        #
        #    open(FASTQ, "$outputID\.fq") or die("not able to open file $outputID\.fq: $!");
        #    open(FASTQ_1, ">$read_1") or die("not able to open file $read_1: $!");
        #    open(FASTQ_2, ">$read_2") or die("not able to open file $read_2: $!");
        #    while(my $head_1  = <FASTQ>) {
        #        my $seq_1     = <FASTQ>;
        #        my $comment_1 = <FASTQ>;
        #        my $quality_1 = <FASTQ>;
        #
        #       my $head_2    = <FASTQ>;
        #        my $seq_2     = <FASTQ>;
        #       my $comment_2 = <FASTQ>;
        #        my $quality_2 = <FASTQ>;
        #
        #        print FASTQ_1 $head_1;
        #        print FASTQ_1 $seq_1;
        #        print FASTQ_1 $comment_1;
        #        print FASTQ_1 $quality_1;
        #
        #        print FASTQ_2 $head_2;
        #        print FASTQ_2 $seq_2;
        #        print FASTQ_2 $comment_2;
        #        print FASTQ_2 $quality_2;
        #
        #    }
        #    close FASTQ;
        #    close FASTQ_1;
        #    close FASTQ_2;
        #
        #    my ($numReadsFile1, $numReadsFile2) = countReads($read_1, $read_2);
        #    print "number of reads pairs after duplicate removal is $numReadsFile1\n";
        #    push @ToBeDeleted, $read_1;
        #    push @ToBeDeleted, $read_2;
        #    push @ToBeDeleted, "$outputID\.fq";
        #
        #}
        
        ### NOW TRIM ADAPTER
       ($read_1, $read_2, $ID) = removeLinker($read_1, $read_2, $ID);
        
        
    }

}




sub alignReads {
    my ($read_1, $read_2, $outputID) = @_;
    
    my $MappingBase="$outputID-to-$referenceID";
    
    my $BAM="$MappingBase.bam";
    my $Stats="$MappingBase.stats.txt";
    my $BAMMapped=$MappingBase.".unsorted.mapped.bam";
    my $BAMaligned="$MappingBase\_onlyAligned.bam";
    my $BAMmarkDup="$MappingBase\_onlyAligned_markDup.bam";
    
    my $Reads1SAI       ="$MappingBase\_1.sai";
    my $Reads2SAI       ="$MappingBase\_2.sai";
    my $Reads1SAI_err   ="$MappingBase\_1.err";
    my $Reads2SAI_err   ="$MappingBase\_2.err";
    
    system("bwa aln -t 8 -o 1 $reference $read_1 > $Reads1SAI 2> $Reads1SAI.err") == 0 or die("error:$!");
    system("bwa aln -t 8 -o 1 $reference $read_2 > $Reads2SAI 2> $Reads2SAI.err") == 0 or die("error:$!");
       
    system("bwa sampe -P -s $reference $Reads1SAI $Reads2SAI $read_1 $read_2 2> $MappingBase\_sampe.err | samtools view -Shb  - > $BAMMapped 2>> $MappingBase\_sampe.err") == 0 or die("error:$!");
    system("samtools sort $BAMMapped $MappingBase  > samtools_sort.std 2> $MappingBase\_samtools_sort.err") == 0 or die("error:$!");
    system("samtools index $BAM") == 0 or die("error:$!");
    system("rm -f $Reads1SAI $Reads2SAI $BAMMapped");
    system("samtools view -b -F 4 $BAM > $BAMaligned 2> $MappingBase\_samtools_view.err") == 0 or die("error:$!");
    system("java -Xmx16g -jar \$PICARD_HOME/CollectInsertSizeMetrics.jar MINIMUM_PCT=0 HISTOGRAM_FILE=\"$MappingBase.pdf\" INPUT=$BAMaligned OUTPUT=$MappingBase\.collectInseSize HISTOGRAM_WIDTH=$histWidth 2> $MappingBase\_collectInsertSize.err") == 0 or die("error:$!");
    system("java -Xmx16g -jar /bubo/sw/apps/bioinfo/picard/1.85/kalkyl/MarkDuplicates.jar INPUT=$BAMaligned OUTPUT=$BAMmarkDup METRICS_FILE=$MappingBase\.markDuplicates ASSUME_SORTED=true 2> $MappingBase\_MarkDuplicates.err") == 0 or die("error:$!");
    
    
    push @ToBeDeleted, "$MappingBase\_sampe.err";
    push @ToBeDeleted, "$Reads1SAI.err";
    push @ToBeDeleted, "$Reads2SAI.err";
    push @ToBeDeleted, "$MappingBase\_samtools_sort.err";
    push @ToBeDeleted, "$MappingBase\_samtools_view.err";
    push @ToBeDeleted, "$MappingBase\_collectInsertSize.err";
    push @ToBeDeleted, "$MappingBase\_MarkDuplicates.err";
}



sub removeLinker {
    my ($read_1, $read_2, $outputID) = @_;
    my $trimTool ="/bubo/home/h18/johann/bin/clc_trimmer/adapter_trim-4.08beta";
    
    print "looking for adapter\n";
    
    my $pairedInterlevedLinker   = "$outputID\.i.linker.fq";
    my $pairedInterlevedNoLinker = "$outputID\.i.nonlinker.fq";
    my $singleLinker             = "$outputID\.se.linker.fq";
    my $singleNoLinker           = "$outputID\.se.nonlinker.fq";
    system("$trimTool  -r -i $read_1 $read_2 -f $pairedInterlevedLinker -g $pairedInterlevedNoLinker -t $singleLinker -u $singleNoLinker -m 30 -a TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG -a CGTAATAACTTCGTATAGCATACATTATACGAAGTTATACGA   > output_adapterTrimming.std 2> output_adapterTrimming.err") == 0 or die("error while running $trimTool: $!");
    
    
    my $pairsWithLinker         = countReadsSingleFile($pairedInterlevedLinker)/2;
    my $pairsWithNoLinker       = countReadsSingleFile($pairedInterlevedNoLinker)/2;
    my $seAfterLinkerRemoval    = countReadsSingleFile($singleLinker) + countReadsSingleFile($singleNoLinker);
    print "pairs with linker: $pairsWithLinker\n";
    print "pairs no linker: $pairsWithNoLinker\n";
    print "se after linker removal: $seAfterLinkerRemoval\n";
    
    
    push @ToBeDeleted, $pairedInterlevedLinker;
    push @ToBeDeleted, $pairedInterlevedNoLinker;
    push @ToBeDeleted, $singleLinker;
    push @ToBeDeleted, $singleNoLinker;
    push @ToBeDeleted, "output_adapterTrimming.std";
    push @ToBeDeleted, "output_adapterTrimming.err";
    
    my $pairedInterleavedTrimmed = "$outputID\.i.trimmed.fq";
    system("cat $pairedInterlevedLinker $pairedInterlevedNoLinker > $pairedInterleavedTrimmed") == 0 or die("error while running cat $pairedInterlevedLinker $pairedInterlevedNoLinker > $pairedInterleavedTrimmed: $!");
    
    my $tailID .= "_trimmed";
    system("$scriptDir/fastqi2fastqp.pl $pairedInterleavedTrimmed $outputID$tailID > output_fastqi2fastqp.pl.std 2> output_fastqi2fastqp.pl.err") == 0 or die("error while running $scriptDir/fastqi2fastqp.pl: $!");
    push @ToBeDeleted, "output_fastqi2fastqp.pl.std";
    push @ToBeDeleted, "output_fastqi2fastqp.pl.err";
    
    $read_1 = "$outputID$tailID\.1\.fq";
    $read_2 = "$outputID$tailID\.2\.fq";
    push @ToBeDeleted, $pairedInterleavedTrimmed;
    
    system("sed -i 's/#0\/1_/#0\/1 /' $read_1");
    system("sed -i 's/#0\/2_/#0\/2 /' $read_2");
    
    my ($numReadsFile1, $numReadsFile2) = countReads($read_1, $read_2);
    print "number of reads pairs after adaptor trimming is $numReadsFile1\n";
    
    return ($read_1, $read_2, "$outputID$tailID");
    
}

sub checkDataConsistency {

    foreach my $reads_1 (@reads_1) {
        my $reads_2 = "";
        my ($ID, $dirname, $suffix);
        my $gz = 0;
        
        if($reads_1 =~ /(.*)\.1.fq$/) {
           ($ID, $dirname, $suffix) = fileparse($reads_1, qr/\.1\.fq/);
        } elsif($reads_1 =~ /(.*)\.1.fq.gz$/) {
            $gz = 1;
           ($ID, $dirname, $suffix) = fileparse($reads_1, qr/\.1\.fq\.gz/);
        } else {
            die("file $reads_1 not in correct format: should look like name.1.fq or name.1.fq.gz: $!");
        }
        
        #print $ID."\n";
        #print $dirname."\n";
        #print $suffix."\n";
       
        if($gz == 0) {
            $reads_1    = $dirname.$ID.".1.fq";
            $reads_2    = $dirname.$ID.".2.fq";
        } else {
            $reads_1    = $dirname.$ID.".1.fq.gz";
            $reads_2    = $dirname.$ID.".2.fq.gz";
        }
        
        if(! -e $reads_1) {
            die("File $reads_1 does not exists: please provide an existing file: $!");
        }
        if(! -e $reads_2) {
            die("File $reads_2 does not exists: check that $reads_1 is a correct mate pair file: $!");
        }
        
        ### Check if I have the rights to write in this folder
        die("Cannot write in current folder: $!") if (! -w "./");
        
        push @reads_2, $reads_2;
        push @IDs, $ID;
    }
    
    if($reference ne "") {
        if(! -e $reference) {
            #die("File $reference does not exists: please provided an existing file: $!");
        }
        my ($referenceDirname, $referenceSuffix);
        ($referenceID, $referenceDirname, $referenceSuffix) = fileparse($reference, qr/\.fasta/);

    }

}



sub countReads {
    my ($firstMate, $secondMate) = @_;
    
    
    
    if ($firstMate =~ /\.gz$/) {
        open(FASTQ, "gunzip -c $firstMate |") || die "not able to open file $firstMate";
    } else {
        open(FASTQ, "$firstMate") or die("not able to open file $firstMate: $!");
    }
    my $readsFile1 = 0;
    while(<FASTQ>) {
        $readsFile1 ++;
    }
    close FASTQ;
    $readsFile1 = $readsFile1/4;
    
    if ($secondMate =~ /\.gz$/) {
        open(FASTQ, "gunzip -c $secondMate |") || die "not able to open file $secondMate";
    } else {
        open(FASTQ, "$secondMate") or die("not able to open file $secondMate: $!");
    }
    my $readsFile2 = 0;
    while(<FASTQ>) {
        $readsFile2 ++;
    }
    $readsFile2 = $readsFile2/4;


    print "$firstMate contains $readsFile2\n";
    print "$secondMate contains $readsFile2\n";

    if($readsFile1 != $readsFile2) {
        print "ERROR: two files must contain same amount of reads\n";
        exit 0;
    }
    
    
    return ($readsFile1,$readsFile2);
    

    
}

sub countReadsSingleFile {
    my ($file) = @_;
    
    
    
    if ($file =~ /\.gz$/) {
        open(FASTQ, "gunzip -c $file |") || die "not able to open file $file";
    } else {
        open(FASTQ, "$file") or die("not able to open file $file: $!");
    }
    my $readsFile1 = 0;
    while(<FASTQ>) {
        $readsFile1 ++;
    }
    close FASTQ;
    $readsFile1 = $readsFile1/4;
    
        
    return $readsFile1;
    
    
    
}


__END__

=head1 SciLifeLab-MP-removeDuplicates.pl
 
=head1 SYNOPSIS
 
 SciLifeLab-MP-removeDuplicates.pl --reads fileToProcess_1.fq [OPTIONS]
 
 Options (Mandatory)
 --reads file containing reads one. Data is assumed to be provided in pairs in two different files with the same base name but ending .1.fq and .2.fq respectively. This option can be repeated sevaral times
 
 Options (Not mandatory)
 --reference
 --hist-width
 --no-adapt-removal
 --no-remove-dupl
 --help displays this help message
 --no-clean do not delete intermediate files
