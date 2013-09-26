#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;


##########################################################################################################
##########																						##########
##########                 MATE-PAIRS-APPLICATION PIPELINE										##########
##########																						##########
##########################################################################################################

### adapter
## -a CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -a AGATGTGTATAAGAGACAG
### references /bubo/nobackup/uppnex/reference/biodata/biodata.stage/genomes/Hsapiens/hg19/bwa/



my @reads_1 			= ();
my @reads_2				= (); 
my @IDs      			= (); ## headers name of all outputs files
my $histWidth			= 10000;
my $noAdaptorRemoval	= 0;
my $noDuplicatesRemoval = 0;
my $noClean 			= 0;
my $aligner             = "bwa";
my $reference 			= "";
my @adaptors            = ();
my $help  				= 0;

my $scriptDir    		= "/proj/b2013064/vezzi/scripts/projects/MP_pipelineProcessing/";

GetOptions(
    'aligner=s'         => \$aligner,
	'reads=s'  			=> \@reads_1,
	'reference=s'		=> \$reference,
    'adaptor=s'         => \@adaptors,
	'hist-width=s'		=> \$histWidth,
	'no-adapt-removal'	=> \$noAdaptorRemoval,
	'no-clean'     		=> \$noClean,
	'help|?'    		=> \$help
) or pod2usage(2);



pod2usage(1) if $help;
pod2usage("--reads cannot be left empty: specify at least one  dataset")  if (@reads_1 == 0);



if(@adaptors == 0) {
    @adaptors = ("CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG", "CTGTCTCTTATACACATCT",  "AGATGTGTATAAGAGACAG");
}

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
        
        ## aling original reads
        
        if($noAdaptorRemoval == 0) { #remove adoptor and re-aling
            ($read_1, $read_2, $ID) = removeLinker($read_1, $read_2, $ID);
            if($aligner eq "bwa") {
                alignReads_bwa($read_1, $read_2, $ID);
            } elsif($aligner eq "bwaMEM") {
                alignReads_bwaMEM($read_1, $read_2, $ID);
            } else {
                alignReads_bowtie($read_1, $read_2, $ID);
            }
        } else {
            if($aligner eq "bwa") {
                alignReads_bwa($read_1, $read_2, $ID);
            } elsif($aligner eq "bwaMEM") {
                alignReads_bwaMEM($read_1, $read_2, $ID);
            } else {
                alignReads_bowtie($read_1, $read_2, $ID);
            }
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
        
        ### NOW TRIM ADAPTER
       ($read_1, $read_2, $ID) = removeLinker($read_1, $read_2, $ID);
        
        
    }

}




sub alignReads_bwa {
    my ($read_1, $read_2, $outputID) = @_;
    my $MappingBase="$outputID-to-$referenceID";
    my $BAM        = "$MappingBase.bam";
    my $BAMMapped  = $MappingBase.".unsorted.mapped.bam";
    my $BAMmarkDup = "$MappingBase\_markDup.bam";
    
    my $Reads1SAI		="$MappingBase\_1.sai";
    my $Reads2SAI		="$MappingBase\_2.sai";
    my $Reads1SAI_err	="$MappingBase\_1.err";
    my $Reads2SAI_err	="$MappingBase\_2.err";
    
    system("bwa aln -t 8 -o 1 $reference $read_1 > $Reads1SAI 2> $Reads1SAI.err") == 0 or die("error:$!");
    system("bwa aln -t 8 -o 1 $reference $read_2 > $Reads2SAI 2> $Reads2SAI.err") == 0 or die("error:$!");
    system("bwa sampe -P -s $reference $Reads1SAI $Reads2SAI $read_1 $read_2 2> $MappingBase\_sampe.err | samtools view -Shb -F 4 - > $BAMMapped 2>> $MappingBase\_sampe.err") == 0 or die("error:$!");
    system("rm -f $Reads1SAI $Reads2SAI");
    push @ToBeDeleted, "$MappingBase\_sampe.err";
    push @ToBeDeleted, "$Reads1SAI.err";
    push @ToBeDeleted, "$Reads2SAI.err";
    
     sort_duplicates_GCcontentCheck($read_1, $read_2, $outputID)
}



sub alignReads_bwaMEM {
    my ($read_1, $read_2, $outputID) = @_;
    my $MappingBase="$outputID-to-$referenceID";
    my $BAM        = "$MappingBase.bam";
    my $SAMMapped  = $MappingBase.".unsorted.mapped.sam";
    my $BAMMapped  = $MappingBase.".unsorted.mapped.bam";
    my $BAMmarkDup = "$MappingBase\_markDup.bam";
    my $BWAmemErr	="$MappingBase\.err";
    
    system("bwa mem -M -t 8  $reference $read_1 $read_2 > $SAMMapped 2> $BWAmemErr") == 0 or die("error:$!");
    system("samtools view -S -F 4 -b $SAMMapped >  $BAMMapped");
    system("rm $SAMMapped");
    
    sort_duplicates_GCcontentCheck($read_1, $read_2, $outputID)
}

sub alignReads_bowtie {
    my ($read_1, $read_2, $outputID) = @_;
    my $MappingBase="$outputID-to-$referenceID";
    my $BAM        = "$MappingBase.bam";
    my $SAMMapped  = $MappingBase.".unsorted.mapped.sam";
    my $BAMMapped  = $MappingBase.".unsorted.mapped.bam";
    my $BAMmarkDup = "$MappingBase\_markDup.bam";
    my $BOWTIEerr	="$MappingBase\.err";
    my $BOWTIEout	="$MappingBase\.out";
    
    print "bowtie2 -p 8 -x $reference -1 $read_1 -2 $read_2 -S $SAMMapped 2> $BOWTIEerr > $BOWTIEout";
    system("bowtie2 -p 8 -x $reference -1 $read_1 -2 $read_2 -S $SAMMapped 2> $BOWTIEerr > $BOWTIEout") == 0 or die("error:$!");
    system("samtools view -S -F 4 -b $SAMMapped >  $BAMMapped");
    system("rm $SAMMapped");
    
    sort_duplicates_GCcontentCheck($read_1, $read_2, $outputID)
}



sub sort_duplicates_GCcontentCheck{
    my ($read_1, $read_2, $outputID) = @_;
    my $MappingBase="$outputID-to-$referenceID";
    my $BAM        = "$MappingBase.bam";
    my $SAMMapped  = $MappingBase.".unsorted.mapped.sam";
    my $BAMMapped  = $MappingBase.".unsorted.mapped.bam";
    my $BAMmarkDup = "$MappingBase\_markDup.bam";
    my $BWAmemErr	="$MappingBase\.err";
    
    system("samtools sort $BAMMapped $MappingBase");
    system("samtools index $BAM") == 0 or die("error:$!");
    
    system("rm $BAMMapped");
    
    print "number of aligned reads ";
    system("samtools view $BAM | wc -l");
    
    system("java -Xmx16g -XX:PermSize=1g -jar \$PICARD_HOME/CollectInsertSizeMetrics.jar MINIMUM_PCT=0 HISTOGRAM_FILE=\"$MappingBase.pdf\" INPUT=$BAM OUTPUT=$MappingBase\.collectInseSize HISTOGRAM_WIDTH=$histWidth 2> $MappingBase\_collectInsertSize.err") == 0 or die("error:$!");
    
    system("java -Xmx16g -XX:PermSize=1g -jar \$PICARD_HOME/MarkDuplicates.jar INPUT=$BAM OUTPUT=$BAMmarkDup METRICS_FILE=$MappingBase\.markDuplicates ASSUME_SORTED=true 2> $MappingBase\_MarkDuplicates.err") == 0 or die("error:$!");
    

    system("java -Xmx16g -XX:PermSize=2g -jar \$PICARD_HOME/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE=$reference  INPUT=$BAM OUTPUT=$MappingBase\.collectGcBias.txt CHART_OUTPUT=$MappingBase\.collectGcBias.pdf ASSUME_SORTED=true 2> $MappingBase\_collectGCbias.std > $MappingBase\_collectGCbias.out  ");

    
    push @ToBeDeleted, "$MappingBase\_collectInsertSize.err";
	push @ToBeDeleted, "$MappingBase\_MarkDuplicates.err";
    push @ToBeDeleted, "$MappingBase\_samtools_sort.err";
    push @ToBeDeleted, "$MappingBase\_collectGCbias.std";
    push @ToBeDeleted, "$MappingBase\_collectGCbias.out";

}



sub removeLinker {
    my ($read_1, $read_2, $outputID) = @_;
    my $trimTool ="/proj/a2012043/software/clc-adapter-trim/adapter_trim-4.08beta";
    
    print "looking for adapter\n";
    
    my $pairedInterlevedLinker   = "$outputID\.i.linker.fastq";
    my $pairedInterlevedNoLinker = "$outputID\.i.nonlinker.fastq";
    my $singleLinker             = "$outputID\.se.linker.fastq";
    my $singleNoLinker           = "$outputID\.se.nonlinker.fastq";
    my $toBeRemoved = "";
    foreach my $adaptor (@adaptors) {
        $toBeRemoved .= " -a $adaptor ";
    }
    
    system("$trimTool  -r -i $read_1 $read_2 -f $pairedInterlevedLinker -g $pairedInterlevedNoLinker -t $singleLinker -u $singleNoLinker -m 30 $toBeRemoved   > output_adapterTrimming.std 2> output_adapterTrimming.err") == 0 or die("error while running $trimTool: $!");
    
    
    my $pairsWithLinker   		= countReadsSingleFile($pairedInterlevedLinker)/2;
    my $pairsWithNoLinker 		= countReadsSingleFile($pairedInterlevedNoLinker)/2;
    my $seAfterLinkerRemoval	= countReadsSingleFile($singleLinker) + countReadsSingleFile($singleNoLinker);
    print "pairs with linker: $pairsWithLinker\n";
    print "pairs no linker: $pairsWithNoLinker\n";
    print "se after linker removal: $seAfterLinkerRemoval\n";
    
    
    push @ToBeDeleted, $pairedInterlevedLinker;
    push @ToBeDeleted, $pairedInterlevedNoLinker;
    push @ToBeDeleted, $singleLinker;
    push @ToBeDeleted, $singleNoLinker;
    push @ToBeDeleted, "output_adapterTrimming.std";
    push @ToBeDeleted, "output_adapterTrimming.err";
    
    my $pairedInterleavedTrimmed = "$outputID\.i.trimmed.fastq";
    system("cat $pairedInterlevedLinker $pairedInterlevedNoLinker > $pairedInterleavedTrimmed") == 0 or die("error while running cat $pairedInterlevedLinker $pairedInterlevedNoLinker > $pairedInterleavedTrimmed: $!");
    
    my $tailID .= "_trimmed";
    system("$scriptDir/fastqi2fastqp.pl $pairedInterleavedTrimmed $outputID$tailID > output_fastqi2fastqp.pl.std 2> output_fastqi2fastqp.pl.err") == 0 or die("error while running $scriptDir/fastqi2fastqp.pl: $!");
    push @ToBeDeleted, "output_fastqi2fastqp.pl.std";
    push @ToBeDeleted, "output_fastqi2fastqp.pl.err";
    
    $read_1 = "$outputID$tailID\_1\.fastq";
    $read_2 = "$outputID$tailID\_2\.fastq";
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
        
        if($reads_1 =~ /(.*)\_1.fastq$/) {
           ($ID, $dirname, $suffix) = fileparse($reads_1, qr/\_1\.fastq/);
        } elsif($reads_1 =~ /(.*)\_1.fastq.gz$/) {
            $gz = 1;
           ($ID, $dirname, $suffix) = fileparse($reads_1, qr/\_1\.fastq\.gz/);
        } else {
            die("file $reads_1 not in correct format: should look like name_1.fastq or name_1.fastq.gz: $!");
        }
    	
        #print $ID."\n";
        #print $dirname."\n";
        #print $suffix."\n";
       
        if($gz == 0) {
	        $reads_1    = $dirname.$ID."_1.fastq";
    		$reads_2 	= $dirname.$ID."_2.fastq";
        } else {
            $reads_1    = $dirname.$ID."_1.fastq.gz";
    		$reads_2 	= $dirname.$ID."_2.fastq.gz";
        }
        
	    if(! -e $reads_1) {
    	    die("File $reads_1 does not exists: please provided an existing file: $!");
    	}
        if(! -e $reads_2) {
            die("File $reads_2 does not exists: check that $reads_1 is a crrectly mate pair file: $!");
        }
        
        ### Check if I have the rights to write inthis folder
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
 
 SciLifeLab-MP-removeDuplicates.pl --reads fileToProcess_1.fastq [OPTIONS]
 
 Options (Mandatory)
 
 --reads file containing reads one. Data is assumed to be provided in pairs in two different files with the header but ending _1.fastq and _2.fastq respectively. This option can be repeated sevaral times
 
  Options (Not mandatory)
--reference
 --hist-width
 --no-adapt-removal
 --no-remove-dupl
 --help displays this help message
 --no-clean do not delete intermediate files








