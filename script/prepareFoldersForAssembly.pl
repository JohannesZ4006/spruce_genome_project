#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;



my @peFiles			= ();
my @peOutputs 		= ();

my @mpFiles			= ();
my @mpOutputs 		= ();

my $merged			= 0;

my $bashScriptFile 	= "";
my $help			= 0;



my $mergedPairs	= ""; # /srv/box2/projects/b2010042/assembly2.0/DATA/LIB_INFO/diploid_pe_180_merged/Z4006c01_lib180_merged.txt
my $notMerged_1 = ""; # /srv/box2/projects/b2010042/assembly2.0/DATA/LIB_INFO/diploid_pe_180_merged/Z4006c01_lib180_unmerged_1.txt
my $notMerged_2 = ""; # /bubo/home/h23/vezzi/spruce/b2010042/nobackup/francesco/assembly_2.0/LIB_INFO/diploid_pe_180_merged/Z4006c01_lib180_unmerged_2.txt

## typical command
## prepareFoldersForAssembly.pl --pe PATH_TO_DATA/LIB_INFO/haploid_466/Z4006s466_raw_libnames_lib180_interleaved.txt --pe PATH_TO_DATA/LIB_INFO/haploid_466/Z4006s466_raw_libnames_lib300_interleaved.txt --pe PATH_TO_DATA/LIB_INFO/haploid_466/Z4006s466_raw_libnames_lib650_interleaved.txt --peOutput lib180 --peOutput lib300 --peOutput lib650 --mp ../../../DATA/LIB_INFO/diploid_mp/libs2Kbp.txt --mp PATH_TO_DATA/LIB_INFO/diploid_mp/libs4Kbp.txt --mp PATH_TO_DATA/LIB_INFO/diploid_mp/libs10Kbp.txt --mpOutput lib2k --mpOutput lib4k --mpOutput lib10k --bash Abyss-run_13_07_12.job


GetOptions(
	'pe=s'              => \@peFiles,
	'peOutput=s'        => \@peOutputs,
	'mp=s'              => \@mpFiles,
	'mpOutput=s'        => \@mpOutputs,
	'merged'            => \$merged,
    'mergedPairs=s'     => \$mergedPairs,
    'merged_paired1=s'  => \$notMerged_1,
    'merged_paired2=s'  => \$notMerged_2,
	'bash=s'            => \$bashScriptFile,
	'help|?'            => \$help
) or pod2usage(2);


pod2usage(1) if $help;
pod2usage("--file cannot be left empty: specify at least one  dataset")  if (@peFiles == 0);
pod2usage("--peOutput same number as --pe")  if (@peFiles != @peOutputs);
pod2usage("--mpOutput same number as --mp")  if (@mpFiles != @mpOutputs);

if($merged != 0) {
    # in this case the files mergedPairs and notmergedpairs1 and 2 must be specified
    if($mergedPairs eq "") {
        die("if --merged is specified also --mergedPairs --merged_paired1 --merged_paired2 must be specified: $!");
    }
    if($notMerged_1 eq "") {
        die("if --merged is specified also --mergedPairs --merged_paired1 --merged_paired2 must be specified: $!");
    }
    if($notMerged_2 eq "") {
        die("if --merged is specified also --mergedPairs --merged_paired1 --merged_paired2 must be specified: $!");
    }
}

pod2usage("--bash")  if ($bashScriptFile eq "");



if(-e $bashScriptFile) {
    print "$bashScriptFile already exists, it is not a good idea to overwrite it, modify the name\n";
    exit 1;
}

open(BASH, ">$bashScriptFile" ) or die("not able to create file $bashScriptFile: $!");

print BASH "#!/bin/bash -l\n";
print BASH "#SBATCH -A b2010042\n";
print BASH "#SBATCH -p node -N 149 -n 149\n";
print BASH "#SBATCH -t 4-00:00:00\n";
print BASH "#SBATCH -J $bashScriptFile\n";
print BASH "#SBATCH -o $bashScriptFile\.out\n";
print BASH "#SBATCH -e $bashScriptFile\.err\n";
print BASH "#SBATCH --mail-type=ALL --mail-user=francesco.vezzi\@scilifelab.se\n";

print BASH "\n";
print BASH "\n";
print BASH "module load abyss/1.3.5\n";


### start by processing pe reads (interleaved)
my $libNum = 1;


my $peLibraryString = "lib='";
my $peReadsString   = "";
my $single 			= "";


if($merged) {
    print "$notMerged_1 test\n";
    open(UNMERGED1, "$notMerged_1");
    open(UNMERGED2, "$notMerged_2");
    open(MERGED, "$mergedPairs");

    $single="se='";
    $peLibraryString.="lib180unmerged ";
    $peReadsString .="lib180unmerged='";
    
    
    while(my $line = <UNMERGED1>) {
        chomp $line;
        system("ln -s $line lib180_$libNum\_1.fastq.gz");
        $line = <UNMERGED2>;
        chomp $line;
        system("ln -s $line lib180_$libNum\_2.fastq.gz");
        $line = <MERGED>;
        chomp $line;
        system("ln -s $line lib180_$libNum\.fastq.gz");
        
        $peReadsString .= "lib180_$libNum\_1.fastq.gz lib180_$libNum\_2.fastq.gz ";
        $single.=" lib180_$libNum\.fastq.gz ";
        $libNum ++;
    }
    $peReadsString .= "' ";
    $single.="' ";
    
}


foreach my $file (@peFiles) {
    my $library = shift @peOutputs;
    
    $peLibraryString.="$library ";
    $peReadsString.="$library='";
    print $library."\n";
	open(FILE, "$file") or die("not able to file $file: $!");
	while(my $line = <FILE>) {
        chomp $line;
    	if($line =~ /^\#/) {
    	    next;
		}
        my $libName = $library."_".$libNum.".i.fq.gz";
    	$peReadsString.="$libName ";
    
    	$libNum++;
   	    #print "$line $libName\n";
    	my $return = system("ln -s $line $libName");
    	if($return > 0) {
	    	die("not able to ln -s $line $libName: $!\n");
    	}

	}
	close FILE;
	$peReadsString.="' ";

}
$peLibraryString.="' ";


my $mpLibraryString = "";
my $mpReadsString   = "";
if(@mpFiles > 0) {
    $mpLibraryString = "mp='";
    foreach my $file (@mpFiles) {
        my $library = shift @mpOutputs;
        
        $mpLibraryString.="$library ";
        $mpReadsString.="$library='";
        print $library."\n";
        
        open(FILE, "$file") or die("not able to file $file: $!");
        while(my $line1 = <FILE>) {
            chomp $line1;
            my $line2 = <FILE>;
            chomp $line2;
      		
            my $libName1 = $library."_".$libNum."_1.fq.gz";
            my $libName2 = $library."_".$libNum."_2.fq.gz";
    		$mpReadsString.= "$libName1 $libName2 ";
            $libNum++;

	    	my $return = system("ln -s $line1 $libName1");
    		if($return > 0) {
	    		die("not able to ln -s $line1 $libName1: $!\n");
    		}
            $return = system("ln -s $line2 $libName2");
    		if($return > 0) {
	    		die("not able to ln -s $line2 $libName2: $!\n");
    		}
        }
        close FILE;
        $mpReadsString.="' ";
    }
     $mpLibraryString.="' ";
}






print BASH "abyss-pe n=10 s=500 v=-v --dry-run j=8 k=50 name=spruce $single $peLibraryString $mpLibraryString $peReadsString $mpReadsString";



print BASH "\n";
print BASH "\n";





__END__










