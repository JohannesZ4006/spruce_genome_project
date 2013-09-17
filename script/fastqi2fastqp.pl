#!/usr/bin/perl -w

use strict;

=head1 FUNCTION

    Convert an interleaved fastq file to two separate (paired) fastqfiles. No sanity checks at this point.

=cut

my ($in, $outID)=@ARGV;

open(IN, $in);
open(OUT1, ">$outID\_1.fastq");
open(OUT2, ">$outID\_2.fastq");

while(my $seq1=<IN>){
    $seq1 .= <IN>;
    $seq1 .= <IN>;
    $seq1 .= <IN>;

    my $seq2 = <IN>;
    $seq2 .= <IN>;
    $seq2 .= <IN>;
    $seq2 .= <IN>;

    print OUT1 $seq1;
    print OUT2 $seq2;
}
