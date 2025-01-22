#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my ( $fasta, $minlen, $maxlen, $prefix, $rm_terminal_ns);

&GetOptions(
    "i=s" => \$fasta,
    "min=i" => \$minlen,
    "max=i" => \$maxlen,
    "prefix=s" => \$prefix,
    "rm_terminal_ns" => \$rm_terminal_ns
);

($fasta)
  || die "Filter sequence by length
    usage: $0 OPTIONS
    where options are:
        -i  <fasta sequence file>
        --min <sequence min length> 
        --max <sequence max length, default is -1>
        --prefix <add prefix to beginning of the sequence id, default is ''>
        --rm_terminal_ns < trim the terminal Ns from the sequences, default is off>
    \n";

$minlen ||= 1;
$maxlen ||= -1;
$prefix ||= "";

open( FASTA, "<$fasta" ) or die "Could not open $fasta to read, $!\n";
$/ = "\n>";
while (<FASTA>) {
    chomp;
    if ( my ( $seq_name, $seq ) = /^>?(\S+.*?)\n(.*)/s ) {
        $seq =~ s/\s+//g;
        #skip the seqs which only contain Ns
        next if $seq =~ /^N+$/i;
        if($rm_terminal_ns){
            $seq =~ s/^N+//i;
            $seq =~ s/N+$//i;
           
        }
        my $seqlen = length($seq);
        next if $seqlen == 0;
        $seq_name = "$prefix-$seq_name" if length($prefix) > 0;

        if ( $seqlen >= $minlen ) {
            if ( $maxlen != -1 ) {
                if ( $seqlen <= $maxlen ) {
                    print ">$seq_name\n$seq\n";
                }
            }
            else {
                print ">$seq_name\n$seq\n";

            }
        }
    }
}

$/ = "\n";

close(FASTA);
