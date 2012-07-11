#!/usr/bin/perl

BEGIN{unshift(@INC,"$ENV{HOME}/work/tsv2sv/src");}
use strict;
use warnings;
use File::Path;
use Utility;

our $opt = Utility::options2hash(\@ARGV,{'w'=>100,'s'=>0.0001,'r'=>1,'q'=>'I','t'=>$Utility::NO_VALUE});

if(@ARGV<5){
    die "Usage: $0 <readlen> <average of width> <sd of width> <depth> <output_base> <fasta ...>\n";
}

our $readlen = shift @ARGV;
our $avg_width = shift @ARGV;
our $sd_width  = shift @ARGV;
our $depth = shift @ARGV;
our $output_base = shift @ARGV;

if($$opt{'r'}>0){
    print STDERR "# Using seed $$opt{r}.\n";
    srand($$opt{'r'});
}else{
    print STDERR "# Using default seed.\n";
}

sub fasta_commit($$$){
    my ($h,$d,$s)=@_;
    unless($d){return;}
    $$h{$d}=$s;
}
sub read_fasta($$){
    my ($hash,$handle)=@_;
    #my $handle = open_data($filename);
    my ($defline,$seq)=('','');
    while(my $l=<$handle>){
        chomp $l;
        if($l=~/^\s*\#/ || $l=~/^\s*$/){next;}
        if($l=~/^>\s*(\S.*)$/){
            print STDERR "# Reading $1\n";
            fasta_commit($hash,$defline,$seq);
            $defline=$1;
            $seq='';
        }else{
            $seq .= $l;
        }
    }
    fasta_commit($hash,$defline,$seq);
    print STDERR "# Reading done. (last: $defline)\n";
}

sub rnorm {
    my ($m, $sigma) = @_;
    my ($r1, $r2) = (rand(), rand());
    while ($r1 == 0) { $r1 = rand(); }
    return ($sigma * sqrt(-2 * log($r1)) * sin(2 * 3.14159265359 * $r2)) + $m;
}

sub generate_pairs($$){
    my ($seq,$output_base)=@_;

    my $len = length($seq);
    print STDERR "# length=$len\n";
    my $total=int($len*$depth);
    my $n_sequences = int($len*$depth/($readlen*2))+1;
    print STDERR "# Generating 2*$readlen*$n_sequences=$total bp.\n";

    my ($handle1,$handle2);
    my $fq1="${output_base}_1.fq";
    my $fq2="${output_base}_2.fq";
    if(-f $fq1){die "$fq1 already exists.\n";}
    if(-f $fq2){die "$fq2 already exists.\n";}
    open($handle1,">$fq1") || die "Cannot create $fq1.\n";
    open($handle2,">$fq2") || die "Cannot create $fq2.\n";

    my $q = $$opt{'q'} x $readlen;
    my $n_hasN=0;
    for(my $i=0;$i<$n_sequences;){
        my $p = int(rand()*$len);
        my $d = 0;
        while($d<100){$d=rnorm($avg_width,$sd_width);}
        $d=int($d);
        if($p+$readlen*2+$d>$len){next;}
        my $s1=substr($seq,$p,$readlen);
        my $s2=substr($seq,$p+$d+$readlen,$readlen);
        #if($s1=~/N/ || $s2=~/N/){next;}
        if($s1=~/N/ || $s2=~/N/){$n_hasN++;}
        my $s2r='';
        for(my $j=0;$j<$readlen;++$j){$s2r.=substr($s2,$readlen-1-$j,1);}
        $s2r=~tr/ATGC/TACG/;
        #print STDERR "$s2 -> $s2r\n";
        print $handle1 "\@P${p}-$d\n";
        print $handle1 "$s1\n";
        print $handle1 "+\n";
        print $handle1 "$q\n";
        print $handle2 "\@P${p}-$d\n";
        print $handle2 "$s2r\n";
        print $handle2 "+\n";
        print $handle2 "$q\n";
        ++$i;
    }
    close($handle2);
    close($handle1);
    print STDERR "# $n_hasN pairs have N.\n";
}

if($$opt{'t'} ne $Utility::NO_VALUE){
    for(my $i=0;$i<$depth;++$i){
        print rnorm($avg_width,$sd_width),"\n";
    }
    exit(0);
}

our %h = ();
read_fasta(\%h,*ARGV);

our @k=keys %h;
if(@k!=1){die "Multiple sequences: @k\n";}
generate_pairs($h{$k[0]},$output_base);

