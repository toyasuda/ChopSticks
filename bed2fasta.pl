#!/usr/bin/perl

BEGIN{unshift(@INC,"$ENV{HOME}/work/tsv2sv/src");}
use strict;
use warnings;
use File::Path;
use Utility;

our $opt = Utility::options2hash(\@ARGV,{'w'=>100,'s'=>0.0001,'r'=>1});

if(@ARGV<1){
    die "Usage: $0 <fasta> <sv bed ...>\n";
}

our $fasta_file = shift @ARGV;

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

sub edit_fasta($$){
    my ($sequences,$handle)=@_;
    my @prev=('',-1,-1,'');
    my $seq='';
    my $newseq='';
    my $cursor=0;
    my $target_chr='';
    #my @k=keys %$sequences;
    #foreach my $c (@k){print STDERR "$c: $$sequences{$c}\n";}
    #print STDERR "# \@k=(@k)\n";
    while(my $l=<$handle>){
        chomp $l;
        if($l=~/^\s*\#/ || $l=~/^\s*$/){next;}
        my ($chr,$begin,$end,$type,$extra)=split(/\s+/,$l);
        #print STDERR "SV: ($chr,$begin,$end,$type,$extra)\n";
        if(($prev[0] && $prev[0] ne $chr) || $prev[1] >= $begin || $begin>$end){
            die "Bad SV data: ($chr,$begin,$end,$type,...) after (@prev)\n";
        }
        if($seq){
            if($target_chr ne $chr){die "Mixed SV data: $target_chr vs. $chr\n";}
        }else{
            #if(! defined $$sequences{$chr}){die "Unknown reference: $chr\n";}
            if(! defined $$sequences{$chr}){die "Unknown reference: $chr ($$sequences{$chr})\n";}
            $target_chr=$chr;
            $seq = $$sequences{$chr};
            $seq = insert_SNPs($chr,$seq,$$opt{'s'});
        }

        #$newseq .= substr($seq,$cursor,$cursor-$begin+1);
        if($begin-$cursor<=0){die "Negative length of substring to be add: $begin-$cursor\n";}
        $newseq .= substr($seq,$cursor,$begin-$cursor);
        if($type eq 'DEL'){
        }
        elsif($type eq 'INS'){
            if($extra && $extra=~/^[ATGCNatgcn]*$/){
                $newseq .= $extra;
            }elsif($extra=~/(\d+)/){
                $newseq .= Utility::random_sequence($1);
            }else{
                die "Inserted sequence must be specified: $l\n";
            }
        }
        else{
            die "Unexpected SV type: $type\n";
        }
        $cursor = $end;
        unless($extra){$extra="''";}
        #if(length($extra)>12){$extra=substr($extra,0,10) . "...";}
        #print STDERR "SV: ($chr,$begin,$end,$type,$extra), newseq:$newseq, cursor:$cursor\n";
        #print STDERR "SV: ($chr,$begin,$end,$type,$extra)\n";
        print STDERR "$l\n";
        @prev = ($chr,$begin,$end,$type,$extra);
    }
    $newseq .= substr($seq,$cursor);
    return ($target_chr,$newseq);
}

sub insert_SNPs($$$){
    my ($chr,$seq,$prob)=@_;
    for(my $i=0;$i<length($seq);++$i){
        if(rand()<$prob){
            my $b0=substr($seq,$i,1);
            my $b=$b0;
            while($b eq $b0){$b=$Utility::bases[ int(rand()*4) ];}
            my $i2=$i+1;
            print STDERR "$chr\t$i\t$i2\tSNP\t$b0\t$b\n";
            substr($seq,$i,1)=$b;
        }
    }
    return $seq;
}


our %sequences = ();
read_fasta(\%sequences,Utility::open_data($fasta_file));
our ($chr,$seq)=edit_fasta(\%sequences,*ARGV);
Utility::write_fasta($chr,$seq,*STDOUT,$$opt{'w'});
