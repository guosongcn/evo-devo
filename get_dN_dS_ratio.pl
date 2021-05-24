#! /usr/bin/perl -w


##$ARGV[0] is the input file of yn
##$ARGV[1] is the gene name from sp
##$ARGV[2] is the gene name from Mm
##$ARGV[3] is the output file

print "$ARGV[0]\n";
open YN,"$ARGV[0]" or die $!;
open SP,"$ARGV[1]" or die $!;
open MM,"$ARGV[2]" or die $!;
open NEW,">$ARGV[3]";

my $line=0;
my $r1;
my $r2;
my $r3;
my $r4;
my $r5;

while(<YN>){
	$line++;
	if($line%114 == 83){
		$r1=(split /\s+/)[1];	
	}
	if($line%114 == 94){
		my $dn=(split /\s+/)[8];
		my $ds=(split /\s+/)[11];
		$r2=$dn/$ds;
	}
	if($line%114 == 111){
		my ($ds,$dn)=(/dS\s=\s+(.*)\sdN\s=\s+(.*)\sw/);
		$r3=$dn/$ds;
	}
	if($line%114 == 112){
		my ($ds,$dn)=(/dS\s=\s+(.*)\sdN\s=\s+(.*)\sw/);
		$r4=$dn/$ds;
	}
	if($line%114 == 113){
		my ($ds,$dn)=(/dS\s=\s+(.*)\sdN\s=\s+(.*)\sw/);
		$r5=$dn/$ds;
		my $sp=<SP>;
		chomp $sp;
		my ($sp_nm,$sp_score)=(split /\s/,$sp);
		my $mm=<MM>;
		chomp $mm;
		print NEW "$sp_nm\t$mm\t$sp_score\t$r1\t$r2\t$r3\t$r4\t$r5\n";
	}
}

close YN;
close SP;
close MM;
close NEW;
