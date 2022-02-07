#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file")
	or die;

while (<FILE>) {
	
	$ID = $_;
	$seq = <FILE>;
	substr($ID,0,1) = "";
	($locus) = split(/\{/,$ID);
        $count{$locus}++;

}
close FILE;

open(FILE, "<$file")
	or die;

while (<FILE>) {
	
	$ID = $_;
	$seq = <FILE>;
	substr($ID,0,1) = "";
	($locus) = split(/\{/,$ID);
        $array{$locus}++;

	if ($count{$locus} >= 1) {
		print ">$locus" . "_" . $array{$locus} . "\n";
		print $seq;
	}
}
close FILE;


