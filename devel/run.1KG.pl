#!/usr/bin/env perl

use warnings "all";
use strict;

my $length;

open(IN, "params") || die "Could not open params file.\n";
while (<IN>) {
    if (m/^length/) {
	$length = $';
	chomp($length);
	$length =~ s/\s+//g;
    }
}
my $cmnd1;
$cmnd1="recomap_hapmap2 /idi/sabeti-scratch/ilya/gsvn/Data/Ilya_Data/genmaps/hm2 $length model.test";
print $cmnd1;
my $resp = `$cmnd1`;
print $resp;
my $resp2 = `coalescent -p params -o out`;
print $resp2;
