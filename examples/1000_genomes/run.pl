#!/util/bin/perl -w
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
my $resp = `../../bin/recomap_hapmap2 ../../recosim/genmap_hapmap_b36 $length model.test`;
print $resp;
$resp = `../../bin/coalescent -p params -o out`;
print $resp;
