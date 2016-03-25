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

my $resp = `../../recosim recParams $length`;
print $resp;
$resp = `../../coalescent -p params -o out`;
print $resp;


