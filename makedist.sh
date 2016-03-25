#!/usr/bin/env bash

#
# Script: makedist.sh
#
# Make the distribution, after saving the git info about the source version from which
# the distribution is being made.  Used by ../Makefile.am .
#
# Syntax: makedist.sh srcdir distdir

set -e -o pipefail

export srcdir=`realpath $1`
export distdir=`realpath $2`

echo srcdir is $srcdir distdir is $distdir

if [[ ! ( -d "$srcdir" ) ]]; then
	echo srcdir does not exist: $srcdir
	exit 1
fi

if [[ ! ( -d "$distdir" ) ]]; then
	echo distdir does not exist: $distdir
	exit 1
fi


prep_file () {
		touch $1
		chmod u+w $1
		echo Creating $1 ...
}

export RI=$distdir/release-info

if [[ -d "$srcdir/.git" ]]; then
	mkdir -p $RI

	prep_file $RI/git-commit.txt
	git --git-dir=$srcdir/.git --work-tree=$srcdir rev-parse HEAD > $RI/git-commit.txt

	prep_file $RI/git-status.txt
	git --git-dir=$srcdir/.git --work-tree=$srcdir status > $RI/git-status.txt
	
	prep_file $RI/git-diff.txt
	git --git-dir=$srcdir/.git --work-tree=$srcdir diff > $RI/git-diff.txt

	echo Saved git info to $RI
fi


