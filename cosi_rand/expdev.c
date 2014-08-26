/* $Id: expdev.c,v 1.3 2011/05/06 23:13:09 sfs Exp $ */

/* calculates exponential deviates. used for interarrival
 * times for poisson processes.
 */

#include <math.h>
#include <cosi_rand/random.h>

double
expdev (void)
{
	double dum = 0;

	while (dum == 0.0)
		dum = (double) 1 - random_double();
	return -log (dum);
}
