
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <clocale>
#include <ctime>
#include <unistd.h>
#include <fcntl.h>
#include <search.h>
#include <sys/resource.h>
#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif
#include <cosi_rand/random.h>
#include <cosi_rand/mtwist.h>
#include <miscutil/utils.h>

namespace miscutil {

const char * DateStr( )
{    static char nowstr[256];
  memset( nowstr, 0, 256 );
  time_t nowbin;
  const struct tm *nowstruct;
  static int locale_has_been_set = 0;
  if( ! locale_has_been_set ) {
		(void)setlocale(LC_ALL, "");
		locale_has_been_set = 1;
  }
  if (time(&nowbin) == (time_t) - 1) return "(date unavailable - time failed)";
  nowstruct = localtime(&nowbin);
  if (strftime(nowstr, 80, "%a %b %d %H:%M:%S %Y", nowstruct) == (size_t) 0)
		 return "(date unavailable - strftime failed)";
  return nowstr;
}

void dbg(const char *fmt, ...) {
  va_list ap;

  printf( "%s: ", DateStr() );
  
  va_start(ap, fmt);
  vprintf (fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  fflush(stdout);
}

/// FuncProt: chk
/// Check that a pointeer is non-NULL.  If it is non-NULL, return the pointer.
/// If it is NULL, cause a fatal error with the specified error message.
void *chk(const void *p, const char *fmt, ...) {
  va_list ap;

  if ( p ) return (void *)p;

  printf( "Error: null ptr - %s: ", DateStr() );
  
  va_start(ap, fmt);
  vprintf (fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  fflush(stdout);
  throw std::logic_error( "error" );
  return NULL;
}

void chkCond_report_error(void) {
  fprintf( stderr, "Error: condition failed - " );
}

/// Function: chk
/// Check that a condition is true; abort with an error if not.
/// Same as assert(), except the check is always done regardless of -DNDEBUG status.
void chkCond(int cond, const char *fmt, ...) {
  va_list ap;

  if ( cond ) return;

  chkCond_report_error();
  fprintf( stderr, "%s: ", DateStr() );
  
  va_start(ap, fmt);
  vfprintf (stderr, fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  //exit(EXIT_FAILURE);
	throw std::logic_error( "error" );
  fflush(stdout);
}

/// Func: fopenChk
/// Open a file, with error-checking.
FILE *fopenChk( const char *fname, const char *mode ) {
  return ( FILE *)chk( fopen( fname, mode ), "could not open file %s in mode %s", fname, mode );
}


long GetMemUsage( )
{   
#ifdef __linux
	long currsize = 0;

	char pszStatFile[1024];
	sprintf(pszStatFile, "/proc/%d/statm", getpid());

	const int maxContentsSize = 1024;
	char pszStatFileContents[ maxContentsSize ];

	int statFileFd = open( pszStatFile, O_RDONLY );
	if ( statFileFd < 0 )
		 perror( "open() in GetMemUsage()" ); 
	else
	{
		long bytesRead = read( statFileFd, pszStatFileContents, maxContentsSize );
		if ( bytesRead < 0 ) 
			 perror( "read() in GetMemUsage()" ); 
		else
		{
			pszStatFileContents[bytesRead] = 0;
			currsize = strtol( pszStatFileContents, 0, 10 );
			currsize *= sysconf( _SC_PAGESIZE );
			currsize /= 1024;
		}

		statFileFd = close( statFileFd );
		if ( statFileFd < 0 )
			 perror( "close() in GetMemUsage()" ); 
	}
	return currsize;
#else
#ifdef HAVE_GETRUSAGE
	struct rusage my_rusage;
	getrusage(RUSAGE_SELF, &my_rusage);
	return my_rusage.ru_maxrss;
#else
	return 0;
#endif	 
#endif
}

#if 0

static char *my_strdup( const char *s );

char *my_strdup( const char *s ) {
  if ( !s ) return NULL;
  {
		char *new_s = (char *)malloc( strlen( s ) + 1 );
		strcpy( new_s, s );
		return new_s;
  }
}

char *get_ith_token( const char *str, int i ) {
  char *saveptr = NULL;
  int j;
  char *strCpy = my_strdup( str );
  char *strCpySave = strCpy;
  assert( i >= 0 );
  for ( j = 0; j < i-1; j++, strCpy = NULL ) {
		const char *tokenHere = strtok_r( strCpy, "\t", &saveptr );
		assert( tokenHere );
  }
  {
		char *result = my_strdup( strtok_r( NULL, "\t", &saveptr ) );
		free( strCpySave );
		return result;
  }
}

int get_ith_token_int( const char *str, int i ) {
  char *token = get_ith_token( str, i );
  int val = 0;
  chkCond( sscanf( token, "%d", &val ) == 1,
					 "error reading token#%d from %s (%s) as int", i, str, token );
  free( token );
  return val;
}

double get_ith_token_double( const char *str, int i ) {
  char *token = get_ith_token( str, i );
  double val = 0.0;
  chkCond( sscanf( token, "%lf", &val ) == 1,
					 "error reading token#%d from %s (%s) as double", i, str, token );
  free( token );
  return val;
}
int index_of_token( const char *str, const char *token ) {
  int i;
  char *saveptr = NULL;
  char *strCpy = my_strdup( str );
  char *strCpySave = strCpy;
  for ( i = 0; ; strCpy = NULL, i++ ) {
		const char *thisToken = strtok_r( strCpy, "\t", &saveptr );
		if ( !thisToken ) goto notfound;
		if ( !strcmp( thisToken, token ) ) {
			free( strCpySave );
			return i;
		}
  }
	notfound:
  fprintf( stderr, "token %s not found in string %s\n", token, str );
  assert(0);
  return -1;
}
#endif
int comparePointers( const void *pp1, const void *pp2 ) {
  const void * p1 = *((const void **)pp1);
  const void * p2 = *((const void **)pp2);
  return p1==p2 ? 0 : ( p1 < p2 ? -1 : 1 );
}

void dontFreePtr( void *p ) { }

#ifndef HAVE_TDESTROY
void tdestroy (void *root, void (*free_node)(void *nodep)) {
  
}
#endif

#if 1
int binarySearch( const cosi_double *a, int from, int to, cosi_double key ) {
  cosi_double midVal;
  to--;
  while (from <= to) {
		int mid = (from + to) >> 1;
		midVal = a[ mid ];
	 
		if (midVal < key) from = mid + 1;
		else if (midVal > key) to = mid - 1;
		else {
			return mid;
		}
  }
  return -( from + 1 );
}


#else
int binarySearch( const cosi_double *a, int from, int to, cosi_double key ) {
  /* actually do a linear search, just for comparison */
  a += from;
  while( from < to ) {
		cosi_double diff = *a - key;
		if ( diff == 0 ) return from;
		if ( diff > 0 ) return -( from + 1 );
		a++; from++; 
  }
  return -( to + 1 );
}

#endif

int *make_random_permutation( size_t n ) {
  int *perm = COSI_CALLOC( int, n );
  int *p = perm;
  size_t i = 0;

  while( i < n )
		 *p++ = i++;

  for (i = n; i > 1; i--)
		 SWAP(int, perm[i-1], perm[(int)(random_double() * n)]);

  return perm;
}

bool_t random_bit(void) {
  return ( mt_lrand() & 0x01 );
}

char *cosi_strdup( const char *s ) {
  size_t len = sizeof( char ) * strlen( s ) + 1;
  char *s_dup = (char *)malloc( len );
  memcpy( s_dup, s, len );
  return s_dup;
}

/* from http://www.freebsd.org/cgi/cvsweb.cgi/~checkout~/src/lib/libc/string/strtok.c?rev=1.9&content-type=text/plain&hideattic=0 */
char *
cosi_strtok_r(char *s, const char *delim, char **last)
{
  char *spanp, *tok;
  int c, sc;

  if (s == NULL && (s = *last) == NULL)
		 return (NULL);

  /*
	 * Skip (span) leading delimiters (s += strspn(s, delim), sort of).
	 */
	cont:
  c = *s++;
  for (spanp = (char *)delim; (sc = *spanp++) != 0;) {
		if (c == sc)
			 goto cont;
  }

  if (c == 0) {/* no non-delimiter characters */
		*last = NULL;
		return (NULL);
  }
  tok = s - 1;

  /*
	 * Scan token (scan for delimiters: s += strcspn(s, delim), sort of).
	 * Note that delim must have one NUL; we stop if we see that, too.
	 */
  for (;;) {
		c = *s++;
		spanp = (char *)delim;
		do {
			if ((sc = *spanp++) == c) {
				if (c == 0)
					 s = NULL;
				else
					 s[-1] = '\0';
				*last = s;
				return (tok);
			}
		} while (sc != 0);
  }
  /* NOTREACHED */
}


void cosi_fwrite_helper(const void *ptr, size_t size, size_t nmemb, FILE *stream, const char *expr, const char *fname, int line) {
  assert( ptr );
  assert( stream );
  size_t nwritten = fwrite( ptr, size, nmemb, stream );
  chkCond( nwritten == nmemb, "fwrite failed writing %s at %s:%d - writing %d, wrote only %d\n", expr, fname, line, nmemb, nwritten );
}
void cosi_fread_helper(void *ptr, size_t size, size_t nmemb, FILE *stream, const char *expr, const char *fname, int line) {
  assert( ptr );
  assert( stream );
  size_t nread = fread( ptr, size, nmemb, stream );
  chkCond( nread == nmemb, "fread failed reading %s at %s:%d - reading %d, read only %d\n", expr, fname, line, nmemb, nread );
}

/* Obtain a backtrace and print it to stdout. */
void print_trace (void)
{
  void *array[60];
#ifdef HAVE_EXECINFO_H     
  int size = backtrace (array, 60);
  backtrace_symbols_fd (array, size, 0);
#endif  
}

typedef const void *const_void_p_t;

static int ptr_dummy;

void print_ptr_ids(void);

int ptr_id( const void *p ) {
#define N_PTRS 1024  
  static const_void_p_t ptrs[ N_PTRS ];
  static int count;

  if ( p == NULL ) return 0;

  if ( p == &ptr_dummy ) {
		printf( "\n---------------\n" );
		printf( "\nptrs:\n" );
		for ( int i = 0; i < count; i++ )
			 printf( "%d %p\n", i, ptrs[ i ] );
		printf( "---------------\n" );
  }

  for ( int i = 0; i < count; i++ )
		 if ( ptrs[ i ] == p ) return (i+1);

  if ( count < N_PTRS-1 ) {
		if ( count == 0 ) atexit( print_ptr_ids );
		ptrs[ count++ ] = p;

		return count;
  }

  return -1;
}

void print_ptr_ids(void) {
  ptr_id( &ptr_dummy );
}

obj_id_t get_obj_id(void) {
  static obj_id_t next_obj_id = 1;
  return next_obj_id++;
}


cosi_double interpolate( const vector<cosi_double>& f, cosi_double i ) {
  int int_i = (int)floor( i );
  cosi_double frac_i = i - int_i;
  return ( fabs( int_i - i ) < 1e-5 ) ? f[ int_i ] :
		 (1-frac_i) * f[int_i] + frac_i*f[int_i+1]; 
}


cosi_double findWhere( const vector<cosi_double>& f, cosi_double c ) {
  for ( int i = 0; i < ((int)f.size())-1; i++ ) {
		if( (f[i] <= c && f[i+1] >= c)
				|| (f[i] >= c && f[i+1] <= c) ) {
			// Interpolate between i and i+1
			cosi_double frac = fabs( (f[i] - c)/(f[i+1]-f[i]) );
			return (1-frac)*i + frac*(i+1);
		}
  }
  // f never crosses c
  return -1;
}

cosi_double integrate( const vector<cosi_double>& f, const vector<cosi_double>& x, cosi_double from, cosi_double to ) {
  assert( 0 <= from && from <= f.size()-1 );
  assert( 0 <= to && to <= f.size()-1 );

  cosi_double firstF = interpolate( f, from );
  cosi_double firstX = interpolate( x, from );
  cosi_double lastF = interpolate( f, to );
  cosi_double lastX = interpolate( x, to );

  if ( ((int)from) == ((int)to) )
		 // The integral doesn't cross any of the f samples,
		 // handle this case separately
		 return 0.5*(firstF+lastF)*(lastX-firstX);

  int ceilFrom = (int)ceil( from );
  int floorTo = (int)floor( to );

  cosi_double I = 0.5*(f[ceilFrom]+firstF)*(x[ceilFrom]-firstX);
  for ( int i = ceilFrom; i < floorTo; i++ ) {
		cosi_double incr = 0.5*(f[i]+f[i+1])*(x[i+1]-x[i]);
		I += incr;
  }

  cosi_double finalAdd = 0.5*(f[floorTo]+lastF)*(lastX-x[floorTo]);
  I += finalAdd;
  return I;
}

std::string Date( )
{    char nowstr[80];
	time_t nowbin;
	const struct tm *nowstruct;
	static bool locale_has_been_set = false;
	if( ! locale_has_been_set ) {
		(void)setlocale(LC_ALL, "");
		locale_has_been_set = true;
	}
	if (time(&nowbin) == (time_t) - 1) return "(date unavailable - time failed)";
	nowstruct = localtime(&nowbin);
	if (strftime(nowstr, 80, "%a %b %d %H:%M:%S %Y", nowstruct) == (size_t) 0)
		 return "(date unavailable - strftime failed)";
	return std::string(nowstr);
}


#if 0
#ifdef _GNU_SOURCE
unsigned long timespec_since(const struct timespec *start)
{
  struct timespec end, temp;
  clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &end );
  
  if ((end.tv_nsec-start->tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start->tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start->tv_nsec;
  } else {
		temp.tv_sec = end.tv_sec-start->tv_sec;
		temp.tv_nsec = end.tv_nsec-start->tv_nsec;
  }
  return temp.tv_sec * 1000000000 + temp.tv_nsec;
}
#endif
#endif

}


