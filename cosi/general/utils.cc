
// * utils: Miscellaneous utils not specific to population genetics

// *** includes
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <clocale>
#include <ctime>
#include <cstdint>
#include <unistd.h>
#include <fcntl.h>
#include <search.h>
#include <string>
#include <utility>
#include <map>
#include <boost/filesystem/fstream.hpp>
#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/throw_exception.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception/errinfo_at_line.hpp>
#include <boost/exception/errinfo_file_name.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/container/vector.hpp>
//#include <sys/resource.h>
#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif
#include <cosi/general/utils.h>

// *** main impl

namespace cosi {
namespace util {

bool noDbgPrint = False;

struct MyDbgPrintInit {
   MyDbgPrintInit() { if ( getenv( "COSI_NO_DBG_PRINT" ) ) noDbgPrint = True; }
} myDbgPrintInit;

// *** Function DateStr
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

// *** Function dbg
void dbg(const char *fmt, ...) {
  va_list ap;

  printf( "%s: ", DateStr() );
  
  va_start(ap, fmt);
  vprintf (fmt, ap);
  va_end(ap);
  
  printf( "\n" );
  fflush(stdout);
}


static void chkCond_report_error(void);

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
  fflush(stdout);
  BOOST_THROW_EXCEPTION( cosi_error() << error_msg( "cosi error" ) );
}

/// Func: fopenChk
/// Open a file, with error-checking.
FILE *fopenChk( const char *fname, const char *mode ) {
  return ( FILE *)chk( fopen( fname, mode ), "could not open file %s in mode %s", fname, mode );
}

// *** Function GetMemUsage
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
  for (spanp = const_cast<char *>(delim); (sc = *spanp++) != 0;) {
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
    spanp = const_cast<char *>(delim);
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
#ifdef HAVE_EXECINFO_H     
  void *array[60];
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

namespace tsv {

// *** Function create_tsv_idx

// Create an index of a tsv file based on one of its columns.
//
// Input params:
//    - tsvFN :: name of the TSV file
//    - colNum :: column within the TSV file containing integer values to index
//
// Output params:
//    - idxFN :: the output file
TSVIdx::TSVIdx( filename_t tsvFN, unsigned colNum, filename_t idxFN ) {
	boost::container::vector< std::pair< index_t, istream::streampos > > idx_streamPos;
	{
  boost::filesystem::ifstream tsvFile( tsvFN );
  tsvFile.exceptions( ios::failbit | ios::badbit );
	int lineNum = 0;
  while ( true ) {
    istream::streampos lineBeg = tsvFile.tellg();
    std::string line;
    try { std::getline( tsvFile, line ); }
    catch( std::ios::failure ) { break; }
		++lineNum;
		if ( lineNum == 1 ) continue;

    // extract the value in the relevant column
		unsigned colsSeen = 0;
    size_t colBeg = 0;
		while ( colsSeen < colNum ) {
			colBeg = line.find_first_of( '\t', colBeg );
			if ( colBeg == std::string::npos )
				 BOOST_THROW_EXCEPTION( cosi_io_error()
																<< boost::errinfo_file_name( tsvFN.string() )
																<< boost::errinfo_at_line( lineNum )
																<< error_msg( "Error reading tsv file: could not find index column" ) );
			++colNum;
		}
		size_t colEnd = line.find_first_of( '\t', colBeg );
		if ( colEnd == std::string::npos ) colEnd = line.size();
		index_t idxVal = boost::lexical_cast<index_t>( line.substr( colBeg, colEnd - colBeg ) );
		idx_streamPos.emplace_back( idxVal, lineBeg );
  }
	}
	boost::filesystem::ofstream idxFile( idxFN );
  idxFile.exceptions( ios::failbit | ios::badbit );

	size_t numEntries = idx_streamPos.size();
	idxFile.write( (const char *)&numEntries, sizeof( numEntries ) );
	idxFile.write( (const char *)idx_streamPos.data(), numEntries * sizeof( index_t ) );
}

}  // namespace tsv

} // namespace util
} // namespace cosi
