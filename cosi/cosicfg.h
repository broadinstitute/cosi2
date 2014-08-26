//
// Header: cosicfg.h
//
// Compile-time configuration options for cosi.
//

#ifndef __INCLUDE_COSICFG_H
#define __INCLUDE_COSICFG_H

#ifdef COSI_SUPPORT_COALAPX
#define IF_COSI_SUPPORT_COALAPX(x) x
#define IF_NOT_COSI_SUPPORT_COALAPX(x)
#define IFELSE_COSI_SUPPORT_COALAPX(x,y) x
#else
#define IF_COSI_SUPPORT_COALAPX(x)
#define IFNOT_COSI_SUPPORT_COALAPX(x) x
#define IFELSE_COSI_SUPPORT_COALAPX(x,y) y
#endif // #ifdef COSI_SUPPORT_COALAPX


#endif  // #ifndef __INCLUDE_COSICFG_H
