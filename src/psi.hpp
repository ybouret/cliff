#ifndef PSI_INCLUDED
#define PSI_INCLUDED 1

#include <cmath>

static inline double PsiExp( const double u )
{
    return u/(1.0-exp(-u));
}

static inline double PsiPol( const double u )
{
    const double u2 = u*u;
    const double u4 = u2*u2;
    const double u6 = u4*u2;
    const double u8 = u4*u4;
    return 1.0 + u/2.0 + u2/12.0 -u4/720.0 + u6/30240.0 - u8/1209600.0;
}

static inline double PsiOf( const double u )
{

    if( fabs(u) <= 1e-4 )
    {
        return PsiPol(u);
    }
    else
    {
        return PsiExp(u);
    }

}

#endif
