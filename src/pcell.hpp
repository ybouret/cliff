
#ifndef PCELL_INCLUDED
#define PCELL_INCLUDED 1

#include "hcell.hpp"

#define LdVAR(ID) ID( vars(atry,#ID) )
#define LdTMX(ID) ID( vars(atry,#ID) * TIME_SCALE )

class PCell : public HCell
{
public:

    

    inline virtual ~PCell() throw() {}

    inline explicit PCell(Lua::VM &_vm) :
    HCell(_vm)
    {
        
    }


    inline void add_leaks(Array           &dYdt,
                          const double      ,
                          ConstArray       &Y,
                          ConstArray       &atry,
                          const Variables  &vars) const
    {
        const double Li6 = Y[Index::Li6];
        const double Li7 = Y[Index::Li7];
        const double Q   = Y[Index::Q];
        const double Vm  = Y[Index::Vm];
        const double A   = Y[Index::A];

        const double zeta = milliVolts2Zeta(Vm);
        const double psi  = PsiOf(zeta);
        const double emz  = exp(-zeta);

        const double LdVAR(Lambda);
        const double LdTMX(k7);
        const double LdVAR(sigma);
        const double k6     = k7*sigma;
        const double Li6out = eps6*Lambda;
        const double Li7out = eps7*Lambda;
        
        const double LdTMX(kQ);
        const double LdVAR(Qout);
        const double LdVAR(Aout);
        const double LdVAR(kA_ratio);

        const double kA    = kQ * kA_ratio;
        const double leak6 = k6 * psi * (Li6out * emz - Li6 );
        const double leak7 = k7 * psi * (Li7out * emz - Li7 );
        const double leakQ = kQ * psi * (Qout   * emz - Q   );
        const double leakA = kA * psi * (Aout   - A * emz   );

        const double LdVAR(gamma);
        dYdt[Index::Vm ] += gamma * (leak6+leak7+leakQ-leakA);
        dYdt[Index::Li6] += leak6;
        dYdt[Index::Li7] += leak7;
        dYdt[Index::Q  ] += leakQ;
        dYdt[Index::A  ] += leakA;

        
        
    }




private:
    Y_DISABLE_COPY_AND_ASSIGN(PCell);
};


#endif

