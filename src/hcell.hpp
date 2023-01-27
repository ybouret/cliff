#ifndef HCELL_INCLUDED
#define HCELL_INCLUDED 1

#include "types.hpp"
#include "y/ptr/arc.hpp"

//! for displaying value
#define DISPLAY_VALUE_ALIGN 12

//! display algined values
template <typename T>
static inline void display_value(const char *text, const T &value, const char *units)
{
    const string s = text;
    std::cerr << s; for(size_t i=s.size();i<=DISPLAY_VALUE_ALIGN;++i) std::cerr << ' ';
    std::cerr << '=' << ' ' << value;
    if(units) std::cerr << ' ' << units;
    std::cerr << std::endl;
}

#define VSHOW(ID)       display_value(#ID,ID,NULL)
#define USHOW(ID,UNITS) display_value(#ID,ID,UNITS)

//! helper macro
#define LdLua(VAR) VAR( vm->get<double>( #VAR ) )

//! all the fixed/experimental parameters
class HCell 
{
public:
    Lua::VM      vm;          //!< virtual machine
    const double A_h;         //!< th parameter, in seconds
    const double B_h;         //!< th parameter, in [mM] !! Warning !!
    const double pH_ref;      //!< reference pH
    const double pH_ini;      //!< pH parameter
    const double Lambda_h;    //!< pH parameter, in [mM] !! Warning
    const double h_ini;       //!< pH parameter

    const double pH_out;      //!< outer pH
    const double h_out;       //!< [H+] out in [mM] !!

    const double rho_s;       //!< reference Lithium ratio
    const double d7out;       //!< external Lithium separation
    const double eps6;        //!< Li6 fraction
    const double eps7;        //!< Li7 fraction
    const double Temperature; //!< in Kelvin
    
    
    inline explicit HCell( const Lua::VM &_vm ) :
    vm(_vm),
    LdLua(A_h),
    LdLua(B_h),
    LdLua(pH_ref),
    LdLua(pH_ini),
    LdLua(Lambda_h),
    h_ini( pow(10.0,-pH_ini) * PROTON_SCALE),

    LdLua(pH_out),
    h_out( pow(10.0,-pH_out) * PROTON_SCALE ),

    LdLua(rho_s),
    LdLua(d7out),
    eps6(1.0/(1.0+rho_s*(1.0+d7out/1000.0))),
    eps7(1.0-eps6),
    LdLua(Temperature)
    {
    }

    inline virtual ~HCell() throw()
    {
    }

    inline void show() const
    {
        std::cerr << "<HCell>" << std::endl;
        std::cerr << "\t<pH recovery>" << std::endl;
        VSHOW(A_h);
        VSHOW(B_h);
        VSHOW(pH_ref);
        VSHOW(pH_ini);
        VSHOW(Lambda_h);
        USHOW(h_ini,"mM");
        std::cerr << "\t<pH recovery/>" << std::endl;

        std::cerr << "\t<environment>" << std::endl;
        VSHOW(pH_out);
        USHOW(h_out,"mM");
        VSHOW(rho_s);
        VSHOW(d7out);
        VSHOW(eps6);
        VSHOW(eps7);
        VSHOW(Temperature);
        std::cerr << "\t<environment/>" << std::endl;

        std::cerr << "<HCell/>" << std::endl;
    }

    //--------------------------------------------------------------------------
    // pH recovery parameters
    //--------------------------------------------------------------------------

    //! compute scaling time from Lambda in mM
    inline double get_th(const double Lambda) const throw()
    {
        return A_h/qerf(B_h*Lambda);
    }

    //! compute h_end from Lambda
    inline double get_pH_end(const double Lambda) const throw()
    {
        return pH_ini + (pH_ref-pH_ini) * Lambda/(Lambda+Lambda_h);
    }

    //--------------------------------------------------------------------------
    // fractionation utils
    //--------------------------------------------------------------------------

#if 0
    //! ratio => delta
    inline double ratio2delta( const double r ) const throw()
    {
        return ( (1.0+d7out/1000.0)*r - 1.0 ) * 1000.0;
    }

    //! delta => ratio
    inline double delta2ratio( const double delta7 ) const throw()
    {
        return (1.0+delta7/1000.0)/(1.0+d7out/1000.0);
    }
#endif

    inline double conc_ratio_to_delta( const double r) const throw()
    {
        return 1000.0 * (r/rho_s-1.0);
    }

    //--------------------------------------------------------------------------
    // potential utils
    //--------------------------------------------------------------------------

    //! mV=>zeta
    inline double milliVolts2Zeta(const double milliVolts) const
    {
        return (Y_FARADAY*1e-3*milliVolts)/(Y_R*Temperature);
    }

    //! mV=>Theta
    inline double milliVolts2Theta(const double milliVolts) const
    {
        return exp( - (Y_FARADAY*milliVolts*1e-3)/(Y_R*Temperature) );
    }



private:
    Y_DISABLE_COPY_AND_ASSIGN(HCell);
};


#endif
