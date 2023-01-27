#ifndef TYPES_INCLUDED
#define TYPES_INCLUDED

#include "y/mkl/fitting/gls.hpp"
#include "y/mkl/fitting/sample/load.hpp"


#include "y/mkl/fcn/functions.hpp"
#include "y/type/physics.hpp"
#include "y/lua++/state.hpp"
#include "psi.hpp"
#include "y/alea.hpp"


using namespace upsylon;
using namespace mkl;
using namespace ODE;

typedef vector<double>            Vector;
typedef addressable<double>       Array;
typedef const accessible<double>  ConstArray;


typedef fitting::variables        Variables;
typedef Variables::const_iterator VariablesIterator;
typedef fitting::gls              GLS;
typedef GLS::sample_type          Sample;
typedef GLS::sample_ptr           SharedSample;
typedef GLS::scheme_type          SchemeType;
typedef GLS::scheme_ptr           SchemePointer;
typedef GLS::explode_type         ExplodeType;
typedef GLS::explode_solver_ptr   SolverPointer;

//! define the time scales
#define TIME_SCALE 1e-3

//! H in M to H in mM
#define PROTON_SCALE    1000.0

// warning: concentrations are in mM
// warning: potential is in       mV



inline void LoadParameters(Lua::VM         &vm,
                           const Variables &vars,
                           Array           &aorg)
{
    std::cerr << "<Loading Parameters>" << std::endl;
    for( VariablesIterator i=vars.begin();i!=vars.end();++i)
    {
        const string &id = (**i).name;
        std::cerr << "\t<" << id << " = ";
        if(vm->exists(id))
        {
            vars(aorg,id) = vm->get<double>( id );
            std::cerr << vars(aorg,id);
        }
        else
        {
            vars(aorg,id) = 0;
            std::cerr << " UNDEFINED, set to 0";
            
        }
        std::cerr << ">" << std::endl;
    }
    std::cerr << "<Loading Parameters/>" << std::endl;

}

#define MAKE_INDX() __LINE__ - mark

struct Index
{
    static const size_t mark = __LINE__;
    static const size_t Vm   = MAKE_INDX();
    static const size_t Li6  = MAKE_INDX();
    static const size_t Li7  = MAKE_INDX();
    static const size_t Q    = MAKE_INDX();
    static const size_t A    = MAKE_INDX();
    static const size_t Max  = MAKE_INDX()-1;
};

 


#include "y/fs/vfs.hpp"

inline void makeNames(const string &fileName,
                      string &baseName,
                      string &saveName,
                      string &funcName,
                      string &logName)
{
    baseName = vfs::get_base_name(fileName);
    saveName = vfs::with_new_extension(baseName,"fit.txt");
    funcName = vfs::with_new_extension(baseName,"fcn.txt");
    logName  = vfs::with_new_extension(baseName,"fit.log");

    std::cerr << "load '" << fileName << "'" << std::endl;
    std::cerr << "save '" << saveName << "'" << std::endl;
    std::cerr << "func '" << funcName << "'" << std::endl;
    std::cerr << "log  '" << logName  << "'" << std::endl;

}


#endif
