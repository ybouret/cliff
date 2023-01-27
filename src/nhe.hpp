
#ifndef NHE_INCLUDED
#define NHE_INCLUDED 1

#include "xcell.hpp"
#include "ko.hpp"
#include "y/jive/pattern/matching.hpp"

static inline
void save_delta7(const string &saveName,
                 const string &funcName,
                 const Sample &sample,
                 ExplodeType  &xDelta7,
                 ConstArray   &aorg,
                 const XCell  &cell)
{

    // save data
    {
        ConstArray &t      = *sample.abscissa;
        ConstArray &delta7 = *sample.ordinate;
        ConstArray &d7fit  = *sample.adjusted;
        ios::ocstream fp(saveName);
        fp("#t delta7 d7fit\n");
        for(size_t i=1;i<=sample.count();++i)
        {
            fp("%.15g %.15g %.15g\n",t[i],delta7[i],d7fit[i]);
        }
    }

    // save as function
    {
        ios::ocstream    fp(funcName);
        cell.saveHeader(fp);

        //const double     lo   = 0;
        const double     hi   = sample.abscissa->back() * 1.1;
        const double     dt   = 5; //(hi-lo)/200;
        const Variables &vars = *sample;

        double      t =  0;
        (void) xDelta7.start(t,aorg,vars);
        ConstArray &Y = *xDelta7;

        cell.saveFrame(fp,t,Y,aorg,vars);
        for(t+=dt;t<=hi;t+=dt)
        {
            (void)xDelta7.reach(t,aorg,vars);
            cell.saveFrame(fp,t,Y,aorg,vars);
        }
    }


}

using namespace Jive;

static inline void fillSequences(sequence<string> &strings,
                                 const string     &name,
                                 Lua::VM          &vm)
{
    strings.free();

    lua_State *L = **vm;
    lua_settop(L,0);
    lua_getglobal(L,*name);
    if(!lua_istable(L,-1))
    {
        throw exception("No table '%s'", *name);
    }
    const unsigned n = unsigned(lua_rawlen(L,-1));
    for(size_t i=1;i<=n;++i)
    {
        lua_rawgeti(L,-1,i);
        if( lua_isstring(L,-1) )
        {
            const string tmp = lua_tostring(L,-1);
            strings.push_back(tmp);
        }
        lua_pop(L,1);
    }
}

static inline
void DoFit_NHE(int             &argc,
               char          **&argv,
               Array           &aorg,
               const Variables &vars,
               ExplodeType     &xIntake,
               ExplodeType     &xDelta7,
               const XCell     &cell)
{
    const string resultsName = "active.lua";

    if(argc<=3) return;

    // prepare files
    bool hasIntake = false;
    bool hasDelta7 = false;

    string intakeFileName;
    string intakeBaseName;
    string intakeSaveName;
    string intakeFuncName;
    string intakeLogName;

    string delta7FileName;
    string delta7BaseName;
    string delta7SaveName;
    string delta7FuncName;
    string delta7LogName;

    {
        Matching  matchIntake = "(intake)&";
        Matching  matchDelta7 = "(delta7)&";

        for(int i=3;i<argc;++i)
        {
            const string fileName = argv[i];
            if( matchIntake.isFoundIn(fileName) )
            {
                std::cerr << "<intake>" << std::endl;
                if(hasIntake) throw exception("multiple intake file");
                hasIntake      = true;
                intakeFileName = fileName;
                makeNames(intakeFileName, intakeBaseName, intakeSaveName, intakeFuncName, intakeLogName);
                std::cerr << "<intake/>" <<std::endl;
                continue;
            }

            if( matchDelta7.isFoundIn(fileName) )
            {
                std::cerr << "<delta7>" << std::endl;
                if(hasDelta7) throw exception("multiple delta7 file");
                hasDelta7      = true;
                delta7FileName = fileName;
                makeNames(delta7FileName, delta7BaseName, delta7SaveName, delta7FuncName, delta7LogName);
                std::cerr << "<delta7s/>" <<std::endl;
                continue;
            }

        }
    }

    SharedSample intakeSample_ = Sample::create("intake", 100);
    Sample      &intakeSample = *intakeSample_;
    *intakeSample = vars;

    SharedSample delta7Sample_ = Sample::create("delta7", 100);
    Sample      &delta7Sample  = *delta7Sample_;
    *delta7Sample = vars;


    // load
    size_t intakeN = 0;
    if( hasIntake )
    {
        intakeN = fitting::load_sample::from_file(intakeFileName,intakeSample,1,2);
        std::cerr << "#intake=" << intakeN << std::endl;
        if(intakeN<=0) throw exception("Not enough data in '%s'", *intakeFileName);
    }

    size_t delta7N = 0;
    if( hasDelta7 )
    {
        delta7N = fitting::load_sample::from_file(delta7FileName,delta7Sample,1,2);
        std::cerr << "#delta7=" << delta7N << std::endl;
        if(delta7N<=0) throw exception("Not enough data in '%s'", *delta7FileName);
    }



    // prepare fit
    GLS LS(false);
    LS.grad().h = 1e-4;
    
    // initialize
    const size_t nvar = aorg.size();
    Vector       aerr(nvar,0);
    vector<bool> used(nvar,false);


    static const char *pfx = "\t(*) ";


#define NHE_INTAKE(EXPR) do {                                                     \
/**/vars.only_on(used,EXPR);                                                      \
/**/if(LS.fit(intakeSample,xIntake,aorg,used,aerr))                               \
/**/{                                                                             \
/**/    fitting::display_variables::errors(std::cerr,pfx,vars, aorg,used,aerr);   \
/**/    save_intake(intakeSaveName,intakeFuncName,intakeSample,xIntake,aorg,cell);\
/**/    ios::ocstream logFile(intakeLogName);                                     \
/**/    fitting::display_sample::results(logFile,intakeSample,aorg,used,aerr);    \
/**/    std::cerr << std::endl;                                                   \
/**/} else {                                                                      \
/**/    throw exception("NHE_INTAKE failure for " #EXPR);                         \
/**/} } while(false)

#define NHE_DELTA7(EXPR) do {                                                     \
/**/vars.only_on(used,EXPR);                                                      \
/**/if(LS.fit(delta7Sample,xDelta7,aorg,used,aerr))                               \
/**/{                                                                             \
/**/    fitting::display_variables::errors(std::cerr,pfx,vars, aorg,used,aerr);   \
/**/    save_delta7(delta7SaveName,delta7FuncName,delta7Sample,xDelta7,aorg,cell);\
/**/    ios::ocstream logFile(delta7LogName);                                     \
/**/    fitting::display_sample::results(logFile,delta7Sample,aorg,used,aerr);    \
/**/    std::cerr << std::endl;                                                   \
/**/} else {                                                                      \
/**/    throw exception("NHE_DELTA7 failure for " #EXPR);                         \
/**/} } while(false)

    const size_t loops = max_of<size_t>(1,aliasing::_(cell.vm)->get<int>("loops"));
    for(size_t cycle=1;cycle<=loops;++cycle)
    {

        if(hasDelta7)
        {
            NHE_DELTA7("");
            vector<string> seqs;
            fillSequences(seqs, "delta7", aliasing::_(cell.vm) );
            for(size_t i=1;i<=seqs.size();++i)
            {
                NHE_DELTA7(seqs[i]);
            }
            std::cerr << "used for Delta7: " << seqs << std::endl;
        }

        if(hasIntake)
        {
            NHE_INTAKE("");
            vector<string> seqs;
            fillSequences(seqs, "intake", aliasing::_(cell.vm) );
            for(size_t i=1;i<=seqs.size();++i)
            {
                NHE_INTAKE(seqs[i]);
            }
            std::cerr << "used for Intake: " << seqs << std::endl;

            if(hasDelta7)
            {
                NHE_DELTA7("");
            }
        }

    }
    

    // save fitted for consistency loop
    {
        ios::ocstream fp( resultsName );
        fp << "-- results from NHE fits, generated @" << string_time::stamp() << '\n';
        for(Variables::const_iterator it=vars.begin();it!=vars.end();++it)
        {
            const fitting::variable &v = **it;
            fp << v.name << " = ";
            fp("%.15g", v(aorg) );
            fp << '\n';
        }
    }

}

#endif

