#ifndef KO_INCLUDED
#define KO_INCLUDED 1

#include "xcell.hpp"
#include "y/ios/ocstream.hpp"
#include "y/mkl/solve/zroot.hpp"
#include "y/string/time.hpp"
#include "y/mkl/fitting/sample/display.hpp"
#include "y/mkl/fitting/variable/display.hpp"
#include "y/mkl/fitting/sample/load.hpp"

struct DeltaFinder
{
    double           t_find; //!< should be 60 seconds
    double           zvalue; //!< value at t_find
    ExplodeType     *pExpl;  //!< handle to delta
    const Variables *pVars;  //!< handle to vars
    Array           *pAorg;  //!< handle to params

    double operator()( double sigma )
    {
        assert(pVars);
        assert(pAorg);
        assert(pExpl);
        const Variables &vars   = *pVars;
        Array           &aorg   = *pAorg;
        double          &target = vars( aorg, "sigma" );
        const double     backup = target;
        target = sigma;
        (void) pExpl->start(0, aorg, vars);
        const double ans = pExpl->reach(t_find, aorg, vars)-zvalue;
        target = backup;
        return ans;
    }

};


static inline
void save_intake(const string &saveName,
                 const string &funcName,
                 const Sample &sample,
                 ExplodeType  &xTotal,
                 ConstArray   &aorg,
                 const XCell  &cell)
{

    // save data
    std::cerr << "<save intake/data to '" << saveName << "'>" << std::endl;
    {
        ConstArray &t      = *sample.abscissa;
        ConstArray &Li     = *sample.ordinate;
        ConstArray &LiFit  = *sample.adjusted;
        ios::ocstream fp(saveName);
        fp("#t [Li] [LiFit]\n");
        for(size_t i=1;i<=sample.count();++i)
        {
            fp("%.15g %.15g %.15g\n",t[i],Li[i],LiFit[i]);
        }
    }

    // save as function
    std::cerr << "<save intake/func to '" << funcName << "'>" << std::endl;
    {
        ios::ocstream    fp(funcName);
        cell.saveHeader(fp);

        //const double     lo   = 0;
        const double     hi   = sample.abscissa->back()*1.1;
        const double     dt   = 5; //(hi-lo)/200;
        const Variables &vars = *sample;

        double      t =  0;
        (void) xTotal.start(t,aorg,*sample);
        ConstArray &Y = *xTotal;

        cell.saveFrame(fp,t,Y,aorg,vars);
        for(t+=dt;t<=hi;t+=dt)
        {
            (void)xTotal.reach(t,aorg,vars);
            cell.saveFrame(fp,t,Y,aorg,vars);
        }
    }
    std::cerr << "..all saved: " << saveName << " AND " << funcName << std::endl;

}

static inline
void writeTo(ios::ostream    &fp,
             const Variables &vars,
             const char      *id,
             ConstArray      &aorg,
             ConstArray      &aerr)
{
    assert(id);
    fp("%s=%.15g -- +/- %.15g\n", id, vars(aorg,id), vars(aerr,id));

}

#include "y/fs/local/fs.hpp"

static inline
void DoFit_KO(int             &argc,
              char **         &argv,
              Array           &aorg,
              const Variables &vars,
              ExplodeType     &xTotal,
              ExplodeType     &xDelta,
              const XCell     &cell)
{
    const string resultsName = "passive.lua";
    vfs &fs = local_fs::instance();
    fs.try_remove_file(resultsName);

    if(argc<=3) return;

    const string fileName = argv[3];
    string       baseName;
    string       saveName;
    string       funcName;
    string       logName ;

    makeNames(fileName, baseName, saveName, funcName, logName);

    // prepare a sample
#if 0
    Series  t     = new Vector();
    Series  Li    = new Vector();
    Series  LiFit = new Vector();
    Sample  sample(t,Li,LiFit);
    sample.variables = vars;
#endif

    SharedSample sample = Sample::create("intake",100);
    **sample = vars;

    // load from file
    const size_t n = fitting::load_sample::from_file(fileName,*sample,1,2);
    std::cerr << "#data=" << n << std::endl;
    if(n<=0) throw exception("Not enough data in '%s'", *fileName);

    // prepare fit
    GLS LS(true);
    LS.grad().h = 1e-4;

    // initialize
    const size_t nvar = aorg.size();
    Vector       aerr(nvar,0);
    vector<bool> used(nvar,false);


    // perform sessions

#define KO_SESSION(EXPR) do {                                                     \
/**/vars.only_on(used,EXPR);                                                      \
/**/if(LS.fit(*sample, xTotal, aorg, used, aerr))                                 \
/**/{                                                                             \
/**/    fitting::display_variables::errors(std::cerr, NULL,vars, aorg,used,aerr); \
/**/    save_intake(saveName,funcName,*sample,xTotal,aorg,cell);                  \
/**/    ios::ocstream logFile(logName);                                           \
/**/    fitting::display_sample::results(logFile, *sample, aorg, used, aerr);     \
/**/    std::cerr << std::endl;                                                   \
/**/} else {                                                                      \
/**/    throw exception("KO session failure for " #EXPR);                         \
/**/} } while(false)



    KO_SESSION("");

    KO_SESSION("k7");
    KO_SESSION("gamma");
    KO_SESSION("k7:gamma");
    KO_SESSION("kQ");
    KO_SESSION("k7:gamma:kQ");
    
#if 0
    KO_SESSION("kA");
    KO_SESSION("kA:kQ");
    KO_SESSION("kA:kQ:k7");
#endif

    DeltaFinder F = { 60.0, 10.8, &xDelta, &vars, &aorg};
    double &sigma = vars(aorg,"sigma");
    std::cerr << "F(" << sigma << ")=" << F(sigma) << std::endl;

    const double    fac = 1.001;
    triplet<double> sig = { sigma, sigma, sigma };
    triplet<double> zfn = { 0,0,0 };

    do {
        zfn.a = F( sig.a /= fac );
        zfn.c = F( sig.c *= fac );
    } while( zfn.a * zfn.c >= 0 );

    zroot<double>::bisection solve;
    //if(!solve(F,sig,zfn)) throw exception("Corrupted Function!!!");
    sigma = solve(F,sig,zfn);
    std::cerr << "F(" << sigma << ")=" << F(sigma) << std::endl;

    {
        std::cerr << "Results for KO fit:" << std::endl;
        ios::ocstream fp( ios::cstderr );
        //sample.writeLog(fp,aorg,used,aerr);
        fitting::display_sample::results(fp, *sample, aorg, used, aerr);
        fp("sigma = %.15g\n", sigma);
    }



    // save fitted for consistency loop
    {
        ios::ocstream fp( resultsName );
        fp << "-- results from KO fits, generated @" << string_time::stamp() << '\n';
        writeTo(fp,vars,"k7",aorg,aerr);
        writeTo(fp,vars,"kQ",aorg,aerr);
        writeTo(fp,vars,"gamma",aorg,aerr);
        writeTo(fp,vars,"sigma",aorg,aerr);
    }





}



#endif

