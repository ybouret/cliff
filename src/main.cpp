
#include "ko.hpp"
#include "nhe.hpp"

#include "y/program.hpp"
#include "y/ios/ocstream.hpp"

#include <iomanip>

#define I_SHOW(VAR) std::cerr << "Index::" << std::setw(4) << #VAR << " = " << Index:: VAR << std::endl


enum Task
{
    KO,
    NHE
};

#define RETURN_TASK(ID) if(#ID==task_id) return ID;
static inline Task GetTask( const string &task_id )
{
    RETURN_TASK(KO);
    RETURN_TASK(NHE);
    throw exception("unknown task '%s'", *task_id);
}




Y_PROGRAM_START()
{
    std::cerr << std::setprecision(15);

    I_SHOW(Vm);
    I_SHOW(Li6);
    I_SHOW(Li7);
    I_SHOW(Q);
    I_SHOW(Max);
    
    //==========================================================================
    //
    // process arguments
    //
    //==========================================================================
    if(argc<=1) throw exception("%s task(=KO|NHE) [ config.lua [delta [total]]]",program);
    const string task_id    = argv[1];
    const Task   task       = GetTask(task_id);

    //==========================================================================
    //
    // declare all variables
    //
    //==========================================================================
    const string   potential_list = "V0:gamma";
    const string   passive_list   = "k7:sigma:kQ:Qout:Qini:kA_ratio:Aout:Aini";
    const string   lithium_list   = "Lambda:E0:K5:f7:decay:eta_f:eta_r:phi:Z5";


    //==========================================================================
    //
    // Lua Virtual Machine for I/O
    //
    //==========================================================================
    Lua::VM        vm = new Lua::State();
    if(argc>2)
    {
        vm->doFile(argv[2]);
    }

    //==========================================================================
    //
    // create variables and read from VM
    //
    //==========================================================================

    Variables    vars = potential_list + ':' + passive_list + ':' + lithium_list;  
    const size_t nvar = vars.size();
    Vector       aorg(nvar,0);
    LoadParameters(vm,vars,aorg);

    //==========================================================================
    //
    // create cells/scheme
    //
    //==========================================================================
    XCell::Pointer cellDelta   = new XCellDelta(vm);
    XCell::Pointer cellTotal   = new XCellTotal(vm);
    SchemePointer  schemeDelta = & *cellDelta;
    SchemePointer  schemeTotal = & *cellTotal;

    //==========================================================================
    //
    // create solver
    //
    //==========================================================================
    SolverPointer  solver      = DriverCK<double>::New();
    solver->eps = vm->get<double>("ftol");

    //==========================================================================
    //
    // create ExplODE
    //
    //==========================================================================
    ExplodeType xDelta(solver,schemeDelta);
    ExplodeType xTotal(solver,schemeTotal);



    cellDelta->show();
    {
        std::cerr << "current variables:" << std::endl;
        //vars.display(std::cerr,aorg);
        fitting::display_variables::values(std::cerr, NULL, vars, aorg, NULL);
    }

    
    switch(task)
    {
        case KO:  DoFit_KO (argc,argv,aorg,vars,xTotal,xDelta,*cellTotal); break;
        case NHE: DoFit_NHE(argc,argv,aorg,vars,xTotal,xDelta,*cellTotal); break;
    }

    //std::cerr << "#terms=" << cellDelta->loop.count << std::endl;



}
Y_PROGRAM_END()

