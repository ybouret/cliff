#ifndef XCELL_INCLUDED
#define XCELL_INCLUDED 1

#include "pcell.hpp"
#include "y/tensor/tensor4d.hpp"
#include "y/tensor/loops.hpp"
#include "y/sort/sorted-sum.hpp"
#include "y/code/utils.hpp"
#include "y/core/ipower.hpp"


#define K5_SCALING 1e-6
#define f7_SCALING 1e3

struct Name
{
    static const size_t mark = __LINE__;
    static const size_t h   = MAKE_INDX();
    static const size_t Li6 = MAKE_INDX();
    static const size_t Li7 = MAKE_INDX();
    static const size_t N   = MAKE_INDX()-1;
};


typedef tensor4d<double>     Tensor4;
typedef tensor3d<double>     Tensor3;
typedef tensor2d<double>     Tensor2;
typedef tensor1d<double>     Tensor1;


struct DimerON
{
    int    i,j,k,l;
    double U,V;
};

class XCell : public PCell, public SchemeType
{
public:
    static const size_t N  = Name::N;
    static const size_t N3 = N*N*N;
    static const size_t N4 = N*N*N*N;

    typedef arc_ptr<XCell> Pointer;

    inline virtual ~XCell() throw() {}

    const double        ctrl;       //!< initial time step
    mutable Vector      Xa;         //!< for tensor evaluation
    mutable Vector      Xb;         //!< for tensor evaluation
    mutable Tensor4     K;
    mutable Tensor4     U;
    mutable Tensor4     V;
    mutable Vector      sum_a;
    mutable Vector      sum_b;
    mutable tensor_loop loop;
    //mutable double      E0;
    //mutable double      K0;
    
    explicit XCell(Lua::VM &_vm) :
    PCell(_vm),
    LdLua(ctrl),
    Xa(N,0),       // inside
    Xb(N,0),       // outside
    K(N,N,N,N),
    U(N,N,N,N),
    V(N,N,N,N),
    sum_a(N3,0),
    sum_b(N3,0),
    loop(U)
    {

        U.ld(0);
        V.ld(0);
        K.ld(0);

    }

    // Scheme interface
    size_t dimensions() const throw() { return Index::Max; }

    //! get the initial coordinate
    virtual double      start() const throw() { return 0; }

    //! set the initial state ( at start() )
    virtual void   setup(Array           &Y,
                         ConstArray      &atry,
                         const Variables &vars) const
    {
        tao::ld(Y,0);
        Y[Index::Vm] = vars(atry,"V0");
        Y[Index::Q]  = vars(atry,"Qini");
        Y[Index::A]  = vars(atry,"Aini");
    }



    static inline
    double sum_of( const Vector &s )
    {
        double ans = 0;
        for(size_t i=s.size();i>0;--i)
        {
            ans += s[i];
        }
        return ans;
    }

    inline void fillRates(const double f7,
                          const double f6,
                          const double r7,
                          const double r6,
                          const double K5) const throw()
    {


        const DimerON dimers[] =
        {
            // forward
            { Name::h, Name::Li6, Name::h, Name::Li6, f6, f6 },
            { Name::h, Name::Li6, Name::h, Name::Li7, f6, f7 },
            { Name::h, Name::Li7, Name::h, Name::Li6, f7, f6 },
            { Name::h, Name::Li7, Name::h, Name::Li7, f7, f7 },

            // reverse
            { Name::Li6, Name::h, Name::Li6, Name::h, r6, r6 },
            { Name::Li6, Name::h, Name::Li7, Name::h, r6, r7 },
            { Name::Li7, Name::h, Name::Li6, Name::h, r7, r6 },
            { Name::Li7, Name::h, Name::Li7, Name::h, r7, r7 }
        };

        // specific flip constants
        for(size_t d=0;d<sizeof(dimers)/sizeof(dimers[0]);++d)
        {
            const DimerON &D = dimers[d];
            const size_t   i = D.i;
            const size_t   j = D.j;
            const size_t   k = D.k;
            const size_t   l = D.l;
            U[i][j][k][l]    = D.U;
            V[i][j][k][l]    = D.V;
            K[i][j][k][l]    = K5;
        }


    }

    inline double compute_Xi() const
    {
        double _[N4] = { 0 };
        lightweight_array<double> arr(_,N4);
        int    indx      = 0;
        for(size_t l=N;l>0;--l)
        {
            const double Xb_l = Xb[l];
            for(size_t k=N;k>0;--k)
            {
                const double Xa_k = Xa[k] * Xb_l;
                for(size_t j=N;j>0;--j)
                {
                    const double Xb_j = Xb[j] * Xa_k;
                    for(size_t i=N;i>0;--i)
                    {
                        arr[++indx] = K[i][j][k][l] * Xa[i] * Xb_j;
                    }
                }
            }
        }
        hsort(arr,comparison::decreasing<double>);
        double ans = 0;
        for(size_t i=arr.size();i>0;--i)
        {
            ans += arr[i];
        }
        return ans;
    }

    inline void computeRatesOverE2(double      &Li6RateOverE2,
                                   double      &Li7RateOverE2,
                                   const double K5 ) const
    {
        double _[4] = {0,0,0,0};
        for(size_t i=Name::Li6;i<=Name::Li7;++i)
        {
            size_t ipos = 0;
            for(size_t j=N;j>0;--j)
            {
                const double Xb_j = Xb[j];
                const double Xa_j = Xa[j];
                for(size_t k=N;k>0;--k)
                {
                    const double Xa_k = Xa[k];
                    for(size_t l=N;l>0;--l)
                    {
                        const double Xb_l = Xb[l];
                        ++ipos;
                        {
                            sum_a[ipos] = (U[i][j][k][l] + V[k][l][i][j]) * Xb_j * Xa_k * Xb_l;
                        }

                        {
                            sum_b[ipos] = (U[j][i][k][l] + V[k][l][i][j]) * Xa_j * Xa_k * Xb_l;
                        }
                    }
                }
            }
            assert(N3==ipos);
            const double sumA = Xa[i] * sorted_sum(sum_a); //std::cerr << "sumA=" << sumA << std::endl;
            const double sumB = Xb[i] * sorted_sum(sum_b); //std::cerr << "sumB=" << sumB << std::endl;

            _[i] = sumB-sumA;
        }

        Li6RateOverE2=_[Name::Li6]*K5;
        Li7RateOverE2=_[Name::Li7]*K5;

        //std::cerr << "RatesOverE2: " << Li6RateOverE2 << " | " << Li7RateOverE2 << " | K5=" << K5 << std::endl;

    }


    //! compute the additional rates
    virtual void   rates(Array          &dYdt,
                         const double     t,
                         ConstArray      &Y,
                         ConstArray      &atry,
                         const Variables &vars) const
    {
        //----------------------------------------------------------------------
        // set no increase
        //----------------------------------------------------------------------
        tao::ld(dYdt,0);

        //----------------------------------------------------------------------
        // append leaks
        //----------------------------------------------------------------------
        add_leaks(dYdt,t,Y,atry,vars);

        //----------------------------------------------------------------------
        // get variables and compute h
        //----------------------------------------------------------------------
        const double LdVAR(Lambda);
        const double th     = get_th(Lambda);
        const double pH_end = get_pH_end(Lambda);
        const double h_end  = pow(10.0,-pH_end)*PROTON_SCALE;
        const double w_end  = t/(t+th);
        const double w_ini  = th/(t+th);

        // three components inside
        const double h   = h_ini * w_ini + h_end * w_end; // in mM
        const double Li6 = Y[Index::Li6];
        const double Li7 = Y[Index::Li7];
        
        Xa[Name::h]   = h;     assert( Xa[Name::h] > 0 );
        Xa[Name::Li6] = max_of<double>(Li6,0);
        Xa[Name::Li7] = max_of<double>(Li7,0);


        // three components outside
        const double Li6out = Lambda * eps6;
        const double Li7out = Lambda * eps7;
        
        Xb[Name::h]   = h_out;  assert(Xb[Name::h]  >=0);
        Xb[Name::Li6] = Li6out; assert(Xb[Name::Li6]>=0);
        Xb[Name::Li7] = Li7out; assert(Xb[Name::Li7]>=0);
        


        const double E0 = vars(atry,"E0");

        if(E0>0)
        {
            const double f7    = vars(atry,"f7") * f7_SCALING;
            const double r7    = vars(atry,"decay") * f7;
            const double K5    = vars(atry,"K5") * K5_SCALING;
            const double eta_f = vars(atry,"eta_f");
            const double f6    = f7 * eta_f;
            const double eta_r = vars(atry,"eta_r");
            const double r6    = r7 * eta_r;

            const double Z5    = vars(atry,"Z5") * K5_SCALING; K.ld(Z5);
            fillRates(f7,f6,r7,r6,K5);



            // compute all rates from the model
            double Li6RateOverE2 = 0;
            double Li7RateOverE2 = 0;

            computeRatesOverE2(Li6RateOverE2,Li7RateOverE2,K5);



            const double Xi = compute_Xi();
            const double W  = 0.5*(1.0 + sqrt( 1.0 + 8.0 * E0 * Xi  ) );
            const double E  = E0/W;
            const double E2 = E*E;

            //const double phi   = vars(atry,"phi"); E2 *= phi;

            dYdt[Index::Li6] += Li6RateOverE2 * E2;
            dYdt[Index::Li7] += Li7RateOverE2 * E2;
        }

    }


    virtual double delta() const throw()
    {
        return ctrl;
    }

    inline double delta7(const double t, ConstArray &Y, ConstArray &aorg, const Variables &vars) const
    {
        if(t<=0)
        {
            vector<double> dYdt(Y.size(),0);
            rates(dYdt,0,Y,aorg,vars);
            return conc_ratio_to_delta(dYdt[Index::Li7]/dYdt[Index::Li6]);
        }
        else
        {
            return conc_ratio_to_delta(Y[Index::Li7]/Y[Index::Li6]);
        }
    }

    inline double Li(ConstArray      &Y,
                     ConstArray      &atry,
                     const Variables &vars
                     ) const throw()
    {
        const double sumLi = (Y[Index::Li6]+Y[Index::Li7]);
        const double phi   = vars(atry,"phi"); return  phi * sumLi;
    }

    inline void saveHeader(ios::ostream &fp) const
    {
        fp("#t Vm Li6 Li7 Q A delta7 Li\n");
    }

    inline void saveFrame(ios::ostream &fp, const double t, ConstArray &Y, ConstArray &aorg, const Variables &vars) const
    {
        fp("%.15g",t);
        for(size_t i=1;i<=Y.size();++i)
        {
            fp(" %.15g", Y[i]);
        }
        fp(" %.15g", delta7(t,Y,aorg,vars));
        fp(" %.15g", Li(Y,aorg,vars) );
        fp << '\n';
    }

    
private:
    Y_DISABLE_COPY_AND_ASSIGN(XCell);
};

class XCellDelta : public XCell
{
public:
    inline explicit XCellDelta(Lua::VM            &_vm) :
    XCell(_vm)
    {
    }

    inline virtual ~XCellDelta() throw()
    {
    }

    virtual double query(const double     t,
                         ConstArray      &Y,
                         ConstArray      &aorg,
                         const Variables &vars) const
    {
        return delta7(t,Y,aorg,vars);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(XCellDelta);
};

class XCellTotal : public XCell
{
public:
    inline explicit XCellTotal(Lua::VM &_vm) :
    XCell(_vm)
    {
    }

    inline virtual ~XCellTotal() throw()
    {
    }

    virtual double query(const double     ,
                         ConstArray      &Y,
                         ConstArray      &atry ,
                         const Variables &vars ) const
    {
        return Li(Y,atry,vars);
    }

private:
    Y_DISABLE_COPY_AND_ASSIGN(XCellTotal);
};





#endif

