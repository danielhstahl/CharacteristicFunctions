#ifndef __CHARACTERISTICFUNCTIONS_H_INCLUDED__
#define __CHARACTERISTICFUNCTIONS_H_INCLUDED__

#include <complex>
#include "FunctionalUtilities.h" 
#include "RungeKutta.h" 


/**
    All these Characteristic Functions are with respect to ui, not u. Hence the Gaussian CF, for example, is exp(u*mu+u*u*sigma*sigma*.5) instead of exp(u*mu*i-u*u*sigma*sigma*.5)

*/
namespace chfunctions { 


    template<typename T, typename Number>
    auto gaussLogCF(const T& u, const Number& mu, const Number& sigma){
        return u*mu+.5*futilities::const_power(sigma*u, 2);
    }
    template<typename T, typename Number>
    auto gaussCF(const T& u, const Number& mu, const Number& sigma){
        return exp(gaussLogCF(u, mu, sigma));
    }
    template<typename T, typename Number>
    auto mertonLogCF(const T& u, const Number& lambda, const Number& muL, const Number& sigL){
        return lambda*(gaussCF(u, muL, sigL)-1.0);
    }
    template<typename T, typename Number>
    auto mertonLogRNCF(const T& u, const Number& lambda, const Number& muL, const Number& sigL, const Number& r,  const Number& sigma){
        return gaussLogCF(u, r-futilities::const_power(sigma, 2)*.5-mertonLogCF(1.0, lambda, muL, sigL), sigma)+mertonLogCF(u, lambda, muL, sigL);
    }

    //see http://finance.martinsewell.com/stylized-facts/distribution/CarrGemanMadanYor2002.pdf pg 10
    template<typename U,  typename Number>
    auto cgmyLogCF(const U& u, const Number& C, const Number& G, const Number& M, const Number& Y){

        return Y==1.0?
            0.0:
            Y==0.0?
            C*log((1.0-u/G)*(1.0+u/M)):
            C*tgamma(-Y)*(pow(M-u, Y)+pow(G+u, Y)-pow(M, Y)-pow(G, Y));
    }
    //see http://finance.martinsewell.com/stylized-facts/distribution/CarrGemanMadanYor2002.pdf pg 12 and 13
    template<typename Number, typename T>
    auto cgmyLogRNCF(const std::complex<T>& u, const Number& C, const Number& G, const Number& M, const Number& Y, const Number& r,  const Number& sigma){
        return gaussLogCF(u, 
            r-futilities::const_power(sigma, 2)*.5-cgmyLogCF(1.0, C, G, M, Y),
            sigma
        )+cgmyLogCF(u, C, G, M, Y);
    }
    //note that this is of the form a-kappa*v NOT a(kappa-v)
    //see https://pdfs.semanticscholar.org/67cd/b553e2624c79a960ff79d0dfe6e6833690a7.pdf pg 14
    template<typename Psi, typename A, typename Kappa, typename Sigma,typename T, typename V>
    auto cirLogMGF(const Psi& psi,  const A& a, const Kappa& kappa, const Sigma& sigma,const T& t, const V& v0){
        const auto delta=sqrt(futilities::const_power(kappa, 2)+2.0*psi*futilities::const_power(sigma, 2));
        const auto expT=exp(-delta*t);
        const auto deltMinusKappa=delta-kappa;
        const auto bT=2.0*psi*(1.0-expT)/(delta+kappa+deltMinusKappa*expT);
        auto cT=sigma>0?(a/futilities::const_power(sigma, 2))*(2.0*log(
            1.0-deltMinusKappa*(1.0-expT)/(2.0*delta)
        )+deltMinusKappa*t):psi*(t-(1.0-expT)/kappa);
        return -bT*v0-cT;
    }
    template<typename Psi, typename A, typename Kappa, typename Sigma,typename T, typename V>
    auto cirMGF(const Psi& psi,  const A& a, const Kappa& kappa, const Sigma& sigma,const T& t, const V& v0){
        return exp(cirLogMGF(psi, a, kappa, sigma, t, v0));
    }

    /**
    CF of a stable distribution.  
    Stable distribution is defined by CF exp(iu*mu-abs(cu)^alpha (1-beta*sign(u)*tan(pi*alpha*.5)))
    @u complex value for phi(u)
    @alpha parameter for stable
    @c parameter for stable
    @mu parameter for stable
    @beta parameter for stable
    @returns CF for stable at u
    @
    */
    template<typename T, typename Number>
    auto stableCF(const std::complex<T>& u, const Number& alpha, const Number& mu, const Number& beta, const Number& c){
        auto phi=tan(alpha*.5*M_PI);
        return exp(u*mu-pow(u*std::complex<T>(0, -1)*c, alpha)*std::complex<T>(1, -beta*phi));
    }
    template<typename Number1, typename Number2>
    auto gammaCF(const Number1& u, const Number2& a, const Number2& b){
        return pow(1-u*b, -a);
    }
    
    template<typename Number1, typename Number2>
    auto inverseGaussianCF(const Number1& u, const Number2& mu, const Number2& lambda){
        return exp((lambda/mu)*(1-sqrt(1-(2*mu*mu*u)/lambda)));
    }

    template<typename Number1, typename Number2>
    auto exponentialCF(const Number1& u, const Number2& lambda){
        return lambda/(lambda-u);
    }


    /**
    Explicit "solution" for Beta'(t, T) or Alpha'(t, T)
    */
    template<typename Number, typename Rho, typename KType, typename HType, typename L, typename CFPart>
    auto explSol(const Number& currVal, const Rho& rho, const KType& K, const HType& H, const L& l, const CFPart& cfPart){
        return rho-K*currVal-.5*currVal*currVal*H-l*cfPart;
    }
    
    /**Helper function to compute ODE series found in http://web.stanford.edu/~duffie/dps.pdf page 10. Because of a "measure change" the addition parameter "u" is introduced.  
    note that this is with respect to T-t not t so the equations have signs switched
    */
    template<typename T, typename Number1, typename Number2, typename CF>
    std::vector<T > duffieODE(
        const Number1& u, 
        const std::vector<T >& currentValues, 
        const Number2& rho0, 
        const Number2& rho1, 
        const Number2& K0, 
        const Number2& K1, 
        const Number2& H0, 
        const Number2& H1, 
        const Number2& l0, 
        const Number2& l1, 
        CF&& cf
    ){ //double alpha, double mu, double beta, double c,
        auto cfPart=cf(u)-1.0;
        return //beta, alpha
        {
            -explSol(currentValues[0], rho1, K1, H1, l1, cfPart),
            -explSol(currentValues[0], rho0, K0, H0, l0, cfPart)
        };
    }

    /**Curried function.  Can be called for either Alpha' or Beta'*/
    template<typename T>
    auto AlphaOrBeta(const T& rho, const T& K, const T& H, const T& l){
        return [&](const auto& val, const auto& cf){
            return -explSol(val, rho, K, H, l, cf);
        };
    }

    /**Curried function for stableCF.*/
    template<typename Number>
    auto augCF(const Number& alpha, const Number& mu, const Number& beta, const Number& c){
        return [&](const auto& u){
            return stableCF(u, alpha, mu, beta, c)-1.0;
        };
    }
    
    /**Helper function to compute ODE series found in http://web.stanford.edu/~duffie/dps.pdf page 10. 
    note that this is with respect to T-t not t so the equations have signs switched
    */
    template<typename T, typename Number1, typename CF>
    std::vector<T > duffieODE(
        const std::vector<T >& currentValues, 
        const Number1& rho0, 
        const Number1& rho1, 
        const Number1& K0, 
        const Number1& K1, 
        const Number1& H0, 
        const Number1& H1, 
        const Number1& l0, 
        const Number1& l1, 
        CF&& cf){ //double alpha, double mu, double beta, double c,
        //auto sig=sigma*sigma*.5;
        return duffieODE(currentValues[0], currentValues, rho0, rho1, K0, K1, H0, H1, l0, l1, cf);
    }


    /**
    Helper function to compute expontential of the Duffie ODE
    */
    template<typename T, typename Number>
    T logAffine(const std::vector<T>& vals, const Number& v0){
        return vals[0]*v0+vals[1];
    }
    /**
    Helper function to compute expontential of the Duffie ODE
    */
    template<typename T, typename Number>
    T expAffine(const std::vector<T>& vals, const Number& v0){
        return exp(logAffine(vals, v0));
    }
}
#endif
