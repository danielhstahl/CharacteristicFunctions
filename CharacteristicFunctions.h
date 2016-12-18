#ifndef __CHARACTERISTICFUNCTIONS_H_INCLUDED__
#define __CHARACTERISTICFUNCTIONS_H_INCLUDED__

#include <complex>
#include "FunctionalUtilities.h" 
#include "RungeKutta.h" 


/**
    All these Characteristic Functions are with respect to ui, not u. Hence the Gaussian CF, for example, is exp(u*mu+u*u*sigma*sigma*.5) instead of exp(u*mu*i-u*u*sigma*sigma*.5)

*/
namespace chfunctions { 
    auto gaussCF(const auto& u, const auto& mu, const auto& sigma){
        return exp(u*mu+u*u*sigma*sigma*.5);
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
    template<typename T>
    auto stableCF(const std::complex<T>& u, const auto& alpha, const auto& mu, const auto& beta, const auto& c){
        auto phi=tan(alpha*.5*M_PI);
        return exp(u*mu-pow(u*std::complex<T>(0, -1)*c, alpha)*std::complex<T>(1, -beta*phi));
    }
    auto gammaCF(const auto& u, const auto& a, const auto& b){
        return pow(1-u*b, -a);
    }
    auto inverseGaussianCF(const auto& u, const auto& mu, const auto& lambda){
        return exp((lambda/mu)*(1-sqrt(1-(2*mu*mu*u)/lambda)));
    }

    /**
    Helper function to compute ODE series found in http://web.stanford.edu/~duffie/dps.pdf page 10. THIS NEEDS TO BE MORE GENERIC!!!  Only a subset of duffie's ODEs
    @u complex value for phi(u)
    @cf characteristic function of the "jump" process
    @currentValues currentValues of the the system of ODEs
    @sigma volatility of diffusion
    @lambda 
    */
   /* template<typename T>
    std::vector<T > duffieODE(const auto& u, auto&& cf, const std::vector<T >& currentValues, double sigma, double lambda, double a, double delta, double b){ //double alpha, double mu, double beta, double c,
        auto sig=sigma*sigma*.5;
        return 
        {
            currentValues[0]*currentValues[0]*sig+cf(u+currentValues[0]*delta)*lambda-lambda-currentValues[0]*a, 
            currentValues[0]*b*a
        };
    }*/


    /**
    Equation for Beta'(t, T) or Alpha'(t, T)
    */
    auto ODE(const auto& currVal, const auto& rho, const auto& K, const auto& H, const auto& l, const auto& cfPart){
        return rho-K*currVal-.5*currVal*currVal*H-l*cfPart;
    }
    /**Helper function to compute ODE series found in http://web.stanford.edu/~duffie/dps.pdf page 10. Because of a "measure change" the addition parameter "u" is introduced.  
    note that this is with respect to T-t not t so the equations have signs switched
    */
    template<typename T>
    std::vector<T > duffieODE(const auto& u, const std::vector<T >& currentValues, const auto& rho0, const auto& rho1, const auto& K0, const auto& K1, const auto& H0, const auto& H1, const auto& l0, const auto& l1, auto&& cf){ //double alpha, double mu, double beta, double c,
        auto cfPart=cf(u)-1.0;
        return //beta, alpha
        {
            -ODE(currentValues[0], rho1, K1, H1, l1, cfPart),
            -ODE(currentValues[0], rho0, K0, H0, l0, cfPart)
        };
    }
    
    /**Helper function to compute ODE series found in http://web.stanford.edu/~duffie/dps.pdf page 10. 
    note that this is with respect to T-t not t so the equations have signs switched
    */
    template<typename T>
    std::vector<T > duffieODE(const std::vector<T >& currentValues, const auto& rho0, const auto& rho1, const auto& K0, const auto& K1, const auto& H0, const auto& H1, const auto& l0, const auto& l1, auto&& cf){ //double alpha, double mu, double beta, double c,
        //auto sig=sigma*sigma*.5;
        return duffieODE(currentValues[0], currentValues, rho0, rho1, K0, K1, H0, H1, l0, l1, cf);
    }


    /**
    Helper function to compute expontential of the Duffie ODE
    */
    template<typename T>
    T expAffine(const std::vector<T>& vals, const auto& v0){
        return exp(vals[0]*v0+vals[1]);
    }
}
#endif