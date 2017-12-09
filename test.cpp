#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "CharacteristicFunctions.h"


TEST_CASE("Test ODE", "[CF]"){
    std::vector<double> currentValues({.5, .5});
    auto sig=.3;
    auto a=.3;
    auto b=1;
    auto lambda=.5;
    auto delta=.9;
    auto cf=[](const auto& u){
        return u;
    };
    auto u=.2;
    auto rho0=0.0;
    auto rho1=0.0;
    auto K0=a*b;
    auto K1=-a;
    auto H0=0.0;
    auto H1=sig*sig;
    auto l0=0.0;
    auto l1=lambda;
    
    auto sigma=sig*sig*.5;
    auto cfPart=cf(u+delta*currentValues[0])-1.0;
    REQUIRE(-chfunctions::explSol(currentValues[0], rho1, K1, H1, l1, cfPart)==currentValues[0]*currentValues[0]*sigma+cf(u+currentValues[0]*delta)*lambda-lambda-currentValues[0]*a);
    REQUIRE(-chfunctions::explSol(currentValues[0], rho0, K0, H0, l0, cfPart)==currentValues[0]*a*b);
}

TEST_CASE("Test CIR", "[CF]"){
    
    auto sig=.3;
    auto a=.3;
    auto b=.05;
    auto r0=.05;
    auto h=sqrt(a*a+2*sig*sig);
    auto T=1.0;
    auto aNum=2*h*exp((a+h)*T*.5);
    auto aDen=2*h+(a+h)*(exp(T*h)-1.0);
    auto AtT=pow(aNum/aDen, (2*a*b)/(sig*sig));
    auto bNum=2*(exp(T*h)-1.0);
    auto bDen=aDen;
    auto BtT=bNum/bDen;
    auto BondPrice=AtT*exp(-BtT*r0);


    auto rho1=1.0;
    auto k0=a*b;
    auto k1=-a;
    auto H1=sig*sig;
    auto approxBondPrice=chfunctions::expAffine(
        rungekutta::computeFunctional(T, 2048, std::vector<double >({0, 0}),
            [&](double t, const std::vector<double>& x){
                return chfunctions::duffieODE(
                    x, //current values
                    0.0, //rho0
                    rho1,
                    k0,
                    k1, 
                    0.0, //H0
                    H1, //H1, 
                    0.0,//l0
                    0.0, //l1
                    [&](const auto& uhat){
                        return 0;
                    }
                );
            }
        ),
        r0
    );
    REQUIRE(approxBondPrice==Approx(BondPrice));
}
TEST_CASE("Test CIR with curried function", "[CF]"){
    
    auto sig=.3;
    auto a=.3;
    auto b=.05;
    auto r0=.05;
    auto h=sqrt(a*a+2*sig*sig);
    auto T=1.0;
    auto aNum=2*h*exp((a+h)*T*.5);
    auto aDen=2*h+(a+h)*(exp(T*h)-1.0);
    auto AtT=pow(aNum/aDen, (2*a*b)/(sig*sig));
    auto bNum=2*(exp(T*h)-1.0);
    auto bDen=aDen;
    auto BtT=bNum/bDen;
    auto BondPrice=AtT*exp(-BtT*r0);


    auto rho1=1.0;
    auto k0=a*b;
    auto k1=-a;
    auto H1=sig*sig;
    auto beta=chfunctions::AlphaOrBeta(rho1, k1, H1, 0.0);
    auto alpha=chfunctions::AlphaOrBeta(0.0, k0, 0.0, 0.0);
    auto approxBondPrice=chfunctions::expAffine(
        rungekutta::computeFunctional(T, 2048, std::vector<double >({0, 0}),
            [&](double t, const std::vector<double>& x){
                return std::vector<double>({
                    beta(x[0], -1.0),
                    alpha(x[0], -1.0)
                });
            }
        ),
        r0
    );
    REQUIRE(approxBondPrice==Approx(BondPrice).epsilon(.0001));
}

/*TEST_CASE("Test CGMY", "[CF]"){

}*/