#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "CharacteristicFunctions.h"

/*
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
*/
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
TEST_CASE("Test CIR with analytical", "[CF]"){
    
    auto sig=.3;
    auto a=.3;
    auto b=.05;
    auto r0=.05;
    auto h=sqrt(a*a+2*sig*sig);
    auto T=.25;
    auto aNum=2*h*exp((a+h)*T*.5);
    auto aDen=2*h+(a+h)*(exp(T*h)-1.0);
    auto AtT=pow(aNum/aDen, (2*a*b)/(sig*sig));
    auto bNum=2*(exp(T*h)-1.0);
    auto bDen=aDen;
    auto BtT=bNum/bDen;
    auto BondPrice=AtT*exp(-BtT*r0);

    auto approxBondPrice=chfunctions::cirMGF(1.0, a*b, a, sig, T, r0);
    REQUIRE(approxBondPrice==Approx(BondPrice));
}
TEST_CASE("Test CIR with Heston", "[CF]"){
    
    auto T=.25;
    auto k=.2;
    //auto a=.2;
    auto v0=.98;
    auto sig=.2;
    auto rho=-.3;
    auto sigTot=.3;
    auto u=std::complex<double>(.5, .5);
    auto negPsi=.5*sigTot*sigTot*(u-u*u);
    auto kStar=k-u*rho*sig*sigTot;
    auto ada=sqrt(kStar*kStar+2*sig*sig*negPsi);
    auto bT=2.0*negPsi*(1.0-exp(-ada*T))/(2.0*ada-(ada-kStar)*(1.0-exp(-ada*T)));
    auto cT=(k/(sig*sig))*(2.0*log(1.0-(1.0-exp(-ada*T))*(ada-kStar)/(2.0*ada))+(ada-kStar)*T);
    auto cfHeston=exp(-bT*v0-cT).real();

    auto approxBondPrice=chfunctions::cirMGF(negPsi, k, kStar, sig, T, v0);
    REQUIRE(approxBondPrice.real()==Approx(cfHeston));
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
TEST_CASE("Test CIR against runge kutta", "[CF]"){
    auto rho1=1.0;
    auto k0=.05; //long run average of .05/.3
    auto k1=.3;
    auto H1=.09;
    auto T=.5;
    auto r0=.15;
    auto beta=chfunctions::AlphaOrBeta(rho1, k1, H1, 0.0);
    auto alpha=chfunctions::AlphaOrBeta(0.0, k0, 0.0, 0.0);
    //auto bPrime=chfunctions::TimeChangeB(rho1, k1, H1, 1.0, 0.0);
    //auto cPrime=chfunctions::TimeChangeC(rho1, k0, 0.0, 0.0, 0.0);

    REQUIRE(chfunctions::logAffine(
        rungekutta::computeFunctional(T, 2048, std::vector<double >({0, 0}),
            [&](double t, const std::vector<double>& x){
                return std::vector<double>({
                    beta(x[0], -1.0), //-1.0 doesnt matter because l is 0
                    alpha(x[0], -1.0)//-1.0 doesnt matter because l is 0
                });
            }
        ),
        r0)==Approx(chfunctions::cirLogMGF(rho1, k0, -k1, sqrt(H1), T, r0)).epsilon(.0001));

}

/*
TEST_CASE("Test TimeChangeC", "[CF]"){
    REQUIRE(chfunctions::TimeChangeC(1.0, .03, 0.0, 0.0, 0.0)(.2, -1.0)==Approx(.2*.03));
}
TEST_CASE("Test AlphaOrBeta", "[CF]"){
    REQUIRE(chfunctions::AlphaOrBeta(0.0, .03, 0.0, 0.0)(.2, -1.0)==Approx(.2*.03));
}
TEST_CASE("Test Beta vs TimeChangeB", "[CF]"){
    std::complex<double> u(1.0, 1.0);
    double lambda=0.0;
    double muJ=0.0;
    double sigJ=0.0;
    double sigma=sqrt(.0398);
    double v0=.0175/.0398;
    double speed=1.5768;
    double adaV=.5751/sigma;
    double rho=-.5711;
    double r=0;
    double T=1.0;
    double q=5;
    double delta=0.0;
    auto beta1=chfunctions::AlphaOrBeta(
        -chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sigma), 
        -(speed+(delta*lambda)/q-u*rho*sigma*adaV),
        -adaV*adaV, //how odd
        -lambda //this feels bizarre
    );

    auto beta2=[&](const auto& val, const auto& cfPart){
        return chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sigma)-(speed+(delta*lambda)/q-u*rho*sigma*adaV)*val-val*val*adaV*adaV*.5-lambda*cfPart;
    };
    auto result1=beta1(.2, -1.0);
    auto result2=beta2(.2, -1.0);

    REQUIRE(result1.real()==Approx(result2.real()));
    REQUIRE(result1.imag()==Approx(result2.imag()));
}*/

TEST_CASE("Test CIR against runge kutta average 1", "[CF]"){
    auto rho1=1.0;
    auto k1=.3; //long run avearge of 1
    auto H1=.09;
    auto T=.5;
    auto r0=.5;
   // auto bPrime=chfunctions::TimeChangeB(rho1, k1, H1, 1.0, 0.0);
    //auto cPrime=chfunctions::TimeChangeC(rho1, k1, 0.0, 0.0, 0.0);

    auto beta=chfunctions::AlphaOrBeta(rho1, -k1, H1, 0.0);
    auto alpha=chfunctions::AlphaOrBeta(0.0, k1, 0.0, 0.0);
    
    REQUIRE(chfunctions::logAffine(
        rungekutta::computeFunctional_move(T, 32, std::vector<double >({0, 0}),
            [&](double t, const std::vector<double>& x){
                return std::vector<double>({
                    beta(x[0], -1.0), //-1.0 doesnt matter because l is 0
                    alpha(x[0], -1.0)//-1.0 doesnt matter because l is 0
                });
                // chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sigma)-(speed+(delta*lambda)/q-u*rho*sigma*adaV)*val-val*val*adaV*adaV*.5-lambda*cfPart;
                //-rho+K*currVal+.5*currVal*currVal*H+l*cfPart;
            }
        ),
        r0)==Approx(chfunctions::cirLogMGF(rho1, k1, k1, sqrt(H1), T, r0)).epsilon(.0001));

}

/*TEST_CASE("Test CGMY", "[CF]"){

}*/