#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fenv.h>
#include <algorithm>
#include "bs.h"

namespace functor {

class f5 {
protected:
   double m_s, m_vol, m_rd, m_rf, m_t, m_T;
public:
   virtual double operator() (double s, double vol,
                  double rd, double rf, double t,
                  bs::types::Greeks type=bs::types::Value) const = 0;
   void get_central(double& s,  double& vol,
                     double& rd, double& rf, double& t) const {
      s=m_s;  vol=m_vol; rd=m_rd; rf=m_rf; t=m_t;
   }
};

class barrier: public f5 {
protected:
   bs::types::PutCall m_pc;
   bs::types::BarrierActive m_bcont;
   bs::types::BarrierKIO m_kio;
   double m_B1, m_B2;
   double m_K;
   double m_rebate;
public:
   barrier(double s, double vol, double rd, double rf, double tau,
                  double K, double B1, double B2, double rebate,
                  bs::types::PutCall pc, bs::types::BarrierKIO kio,
                  bs::types::BarrierActive bcont) {
      m_s=s; m_vol=vol; m_rd=rd; m_rf=rf; m_t=0.0; m_T=tau;
      m_K=K; m_B1=B1; m_B2=B2; m_rebate=rebate;
      m_pc=pc; m_kio=kio; m_bcont=bcont;
   }
               
   double operator()(double s, double vol,
                  double rd, double rf, double t,
                  bs::types::Greeks greek) const {
      return bs::barrier(s,vol,rd,rf,m_T-t,
                        m_K,m_B1,m_B2,m_rebate,m_pc,m_kio,m_bcont,greek);
   }
};

class touch: public f5 {
protected:
   bs::types::BarrierActive m_bcont;
   bs::types::BarrierKIO m_kio;
   bs::types::ForDom m_fd;
   double m_B1, m_B2;
public:
   touch(double s, double vol, double rd, double rf, double tau,
                  double B1, double B2,
                  bs::types::ForDom fd, bs::types::BarrierKIO kio,
                  bs::types::BarrierActive bcont) {
      m_s=s; m_vol=vol; m_rd=rd; m_rf=rf; m_t=0.0; m_T=tau;
      m_B1=B1; m_B2=B2;
      m_fd=fd; m_kio=kio; m_bcont=bcont;
   }
               
   double operator()(double s, double vol,
                  double rd, double rf, double t,
                  bs::types::Greeks greek) const {
      return bs::touch(s,vol,rd,rf,m_T-t,m_B1,m_B2,m_fd,m_kio,m_bcont,greek);
   }
};

class putcall: public f5 {
protected:
   bs::types::PutCall m_pc;
   double m_K;
public:
   putcall(double s, double vol, double rd, double rf, double tau,
            double K, bs::types::PutCall pc) {
      m_s=s; m_vol=vol; m_rd=rd; m_rf=rf; m_t=0.0; m_T=tau;
      m_K=K; m_pc=pc;
   }
               
   double operator()(double s, double vol,
                  double rd, double rf, double t,
                  bs::types::Greeks greek) const {
      return bs::putcall(s,vol,rd,rf,m_T-t,m_K,m_pc,greek);
   }
};



} // namespace



// returns the number of identical digits
int digits(double x, double y){
   double d = -log10(fabs((x-y)/x));
   return std::max( 0, (int) floor(d) );
}


void test_greeks(const functor::f5& f){

   double s, vol, rd, rf, t;
   f.get_central(s,vol,rd,rf,t);

   double ds=s*1e-6;
   double ds2=s*1e-4;
   double dt=1e-6;
   double dvol=1e-6;
   double dvol2=1e-4;
   double dr=1e-7;

   // analytic calculation of price sensitivities (greeks)
   double v = f(s, vol, rd, rf, t);
   double delta = f(s, vol, rd, rf, t, bs::types::Delta);
   double gamma = f(s, vol, rd, rf, t, bs::types::Gamma);
   double vega = f(s, vol, rd, rf, t, bs::types::Vega);
   double volga = f(s, vol, rd, rf, t, bs::types::Volga);
   double vanna = f(s, vol, rd, rf, t, bs::types::Vanna);
   double rho = f(s, vol, rd, rf, t, bs::types::Rho_d);
   double rhof = f(s, vol, rd, rf, t, bs::types::Rho_f);
   double theta = f(s, vol, rd, rf, t, bs::types::Theta);

   // calculation via finite differences
   double v1, v2, v11, v12, v21, v22;

   // value

   // delta
   v1 = f(s-ds, vol, rd, rf, t);
   v2 = f(s+ds, vol, rd, rf, t);
   double delta_1 = (v2-v1)/(2.0*ds);

   // gamma
   v1 = f(s-ds2, vol, rd, rf, t);
   v2 = f(s+ds2, vol, rd, rf, t);
   double gamma_1 = (v1-2.0*v+v2)/(ds2*ds2);
   // gamma (calculated via delta)
   v1 = f(s-ds, vol, rd, rf, t, bs::types::Delta);
   v2 = f(s+ds, vol, rd, rf, t, bs::types::Delta);
   double gamma_2 = (v2-v1)/(2.0*ds);

   // vega
   v1 = f(s, vol-dvol, rd, rf, t);
   v2 = f(s, vol+dvol, rd, rf, t);
   double vega_1 = (v2-v1)/(2.0*dvol);

   // volga (vol gamma)
   v1 = f(s, vol-dvol2, rd, rf, t);
   v2 = f(s, vol+dvol2, rd, rf, t);
   double volga_1 = (v1-2.0*v+v2)/(dvol2*dvol2);
   // volga (calculated via vega)
   v1 = f(s, vol-dvol, rd, rf, t, bs::types::Vega);
   v2 = f(s, vol+dvol, rd, rf, t, bs::types::Vega);
   double volga_2 = (v2-v1)/(2.0*dvol);

   // vanna
   v11 = f(s-ds2, vol-dvol2, rd, rf, t);
   v12 = f(s-ds2, vol+dvol2, rd, rf, t);
   v21 = f(s+ds2, vol-dvol2, rd, rf, t);
   v22 = f(s+ds2, vol+dvol2, rd, rf, t);
   double vanna_1 = (v22-v12-v21+v11) / (4.0*dvol2*ds2);
   // vanna (calculated via delta)
   v1 = f(s, vol-dvol, rd, rf, t, bs::types::Delta);
   v2 = f(s, vol+dvol, rd, rf, t, bs::types::Delta);
   double vanna_2 = (v2-v1)/(2.0*dvol);
   // vanna (calculated via vega)
   v1 = f(s-ds, vol, rd, rf, t, bs::types::Vega);
   v2 = f(s+ds, vol, rd, rf, t, bs::types::Vega);
   double vanna_3 = (v2-v1)/(2.0*ds);

   // rho
   v1 = f(s, vol, rd-dr, rf, t);
   v2 = f(s, vol, rd+dr, rf, t);
   double rho_1 = (v2-v1)/(2.0*dr);

   // rhof
   v1 = f(s, vol, rd, rf-dr, t);
   v2 = f(s, vol, rd, rf+dr, t);
   double rhof_1 = (v2-v1)/(2.0*dr);

   // theta
   v1 = f(s, vol, rd, rf, t-dt);
   v2 = f(s, vol, rd, rf, t+dt);
   double theta_1 = (v2-v1)/(2.0*dt);

   printf("        analytic       finite diff        ");
   printf("finite diff using analytic first order\n");
   printf("val   %14.10f\n",v);
   printf("delta %14.10f %14.10f [%i]\n",delta,
            delta_1, digits(delta, delta_1) );
   printf("gamma %14.10f %14.10f [%i] %14.10f [%i]\n", gamma,
            gamma_1, digits(gamma,gamma_1),
            gamma_2, digits(gamma,gamma_2) );
   printf("vega  %14.10f %14.10f [%i]\n", vega,
            vega_1, digits(vega,vega_1) );
   printf("volga %14.10f %14.10f [%i] %14.10f [%i]\n", volga,
            volga_1, digits(volga, volga_1),
            volga_2, digits(volga, volga_2) );
   printf("vanna %14.10f %14.10f [%i] %14.10f [%i] %14.10f [%i]\n", vanna,
            vanna_1, digits(vanna, vanna_1),
            vanna_2, digits(vanna, vanna_2),
            vanna_3, digits(vanna, vanna_3) );
   printf("rhod  %14.10f %14.10f [%i]\n", rho,
            rho_1, digits(rho, rho_1) );
   printf("rhof  %14.10f %14.10f [%i]\n", rhof,
            rhof_1, digits(rhof, rhof_1) );
   printf("theta %14.10f %14.10f [%i]\n", theta,
            theta_1, digits(theta, theta_1) );

}


// ======================================================================
//         Name:  main
//  Description:  main function, executed at runtime
// ======================================================================
int main(int argc, char** argv) {
   if(argc<1) {
      printf("usage: %s <>\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);


   // barrier option specification
   bs::types::PutCall pc = bs::types::Call;
   bs::types::BarrierActive bcont = bs::types::Continuous;
   bs::types::BarrierKIO kio = bs::types::KnockOut;
   double tau=1.0;
   double B1=0.95;
   double B2=1.4;
   double K=1.0;
   double rebate=0.01;
   double s=1.0;
   double vol=0.15;
   double rd=log(1.01);
   double rf=log(1.02);

   // additional for touch option
   bs::types::ForDom  fd=bs::types::Domestic;

   // declare test functions
   functor::putcall  f1(s,vol,rd,rf,tau,K,pc);
   functor::barrier  f2(s,vol,rd,rf,tau,K,B1,B2,rebate,pc,kio,bcont);
   functor::touch    f3(s,vol,rd,rf,tau,B1,B2,fd,kio,bcont);

   // run the test
   printf("\n");
   printf("comparison of analytical derivatives vs finite differences\n");
   printf("----------------------------------------------------------\n");

   printf("market: S=%f, vol=%f, rd=%f, rf=%f\n\n",s,vol,rd,rf);

   printf("%s option:\n   K=%f, T=%f\n\n",
         pc==bs::types::Call ? "call" : "put", K, tau);
   test_greeks(f1);

   printf("\nknock-%s barrier option:\n",
         kio==bs::types::KnockIn ? "in" : "out");
   printf("   barriers={%f,%f}, rebate=%f\n",  B1, B2, rebate);
   printf("   barriers observed %s\n\n",
         bcont==bs::types::Continuous ? "continuously" : "at maturity only");
   test_greeks(f2);

   printf("\ntouch option:\n   barriers={%f,%f}, payment in %s\n\n",
         B1, B2, fd==bs::types::Domestic ? "domestic" : "foreign");
   test_greeks(f3);

   printf("\nNote, values in [] indicate the number of matching digits\n\n");

   return EXIT_SUCCESS;
}


