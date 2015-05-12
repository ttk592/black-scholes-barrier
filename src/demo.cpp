#include <cstdio>
#include <cstdlib>
#include "bs.h"

int main(int argc, char** argv) {
   // define market data
   double S=1.0, vol=0.15, r=0.03, rf=0.01;
   // define call option
   double T=1.0, K=1.0;
   bs::types::PutCall pc = bs::types::Call;
   // define barrier option (up-and-out)
   double B1=0.0, B2=1.2, rebate=0.0;
   bs::types::BarrierKIO kio = bs::types::KnockOut;
   bs::types::BarrierActive bcont = bs::types::Continuous;

   double price, delta, prob, bhit;

   // calculate call option price
   price=bs::putcall(S,vol,r,rf,T,K,pc);
   delta=bs::putcall(S,vol,r,rf,T,K,pc,bs::types::Delta);
   prob=bs::prob_in_money(S,vol,r-rf,T,K,0.0,0.0,pc);
   printf("call option: price = %f, delta = %f\n", price, delta);
   printf("   probability of S_T in the money    = %.2f%%\n\n", prob*100.0);

   // calculate barrier option price
   price=bs::barrier(S,vol,r,rf,T,K,B1,B2,rebate,pc,kio,bcont);
   delta=bs::barrier(S,vol,r,rf,T,K,B1,B2,rebate,pc,kio,bcont,bs::types::Delta);
   prob=bs::prob_in_money(S,vol,r-rf,T,K,B1,B2,pc);
   bhit=bs::prob_hit(S,vol,r-rf,T,B1,B2);
   printf("barrier option: price = %f, delta = %f\n", price, delta);
   printf("   probability of S_T in the money    = %.2f%%\n", prob*100.0);
   printf("   probability of hitting the barrier = %.2f%%\n\n", bhit*100.0);

   return EXIT_SUCCESS;
}
