#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cassert>
#include <time.h>
#include <sys/time.h>
#include <fenv.h>
#include "bs.h"


double stoptime(void) {
   struct timeval	t;
   gettimeofday(&t,NULL);
   return (double) t.tv_sec + t.tv_usec/1000000.0;
}

double pnorm(double x) {
   return 0.5 * erfc(-M_SQRT1_2*x);
}


// wrap market, contracts and pricing into classes
// -----------------------------------------------
class MarketData {
public:
   double S, vol, rd, rf;
};

class Option {
public:
   MarketData m;
   double T;
   virtual double price(bs::types::Greeks) = 0;
};

class OptionDummy : public Option {
public:
   double price(bs::types::Greeks greek) {
      return 0.0;
   }
};
class OptionAdd : public Option {
public:
   double price(bs::types::Greeks greek) {
      return T+T;
   }
};
class OptionNormCdf : public Option {
public:
   double price(bs::types::Greeks greek) {
      return pnorm(T);
   }
};


class OptionBinCash : public Option {
public:
   double B;
   bs::types::PutCall pc;
   double price(bs::types::Greeks greek) {
      return bs::bincash(m.S,m.vol,m.rd,m.rf,T,B,pc,greek);
   }
};
class OptionBinary: public Option {
public:
   double B1, B2;
   bs::types::ForDom fd;
   double price(bs::types::Greeks greek) {
      return bs::binary(m.S,m.vol,m.rd,m.rf,T,B1,B2,fd,greek);
   }
};
class OptionTouch: public Option {
public:
   double B1, B2;
   bs::types::ForDom fd;
   bs::types::BarrierKIO kio;
   bs::types::BarrierActive bcont;
   double price(bs::types::Greeks greek) {
      return bs::touch(m.S,m.vol,m.rd,m.rf,T,B1,B2,fd,kio,bcont,greek);
   }
};
class OptionPutCall: public Option {
public:
   double K;
   bs::types::PutCall pc;
   double price(bs::types::Greeks greek) {
      return bs::putcall(m.S,m.vol,m.rd,m.rf,T,K,pc,greek);
   }
};
class OptionBarrier: public Option {
public:
   double K,B1,B2,rebate;
   bs::types::PutCall pc;
   bs::types::BarrierKIO kio;
   bs::types::BarrierActive bcont;
   double price(bs::types::Greeks greek) {
      return bs::barrier(m.S,m.vol,m.rd,m.rf,T,K,B1,B2,rebate,pc,kio,bcont,greek);
   }
};



void speed_test(const char* message, Option& o, int iter, double mhz) {
   const int maxgreeks=9;
   double T=o.T;
   //printf("%s: V=%.3f, ", message, o.price(bs::types::Value));
   printf("%s", message);
   double x=0.0;

   // benchmark all price sensitivities individually
   for(int greek=0; greek<maxgreeks; greek++) {
      x=0.0;
      double t=stoptime();
      for(int i=0; i<iter; i++) {
         x+=o.price((bs::types::Greeks)greek);
         o.T+=1e-10;       // so gcc doesn't optimise loop away
      }
      t=stoptime()-t;
      assert(x!=-1.23456); // so gcc doesn't optimise loop away
      printf("%6i ", (int) round(t*mhz*1e6/iter));
   }

   // benchmark computation of all price sensitivities at once
   /*
   x=0.0;
   double t=stoptime();
   for(int i=0; i<iter; i++){
      for(int greek=0; greek<maxgreeks; greek++) {
         x+=o.price((bs::types::Greeks)greek);
      }
      o.T+=1e-10;       // so gcc doesn't optimise loop away
   }
   t=stoptime()-t;
   assert(x!=-1.23456); // so gcc doesn't optimise loop away
   printf("%6i", (int) round(t*mhz*1e6/iter));
   */
   printf("\n");
   o.T=T;
}



// ======================================================================
//         Name:  main
//  Description:  main function, executed at runtime
// ======================================================================
int main(int argc, char** argv) {
   // we don't want performance degradation due to nan's and inf's
   // so if anything occurs we exit
   feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);

   if(argc!=3) {
      printf("usage: %s <cpu MHz> <thousand ops>\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   double mhz=atof(argv[1]);
   int      n=atoi(argv[2])*1000;
   if(n<=0) n=1000;

   MarketData m;
   m.S=1.2;
   m.rd=log(1.01);
   m.rf=log(1.02);
   m.vol=0.13;

   double T=1.5;

   OptionDummy    o_dummy;
   OptionAdd      o_add;
   OptionNormCdf  o_pnorm;
   OptionBinCash  o_bincash;
   OptionBinary   o_binary;
   OptionTouch    o_touch;
   OptionPutCall  o_putcall;
   OptionBarrier  o_barrier;


   o_add.T=T;
   o_pnorm.T=T;

   o_bincash.m=m;
   o_bincash.T=T;
   o_bincash.B=1.0;
   o_bincash.pc=bs::types::Call;

   o_binary.m=m;
   o_binary.T=T;
   o_binary.B1=1.0;
   o_binary.B2=0.0;
   o_binary.fd=bs::types::Domestic;

   o_touch.m=m;
   o_touch.T=T;
   o_touch.B1=1.0;
   o_touch.B2=0.0;
   o_touch.fd=bs::types::Domestic;
   o_touch.kio=bs::types::KnockOut;
   o_touch.bcont=bs::types::Continuous;

   o_putcall.m=m;
   o_putcall.T=T;
   o_putcall.K=1.2;
   o_putcall.pc=bs::types::Call;
   o_barrier.bcont=bs::types::Maturity;


   o_barrier.m=m;
   o_barrier.T=T;
   o_barrier.B1=0.0;
   o_barrier.B2=1.4;
   o_barrier.K=1.2;
   o_barrier.rebate=0.0;
   o_barrier.pc=bs::types::Call;
   o_barrier.kio=bs::types::KnockOut;

   printf("\nnumber of cpu-cycles per function:\n\n");
   printf("            price  delta  gamma  theta   vega  volga  vanna");
   printf("    rho  rho_f\n");
   speed_test("dummy      ", o_dummy, n, mhz);
   speed_test("add        ", o_add, n, mhz);
   speed_test("norm cdf   ", o_pnorm, n, mhz);

   printf("\n");

   speed_test("bincash()  ", o_bincash, n, mhz);
   speed_test("binary()   ", o_binary, n, mhz);
   speed_test("1 touch()  ", o_touch, n/4, mhz);
   o_touch.B2=1.4;
   speed_test("2 touch()  ", o_touch, n/40, mhz);
   printf("\n");
   speed_test("putcall()  ", o_putcall, n, mhz);
   speed_test("T barrier()", o_barrier, n/2, mhz);  // barrier at T only
   o_barrier.bcont=bs::types::Continuous;
   speed_test("1 barrier()", o_barrier, n/8, mhz);
   o_barrier.B1=1.0;
   speed_test("2 barrier()", o_barrier, n/50, mhz);
   printf("\n");


   return EXIT_SUCCESS;
}

