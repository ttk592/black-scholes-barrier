/*
 * bs.h
 *
 * vanilla and barrier option pricing within the Black-Scholes model,
 * i.e. TV or theoretical values of options,
 * also calculates first order and a few second order price
 * sensitivities analytically
 *
 * Model: dS/S = (rd-rf) dt + sigma dW
 *
 * - S     ... underlying spot price
 * - rd    ... (domestic) interest rate
 * - rf    ... FX: foreign interest rate
 *             Equity: continuous dividend yield / borrow cost
 *             Commodity: convenience yield
 * - sigma ... volatility
 * - W     ... Brownian motion
 *
 * This library does not support discrete dividends.
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2012, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */

#ifndef _bs_h
#define _bs_h

namespace bs {

namespace types {

// price sensitivities
enum Greeks {
   Value    = 0,        //
   Delta    = 1,        // d/dS
   Gamma    = 2,        // d^2/dS^2
   Theta    = 3,        // d/dt
   Vega     = 4,        // d/dsigma
   Volga    = 5,        // d^2/dsigma^2
   Vanna    = 6,        // d^2/dsigma dS
   Rho_d    = 7,        // d/dr_d
   Rho_f    = 8         // d/dr_f
};

// put or call option
enum PutCall {
   Call     = 1,
   Put      = -1
};

// regular/reverse barrier
enum KOType {
   Regular = 0,
   Reverse = 1
};

// knock-in or knock-out
enum BarrierKIO {
   KnockIn = -1,
   KnockOut = 1
};

// barrier observed continuously or just at maturity (truncated payoff)
enum BarrierActive {
   Continuous = 0,
   Maturity   = 1
};


enum ForDom {
   Domestic = 0,     // cash or nothing
   Foreign  = 1      // asset or nothing
};

} // namespace types




// ---------------------------------------------------------------------
// main wrapper functions for barrier option pricing: barrier(), touch()
//
// common features:
// - flag types::BarrierActive, determines if barriers are monitored
//   continuously or only at maturity (ie truncated payoff)
// - flag types::BarrierKIO, determines if the barrier is of knock-in
//   or knock-out type
// - B1 is the lower barrier, B2 is the upper barrier (B1<B2)
// - if barriers B1 or B2 are set to 0 or negative, they will be
//   ignored, ie these wrapper functions can price double-barrier,
//   single-barrier, or options with no barriers at all
//
// the differences only lie in the payoff at maturity:
// - barrier() has a put/call payoff (typically called a barrier option)
// - touch() has a cash/asset-or-nothing payoff (called a touch option)
// ---------------------------------------------------------------------

// barrier option (put/call payoff profile)
// rebate is assumed to be paid immediately when a barrier is hit
double barrier(double S, double vol, double rd, double rf,
               double tau, double K, double B1, double B2,
               double rebate,
               types::PutCall pc, types::BarrierKIO kio,
               types::BarrierActive bcont,
               types::Greeks greek=types::Value);

// touch/no-touch options (cash/asset or nothing payoff profile)
double touch(double S, double vol, double rd, double rf,
             double tau, double B1, double B2, types::ForDom fd,
             types::BarrierKIO kio, types::BarrierActive bcont,
             types::Greeks greek=types::Value);

// probability of hitting a barrier
// this is almost the same as the price of a touch option (domestic)
// as it pays one if a barrier is hit; we only have to offset the
// discounting and we get the probability
double prob_hit(double S, double vol, double mu,
                double tau, double B1, double B2);

// probability of being in-the-money, ie payoff is greater zero,
// assuming payoff(S_T) > 0 iff S_T in [B1, B2]
// this the same as the price of a cash or nothing option
// with no discounting
double prob_in_money(double S, double vol, double mu,
                     double tau, double K, double B1, double B2,
                     types::PutCall pc);



// ---------------------------------------------------------------------
// direct functions for options without barriers, ie vanilla options or
// options with truncated payoffs (barrier observed at maturity)
//
// although all functionality can be accessed via the wrapper functions
// above, it may be convenient to call the pricing functions of
// vanilla options directly:
//
//  - bincash(): cash or nothing option
//  - binasset(): asset or noting option
//  - binary(): wrapper for both cash/asset or nothing option
//  - putcall(): vanilla put or call option
//  - putcalltrunc(): put/call option with truncated payoff
// ---------------------------------------------------------------------

// binary option: cash or nothing (domestic)
//   call - pays 1 if S_T is above strike K
//   put  - pays 1 if S_T is below strike K
double bincash(double S, double vol, double rd, double rf,
               double tau, double K,
               types::PutCall pc, types::Greeks greeks=types::Value);

// binary option: asset or nothing (foreign)
//   call - pays S_T if S_T is above strike K
//   put  - pays S_T if S_T is below strike K
double binasset(double S, double vol, double rd, double rf,
                double tau, double K,
                types::PutCall pc, types::Greeks greeks=types::Value);

// just for convenience we can combine bincash and binasset into
// one function binary
// using bincash()  if fd==types::Domestic
// using binasset() if fd==types::Foreign
double binary(double S, double vol, double rd, double rf,
              double tau, double K,
              types::PutCall pc, types::ForDom fd,
              types::Greeks greek=types::Value);

// further wrapper to combine single/double barrier binary options
// into one function
// B1<=0 - it is assumed lower barrier not set
// B2<=0 - it is assumed upper barrier not set
double binary(double S, double vol, double rd, double rf,
              double tau, double B1, double B2,
              types::ForDom fd, types::Greeks greek=types::Value);

// vanilla put/call option
//   call pays (S_T-K)^+
//   put  pays (K-S_T)^+
// this is the same as: +/- (binasset - K*bincash)
double putcall(double S, double vol, double rd, double rf,
               double tau, double K,
               types::PutCall putcall, types::Greeks greeks=types::Value);

// truncated put/call option, single barrier
// need to specify whether it's down-and-out or up-and-out
// regular (keeps monotonicity): down-and-out for call, up-and-out for put
// reverse (destroys monoton):   up-and-out for call, down-and-out for put
//   call pays (S_T-K)^+
//   put  pays (K-S_T)^+
double putcalltrunc(double S, double vol, double rd, double rf,
                    double tau, double K, double B,
                    types::PutCall pc, types::KOType kotype,
                    types::Greeks greeks=types::Value);

// wrapper function for put/call option which combines
// double/single/no truncation barrier
// B1<=0 - assume no lower barrier
// B2<=0 - assume no upper barrier
double putcalltrunc(double S, double vol, double rd, double rf,
                    double tau, double K, double B1, double B2,
                    types::PutCall pc, types::Greeks greek=types::Value);


} // namespace bs



#endif   // ifdef _bs_h
