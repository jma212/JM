#include <RActionCounter.h>

/**Create a random action counter, the probability of increase is a decreasing function of v, (exp(1-Beta*v))
@L label
@Act is the Target
@Beta, curvature parameters
**/

RActionCounter::RActionCounter(L, N, Act, ToTrack, Beta){
   this.Beta=Beta;
   this.ToTrack=ToTrack;
   this.Target=Act;
   this.N=N;
   StateVariable(L,N);
   }

RActionCounter::Transit(const FeasA){
   decl F, P, inda;
   if (v==N-1) return {matrix(v), ones(rows(FeasA),1)};
   inda=(FeasA[][Target.pos].==ToTrack);
   P=exp(-Beta*v)*inda;
   return{v~(v+1), (1-P)~P};  
   }

