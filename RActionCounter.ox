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
   F=setbounds(v~v+1, 0 , N);
   inda=sumr(FeasA[][Target.pos].==ToTrack);
   P=(1-inda*exp(-Beta*v))~(inda*exp(-Beta*v));
   return {F,P};  
   }
