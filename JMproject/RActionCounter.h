# import "DDP"
/**Create a random action counter, the probability of increase is a decreasing function of v, (exp(1-Beta*v))**/

struct RActionCounter:StateVariable{
   decl Beta, Target, ToTrack, N;  // the paramter that capture the curvature of time variant probability
   RActionCounter(L, N, Act, ToTrack, Beta);
   Transit(const FeasA);
   }