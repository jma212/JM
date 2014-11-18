# import "DDP"
# include "GaussianHermite.ox"

/** Create RandomEffect, discretized bivariate normal distribution**/
struct BiNormalRandomEffect:RandomEffect{
   decl mu, sigma, chol, order;
   decl nodes;
   BiNormalRandomEffect(L, N, mu, sigma);
   }