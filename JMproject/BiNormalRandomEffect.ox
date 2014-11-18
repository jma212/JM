#include <BiNormalRandomEffect.h>
/** Create RandomEffect, discretized bivariate normal distribution
@N is number of points in each dimension
@mu 2*1 vector
@sigma 2*2 variance covariance matrix
**/

BiNormalRandomEffect::BiNormalRandomEffect(L,N, mu, sigma){
   decl M=N^2, wights=<>,i;
   RandomEffect(L,M);
   this.order=N;
   this.mu=mu;
   this.sigma=sigma;
   this.chol=choleski(sigma);
   nodes=<>;
   GaussianHermite::Initialize(order);
   for (i=0;i<order; ++i){
       nodes|=ones(order,1)*GaussianHermite::nodes[i]~GaussianHermite::nodes';
       wights|=ones(order,1)*GaussianHermite::wght[i]~GaussianHermite::wght';
       }
   nodes=mu+chol*nodes';
   pdf[]=wights[][0].*wights[][1]/(M_PI);
   }

/**
BiNormalRandomEffect::Distribution{

}
**/
	
