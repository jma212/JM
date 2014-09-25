# include <oxstd.h>
#include <oxfloat.oxh>
/**Create discretized univariate normal GaussianHermite quadrature**/

struct GaussianHermite{
   static decl nodes, wght, n, polycoef;  
   static Initialize(order);
   static coef();
   }