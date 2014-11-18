# include<GaussianHermite.h>
/**Create discretized univariate normal GaussianHermite quadrature
@order, order of quadrature**/

GaussianHermite::coef(){
   decl i;
   polycoef=new array[max(this.n+1,3)];
   polycoef[0]=<1>;
   polycoef[1]=<2,0>;
   polycoef[2]=<4,0,-2>;
   if (this.n>2) {
      for (i=3; i<this.n+1;++i) polycoef[i]=(2*polycoef[i-1]~0)+(0~0~-2*(i-1)*polycoef[i-2]);	  // compute the coefficients recursively
      }
   return polycoef;
   }

GaussianHermite::Initialize(order){
   this.n=order;
   if (this.n<1) oxrunerror("GaussianHermite("+sprint("%d",this.n)+") not supported");
   coef();
   decl i, pcoef=reverser(polycoef[this.n-1]), pr, num=2^(this.n-1)*factorial(this.n)*sqrt(M_PI);
   polyroots(polycoef[this.n], &pr);
   nodes=sortr(pr[0][]);
   wght=zeros(1, this.n);
   for (i=0; i<this.n;++i) wght[i]=num/sqr(n*polyeval(pcoef, nodes[i]));
   return TRUE;
   }



