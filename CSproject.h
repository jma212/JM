//# include "NAging.ox"
#import "niqlow"
# include "RActionCounter.ox"
# include "BiNormalRandomEffect.ox"
#include "Wage.ox"

struct CSproject:ExPostSmoothing{
       enum {cs, ncs, Msector};
       enum {A0=22,
             A1=65,
	         Ntypes=2};
       static const decl
       /**log efficiency wage  **/			  omiga=<2.0;2.9>,
       /**wage parameters**/                  alpha={<1.01; 0.8>,
                                                     <0.73; 0.97>},
       /**experience profile**/				  beta=<0.01;0.02>,						  
       /**dynamic roy **/					  mu=<0;0>,
                                              sigma=<3.0,-1;-1,10.0>,
       /**discount rate**/					  delta=0.9,
       /**ExPostSmoothing **/                 rho=1/4;
                                     
									  									 
      static decl
       /**choices of occupation**/            occupations,
       /**experience**/                       xper,
       /**random effect **/			          types,
       /**auxiliary wage**/                   AuxiWage,
       /**auxiliary efficient unit**/	      realizedutility,
       /**solving method**/                   Method,                                  
	   /**random effect nodes **/			  nodes;
	   
             static Run();
             static Reachable();
             static fPi0(FeasA);
             static fPi1(FeasA);
                    Auxi();
                    Utility();
	   
}



	  