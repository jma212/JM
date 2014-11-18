/**estimate of CSproject **/
# include "Clock.h"
# include "RActionCounter.ox"
# include "BiNormalRandomEffect.ox"
# import "FiveO"

/**put (estimation with simulated method of moments) and (solving the DDP & market equilibirum) into a catchall class**/
struct CSprojectEstimates {
	   static decl
	   /**smm object**/                                 ssmobj,
	   /**method of solving smm**/						ssmMethod;
	   static DoAll();
       }

/**The CS project model with estimated parameters**/
struct CSprojectHat:ExPostSmoothing {
	   static decl hat1, hat2;
	   enum {cs, ncs, Msector};
	   enum {A0=22,
	         A1=40,
			 Ntypes=2};
	   enum {alpha, beta, mu, sigma, delta, rho, tfp, share, elas, Nparams};	  // tages for parameters to estimate omiga put separatedly
	   static const decl
	          DGP1=<2.0;2.9>,
			  DGP2={<1.01; 0.35; 0.4; 0.97>, <0.1;0.2>, <0;0>, <3.0,-2.5;-2.5,10>, 0.9, 1/20, <1.1;0.98>, <0.2;0.3>, <0.3;0.7> };
	   static decl
	     /**choices of occupation**/            occupations,
         /**experience**/                       xper,
         /**random effect **/			        types,
         /**auxiliary wage**/                   AuxiWage,
         /**auxiliary efficient unit**/	        realizedutility,
         /**solving method**/                   Method,                                  
	     /**random effect nodes **/			    nodes,
		 /**aggregate labor supply **/          laborsupply,
		 /**aggregate labor demand **/          labordemand;
		static const decl  
	   	 /**immigrant shock **/				    M=<11463,48794>;
        static  SetUp();
		static  Run();
		static  SimulateMoments();
		static 	FirstStage();	  
		static 	Reachable();
		static  ldemand(const vP, const adfunc, const avScore, const amHess);
				Utility();
                Auxi();
		}		
/**the FiveO objective to solve  for market equilibirum price**/

struct MarketClear:BlackBox{	  // this is the black box that compute the inner market clear wage
       MarketClear(L, hat1, hat2);
	   vfunc();
       }

struct SMM:BlackBox{              // this is the black box that compute the outter ssm
       static decl w, datamoments;          // data holds the data moments
	   SMM(L,data, hat1, hat2);
	   vfunc();
       }
	   