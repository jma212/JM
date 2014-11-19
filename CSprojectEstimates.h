/**estimate of CSproject **/
# import "niqlow"            //CF:  just import all of niqlow if doing DP and optimization.
# include "RActionCounter.ox"
# include "BiNormalRandomEffect.ox"
#include "Wage.ox"

/** put estimation and solving the DDP in to a catchall class **/

struct CSprojectEstimates {
	static decl
    /** Value function iteration method. **/                     EMax,
    /** Solving system of equations to find efficency wages**/   Roots,
    /** Simulated Method of Moments**/                           SMM;
	static DoAll();
	}

/** The CS project model with estimated parameters **/
struct CSprojectHat : ExPostSmoothing	{
        static 	decl omigahat, alphahat, betahat, muhat, sigmahat, deltahat, rhohat, tfphat, sharehat, elashat;
        enum {cs, ncs, Msector};
        enum {A0=22,
              A1=40,
	          Ntypes=2};
		enum {omiga, alpha, beta, mu, sigma, delta, rho, tfp, share, elas, Nparams};   // tags for paramters to estimate
		static const  decl
		        pars={<2.0;2.9>, <1.01, 0.35, 0.4, 0.97>, <0.1;0.2>, <0;0>, <3.0,10.0,-0.5>, 0.9, 1/20, <1.1;0.98>, <0.2;0.3>, <0.3;0.7> };
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
		
		static 	decl hat;// hat store the paramter values
		static  Run();          
		static 	FirstStage();	  
        static  SecondStage();
		static 	Reachable();
		static  ldemand(const vP, const adfunc, const avScore, const amHess);
				Utility();
                Auxi();
		}

/**the FiveO objective to solve  for market equilibirum price**/

struct MarketClear:BlackBox{
       MarketClear(L,...);
	   vfunc();
       }
