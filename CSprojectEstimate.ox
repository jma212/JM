#include "CSprojectEstimates.h"
#import <maximize>

/**find the market clear wage**/
CSprojectHat::FirstStage(){
    decl i;
	omigahat=new array[Msector];
	omigahat[0]=new Positive("logwagecs", pars[omiga][0]);
	omigahat[1]=new Positive("logwagencs", pars[omiga][1]);
	alphahat=new array[4];
	for (i=0;i<4;++i) alphahat[i]=new Free("ParaWage"+sprint(i), pars[alpha][i]);
	betahat= new array[Msector];
	betahat[0]=new Positive("bcs", pars[beta][0]);
	betahat[1]=new Positive("bncs", pars[beta][1]);
//	muhat=new array[Msector];
//	muhat[0]=new Determined("mucs",pars[mu][0]);
//	muhat[1]=new Determined("muncs", pars[mu][1]);
	sigmahat=new array[3];	 // sigma should be code as a block of parameters that satisfied the positive definite 
	sigmahat[0]=new Positive("vcs", pars[sigma][0]);
	sigmahat[1]=new Positive("vncs", pars[sigma][1]);
	sigmahat[2]=new Correlation("corr", pars[sigma][2]);
	deltahat= new Determined("disc", pars[delta]);
	rhohat=new Determined("rho", pars[rho]);
	tfphat=new array[Msector];
	tfphat[0]=new Positive("tfpcs", pars[tfp][0]);
	tfphat[1]=new Positive("tfpncs", pars[tfp][1]);
	sharehat=new array[Msector];
	sharehat[0]=new Probability("sharecs",pars[share][0]);
	sharehat[1]=new Probability("sharencs",pars[share][1]);
	elashat=new array[Msector];
	elashat[0]=new Free("elas", pars[elas][0]);
	elashat[1]=new Free("elas", pars[elas][1]);
	decl obj, equimethod;
	obj=new MarketClear("Equilibrium", omigahat, alphahat, betahat, sigmahat, deltahat, rhohat, tfphat, sharehat, elashat);
	ToggleParams(alphahat, betahat, sigmahat, tfphat, sharehat, elashat);
	equimethod=new NelderMead(obj);	   
	equimethod.Volume=LOUD;
	equimethod-> Iterate(0);

    }


/** Setup the DP model and first stage estimation.
@param row 0 or 1, row of Rust table to replicate (&delta;=0.0 or &delta;=0.999)
**/
/**set up the DP model**/
CSprojectHat::Run(){
    decl Ttotal=A1-A0, i, AMat, ps, simdata, acs, ancs, cho, occw, wcs, wncs, occl, rls, Nsim=1000, adFunc;
	decl tempcov, tempsigma;
	Initialize(Reachable, TRUE, 0);
	SetClock(NormalAging, Ttotal-1);
	Actions(occupations=new ActionVariable("Occupation Choice", Msector));
	xper= new array[Msector];
//	xper[0]=new RActionCounter("ExCS", Ttotal-1, occupations, 0, betahat[0]);
    xper[0]=new ActionCounter("ExCS", Ttotal-1, occupations, 0, 0);
//	xper[1]=new RActionCounter("ExNcs", Ttotal-1, occupations, 1, betahat[1]);
    xper[1]=new ActionCounter("ExNCS", Ttotal-1, occupations, 1, 0);
	EndogenousStates(xper[0], xper[1]);
	tempcov=CV(sigmahat[2])*sqrt(CV(sigmahat[0])*CV(sigmahat[1]));
	tempsigma=(CV(sigmahat[0])~tempcov)|(tempcov~CV(sigmahat[1]));
	types=new BiNormalRandomEffect("Heterogeniety", Ntypes, pars[mu], tempsigma);
	nodes=<>;
	GaussianHermite::Initialize(Ntypes);
	for (i=0;i<Ntypes;++i){
	    this.nodes|=ones(Ntypes,1)*GaussianHermite::nodes[i]~GaussianHermite::nodes';
		}
	nodes=choleski(tempsigma)*nodes';
	GroupVariables(types);
	AuxiWage=new Wage();
	AuxiliaryOutcomes(AuxiWage);
	SetDelta(deltahat);
	CreateSpaces(LogitKernel,rhohat);
	Method= new ValueIteration(0);
	// Volume=LOUD;
	Method->Solve();
	//DPDebug::outV(FALSE, &AMat);
	ps= new Panel(0, zeros(8,1));
	ps->Simulate(Nsim,Ttotal, zeros(8,1),TRUE);
	ps->Flat();
	simdata=ps.flat;
	acs=reshape(simdata[][columns(simdata)-3], Ttotal, Nsim);
	ancs=1-acs;
	occw=reshape(simdata[][columns(simdata)-2],Ttotal, Nsim);
	occl=reshape(simdata[][columns(simdata)-1],Ttotal, Nsim);
    laborsupply=sumc(sumr(occl.*acs))/(Ttotal*Nsim)|sumc(sumr(occl.*ancs))/(Ttotal*Nsim);
	/**labor demand side **/
	labordemand=<50;100>;
	MaxBFGS(ldemand, &labordemand, &adFunc, 0, 1);
	return -(laborsupply-labordemand)'*(laborsupply-labordemand);
	//println(laborsupply~labordemand);
    }

CSprojectHat::Reachable(){
   decl i,totexp;
   for (i=0,totexp=0;i<Msector;++i) totexp +=xper[i].v;
   return curt<totexp ? 0: new CSprojectHat();
   }

CSprojectHat::Utility(){
   decl xcs=xper[cs].v, xncs=xper[ncs].v, R;
   decl rr=((xcs*CV(alphahat[0])+xncs*CV(alphahat[1]))|(xcs*CV(alphahat[2])+xncs*CV(alphahat[3])))+this.nodes[][types.v];
   //decl rr=(x*<CV(alphahat[0]);CV(alphahat[1])>|x*<CV(alphahat[2]);CV(alphahat[3])>)+this.nodes[][types.v];
   R=(CV(omigahat[0])|CV(omigahat[1]))+rr;
   return R[A[Aind]];
   }   

CSprojectHat::Auxi(){
   decl xcs=xper[cs].v	, xncs=xper[ncs].v, R;
   decl rr=((xcs*CV(alphahat[0])+xncs*CV(alphahat[1]))|(xcs*CV(alphahat[2])+xncs*CV(alphahat[3])))+this.nodes[][types.v];
   //decl rr=(x*<CV(alphahat[0]);CV(alphahat[1])>|x*<CV(alphahat[2]);CV(alphahat[3])>)+this.nodes[][types.v];
   R=exp(rr);
   return R[A[Aind]];
   }

//AuxiliaryVariable::Realize(q,y){
//   decl v1=q->Utility()[q.ialpha][];
//   decl v2=q->Auxi()[q.ialpha][];
//   v=v1~v2;
//   }

CSprojectHat::ldemand(const vP, const adfunc, const avScore, const amHess){
   decl dif1, dif2;
   dif1=CV(tfphat[0])*(CV(sharehat[0])*vP[0]^CV(elashat[0])+(1-CV(sharehat[0]))*M[0]^CV(elashat[0]))^(1/CV(elashat[0])-1)
       *CV(sharehat[0])*vP[0]^(CV(elashat[0])-1)-exp(CV(omigahat[0]));
   dif2=CV(tfphat[1])*(CV(sharehat[1])*vP[1]^CV(elashat[1])+(1-CV(sharehat[1]))*M[1]^CV(elashat[1]))^(1/CV(elashat[1])-1)
       *CV(sharehat[1])*vP[1]^(CV(elashat[1])-1)-exp(CV(omigahat[1]));
   adfunc[0]=-dif1^2-dif2^2;
   return 1;
   }   
   
MarketClear::MarketClear(L,...){
   BlackBox(L);
   decl va=va_arglist(), i;
   if (sizeof(va)){
      for(i=0;i<sizeof(va); ++i) Parameters(va[i]);	   // adding parameters to be estimated
	  Encode(0);
      }
   else oxwarning("No estimated parameters added to "+L+"estimation");
   }

MarketClear::vfunc(){
   decl r=CSprojectHat::Run();
   return r;

   }



	