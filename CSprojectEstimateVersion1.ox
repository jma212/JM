#include "CSprojectEstimateVersion1.h"
#import <maximize>

/**find the market clear wage**/
CSprojectHat::FirstStage(){
	decl obj, objsim, method, methodsim, data;
	CSprojectHat::SetUp();
//	obj=new MarketClear("MarketClear", hat1, hat2);
//	ToggleParams(hat2);
//	method=new NelderMead(obj);
//	method.Volume=LOUD;
//	method->Iterate(0);		  // now the true parameter that generates the model is find
	data=CSprojectHat::SimulateMoments();
	data=data[1:][]~150;
//	ToggleParams(hat2);
	objsim= new SMM("SMM", data, hat1, hat2);
	methodsim= new NelderMead(objsim);
	methodsim.Volume=QUIET;
	methodsim->Iterate(0);
	}

/**set up the DP model**/
CSprojectHat::SetUp(){
    decl Ttotal=A1-A0, i, crtbeta;
//	hat1=new StDeviations("logw",DGP1);
    hat1=new StDeviations("logw", <2.418;2.711>);
	hat2=new array[Nparams];
	hat2[alpha]=new Coefficients("ParaWage", DGP2[alpha]);
	hat2[beta]=new ParameterBlock("b", new Positive("bcs",DGP2[beta][0]),new Positive("bncs", DGP2[beta][1]));
	hat2[mu]=new FixedBlock("mu",DGP2[mu]);
	hat2[sigma]=new Coefficients("cholesky", vech(choleski(DGP2[sigma])));
	hat2[delta]=new Determined("disc",DGP2[delta]);
	hat2[rho]=new Determined("rho",DGP2[rho]);
	hat2[tfp]=new StDeviations("tfp",DGP2[tfp]);
	hat2[share]=new ParameterBlock("share", new Probability("sharecs", DGP2[share][0]), new Probability("sharencs", DGP2[share][1]));
	hat2[elas]=new Coefficients("elas",DGP2[elas]);
	decl chol, tempsigma;
	Initialize(Reachable, TRUE, 0);
	SetClock(NormalAging, Ttotal-1);
	Actions(occupations=new ActionVariable("Occupation Choice", Msector));
	xper=new array[Msector];
	crtbeta=CV(hat2[beta]);
	xper[0]=new RActionCounter("ExCS", Ttotal-1, occupations, 0 , crtbeta[0]);
	xper[1]=new RActionCounter("ExNCS", Ttotal-1, occupations, 1, crtbeta[1]);
	EndogenousStates(xper[0], xper[1]);
	chol=lower(unvech(CV(hat2[sigma])));
	tempsigma=chol*chol';
	types=new BiNormalRandomEffect("Heterogeniety", Ntypes, CV(hat2[mu]), tempsigma);
	nodes=<>;
	GaussianHermite::Initialize(Ntypes);
	for (i=0;i<Ntypes;++i){
	    this.nodes|=ones(Ntypes,1)*GaussianHermite::nodes[i]~GaussianHermite::nodes';
		}
	nodes=chol*nodes';
	GroupVariables(types);
	AuxiWage= new AuxiliaryVariable("Wage");
	AuxiliaryOutcomes(AuxiWage);
	SetDelta(CV(hat2[delta]));
	CreateSpaces(LogitKernel, CV(hat2[rho]));
    }


CSprojectHat::Run(){
	decl Ttotal=A1-A0, i, AMat, ps, simdata, acs, ancs, occw, occl, rls, Nsim=50, adFunc, crtbeta;
	Method=new ValueIteration(0);
	Method->Solve();
	ps=new Panel(0,zeros(8,1));
	ps->Simulate(Nsim, Ttotal, zeros(8,1),TRUE);
	ps->Flat();
	simdata=ps.flat;
    acs=reshape(simdata[][columns(simdata)-3],Ttotal,Nsim);
	ancs=1-acs;
	occw=reshape(simdata[][columns(simdata)-2],Ttotal, Nsim);
	occl=reshape(simdata[][columns(simdata)-1],Ttotal, Nsim);
	laborsupply=sumc(sumr(occl.*acs))/(Ttotal*Nsim)|sumc(sumr(occl.*ancs))/(Ttotal*Nsim);
	labordemand=<50;100>;
	MaxBFGS(ldemand, &labordemand, &adFunc, 0, 1);
	return -(laborsupply-labordemand)'*(laborsupply-labordemand);
	delete Method, ps;
	}

CSprojectHat::SimulateMoments(){
    decl Ttotal=A1-A0, i, simdata, ps, acs, ancs, occw, occl, rls, Nsim=200, adFunc, momentscs, momentsncs, R;
	Method=new ValueIteration(0);
	Method.vtoler = 1E-01;
	Method->Solve();
	ps=new Panel(0,zeros(8,1));
	ps->Simulate(Nsim, Ttotal, zeros(8,1),TRUE);
	ps->Flat();
	simdata=ps.flat;
	acs=reshape(simdata[][columns(simdata)-3],Ttotal,Nsim);
	ancs=1-acs;
	occw=reshape(simdata[][columns(simdata)-2],Ttotal, Nsim);
	occl=reshape(simdata[][columns(simdata)-1],Ttotal, Nsim);
	laborsupply=sumc(sumr(occl.*acs))/(Ttotal*Nsim)|sumc(sumr(occl.*ancs))/(Ttotal*Nsim);
	labordemand=<50;100>;
	MaxBFGS(ldemand, &labordemand, &adFunc, 0, 1);
	momentscs=moments((replace(occw.*acs, 0, .NaN))',2);     // store wage in cs
	momentsncs=moments((replace(occw.*ancs, 0, .NaN))',2);	 // store wage in ncs
	momentscs[0][]/=Nsim;
	momentsncs=momentsncs[1:][];
	return R=(laborsupply-labordemand)|vecr(momentscs)|vecr(momentsncs); 
	delete Method, ps;
    }
	
CSprojectHat::Reachable(){
	decl i,totexp;
	for (i=0,totexp=0;i<Msector; ++i) totexp +=xper[i].v;
	return curt<totexp ? 0: new CSprojectHat();
    }

CSprojectHat::Utility(){
	decl  xcs=xper[cs].v,xncs=xper[ncs].v, R, crtalpha, crtomiga;
	crtalpha=CV(hat2[alpha]);
	crtomiga=CV(hat1);
	decl rr=((xcs*crtalpha[0]+xncs*crtalpha[1])|(xcs*crtalpha[2]+xncs*crtalpha[3]))+this.nodes[][types.v];
	R=crtomiga+rr;
    return R[A[Aind]];
    }

	
CSprojectHat::Auxi(){
	decl xcs=xper[cs].v, xncs=xper[ncs].v, R, crtalpha;
	crtalpha=CV(hat2[alpha]);
	decl rr=((xcs*crtalpha[0]+xncs*crtalpha[1])|(xcs*crtalpha[2]+xncs*crtalpha[3]))+this.nodes[][types.v];
	R=exp(rr);
	return R[A[Aind]];
	}

AuxiliaryVariable::Realize(q,y){
	decl v1=q->Utility()[q.ialpha][];
	decl v2=q->Auxi()[q.ialpha][];
	v=v1~v2;
	}

CSprojectHat::ldemand(const vP, const adfunc, const avScore, const amHess){
	 decl dif1, dif2, crttfp, crtomiga, crtelas, crtshare;
	 crttfp=CV(hat2[tfp]);
	 crtomiga=CV(hat1);
	 crtshare=CV(hat2[share]);
	 crtelas=CV(hat2[elas]);
	 dif1=crttfp[0]*(crtshare[0]*vP[0]^crtelas[0]+(1-crtshare[0])*M[0]^crtelas[0])^(1/crtelas[0]-1)
       *crtshare[0]*vP[0]^(crtelas[0]-1)-exp(crtomiga[0]);
     dif2=crttfp[1]*(crtshare[1]*vP[1]^crtelas[1]+(1-crtshare[1])*M[1]^crtelas[1])^(1/crtelas[1]-1)
       *crtshare[1]*vP[1]^(crtelas[1]-1)-exp(crtomiga[1]);
     adfunc[0]=-dif1^2-dif2^2;
     return 1;
    }
	
MarketClear::MarketClear(L, hat1, hat2){
	 BlackBox(L);
	 Parameters(hat1, hat2);
	 Encode(0);
	 }

MarketClear::vfunc(){
	 decl r=CSprojectHat::Run();
	 return r;
     }

SMM::SMM(L, data, hat1, hat2){  // the first column is the data moment and second column is the corresponding sample size
	 this.w=diag(max(data[][1])|data[][1]);
	 this.datamoments=(0|data[][0]);
	 BlackBox(L);
	 Parameters(hat1,hat2);
	 Encode(0);
     }

SMM::vfunc(){
     decl simulatemoments=CSprojectHat::SimulateMoments();
	 decl diff=(this.datamoments-simulatemoments);
	 decl minidis=diff'*this.w*diff;
	 return minidis;	 
     }


