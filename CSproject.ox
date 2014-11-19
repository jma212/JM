#include <CSproject.h>
/** Setup and solve the model.**/

CSproject::Run(){
   decl Ttotal=A1-A0, i, AMat, ps, simdata, cho;
   Initialize(Reachable,TRUE,0);
   SetClock(NormalAging, Ttotal-1);
   Actions(occupations= new ActionVariable("Occupation Choice", Msector));
   xper=new array[Msector];
   //xper[0]=new RandomUpDown("ExCS", Ttotal,  fPi0);
   //xper[1]=new RandomUpDown("ExNCS", Ttotal, fPi1);
   xper[0]=new RActionCounter("ExCS", Ttotal-1, occupations, 0, beta[0]);
   xper[1]=new RActionCounter("ExNCS", Ttotal, occupations, 1, beta[1]);
   EndogenousStates(xper[0],xper[1]);
   types=new BiNormalRandomEffect("Heterogeniety", Ntypes, mu, sigma);
   nodes=<>;
   GaussianHermite::Initialize(Ntypes);
   for (i=0;i<Ntypes; ++i){
       this.nodes|=ones(Ntypes,1)*GaussianHermite::nodes[i]~GaussianHermite::nodes';
       }
   nodes=mu+choleski(sigma)*nodes';
   GroupVariables(types);
   AuxiWage= new Wage();
   AuxiliaryOutcomes(AuxiWage);
   SetDelta(delta);
   CreateSpaces(LogitKernel,rho);
   
   Method= new ValueIteration(0);
   //Volume=LOUD;
   Method->Solve();
   //DPDebug::outV(FALSE,&AMat);
   ps=new Panel(0,zeros(8,1));
   ps->Simulate(100,Ttotal,zeros(8,1),TRUE);	 
   ps->Flat();
   simdata=ps.flat;
   cho=sumr(reshape(simdata[][columns(simdata)-3], Ttotal, 100))/100;
   println(cho);   // simulated percentage of people choose to be CS worker
   
   }

CSproject::Reachable(){
   decl i,totexp;
   for (i=0,totexp=0;i<Msector;++i) totexp +=xper[i].v;
   return curt<totexp ? 0: new CSproject();
   }

CSproject::Utility(){
   decl xcs=xper[cs].v, xncs=xper[ncs].v, x=xcs~xncs, R;
   decl rr=(x*alpha[0]|x*alpha[1])+this.nodes[][types.v];
   R=omiga+rr;
   return R[A[Aind]];
   }

CSproject::Auxi(){
   decl xcs=xper[cs].v	, xncs=xper[ncs].v, x=xcs~xncs, R;
   decl rr=(x*alpha[0]|x*alpha[1])+this.nodes[][types.v];
   R=exp(rr);
   return R[A[Aind]];
   }


CSproject::fPi0(FeasA){
   decl inda=sumr(FeasA[][occupations.pos].==cs);
   return zeros(rows(FeasA),1)~(1-inda*exp(-beta[0]*curt))~(inda*exp(-beta[0]*curt));
   }

CSproject::fPi1(FeasA){
   decl inda=sumr(FeasA[][occupations.pos].==ncs);
   return zeros(rows(FeasA),1)~(1-inda*exp(-beta[1]*curt))~(inda*exp(-beta[1]*curt));
   }

