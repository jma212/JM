#include <KW1997.h>
KW1997::Run(){
decl i, sigma, Ttotal=A1-A0, BF, KW;
decl data, newd, ps;
Initialize(Reachable, TRUE, 0);
SetClock(NormalAging, Ttotal);
Actions(accept= new ActionVariable("Accept", Msector));
sigma=vech(choleski(sd.*unvech(rho).*sd'));
ExogenousStates(offers= new MVNormal("eps", Msector, Noffers, zeros(Msector,1), sigma));
GroupVariables(types= new RandomEffect("Type", Ntypes));
xper= new array[Msector-1];
for (i=0;i<Msector-1;++i){
EndogenousStates(xper[i]=new ActionCounter("X"+sprint(i), Ttotal, accept, i, 0));
}
SetDelta(0.7870);
CreateSpaces(LogitKernel, 1/40000.0);
Volume = LOUD;
//BF = new ValueIteration();
//BF ->Solve();
KW = new KeaneWolpin(ones(1,5)~constant(0.8,1,5),0);
KW->Solve();
//Vprint(TRUE);
ps=new FPanel(0,0,FALSE);
ps->Simulate(10, 3, ones(14,1), TRUE);
//newd=ps->Print(0);
//DPDebug::outV(TRUE,0);
}

KW1997::Reachable(){
decl i,totexp;
for (i=0, totexp=0; i<Msector-1 ; ++i ) totexp +=xper[i].v;
return curt<totexp ? 0 : new KW1997();
}

KW1997::Utility(){
//decl chol=choleski(sd.*unvech(rho).*sd');
//decl eps=sqrt(2)*chol*(selectrc(offers.Grid, accept.vals, offers.v))';
decl eps=selectrc(offers.Grid, accept.vals, offers.v)';
decl xs=xper[school].v, xw=xper[white].v, xb=xper[blue].v, xm=xper[military].v,
     xbw=(xs~xw~xb~xm~sqr(xw))*alpha[white],
	 xbb=(xs~xw~xb~xm~sqr(xb))*alpha[blue],
	 xbm=(xs~xw~xb~xm~sqr(xm))*alpha[military],
	R=xbw
	 |xbb
	 |xbm
	 |tc[0]*(xs+School0>=HSGrad)+tc[1]*(xs+School0>=CGrad)
	 |0;
decl betaconst=beta[0]~(beta[0]+beta[1])~(beta[0]+beta[2])~(beta[0]+beta[3]);
R[:military]=exp(R[:military]+eps[:military]+betaconst[:military][types.v]);
R[school:]=R[school:]+eps[school:]+betaconst[school:][types.v];
return R[A[Aind]];
}