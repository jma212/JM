 // replicate Keane and Wolpin 1997 basic model 
 // include title
# import "DDP"
# import "DP"
struct KW1997:ExPostSmoothing{
	decl InSubSample; // the public data used in the Keane and Wolpin approximation method
	enum {white, blue, military, school, home, Msector};
	enum {A1=36,
      A0=16,
	  MxExper = 10,
	  Ttotal=A1-A0, 
	  HSGrad=12,
	  CGrad=16,
	  Noffers=5,
	  Ntypes=1,
	  School0=10};
	static const decl
	/**  wages   **/       alpha={<0.0938;0.1170;0.0748;0.0077;-0.0461>,
	                          <0.0189;0.0674;0.1424;0.1021;-0.1774>,
							  <0.0443;0.0;0.0;0.3391;-29900>},
	/** e(16) types **/     beta={<8.8043;8.9156;8.4704;43948;16887>,
	                          <-0.0668;0.2996;0;-26352;215>,
							  <-0.4221;-0.1223;0;-30541;-16966>,
							  <-0.4998;0.0756;0;226;-13128>},
	/**schooling cost**/      tc=<-2983;-26357>,
	/**correlation matrix**/ rho=<1.001;-0.3806;-0.3688;0;0;1;0.4120;0;0;1;0;0;1;0;1>,
	/**  std  **/             sd=<0.3301;0.3329;0.3308;2312;13394>,
	/**  discount  **/	   delta=0.7870;

	static decl
	/** index accepted offer **/       accept,
	/** offer block **/                offers,
	/**occupation experience array**/  xper,
	/** type indicator **/             types;
	static Run();
	static Reachable();
       Utility();	   							  
}