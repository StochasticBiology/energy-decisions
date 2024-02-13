#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  // parameters
  double gamma2, gamma3, gamma4;
  double cd1, ca1, cd2, ca2;
  double lambda2, lambda3;
} Params;

int main(int argc, char *argv[])
{
  Params Pdef, P;
  
  // state variables
  double pro_1, pro_2;
  double prooff_1, prooff_2;
  double rna_1, rna_2;
  double p_1, p_2;
  double pp_1, pp_2;
  double ppp_1, ppp_2;
  double pppp_1, pppp_2;
  double z_1, z_2;

  // derivatives
  double dpro_1dt, dpro_2dt;
  double dprooff_1dt, dprooff_2dt;
  double drna_1dt, drna_2dt;
  double dp_1dt, dp_2dt;
  double dpp_1dt, dpp_2dt;
  double dppp_1dt, dppp_2dt;
  double dpppp_1dt, dpppp_2dt;
  
  double ip1, ip2;
  FILE *fp;
  double t;
  double dt = 0.005;
  double eps = 1e-6;
  int converge;
  int param, param2, scale;
  int i;
  int n;
  double s1, s2;

  Pdef.gamma2 = 0.00023;
  Pdef.gamma3 = 0.00077;
  Pdef.gamma4 = 0.00058;
  Pdef.cd1 = 0.00023;
  Pdef.cd2 = 1.04;
  Pdef.ca1 = 0.0023;
  Pdef.ca2 = 0.1038; 
  Pdef.lambda2 = 0.0067;
  Pdef.lambda3 = 0.0474;
  
  if(argc != 3)
    {
      printf("I need a cooperativity value (1-4) and a scale perturbation (0-3)\n");
      exit(0);
    }

  n = atoi(argv[1]);
  scale = atoi(argv[2]);

  // for(n = 1; n <= 4; n++)
  {
    
    for(scale = 0; scale <= 7; scale++)
      {
	switch(scale)
	  {
	  case 0: s1 = 0.5; s2 = 0.5; break;
	  case 1: s1 = 2; s2 = 0.5; break;
	  case 2: s1 = 0.5; s2 = 2; break;
	  case 3: s1 = 2; s2 = 2; break;
	  case 4: s1 = 0.1; s2 = 0.1; break;
	  case 5: s1 = 10; s2 = 0.1; break;
	  case 6: s1 = 0.1; s2 = 10; break;
	  case 7: s1 = 10; s2 = 10; break;
	  }
					
	for(param = 0; param < 9; param++)
	  {
	    for(param2 = param+1; param2 < 9; param2++)
	      {
		P = Pdef;

		switch(param)
		  {
		  case 0: P.gamma2 *= s1; break;
		  case 1: P.gamma3 *= s1; break;
		  case 2: P.gamma4 *= s1; break;
		  case 3: P.cd1 *= s1; break;
		  case 4: P.cd2 *= s1; break;
		  case 5: P.ca1 *= s1; break;
		  case 6: P.ca2 *= s1; break;
		  case 7: P.lambda2 *= s1; break;
		  case 8: P.lambda3 *= s1; break;
		  }

		switch(param2)
		  {
		  case 0: P.gamma2 *= s2; break;
		  case 1: P.gamma3 *= s2; break;
		  case 2: P.gamma4 *= s2; break;
		  case 3: P.cd1 *= s2; break;
		  case 4: P.cd2 *= s2; break;
		  case 5: P.ca1 *= s2; break;
		  case 6: P.ca2 *= s2; break;
		  case 7: P.lambda2 *= s2; break;
		  case 8: P.lambda3 *= s2; break;
		  }

		// loop through initial protein concentrations
		for(ip1 = 0; ip1 < 30; ip1 += 3)
		  {
		    for(ip2 = 0; ip2 < 30; ip2 += 3)
		      {
			printf("%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", n,param,param2,s1,s2,P.gamma2,P.gamma3,P.gamma4,P.cd1,P.cd2,P.ca1,P.ca2,P.lambda2,P.lambda3, ip1, ip2);
		      
			// initialise system
			p_1 = ip1; p_2 = ip2;
			pro_1 = 1; pro_2 = 1; z_1 = 1; z_2 = 1;
			rna_1 = rna_2 = prooff_1 = prooff_2 = pp_1 = pp_2 = ppp_1 = ppp_2 = pppp_1 = pppp_2 = 0;
			converge = 0;

			// simple euler solver
			for(t = 0; converge < 500./dt; t += dt)
			  {
			    // time derivatives from Rajneesh's doc
			    drna_1dt   = P.lambda2*pro_1*z_1-P.gamma2*rna_1;
			    drna_2dt   = P.lambda2*pro_2*z_2-P.gamma2*rna_2;

			    if(n == 1)
			      {
				//           translate     degrade    DNAbind       DNAunbind
				dp_1dt     = P.lambda3*rna_1-P.gamma3*p_1-P.ca2*pro_2*p_1+P.cd2*prooff_2;
				dp_2dt     = P.lambda3*rna_2-P.gamma3*p_2-P.ca2*pro_1*p_2+P.cd2*prooff_1;
				//           DNAunbind    DNAbind
				dpro_1dt   = P.cd2*prooff_1-P.ca2*pro_1*p_2;            
				dpro_2dt   = P.cd2*prooff_2-P.ca2*pro_2*p_1;
				//             DNAunbind    DNAbind
				dprooff_1dt = -P.cd2*prooff_1+P.ca2*pro_1*p_2;
				dprooff_2dt = -P.cd2*prooff_2+P.ca2*pro_2*p_1;
				dpp_1dt    = 0;
				dpp_2dt    = 0;
				dppp_1dt    = 0;
				dppp_2dt    = 0;
				dpppp_1dt    = 0;
				dpppp_2dt    = 0;
			      }
		      
			    if(n == 2)
			      {
				//           translate     dimer         dedimer    degrade
				dp_1dt     = P.lambda3*rna_1-2*P.ca1*p_1*p_1+2*P.cd1*pp_1-P.gamma3*p_1;
				dp_2dt     = P.lambda3*rna_2-2*P.ca1*p_2*p_2+2*P.cd1*pp_2-P.gamma3*p_2;
				//           DNAunbind    DNAbind
				dpro_1dt   = P.cd2*prooff_1-P.ca2*pro_1*pp_2;            
				dpro_2dt   = P.cd2*prooff_2-P.ca2*pro_2*pp_1;
				//             DNAunbind    DNAbind
				dprooff_1dt = -P.cd2*prooff_1+P.ca2*pro_1*pp_2;
				dprooff_2dt = -P.cd2*prooff_2+P.ca2*pro_2*pp_1;
				//           dimer       dedimer  DNAbind        degrade     DNAunbind
				dpp_1dt    = P.ca1*p_1*p_1-P.cd1*pp_1-P.ca2*pp_1*pro_2-P.gamma4*pp_1+P.cd2*prooff_2;
				dpp_2dt    = P.ca1*p_2*p_2-P.cd1*pp_2-P.ca2*pp_2*pro_1-P.gamma4*pp_2+P.cd2*prooff_1;
				dppp_1dt    = 0;
				dppp_2dt    = 0;
				dpppp_1dt    = 0;
				dpppp_2dt    = 0;

			      }
			    if(n == 3)
			      {
				//           translate     dimer         dedimer    degrade    trimer       detrimer
				dp_1dt     = P.lambda3*rna_1-2*P.ca1*p_1*p_1+2*P.cd1*pp_1-P.gamma3*p_1-P.ca1*p_1*pp_1+P.cd1*ppp_1;
				dp_2dt     = P.lambda3*rna_2-2*P.ca1*p_2*p_2+2*P.cd1*pp_2-P.gamma3*p_2-P.ca1*p_2*pp_2+P.cd1*ppp_2;
				//           DNAunbind    DNAbind
				dpro_1dt   = P.cd2*prooff_1-P.ca2*pro_1*ppp_2;            
				dpro_2dt   = P.cd2*prooff_2-P.ca2*pro_2*ppp_1;
				//             DNAunbind    DNAbind
				dprooff_1dt = -P.cd2*prooff_1+P.ca2*pro_1*ppp_2;
				dprooff_2dt = -P.cd2*prooff_2+P.ca2*pro_2*ppp_1;
				//           dimer       dedimer  degrade     trimer       detrimer
				dpp_1dt    = P.ca1*p_1*p_1-P.cd1*pp_1-P.gamma4*pp_1-P.ca1*p_1*pp_1+P.cd1*ppp_1;
				dpp_2dt    = P.ca1*p_2*p_2-P.cd1*pp_2-P.gamma4*pp_2-P.ca1*p_2*pp_2+P.cd1*ppp_2;
				//            trimer       detrimer  DNAunbind       degrade      DNAbind      
						dppp_1dt    = P.ca1*pp_1*p_1-P.cd1*ppp_1-P.ca2*ppp_1*pro_2-P.gamma4*ppp_1+P.cd2*prooff_2;
				dppp_2dt    = P.ca1*pp_2*p_2-P.cd1*ppp_2-P.ca2*ppp_2*pro_1-P.gamma4*ppp_2+P.cd2*prooff_1;
				dpppp_1dt    = 0;
				dpppp_2dt    = 0;

			      }
			    if(n == 4)
			      {
				//           translate     dimer         dedimer    degrade    trimer       detrimer  tetramer      detetramer
				dp_1dt     = P.lambda3*rna_1-2*P.ca1*p_1*p_1+2*P.cd1*pp_1-P.gamma3*p_1-P.ca1*p_1*pp_1+P.cd1*ppp_1-P.ca1*p_1*ppp_1+P.cd1*pppp_1;
				dp_2dt     = P.lambda3*rna_2-2*P.ca1*p_2*p_2+2*P.cd1*pp_2-P.gamma3*p_2-P.ca1*p_2*pp_2+P.cd1*ppp_2-P.ca1*p_2*ppp_2+P.cd1*pppp_2;
				//           DNAunbind    DNAbind
				dpro_1dt   = P.cd2*prooff_1-P.ca2*pro_1*pppp_2;            
				dpro_2dt   = P.cd2*prooff_2-P.ca2*pro_2*pppp_1;
				//             DNAunbind    DNAbind
				dprooff_1dt = -P.cd2*prooff_1+P.ca2*pro_1*pppp_2;
				dprooff_2dt = -P.cd2*prooff_2+P.ca2*pro_2*pppp_1;
				//           dimer       dedimer  degrade     trimer       detrimer
				dpp_1dt    = P.ca1*p_1*p_1-P.cd1*pp_1-P.gamma4*pp_1-P.ca1*p_1*pp_1+P.cd1*ppp_1;
				dpp_2dt    = P.ca1*p_2*p_2-P.cd1*pp_2-P.gamma4*pp_2-P.ca1*p_2*pp_2+P.cd1*ppp_2;
				//            trimer       detrimer  degrade      tetramer      detetramer
				dppp_1dt    = P.ca1*pp_1*p_1-P.cd1*ppp_1-P.gamma4*ppp_1-P.ca1*ppp_1*p_1+P.cd1*pppp_2;
				dppp_2dt    = P.ca1*pp_2*p_2-P.cd1*ppp_2-P.gamma4*ppp_2-P.ca1*ppp_2*p_2+P.cd1*pppp_2;
				//             tetramer      detetramer DNAbind          degrade       DNAunbind
				dpppp_1dt    = P.ca1*ppp_1*p_1-P.cd1*pppp_1-P.ca2*pppp_1*pro_2-P.gamma4*pppp_1+P.cd2*prooff_2;
				dpppp_2dt    = P.ca1*ppp_2*p_2-P.cd1*pppp_2-P.ca2*pppp_2*pro_1-P.gamma4*pppp_2+P.cd2*prooff_1;
			      }
		      
			    // test for convergence
			    if(!(dpro_1dt > eps*pro_1 || dpro_2dt > eps*pro_2 || dprooff_1dt > eps*prooff_1 || dprooff_2dt > eps*prooff_2 || drna_1dt > eps*rna_1 || drna_2dt > eps*rna_2 || dp_1dt > eps*p_1 || dp_2dt > eps*p_2 || dpp_1dt > eps*pp_1 || dpp_2dt > eps*pp_2 || dppp_1dt > eps*ppp_1 || dppp_2dt > eps*ppp_2 || dpppp_1dt > eps*pppp_1 || dpppp_2dt > eps*pppp_2))
			      converge++;
			    else converge = 0;

			    // update state of system
			    pro_1      += dt*dpro_1dt;
			    pro_2      += dt*dpro_2dt;
			    prooff_1    += dt*dprooff_1dt;
			    prooff_2    += dt*dprooff_2dt;
			    rna_1      += dt*drna_1dt;
			    rna_2      += dt*drna_2dt;
			    p_1        += dt*dp_1dt;
			    p_2        += dt*dp_2dt;
			    pp_1       += dt*dpp_1dt;
			    pp_2       += dt*dpp_2dt;
			    ppp_1       += dt*dppp_1dt;
			    ppp_2       += dt*dppp_2dt;
			    pppp_1       += dt*dpppp_1dt;
			    pppp_2       += dt*dpppp_2dt;
			  }

			// output start and end points for this sim
		        fp = fopen("param-matrix-scan.csv", "a");
			fprintf(fp, "%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", n,param,param2,s1,s2,P.gamma2,P.gamma3,P.gamma4,P.cd1,P.cd2,P.ca1,P.ca2,P.lambda2,P.lambda3, ip1, ip2, p_1, p_2);
			fclose(fp);
		      }
		  }
	      }

	  }
      }
  }    
  
  return 0;
}
