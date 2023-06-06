#include <stdio.h>
#include <math.h>

int main(void)
{
  // parameters
  double gamma2, gamma3, gamma4;
  double cd1, ca1, cd2, ca2;
  double lambda2, lambda3;

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
  int param, scale;
  int i;
  int n;
  double ATP;
  
  gamma2 = 0.00023;
  gamma3 = 0.00077;
  gamma4 = 0.00058;
  cd1 = 0.00023;
  cd2 = 1.04;
  ca1 = 0.0023;
  ca2 = 0.1038; 
  lambda2 = 0.0067;
  lambda3 = 0.0474;
  
  // open file for output
  fp = fopen("zoom-scan.csv", "w");
  fprintf(fp, "n,param,scale,ATP,gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3,ip1,ip2,p1,p2\n");
  fclose(fp);

  n = 2;
  
  for(param = 0; param < 2; param++)
    {
      switch(param)
	{
	case 0: gamma3 = 0.00077/4; cd2 = 1.04; break;
	case 1: gamma3 = 0.00077; cd2 = 0.01/4; break;
	}
      
      for(scale = 0; scale <= 5; scale++)
	{
	  switch(param)
	    {
	    case 0: gamma3 *= 4; break;
	    case 1: cd2 *= 2; break;
	    }

	  for(ATP = 0.25; ATP <= 2; ATP *= 2)
	    {
	      lambda2 = 0.0067*ATP;
	      lambda3 = 0.0474*ATP;
 
	      // loop through initial protein concentrations
	      for(ip1 = 0; ip1 < 30; ip1 += 3)
		{
		  for(ip2 = 0; ip2 < 30; ip2 += 3)
		    {
		      printf("%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", n,param,scale,ATP,gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3, ip1, ip2);

		      // initialise system
		      p_1 = ip1; p_2 = ip2;
		      pro_1 = 1; pro_2 = 1; z_1 = 1; z_2 = 1;
		      rna_1 = rna_2 = prooff_1 = prooff_2 = pp_1 = pp_2 = ppp_1 = ppp_2 = pppp_1 = pppp_2 = 0;
		      converge = 0;

		      // simple euler solver
		      for(t = 0; converge < 500./dt; t += dt)
			{
			  // time derivatives from Rajneesh's doc
			  drna_1dt   = lambda2*pro_1*z_1-gamma2*rna_1;
			  drna_2dt   = lambda2*pro_2*z_2-gamma2*rna_2;

			  //           translate     dimer         dedimer    degrade
			  dp_1dt     = lambda3*rna_1-2*ca1*p_1*p_1+2*cd1*pp_1-gamma3*p_1;
			  dp_2dt     = lambda3*rna_2-2*ca1*p_2*p_2+2*cd1*pp_2-gamma3*p_2;
			  //           DNAunbind    DNAbind
			  dpro_1dt   = cd2*prooff_1-ca2*pro_1*pp_2;            
			  dpro_2dt   = cd2*prooff_2-ca2*pro_2*pp_1;
			  //             DNAunbind    DNAbind
			  dprooff_1dt = -cd2*prooff_1+ca2*pro_1*pp_2;
			  dprooff_2dt = -cd2*prooff_2+ca2*pro_2*pp_1;
			  //           dimer       dedimer  DNAbind        degrade     DNAunbind
			  dpp_1dt    = ca1*p_1*p_1-cd1*pp_1-ca2*pp_1*pro_2-gamma4*pp_1+cd2*prooff_2;
			  dpp_2dt    = ca1*p_2*p_2-cd1*pp_2-ca2*pp_2*pro_1-gamma4*pp_2+cd2*prooff_1;
			  dppp_1dt    = 0;
			  dppp_2dt    = 0;
			  dpppp_1dt    = 0;
			  dpppp_2dt    = 0;

		      
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
		      fp = fopen("zoom-scan.csv", "a");
		      fprintf(fp, "%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", n,param,scale,ATP,gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3, ip1, ip2, p_1, p_2);
		      fclose(fp);
		    }
		}
	    }

	}
    }

    
  
  return 0;
}
