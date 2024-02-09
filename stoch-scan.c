#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()
#define MAXT 1e8      // max time for simulation
#define MAXSWITCH 100000

int main(int argc, char *argv[])
{
  // parameters
  double gamma2, gamma3, gamma4;
  double cd1, ca1, cd2, ca2;
  double lambda2, lambda3;

  // state variables
  int pro_1, pro_2;
  int prooff_1, prooff_2;
  int rna_1, rna_2;
  int p_1, p_2;
  int pp_1, pp_2;
  int ppp_1, ppp_2;
  int pppp_1, pppp_2;
  int z_1, z_2;

  FILE *fp, *fp1;
  int i, reaction;
  
  double rates[18], cumsum[18];
  double r;
  double t, lastt;
  int ip1, ip2;
  double tau;
  int expt;
  char fstr[100];

  double *switchlength;
  int nswitch;
  double switcht;
  int currstate;
  int param, scale;

  if(argc != 2)
    {
      printf("Which parameter to vary? (0-8)\n");
      exit(0);
    }

  param = atoi(argv[1]);
  
  switchlength = (double*)malloc(sizeof(double)*MAXSWITCH);

  gamma2 = 0.00023;
  gamma3 = 0.00077;
  gamma4 = 0.00058;
  cd1 = 0.00023;
  cd2 = 1.04;
  ca1 = 0.0023;
  ca2 = 0.1038; 
  lambda2 = 0.0067;
  lambda3 = 0.0474;
  

  switch(param)
    {
    case 0: gamma2 /= 1000; break;
    case 1: gamma3 /= 1000; break;
    case 2: gamma4 /= 1000; break;
    case 3: cd1 /= 1000; break;
    case 4: cd2 /= 1000; break;
    case 5: ca1 /= 1000; break;
    case 6: ca2 /= 1000; break;
    case 7: lambda2 /= 1000; break;
    case 8: lambda3 /= 1000; break;
    }
      
  for(scale = 0; scale <= 4; scale++)
    {
      // 0 /100; 1 /10; 2 1; 3 *10; 4 *100
      switch(param)
	{
	case 0: gamma2 *= 10; break;
	case 1: gamma3 *= 10; break;
	case 2: gamma4 *= 10; break;
	case 3: cd1 *= 10; break;
	case 4: cd2 *= 10; break;
	case 5: ca1 *= 10; break;
	case 6: ca2 *= 10; break;
	case 7: lambda2 *= 10; break;
	case 8: lambda3 *= 10; break;
	}

      expt = 10*param+scale;
	      
      fp = fopen("gillespie-scan-params.csv", "a");
      fprintf(fp, "%i,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", expt, gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3);
      printf("%i,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", expt, gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3);
      fclose(fp);

      
      sprintf(fstr, "gillespie-scan-series-%i.csv", expt);
      fp = fopen(fstr, "w");
      fprintf(fp, "t, pro_1, pro_2, prooff_1, prooff_2, rna_1, rna_2, p_1, p_2, pp_1, pp_2\n");

      sprintf(fstr, "gillespie-scan-switches-%i.csv", expt);
      fp1 = fopen(fstr, "w");
      fprintf(fp1, "i, dt\n");

      ip1 = ip2 = 0;
  
      // initialise system
      p_1 = ip1; p_2 = ip2;
      pro_1 = 1; pro_2 = 1; z_1 = 1; z_2 = 1;
      rna_1 = rna_2 = prooff_1 = prooff_2 = pp_1 = pp_2 = ppp_1 = ppp_2 = pppp_1 = pppp_2 = 0;
      tau = 0;

      nswitch = 0;
      for(i = 0; i < MAXSWITCH; i++)
	switchlength[i] = 0;
      
      // run simulation until time MAXT
      for(t = 0; t < MAXT; )
	{
	  // compute rates for each reaction
	  rates[0] = lambda3*rna_1;
	  rates[1] = ca1*p_1*(p_1-1);
	  rates[2] = cd1*pp_1;
	  rates[3] = gamma3*p_1;

	  rates[4] = lambda3*rna_2;
	  rates[5] = ca1*p_2*(p_2-1);
	  rates[6] = cd1*pp_2;
	  rates[7] = gamma3*p_2;

	  rates[8] = cd2*prooff_1;
	  rates[9] = ca2*pro_1*pp_2;
	  rates[10] = cd2*prooff_2;
	  rates[11] = ca2*pro_2*pp_1;

	  rates[12] = gamma4*pp_1;
	  rates[13] = gamma4*pp_2;

	  rates[14] = lambda2*pro_1*z_1;
	  rates[15] = gamma2*rna_1;
	  rates[16] = lambda2*pro_2*z_2;
	  rates[17] = gamma2*rna_2;

	  // build cumulative rate vector
	  for(i = 0; i <= 17; i++)
	    cumsum[i] = (i == 0 ? rates[0] : cumsum[i-1]+rates[i]);

	  // pick next reaction
	  r = RND*cumsum[17];
	  //      printf("%f %f\n", r, cumsum[17]);
	  for(reaction = 0; cumsum[reaction] <= r; reaction++);

	  // apply stoichiometric change for chosen reaction
	  switch(reaction) {
	  case 0: p_1++; break;
	  case 1: p_1 -= 2; pp_1++; break;
	  case 2: p_1 += 2; pp_1--; break;
	  case 3: p_1--; break;
	  case 4: p_2++; break;
	  case 5: p_2 -= 2; pp_2++; break;
	  case 6: p_2 += 2; pp_2--; break;
	  case 7: p_2--; break;
	  case 8: pro_1++; pp_2++; prooff_1--; break;
	  case 9: pro_1--; pp_2--; prooff_1++; break;
	  case 10: pro_2++; pp_1++; prooff_2--; break;
	  case 11: pro_2--; pp_1--; prooff_2++; break;
	  case 12: pp_1--; break;
	  case 13: pp_2--; break;
	  case 14: rna_1++; break;
	  case 15: rna_1--; break;
	  case 16: rna_2++; break;
	  case 17: rna_2--; break;
	  default: printf("bad reaction %i\n", reaction); exit(0);
	  }
      
	  if(p_1 < 0 || pp_1 < 0 || rna_1 < 0 || prooff_1 < 0 || p_2 < 0 || pp_2 < 0 || rna_2 < 0 || prooff_2 < 0)
	    { printf("negative value after reaction %i\n", reaction); exit(0); }

	  // output system state periodically
	  if((int)t % 10000 == 0 && (int)(t-tau) % 10000 != 0)
	    {
	      fprintf(fp, "%.3e,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n", t, pro_1, pro_2, prooff_1, prooff_2, rna_1, rna_2, p_1, p_2, pp_1, pp_2);
	    }

	  if((int)t % 1000000 == 0 && (int)(t-tau) % 1000000 != 0)
	    {
	      printf("%i\n", (int)t);
	    }

	  // if we're burnt in, initialise our counter for state switches
	  if(t+tau > MAXT/10 && t < MAXT/10)
	    {
	      currstate = (p_1 > p_2 ? 1 : 2);
	      switcht = t;
	    }

	  // catch state switches and record time since last one
	  if(t > MAXT/10 && p_1 > p_2 && currstate == 2 && nswitch < MAXSWITCH)
	    {
	      currstate = 1;
	      switchlength[nswitch++] = t-switcht;
	      switcht = t;
	    }
	  if(t > MAXT/10 && p_2 > p_1 && currstate == 1 && nswitch < MAXSWITCH)
	    {
	      currstate = 2;
	      switchlength[nswitch++] = t-switcht;
	      switcht = t;
	    }

	  // sample time increment and apply
	  lastt = t;
	  tau = 1./cumsum[17]*log(1./RND);
	  t += tau;      
	}
      for(i = 0; i < nswitch; i++)
	fprintf(fp1, "%i, %e\n", i, switchlength[i]);

      fclose(fp);
      fclose(fp1);
    }
      
  return 0;
}
