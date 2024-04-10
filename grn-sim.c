#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// parameters of model
typedef struct {
  double gamma2, gamma3, gamma4;
  double cd1, ca1, cd2, ca2;
  double lambda2, lambda3;

  // for bookkeeping
  double ATP;
  int param, param2, scale;
} Params;

// state variables
typedef struct {
  double z_1, z_2;

  double pro_1, pro_2;
  double prooff_1, prooff_2;
  double rna_1, rna_2;
  double p_1, p_2;
  double pp_1, pp_2;

  // used for higher-order models
  double ppp_1, ppp_2;
  double pppp_1, pppp_2;

  // used for alternative dimerisation model
  double propoff_1, propoff_2;
} State;

// default parameterisation
void defaultParams(Params *P)
{
  P->gamma2 = 0.00023;
  P->gamma3 = 0.00077;
  P->gamma4 = 0.00058;
  P->cd1 = 0.00023;
  P->cd2 = 1.04;
  P->ca1 = 0.0023;
  P->ca2 = 0.1038; 
  P->lambda2 = 0.0067;
  P->lambda3 = 0.0474;

  // for bookkeeping
  P->ATP = 1;
  P->param = P->param2 = -1;
  P->scale = -1;
}

// add two state vectors
State addState(State S1, State S2)
{
  State nS;
  nS.pro_1 = S1.pro_1+S2.pro_1;
  nS.pro_2 = S1.pro_2+S2.pro_2;
  nS.prooff_1 = S1.prooff_1+S2.prooff_1;
  nS.prooff_2 = S1.prooff_2+S2.prooff_2;
  nS.rna_1 = S1.rna_1+S2.rna_1;
  nS.rna_2 = S1.rna_2+S2.rna_2;
  nS.p_1 = S1.p_1+S2.p_1;
  nS.p_2 = S1.p_2+S2.p_2;
  nS.pp_1 = S1.pp_1+S2.pp_1;
  nS.pp_2 = S1.pp_2+S2.pp_2;
  nS.ppp_1 = S1.ppp_1+S2.ppp_1;
  nS.ppp_2 = S1.ppp_2+S2.ppp_2;
  nS.pppp_1 = S1.pppp_1+S2.pppp_1;
  nS.pppp_2 = S1.pppp_2+S2.pppp_2;
  nS.z_1 = S1.z_1+S2.z_1;
  nS.z_2 = S1.z_2+S2.z_2;
  
  return nS;
}

// scalar times state vector
State scaleState(State S1, double alpha)
{
  State nS;
  nS.pro_1 = S1.pro_1*alpha;
  nS.pro_2 = S1.pro_2*alpha;
  nS.prooff_1 = S1.prooff_1*alpha;
  nS.prooff_2 = S1.prooff_2*alpha;
  nS.rna_1 = S1.rna_1*alpha;
  nS.rna_2 = S1.rna_2*alpha;
  nS.p_1 = S1.p_1*alpha;
  nS.p_2 = S1.p_2*alpha;
  nS.pp_1 = S1.pp_1*alpha;
  nS.pp_2 = S1.pp_2*alpha;
  nS.ppp_1 = S1.ppp_1*alpha;
  nS.ppp_2 = S1.ppp_2*alpha;
  nS.pppp_1 = S1.pppp_1*alpha;
  nS.pppp_2 = S1.pppp_2*alpha;
  nS.z_1 = S1.z_1*alpha;
  nS.z_2 = S1.z_2*alpha;
  
  return nS;
}

// time derivatives for various model structures
void derivatives(int model, int n, Params P, State S, State *dS)
{
  dS->rna_1   = P.lambda2*S.pro_1*S.z_1-P.gamma2*S.rna_1;
  dS->rna_2   = P.lambda2*S.pro_2*S.z_2-P.gamma2*S.rna_2;
  dS->z_1 = dS->z_2 = 0;

  if(model == 1)
    {
      if(n == 1)
	{
	  //           translate     degrade    DNAbind       DNAunbind
	  dS->p_1     = P.lambda3*S.rna_1-P.gamma3*S.p_1-P.ca2*S.pro_2*S.p_1+P.cd2*S.prooff_2;
	  dS->p_2     = P.lambda3*S.rna_2-P.gamma3*S.p_2-P.ca2*S.pro_1*S.p_2+P.cd2*S.prooff_1;
	  //           DNAunbind    DNAbind
	  dS->pro_1   = P.cd2*S.prooff_1-P.ca2*S.pro_1*S.p_2;            
	  dS->pro_2   = P.cd2*S.prooff_2-P.ca2*S.pro_2*S.p_1;
	  //             DNAunbind    DNAbind
	  dS->prooff_1 = -P.cd2*S.prooff_1+P.ca2*S.pro_1*S.p_2;
	  dS->prooff_2 = -P.cd2*S.prooff_2+P.ca2*S.pro_2*S.p_1;
	  dS->pp_1    = 0;
	  dS->pp_2    = 0;
	  dS->ppp_1    = 0;
	  dS->ppp_2    = 0;
	  dS->pppp_1    = 0;
	  dS->pppp_2    = 0;
	}
		      
      if(n == 2)
	{
	  //           translate     dimer         dedimer    degrade
	  dS->p_1     = P.lambda3*S.rna_1-2*P.ca1*S.p_1*S.p_1+2*P.cd1*S.pp_1-P.gamma3*S.p_1;
	  dS->p_2     = P.lambda3*S.rna_2-2*P.ca1*S.p_2*S.p_2+2*P.cd1*S.pp_2-P.gamma3*S.p_2;
	  //           DNAunbind    DNAbind
	  dS->pro_1   = P.cd2*S.prooff_1-P.ca2*S.pro_1*S.pp_2;            
	  dS->pro_2   = P.cd2*S.prooff_2-P.ca2*S.pro_2*S.pp_1;
	  //             DNAunbind    DNAbind
	  dS->prooff_1 = -P.cd2*S.prooff_1+P.ca2*S.pro_1*S.pp_2;
	  dS->prooff_2 = -P.cd2*S.prooff_2+P.ca2*S.pro_2*S.pp_1;
	  //           dimer       dedimer  DNAbind        degrade     DNAunbind
	  dS->pp_1    = P.ca1*S.p_1*S.p_1-P.cd1*S.pp_1-P.ca2*S.pp_1*S.pro_2-P.gamma4*S.pp_1+P.cd2*S.prooff_2;
	  dS->pp_2    = P.ca1*S.p_2*S.p_2-P.cd1*S.pp_2-P.ca2*S.pp_2*S.pro_1-P.gamma4*S.pp_2+P.cd2*S.prooff_1;
	  dS->ppp_1    = 0;
	  dS->ppp_2    = 0;
	  dS->pppp_1    = 0;
	  dS->pppp_2    = 0;

	}
      if(n == 3)
	{
	  //           translate     dimer         dedimer    degrade    trimer       detrimer
	  dS->p_1     = P.lambda3*S.rna_1-2*P.ca1*S.p_1*S.p_1+2*P.cd1*S.pp_1-P.gamma3*S.p_1-P.ca1*S.p_1*S.pp_1+P.cd1*S.ppp_1;
	  dS->p_2     = P.lambda3*S.rna_2-2*P.ca1*S.p_2*S.p_2+2*P.cd1*S.pp_2-P.gamma3*S.p_2-P.ca1*S.p_2*S.pp_2+P.cd1*S.ppp_2;
	  //           DNAunbind    DNAbind
	  dS->pro_1   = P.cd2*S.prooff_1-P.ca2*S.pro_1*S.ppp_2;            
	  dS->pro_2   = P.cd2*S.prooff_2-P.ca2*S.pro_2*S.ppp_1;
	  //             DNAunbind    DNAbind
	  dS->prooff_1 = -P.cd2*S.prooff_1+P.ca2*S.pro_1*S.ppp_2;
	  dS->prooff_2 = -P.cd2*S.prooff_2+P.ca2*S.pro_2*S.ppp_1;
	  //           dimer       dedimer  degrade     trimer       detrimer
	  dS->pp_1    = P.ca1*S.p_1*S.p_1-P.cd1*S.pp_1-P.gamma4*S.pp_1-P.ca1*S.p_1*S.pp_1+P.cd1*S.ppp_1;
	  dS->pp_2    = P.ca1*S.p_2*S.p_2-P.cd1*S.pp_2-P.gamma4*S.pp_2-P.ca1*S.p_2*S.pp_2+P.cd1*S.ppp_2;
	  //            trimer       detrimer  DNAunbind       degrade      DNAbind      
	  dS->ppp_1    = P.ca1*S.pp_1*S.p_1-P.cd1*S.ppp_1-P.ca2*S.ppp_1*S.pro_2-P.gamma4*S.ppp_1+P.cd2*S.prooff_2;
	  dS->ppp_2    = P.ca1*S.pp_2*S.p_2-P.cd1*S.ppp_2-P.ca2*S.ppp_2*S.pro_1-P.gamma4*S.ppp_2+P.cd2*S.prooff_1;
	  dS->pppp_1    = 0;
	  dS->pppp_2    = 0;

	}
      if(n == 4)
	{
	  //           translate     dimer         dedimer    degrade    trimer       detrimer  tetramer      detetramer
	  dS->p_1     = P.lambda3*S.rna_1-2*P.ca1*S.p_1*S.p_1+2*P.cd1*S.pp_1-P.gamma3*S.p_1-P.ca1*S.p_1*S.pp_1+P.cd1*S.ppp_1-P.ca1*S.p_1*S.ppp_1+P.cd1*S.pppp_1;
	  dS->p_2     = P.lambda3*S.rna_2-2*P.ca1*S.p_2*S.p_2+2*P.cd1*S.pp_2-P.gamma3*S.p_2-P.ca1*S.p_2*S.pp_2+P.cd1*S.ppp_2-P.ca1*S.p_2*S.ppp_2+P.cd1*S.pppp_2;
	  //           DNAunbind    DNAbind
	  dS->pro_1   = P.cd2*S.prooff_1-P.ca2*S.pro_1*S.pppp_2;            
	  dS->pro_2   = P.cd2*S.prooff_2-P.ca2*S.pro_2*S.pppp_1;
	  //             DNAunbind    DNAbind
	  dS->prooff_1 = -P.cd2*S.prooff_1+P.ca2*S.pro_1*S.pppp_2;
	  dS->prooff_2 = -P.cd2*S.prooff_2+P.ca2*S.pro_2*S.pppp_1;
	  //           dimer       dedimer  degrade     trimer       detrimer
	  dS->pp_1    = P.ca1*S.p_1*S.p_1-P.cd1*S.pp_1-P.gamma4*S.pp_1-P.ca1*S.p_1*S.pp_1+P.cd1*S.ppp_1;
	  dS->pp_2    = P.ca1*S.p_2*S.p_2-P.cd1*S.pp_2-P.gamma4*S.pp_2-P.ca1*S.p_2*S.pp_2+P.cd1*S.ppp_2;
	  //            trimer       detrimer  degrade      tetramer      detetramer
	  dS->ppp_1    = P.ca1*S.pp_1*S.p_1-P.cd1*S.ppp_1-P.gamma4*S.ppp_1-P.ca1*S.ppp_1*S.p_1+P.cd1*S.pppp_2;
	  dS->ppp_2    = P.ca1*S.pp_2*S.p_2-P.cd1*S.ppp_2-P.gamma4*S.ppp_2-P.ca1*S.ppp_2*S.p_2+P.cd1*S.pppp_2;
	  //             tetramer      detetramer DNAbind          degrade       DNAunbind
	  dS->pppp_1    = P.ca1*S.ppp_1*S.p_1-P.cd1*S.pppp_1-P.ca2*S.pppp_1*S.pro_2-P.gamma4*S.pppp_1+P.cd2*S.prooff_2;
	  dS->pppp_2    = P.ca1*S.ppp_2*S.p_2-P.cd1*S.pppp_2-P.ca2*S.pppp_2*S.pro_1-P.gamma4*S.pppp_2+P.cd2*S.prooff_1;
	}
    }
  
  if(model == 2)
    {
      double lambda2p = P.lambda2/2;
      double ca1p = 10*P.ca1;
      double cd1p = P.cd1;
      double ca3 = P.ca2;
      double cd3 = P.cd2;

      dS->p_1     = P.lambda3*S.rna_1-2*P.ca1*S.p_1*S.p_1+2*P.cd1*S.pp_1-P.gamma3*S.p_1   + cd3*S.propoff_2-ca3*S.pro_2*S.p_1+cd1p*S.prooff_2-ca1p*S.propoff_2*S.p_1;
      dS->p_2     = P.lambda3*S.rna_2-2*P.ca1*S.p_2*S.p_2+2*P.cd1*S.pp_2-P.gamma3*S.p_2   + cd3*S.propoff_1-ca3*S.pro_1*S.p_2+cd1p*S.prooff_1-ca1p*S.propoff_1*S.p_2;
      //           DNAunbind    DNAbind           DNA unbind m  DNA bind m
      dS->pro_1   = P.cd2*S.prooff_1-P.ca2*S.pro_1*S.pp_2  + cd3*S.propoff_1-ca3*S.pro_1*S.p_2;            
      dS->pro_2   = P.cd2*S.prooff_2-P.ca2*S.pro_2*S.pp_1  + cd3*S.propoff_2-ca3*S.pro_2*S.p_1;
      //             DNAunbind    DNAbind           DNA dimer          DNA dedimer
      dS->prooff_1 = -P.cd2*S.prooff_1+P.ca2*S.pro_1*S.pp_2  + ca1p*S.propoff_1*S.p_2-cd1p*S.prooff_1;
      dS->prooff_2 = -P.cd2*S.prooff_2+P.ca2*S.pro_2*S.pp_1  + ca1p*S.propoff_2*S.p_1-cd1p*S.prooff_2;
      //             DNA bind m    DNA unbind m  DNA dedimer   DNA dimer
      dS->propoff_1 = ca3*S.pro_1*S.p_2-cd3*S.propoff_1+cd1p*S.prooff_1-ca1p*S.propoff_1*S.p_2;
      dS->propoff_2 = ca3*S.pro_2*S.p_1-cd3*S.propoff_2+cd1p*S.prooff_2-ca1p*S.propoff_2*S.p_1;			      
      //           dimer       dedimer  DNAbind        degrade     DNAunbind
      dS->pp_1    = P.ca1*S.p_1*S.p_1-P.cd1*S.pp_1-P.ca2*S.pp_1*S.pro_2-P.gamma4*S.pp_1+P.cd2*S.prooff_2;
      dS->pp_2    = P.ca1*S.p_2*S.p_2-P.cd1*S.pp_2-P.ca2*S.pp_2*S.pro_1-P.gamma4*S.pp_2+P.cd2*S.prooff_1;
      dS->ppp_1    = 0;
      dS->ppp_2    = 0;
      dS->pppp_1    = 0;
      dS->pppp_2    = 0;
    }
  
}

void initialState(State *S, double ip1, double ip2, int ICs)
{
  S->p_1 = ip1; S->p_2 = ip2;
  S->pro_1 = 1; S->pro_2 = 1; S->z_1 = 1; S->z_2 = 1;
  S->rna_1 = S->rna_2 = S->prooff_1 = S->prooff_2 = S->pp_1 = S->pp_2 = S->ppp_1 = S->ppp_2 = S->pppp_1 = S->pppp_2 = 0;

  switch(ICs) {
  case 0: break;
  case 1: S->rna_1 = 10; break;
  case 2: S->rna_2 = 10; break;
  case 3: S->rna_1 = S->rna_2 = 10; break;
  }
    
}
       
// test for convergence in terms of proportional state change
int convergence(State S, State dS, double eps)
{
  if(!(dS.pro_1 > eps*S.pro_1 || dS.pro_2 > eps*S.pro_2 || dS.prooff_1 > eps*S.prooff_1 || dS.prooff_2 > eps*S.prooff_2 || dS.rna_1 > eps*S.rna_1 || dS.rna_2 > eps*S.rna_2 || dS.p_1 > eps*S.p_1 || dS.p_2 > eps*S.p_2 || dS.pp_1 > eps*S.pp_1 || dS.pp_2 > eps*S.pp_2 || dS.ppp_1 > eps*S.ppp_1 || dS.ppp_2 > eps*S.ppp_2 || dS.pppp_1 > eps*S.pppp_1 || dS.pppp_2 > eps*S.pppp_2))
    return 0;
  else return 1;
}

void simulate(Params P, int model, int euler, int n, int ICs, char *filename, char *tfilename, int TIMEOUT, double dt)
{
  double ip1, ip2;
  State S, k1, k2, k3, k4, dS;
  int converge;
  FILE *fp;
  double t;
  double eps = 1e-6;
  
  // loop through initial protein concentrations
  for(ip1 = 0; ip1 < 30; ip1 += 3)
    {
      for(ip2 = 0; ip2 < 30; ip2 += 3)
	{
	  printf("%i,%i,%i,%i,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", model,n,euler,ICs,P.ATP,P.param,P.param2,P.scale,P.gamma2,P.gamma3,P.gamma4,P.cd1,P.cd2,P.ca1,P.ca2,P.lambda2,P.lambda3, ip1, ip2);

	  // initialise system
	  converge = 0;
	  initialState(&S, ip1, ip2, ICs);

	  for(t = 0; converge < 500./dt; t += dt)
	    {
	      if(!euler)
		{
		  // RK4 solver
		  derivatives(model, n, P, S, &k1);
		  derivatives(model, n, P, addState(S, scaleState(k1, dt/2.)), &k2);
		  derivatives(model, n, P, addState(S, scaleState(k2, dt/2.)), &k3);
		  derivatives(model, n, P, addState(S, scaleState(k3, dt)), &k4);

		  dS = scaleState( addState(k1, addState(scaleState(k2, 2), addState(scaleState(k3, 2), k4))), dt/6.);
		}
	      else
		{
		  // simple euler
		  derivatives(model, n, P, S, &k1);
		  dS = scaleState(k1, dt);
		}
	      
	      if(!convergence(S, dS, eps*dt))
		converge++;
	      else converge = 0;
		      
	      S = addState(S, dS);
		    
	      if(TIMEOUT && (int)(0.01*t) > (int)(0.01*(t-dt)))
		{
		  fp = fopen(tfilename, "a");
		  fprintf(fp, "%i,%i,%i,%i,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", model,n,euler,ICs,P.ATP,P.param,P.param2,P.scale,P.gamma2,P.gamma3,P.gamma4,P.cd1,P.cd2,P.ca1,P.ca2,P.lambda2,P.lambda3, ip1, ip2, t, S.p_1, S.p_2);
		  fclose(fp);

		}
	    }
		      
	  // output start and end points for this sim
	  fp = fopen(filename, "a");
	  fprintf(fp, "%i,%i,%i,%i,%e,%i,%i,%i,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", model,n,euler,ICs,P.ATP,P.param,P.param2,P.scale,P.gamma2,P.gamma3,P.gamma4,P.cd1,P.cd2,P.ca1,P.ca2,P.lambda2,P.lambda3, ip1, ip2, S.p_1, S.p_2);
	  fclose(fp);
	}
    }
}

void initialiseFile(char *fname, int t)
{
  FILE *fp;
  
  fp = fopen(fname, "w");
  if(t == 0)
    fprintf(fp, "model,n,euler,ICs,ATP,param,param2,scale,gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3,ip1,ip2,p1,p2\n");
  else
    fprintf(fp, "model,n,euler,ICs,ATP,param,param2,scale,gamma2,gamma3,gamma4,cd1,cd2,ca1,ca2,lambda2,lambda3,ip1,ip2,t,p1,p2\n");
  fclose(fp);
}

int main(int argc, char *argv[])
{
  // parameters
  Params P;
  
  double dt = 0.01;
  int n;
  int expt;
  int ICs;
  char fname[200];
  
  if(argc < 2)
    {
      printf("Which experiment to run?\n  0 -- Euler vs RK4\n  1 -- param sweep with default model\n  2 -- param sweep with on-DNA dimerisation\n  3 -- ATP x gamma3 influence\n  4 -- ATP x cd2 influence\n  5 -- different ICs\n  6 -- different ATP influences\n  7 -- bifurcation plots\n\n");
      exit(0);
    }
  expt = atoi(argv[1]);
  
  // void simulate(Params P, int model, int euler, int n, int ICs, char *filename, char *tfilename, int TIMEOUT, double dt)

  // default setup
  defaultParams(&P);
  n = 2;
  
  // show that RK4 matches Euler
  if(expt == 0)
    {
      initialiseFile("grn-sim-0.csv", 0);
      initialiseFile("grn-sim-0-t.csv", 1);
     
      simulate(P, 1, 0, n, 0, "grn-sim-0.csv", "grn-sim-0-t.csv", 1, dt);
      simulate(P, 1, 1, n, 0, "grn-sim-0.csv", "grn-sim-0-t.csv", 1, dt);
    }

  // loop through param changes for default model (1) or on-DNA dimerisation (2)
  if(expt == 1 || expt == 2)
    {
      if(argc < 3)
	{
	  printf("For parameter sweeps I need another argument 1-4 specifying the cooperativity structure\n");
	  return 0;
	}

      n = atoi(argv[2]);
      
      sprintf(fname, "grn-sim-%i.%i.csv", expt, n);
      initialiseFile(fname, 0);
  
      for(P.param = 0; P.param < 9; P.param++)
	{
	  switch(P.param)
	    {
	    case 0: P.gamma2 /= 1000; break;
	    case 1: P.gamma3 /= 1000; break;
	    case 2: P.gamma4 /= 1000; break;
	    case 3: P.cd1 /= 1000; break;
	    case 4: P.cd2 /= 1000; break;
	    case 5: P.ca1 /= 1000; break;
	    case 6: P.ca2 /= 1000; break;
	    case 7: P.lambda2 /= 1000; break;
	    case 8: P.lambda3 /= 1000; break;
	    }
      
	  for(P.scale = 0; P.scale <= 4; P.scale++)
	    {
	      // 0 /100; 1 /10; 2 1; 3 *10; 4 *100
	      switch(P.param)
		{
		case 0: P.gamma2 *= 10; break;
		case 1: P.gamma3 *= 10; break;
		case 2: P.gamma4 *= 10; break;
		case 3: P.cd1 *= 10; break;
		case 4: P.cd2 *= 10; break;
		case 5: P.ca1 *= 10; break;
		case 6: P.ca2 *= 10; break;
		case 7: P.lambda2 *= 10; break;
		case 8: P.lambda3 *= 10; break;
		}
		  
	      // the particular case of 100x ca2 needs a smaller timestep for stability
	      if(P.param == 6 && P.scale == 4) dt = 0.0001;
	      else dt = 0.005;

	      if(expt == 1) simulate(P, 1, 1, n, 0, fname, "", 0, dt);
	      else simulate(P, 2, 1, n, 0, fname, "", 0, dt);
			  
	    }
	  switch(P.param)
	    {
	    case 0: P.gamma2 /= 100; break;
	    case 1: P.gamma3 /= 100; break;
	    case 2: P.gamma4 /= 100; break;
	    case 3: P.cd1 /= 100; break;
	    case 4: P.cd2 /= 100; break;
	    case 5: P.ca1 /= 100; break;
	    case 6: P.ca2 /= 100; break;
	    case 7: P.lambda2 /= 100; break;
	    case 8: P.lambda3 /= 100; break;
	    }

	}
    }
    

  // ATP influence through gamma3 (3) and cd2 (4)
  if(expt == 3 || expt == 4)
    {
      // open file for output
      if(expt == 3) initialiseFile("grn-sim-3.csv", 0);
      else initialiseFile("grn-sim-4.csv", 0);
  
      n = 2;

      if(expt == 3) P.gamma3 = 0.00077/4;
      if(expt == 4) P.cd2 = 0.005/2;
      for(P.scale = 0; P.scale <= 4; P.scale++)
	{
	  if(expt == 3) P.gamma3 *= 4;
	  if(expt == 4) P.cd2 *= 2;
	  for(P.ATP = 0.5; P.ATP <= 2.1; P.ATP *= 2)
	    {
	      P.lambda2 *= P.ATP; P.lambda3 *= P.ATP;
	      if(expt == 3) simulate(P, 1, 1, n, 0, "grn-sim-3.csv", "", 0, dt);
	      if(expt == 4) simulate(P, 1, 1, n, 0, "grn-sim-4.csv", "", 0, dt);
	      P.lambda2 /= P.ATP; P.lambda3 /= P.ATP;
	    }
	}
    }

  // different initial conditions
  if(expt == 5)
    {
      // open file for output
      initialiseFile("grn-sim-5.csv", 0);
  
      n = 2;

      for(ICs = 0; ICs <= 4; ICs++)
	{
	  for(P.ATP = 0.5; P.ATP <= 2.1; P.ATP *= 2)
	    {
	      P.lambda2 *= P.ATP; P.lambda3 *= P.ATP;
	      simulate(P, 1, 1, n, ICs, "grn-sim-5.csv", "", 0, 0.01);
	      P.lambda2 /= P.ATP; P.lambda3 /= P.ATP;
	    }
	}
    }

  // different modes of ATP influence
  if(expt == 6)
    {
      // open file for output
      initialiseFile("grn-sim-6a.csv", 0);
      initialiseFile("grn-sim-6b.csv", 0);
      initialiseFile("grn-sim-6c.csv", 0);
  
      n = 2;

      P.gamma4 /= 4;
      for(P.scale = 0; P.scale <= 4; P.scale++)
	{
	  P.gamma3 *= 4;

      	  for(P.ATP = 0.5; P.ATP <= 2.1; P.ATP *= 2)
	    {
	      P.lambda2 *= P.ATP;
	      simulate(P, 1, 1, n, 0, "grn-sim-6a.csv", "", 0, dt);
	      P.lambda2 /= P.ATP;

	      P.lambda3 *= P.ATP;
	      simulate(P, 1, 1, n, 0, "grn-sim-6b.csv", "", 0, dt);
	      P.lambda3 /= P.ATP;

     	      P.lambda2 *= P.ATP; P.lambda3 *= P.ATP; P.gamma2 *= P.ATP; P.gamma3 *= P.ATP;
	      simulate(P, 1, 1, n, 0, "grn-sim-6c.csv", "", 0, dt);
     	      P.lambda2 /= P.ATP; P.lambda3 /= P.ATP; P.gamma2 /= P.ATP; P.gamma3 /= P.ATP;
	    }
	}
    }

  // bifurcation plots
  if(expt == 7)
    {
      initialiseFile("grn-sim-7.csv", 0);
  
      n = 2;

      for(P.ATP = 0.1; P.ATP <= 2; P.ATP += 0.1)
	{
	  P.lambda2 *= P.ATP; P.lambda3 *= P.ATP;
	  P.gamma3 *= 4;
	  simulate(P, 1, 1, n, 0, "grn-sim-7.csv", "", 0, dt);
	  P.gamma3 /= 4;
	  P.cd2 /= 4;
	  simulate(P, 1, 1, n, 0, "grn-sim-7.csv", "", 0, dt);
	  P.cd2 *= 4;
	  P.lambda2 /= P.ATP; P.lambda3 /= P.ATP;
	}
    }

  // matrix of pairwise parameter changes
  if(expt == 8)
    {
      //      initialiseFile("grn-sim-8.csv", 0);
	    
      if(argc < 3)
	{
	  printf("For pairwise param changes, I need another argument 0-7 specifying the scale of changes to apply\n");
	  return 0;
	}
      P.scale = atoi(argv[2]);
      double s1, s2;
      
      for(P.param = 0; P.param < 9; P.param++)
	{
	  for(P.param2 = P.param+1; P.param2 < 9; P.param2++)
	    {
	      switch(P.scale)
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
		
	      switch(P.param)
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

	      switch(P.param2)
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

	      simulate(P, 1, 1, n, 0, "grn-sim-8.csv", "", 0, dt);
	      switch(P.param)
		{
		case 0: P.gamma2 /= s1; break;
		case 1: P.gamma3 /= s1; break;
		case 2: P.gamma4 /= s1; break;
		case 3: P.cd1 /= s1; break;
		case 4: P.cd2 /= s1; break;
		case 5: P.ca1 /= s1; break;
		case 6: P.ca2 /= s1; break;
		case 7: P.lambda2 /= s1; break;
		case 8: P.lambda3 /= s1; break;
		}

	      switch(P.param2)
		{
		case 0: P.gamma2 /= s2; break;
		case 1: P.gamma3 /= s2; break;
		case 2: P.gamma4 /= s2; break;
		case 3: P.cd1 /= s2; break;
		case 4: P.cd2 /= s2; break;
		case 5: P.ca1 /= s2; break;
		case 6: P.ca2 /= s2; break;
		case 7: P.lambda2 /= s2; break;
		case 8: P.lambda3 /= s2; break;
		}
	    }
	}
    }
      
  return 0;
}
