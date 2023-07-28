#include <stdio.h>

int main(void)
{
  double x[4], dx[4];
  double b = 0.05, u = 1;
  double t, dt = 0.001;
  FILE *fp;
  double x0;
  double p0, p1, dp0, dp1;
  double lambda = 0.001;
  int n;
  double p11, p12, p14;
  
  fp = fopen("hill.csv", "w");
  fprintf(fp, "x0,p11,p12,p14\n");
  
  // x[n] is concentration of (n+1)-mer
  // p0 is inactive proportion of DNA; p1 active proportion

  // b is monomer binding; u unbinding
  // lambda is deactivation rate?
  
  // loop through initial monomer concentration
  for(x0 = 0; x0 < 15; x0 += 0.05)
    {
      //// monomer x[0] is active regulator
      // p0 + x0 -> p1 , p1 -> p0 + x0
      // initial conditions
      x[0] = x0; 
      p0 = 1; p1 = 0;
      for(t = 0; t < 100; t += dt)
	{
	  dx[0] = -p0*x[0] + p1*lambda;
	  dp1 = p0*x[0] - p1*lambda;
	  dp0 = p1*lambda - p0*x[0];
      
	  x[0] += dx[0]*dt;
	  p0 += dp0*dt;
	  p1 += dp1*dt;
	}
      p11 = p1;

      //// dimer x[1] is active regulator
      x[0] = x0; x[1] = 0;
      p0 = 1; p1 = 0;
      for(t = 0; t < 100; t += dt)
	{
	  dx[0] = -2*b*x[0]*x[0] + 2*u*x[1];
	  dx[1] = b*x[0]*x[0] - u*x[1] - p0*x[1] + p1*lambda;
	  dp1 = p0*x[1] - p1*lambda;
	  dp0 = p1*lambda - p0*x[1];

	  x[0] += dx[0]*dt;
	  x[1] += dx[1]*dt;
	  p0 += dp0*dt;
	  p1 += dp1*dt;
	}
      p12 = p1;

      ///// tetramer x[3] is active regulator
      x[0] = x0; x[1] = x[2] = x[3] = 0;
      p0 = 1; p1 = 0;
      for(t = 0; t < 100; t += dt)
	{
	  dx[0] = -b*x[0]*(2*x[0]+x[1]+x[2]) + u*(2*x[1]+x[2]+x[3]);
	  dx[1] = b*(x[0]*x[0] - x[0]*x[1]) + u*(x[2]-x[1]);
	  dx[2] = b*(x[1]*x[0] - x[0]*x[2]) + u*(x[3]-x[2]);
	  dx[3] = b*(x[2]*x[0]) - u*x[3] - p0*x[3] + p1*lambda;
	  dp1 = p0*x[3] - p1*lambda;
	  dp0 = p1*lambda - p0*x[3];
      
	  x[0] += dx[0]*dt;
	  x[1] += dx[1]*dt;
	  x[2] += dx[2]*dt;
	  x[3] += dx[3]*dt;
	  p0 += dp0*dt;
	  p1 += dp1*dt;
	}
      p14 = p1;
      
      fprintf(fp, "%e,%e,%e,%e\n", x0, p11, p12, p14);
      //      fprintf(fp, "%e %e %e %e %e %e\n", x0, x[3], p10, p12, p14, dx[3]);
    }
  
  return 0;
}
