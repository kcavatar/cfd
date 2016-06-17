#include<iostream>
#include<math.h>
#include<cmath>
#include<algorithm>
#define N 5  //Number of species
#define tol 0.00001 //
#define R 8.315   //Ideal gas constant

using namespace std;
void eqmCons(double Kp[], float T);
void swapRow(double a[N][N+1], int con, int j);
void delCalc(double x[], double del[], double Kp[], float P, float Z);
void gibbs(float a[N][7], float G[], float T);

int main()
{
 //Rxns are : 1. O2 -> 2O    2. O2 + N2 -> 2NO  3. N2 -> 2N
 //Species are : 1(0) -> O2  2(1) -> N2  3(2) -> NO  4(3) -> O  5(4) -> N
 //Input initial parameters
  cout<<"MOLE FRACTION-BASED EQUILIBRIUM COMPOSITION SOLVER"<<endl<<endl;
  cout<<"Input the  values of temperature  (T), pressure (P) and N-O ratio (Z): \n\n";
  float T, P, Z;
  cin>>T>>P>>Z;
 //Count variables 
  int i=0, j=0;
 //Calculate equilibrium constants via external functions
  double Kp[3];
  eqmCons(Kp, T);
  for(i=0;i<3;i++)
   cout<<"\nThe value of Kp["<<i<<"] : "<<Kp[i]<<endl;
 //Conserving nitrogen and oxygen atoms using N/O ratio and sum of mole fractions
 // double x[N] = {0.5, 0.5, 0.0, 0.0, 0.0};       // Initial guess for unscaled coeff matrix
  double x[N] = {0.21, 0.79, 0.0, 0.0, 0.0};        //  Initial guess for scaled coeff. matrix by 100
  double del[N];
  
  delCalc(x, del, Kp, P, Z);
   
 //solve using Newton-Raphson's multivariate method
 double xNew[N];
 while(fabs(*max_element(del, del+5))>tol)
 {
   for(i=0;i<N;i++)
   {
     xNew[i] = x[i] + del[i];
     x[i] = xNew[i];
   }
   delCalc(x, del, Kp, P, Z);
 }
 
 //Display results obtained as such if required
 cout<<"\n\nThe results are as : \n";
 for(i=0;i<N;i++)
   cout<<"Species "<<i+1<<" : "<<x[i]<<endl;
  return true;
}

void eqmCons(double Kp[], float T)
{
                        /* COEFFICIENT TABLE FOR VARIOUS TEMPERATURE RANGES */	
  //Temperature : 200 - 800 K
  float a[N][7] = {{0.377037e1, -0.289522e-2, 0.953322e-5, -0.924699e-8, 0.301919e-11, -0.188598e2, 0.369335e1},
                   {0.346226e1,  0.582024e-3,-0.305255e-5,  0.622801e-8,-0.337560e-11, 0.879517e0,  0.321926e1},
                   {0.420644e1, -0.450984e-2, 0.105574e-4, -0.859194e-8, 0.240471e-11, 0.108890e5, 0.231379e1},
                   {0.321671e1, -0.378227e-2, 0.847468e-5, -0.886582e-8, 0.353656e-11, 0.296405e5, 0.185264e1},
                   {0.250000e1,          0.0,         0.0,          0.0,          0.0, 0.566267e5, 0.418073e1}
                  };
  //Temperature : 800 K - 3000 K
  float b[N][7] = {{0.289692e1, 0.237365e-2, -0.149171e-5, 0.466034e-9, -0.539452e-13, 0.822404e2, 0.750194e1},
                   {0.270224e1, 0.194439e-2, -0.893000e-6, 0.197392e-9, -0.169678e-13, 0.201802e3, 0.720408e1},
                   {0.275438e1, 0.230933e-2, -0.128234e-5, 0.340435e-9, -0.348075e-13, 0.111343e5, 0.907897e1},
                   {0.260454e1,-0.172355e-3, 0.115741e-6,-0.364179e-10, 0.460115e-14, 0.297290e5, 0.458653e1},
				   {0.250751e1,-0.247979e-4,  0.296415e-7,-0.152881e-10, 0.289137e-14, 0.566249e5, 0.414319e1}  
                  };
  //Temperature : 3000K - 6000K
  float c[N][7] = {{0.284211e1, 0.133206e-2, -0.339159e-6, 0.446522e-10, -0.229148e-14, 0.583504e3, 0.862550e1},
                   {0.391435e1, 0.315371e-3, -0.564810e-7, 0.360123e-11, 0.523594e-16, -0.535520e3, 0.216720e-1},
  	               {0.380154e1, 0.498575e-3, -0.125313e-6, 0.140939e-10, -0.448202e-15, 0.106666e5, 0.316194e1},
  	               {0.281017e1,-0.290399e-3,  0.908333e-7,-0.994278e-11, 0.377044e-15,  0.295366e5, 0.325365e1},
  	               {0.263761e1,-0.873733e-5, -0.647727e-7, 0.234734e-10, -0.173962e-14, 0.564523e5, 0.322321e1} 
                  };
  //Temperature : 6000K - 10000K
  float d[N][7] = {{0.568211e1, -0.753006e-3, 0.229801e-6, -0.239559e-10, 0.800485e-15, -0.247004e4, -0.987401e1},
  	               {0.126575e1,  0.182698e-2,-0.375839e-6, 0.330337e-10, -0.936511e-15,  0.314270e4,  0.179433e2 },
                   {0.491332e1, 0.617553e-6, -0.552224e-7, 0.116865e-10, -0.546425e-15,  0.884538e4, -0.457869e1},
			       {0.192093e1, 0.217766e-3, -0.182884e-7, 0.490510e-12,  0.285073e-17,  0.307834e5,  0.927486e1 },	    
				   {0.337206e1, -0.885546e-3, 0.252933e-6, -0.231879e-10, 0.704714e-15,  0.562702e5, -0.105640e1 }
				  }; 
  //Temperature : 10000K - 15000K
                
  float Gibbs[N];
  if(T>=200&&T<800)
    gibbs(a, Gibbs, T);
  if(T>=800&&T<3000)
    gibbs(b, Gibbs, T);
  if(T>=3000&&T<6000)
    gibbs(c, Gibbs, T);
  if(T>=6000&&T<=10000)
    gibbs(d, Gibbs, T);
     
  //Crude Kp; Kp acting as Gibbs free energy change for the reaction
  Kp[0] = 2*Gibbs[3] - Gibbs[0]; //Reaction : O2 -> 2O : Gibbs for O2 and N2 = 0
  Kp[1] = 2*Gibbs[2] - Gibbs[0] - Gibbs[1]; //Reaction : O2 + N2 -> 2NO :
  Kp[2] = 2*Gibbs[4] - Gibbs[1]; //Reaction : N2 -> 2N
  //Redefining Kp as Kp = exp(-delG/RT)
  for(int i=0;i<3;i++)
  {
    Kp[i] = exp(-Kp[i]/(R*T));
  }
}
void gibbs(float a[][7], float Gibbs[], float T)
{
  float H[N], S[N]; //defined as : H = H/RT, S = S/R; Gibbs = Gibbs/RT
  //Calculating the values of H, S and Gibbs
  /* H[i]/RT = Summ:j=0->4:(a[i][j]*T^j/(j+1)) + a[i][5]/T
     S[i]/R = Summ:j=1->4:(a[i][j]*T^j/j) + a[i][5]*log(T)  +  a[i][6]
     Gibbs[i]/RT = H[i]/RT - S[i]/R
  */
  for(int i=0;i<N;i++)
  {
    H[i] = 0; S[i] = 0;
    for(int j=0;j<=4;j++)
     H[i] = H[i] + a[i][j]*pow(T,j)/(j+1);

    H[i]=H[i]+a[i][5]/T;
    for(int j=1;j<=4;j++)
     S[i] = S[i] + a[i][j]*pow(T,j)/j;
    S[i] = S[i] + a[i][0]*log(T) + a[i][6];

    Gibbs[i] = H[i] - S[i];
    //Re-defining Gibbs
    Gibbs[i] = Gibbs[i]*R*T;
  }
}

void delCalc(double x[], double del[], double Kp[], float P, float Z)
{
 //Po = 1 atmosphere
  //Forming Jacobian Matrix

 double a[N][N+1] = {{ (1/P)*Kp[0],         0.0,     0.0, -2*x[3],     0.0,                    x[3]*x[3] - x[0]*(1/P)*Kp[0]},
                    {  x[1]*Kp[1],  x[0]*Kp[1], -2*x[2],     0.0,     0.0,                     x[2]*x[2] - x[0]*x[1]*Kp[1]},
                    {         0.0, (1/P)*Kp[2],     0.0,     0.0, -2*x[4],                    x[4]*x[4] - x[1]*(1/P)*Kp[2]},
                    {       2.0*Z,        -2.0,     Z-1,       Z,    -1.0, -2*Z*x[0] + 2*x[1] - (Z-1)*x[2] - Z*x[3] + x[4]},
                    {         1.0,         1.0,     1.0,     1.0,     1.0,            1 - x[0] - x[1] - x[2] - x[3] - x[4]}
                   };

 //GAUSS ELIMINATION Implementation
   int i,j,con;
   float alpha;
  //Generating upper triangular matrix
  //con controls the elimination process
  for(con=0;con<N;con++)
  {
  	if(a[con][con]!=0.0)
  	{
  	  for(i=con;i<N-1;i++)
  	  {
  	  	alpha = a[i+1][con]/a[con][con];
  	  	for(j=con;j<=N;j++)
        a[i+1][j] = a[i+1][j] - alpha*a[con][j];
  	  }     
  	}	
  	else 
	{
	   for(i=con+1;i<N;i++)
  	   {
  	       if(a[i][con]!=0)
  	       {
  	       	j = i;
  	       	break;
  	       }
  	   }
  	   swapRow(a, con, j);    //Initiating swap function
  	   for(i=con;i<N-1;i++)
  	   {
  	   	 alpha = a[i+1][con]/a[con][con];
  	   	 for(j=con;j<=N;j++)
         a[i+1][j] = a[i+1][j] - alpha*a[con][j];
  	   }	   
    }
  }  
 //Back-substitution
  del[N-1] = a[N-1][N]/a[N-1][N-1];
  float sum;
  for(i=N-2;i>=0;i--)
  {
    sum=0;
    for(j=i+1;j<N;j++)
      sum+=a[i][j]*del[j];
    del[i] = (a[i][N] - sum)/a[i][i];
  }	
}
//Swap function to incorporate partial pivoting 
void swapRow(double a[N][N+1], int con, int j)
{
   	float temp;
   	for(int i=con;i<=N; i++)
   	 {
   	   temp = a[con][i];
	   a[con][i] = a[j][i];
	   a[j][i] = temp;	 	
   	 }
}
