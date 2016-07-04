#include<iostream>
#include<math.h>
#define R 8.314  // J/mol/K
#define N 5      //No of species = 5
double tol = 0.001;
using namespace std;
void eqmCons( double Kp[], double T, double rho);
void gibbs(double a[N][7], double Gibbs[], double T);
   /*      Various models to solve for the composition of the species 
            Model 1 : for low temperatures                                                 */
void model_1(double mol[], double K[], double mol_O, double mol_N);
void model_2(double mol[], double K[], double mol_O, double mol_N);
void model_3(double mol[], double K[], double mol_O, double mol_N);
   
int main()
{
  cout<<"*********** MOLE NUMBER BASED EQUILIBRIUM COMPOSITION  SOLVER  ***************\n\n";
  /*Rxns are : 1. O2 -> 2O    2. O2 + N2 -> 2NO  3. N2 -> 2N
    Species are : 1(0) -> N2  2(1) -> N  3(2) -> O  4(3) -> O2  5(4) -> NO
    Input initial parameters                                 */
  cout<<"Input the  values of temperature  (T), density (rho): \n\n";
  double T, rho;
  cin>>T>>rho;
  //Count variables
  int i=0, j=0;
 //Calculate equilibrium constants 
  double K[3];
  eqmCons(K, T, rho);
  for(i=0;i<3;i++)
   cout<<"\nThe value of Kp["<<i<<"] : "<<K[i]<<endl;
  //Conserving nitrogen and oxygen atoms using N/O ratio and sum of mole fractions 
  double mol[N] = {0.0, 0.0, 0.0, 0.0, 0.0};        // "Initial" composition
  double mol_O = 42.0, mol_N = 158.0;     //Conserving moles of Oxygen and Nitrogen mole percents
  if(T>=200&&T<=1700)
    model_1(mol, K, mol_O, mol_N);
  if(T>1700&&T<=3500)
    model_2(mol, K, mol_O, mol_N);
  if(T>3500&&T<=9000)
    model_3(mol, K, mol_O, mol_N);
  
  //Display results obtained as such if required
  double total = mol[0]+mol[1]+mol[2]+mol[3]+mol[4];
  cout<<"\n\n**** The results are as : ***** \n";
  cout<<"\nTotal number of moles :  "<<total<<endl<<endl;
  
  for(i=0;i<N;i++)
   cout<<"Species "<<i+1<<" : "<<mol[i]<<"       "<<log10(mol[i]/total)<<endl;
  return 1;
}

void eqmCons(double Kp[], double T, double rho)
{
                        /* COEFFICIENT TABLE FOR VARIOUS TEMPERATURE RANGES */	
  //Temperature : 200 - 800 K
       double a[N][7] = {{0.346226e1,  0.582024e-3,-0.305255e-5,  0.622801e-8,-0.337560e-11, 0.879517e0,  0.321926e1}, //N2
                         {0.250000e1,          0.0,         0.0,          0.0,          0.0, 0.566267e5, 0.418073e1},  //N
                         {0.321671e1, -0.378227e-2, 0.847468e-5, -0.886582e-8, 0.353656e-11, 0.296405e5, 0.185264e1},  //O
	                     {0.377037e1, -0.289522e-2, 0.953322e-5, -0.924699e-8, 0.301919e-11, -0.188598e2, 0.369335e1}, //O2
                         {0.420644e1, -0.450984e-2, 0.105574e-4, -0.859194e-8, 0.240471e-11, 0.108890e5, 0.231379e1}   //NO
                        };
  //Temperature : 800 K - 3000 K
       double b[N][7] = {{0.270224e1, 0.194439e-2, -0.893000e-6, 0.197392e-9, -0.169678e-13, 0.201802e3, 0.720408e1},  //N2
                         {0.250751e1,-0.247979e-4,  0.296415e-7,-0.152881e-10, 0.289137e-14, 0.566249e5, 0.414319e1},  //N
                         {0.260454e1,-0.172355e-3, 0.115741e-6,-0.364179e-10, 0.460115e-14, 0.297290e5, 0.458653e1},   //O
                         {0.289692e1, 0.237365e-2, -0.149171e-5, 0.466034e-9, -0.539452e-13, 0.822404e2, 0.750194e1},  //O2
                         {0.275438e1, 0.230933e-2, -0.128234e-5, 0.340435e-9, -0.348075e-13, 0.111343e5, 0.907897e1}   //NO
				         
                        };
  //Temperature : 3000K - 6000K
       double c[N][7] = {{0.391435e1, 0.315371e-3, -0.564810e-7, 0.360123e-11, 0.523594e-16, -0.535520e3, 0.216720e-1}, //N2
                         {0.263761e1,-0.873733e-5, -0.647727e-7, 0.234734e-10, -0.173962e-14, 0.564523e5, 0.322321e1},  //N
                         {0.281017e1,-0.290399e-3,  0.908333e-7,-0.994278e-11, 0.377044e-15,  0.295366e5, 0.325365e1},  //O
	                     {0.284211e1, 0.133206e-2, -0.339159e-6, 0.446522e-10, -0.229148e-14, 0.583504e3, 0.862550e1},  //O2 
  	                     {0.380154e1, 0.498575e-3, -0.125313e-6, 0.140939e-10, -0.448202e-15, 0.106666e5, 0.316194e1}   //NO
                        };
  //Temperature : 6000K - 10000K
       double d[N][7] = {{0.126575e1,  0.182698e-2,-0.375839e-6, 0.330337e-10, -0.936511e-15,  0.314270e4,  0.179433e2 }, //N2
                         {0.337206e1, -0.885546e-3, 0.252933e-6, -0.231879e-10, 0.704714e-15,  0.562702e5, -0.105640e1 }, //N
		            	 {0.192093e1, 0.217766e-3, -0.182884e-7, 0.490510e-12,  0.285073e-17,  0.307834e5,  0.927486e1 }, //O
		             	 {0.568211e1, -0.753006e-3, 0.229801e-6, -0.239559e-10, 0.800485e-15, -0.247004e4, -0.987401e1},  //O2
		                 {0.491332e1, 0.617553e-6, -0.552224e-7, 0.116865e-10, -0.546425e-15,  0.884538e4, -0.457869e1}   //NO	    
                         }; 
  //Temperature : 10000K - 15000K
                
  double Gibbs[N];
  if(T>=200&&T<800)
    gibbs(a, Gibbs, T);
  if(T>=800&&T<3000)
    gibbs(b, Gibbs, T);
  if(T>=3000&&T<6000)
    gibbs(c, Gibbs, T);
  if(T>=6000&&T<=10000)
    gibbs(d, Gibbs, T);
 // Species are : 1(0) -> N2  2(1) -> N  3(2) -> O  4(3) -> O2  5(4) -> NO
  //Crude Kp; Kp acting as Gibbs free energy change for the reaction
  Kp[0] = 2*Gibbs[2] - Gibbs[3]; //Reaction : O2 -> 2O 
  Kp[1] = 2*Gibbs[4] - Gibbs[0] - Gibbs[3]; //Reaction : O2 + N2 -> 2NO :
  Kp[2] = 2*Gibbs[1] - Gibbs[0]; //Reaction : N2 -> 2N
  //Redefining Kp as Kp = exp(-delG/RT)
  for(int i=0;i<3;i++)
  {
    Kp[i] = exp(-Kp[i]/(R*T));
  }
  //Redefine Kp in terms of mole numbers : P0 = 1 atm  -> Kp[i] = Kp[i]*(P0/(rho*R*T))^delta_N   ; i=0,2
  double p0 = 101325; //Unit : Pascals
  Kp[0] = Kp[0]*p0/(rho*R*T);
  Kp[2] = Kp[2]*p0/(rho*R*T);
}
void gibbs(double a[N][7], double Gibbs[], double T)
{
  double H[N], S[N]; //defined as : H = H/RT, S = S/R; Gibbs = Gibbs/RT
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
            /*                  MODEL - 1                */
void model_1(double mol[], double K[], double mol_O, double mol_N)
{
  	double molNew_O2 = 0.0, molNew_N2 = 0.0;     // Required only for two species : O2 and N2
  	mol[3] = mol_O/2;
  	mol[0] = mol_N/2;
  	while(fabs(molNew_N2-mol[0])>tol||fabs(molNew_O2-mol[3])>tol)
  	{
  	  molNew_O2 = mol[3];
	  molNew_N2 = mol[0];
  	  mol[2] = sqrt(K[0]*molNew_O2);
	  mol[1] = sqrt(K[2]*molNew_N2);
	  mol[4] = sqrt(K[1]/(K[0]*K[2]))*mol[2]*mol[1];
	  mol[3] = (mol_O - mol[2]- mol[4])/2;
	  mol[0] = (mol_N - mol[1]- mol[4])/2;
  	}  	
}

            /*                  MODEL - 2                */ 
void model_2(double mol[], double K[], double mol_O, double mol_N)
{
    double b0, b1, b2, b3, b4, F, dF;
    double Ke = sqrt(K[1]/(K[0]*K[2])); 
    b0 = 2*K[0]*mol_O*mol_O;
    b1 = (-4 + K[2]*Ke)*K[0]*mol_O;
    b2 = -8*mol_O + 2*K[0] + (mol_O - mol_N)*K[1] - K[0]*K[2]*Ke;
    b3 = 8 - K[1] - 2*K[2]*Ke;
    b4 = (8 - 2*K[1])/K[0];
   //Initial guess for species : species : O (nascent oxygen)  from quadratic eqn
    double t = 0.5*K[0]*(1+sqrt(0.5*mol_N*K[1]/K[0]));
	double molNew_sp_O = 0.5*(-t + sqrt(t*t + 2*mol_O*K[0]));
    //Iterative computation
    while(fabs(molNew_sp_O - mol[2])>tol)
    {
      mol[2] = molNew_sp_O;
      F = b0 + b1*mol[2] + b2*mol[2]*mol[2] + b3*pow(mol[2], 3) + b4*pow(mol[2], 4);
      dF = b1 + 2*b2*mol[2] + 3*b3*mol[2]*mol[2] + 4*b4*pow(mol[2], 3);
	  molNew_sp_O = mol[2] - F/dF;    //Newton-Raphson iteration
	  mol[1] = (mol_O - 2*molNew_sp_O*molNew_sp_O/K[0] - molNew_sp_O)/(Ke*molNew_sp_O);
	  mol[4] = Ke*mol[1]*molNew_sp_O;
	  mol[0] = (mol_N - mol[1] - mol[4])*0.5;
	  mol[3] = (mol_O - molNew_sp_O - mol[4])*0.5;	
    } 	
}   
            /*                  MODEL - 3                */ 
void model_3(double mol[], double K[], double mol_O, double mol_N)
{
	double molNew_sp_O = mol_O;   // mole number of species : nascent O
	double molNew_N;         // mol_N= mol_N - mol(Species : NO)
	double Ke = sqrt(K[1]/(K[0]*K[2]));
	while(fabs(molNew_sp_O - mol[2])>tol)
	{
	  mol[2] = molNew_sp_O;
	  molNew_N = mol_N - mol[4];
	  //solution of Quadratic equation for species : nascent N
	  mol[1] = (-K[2]*0.5 + sqrt(K[2]*K[2]*0.25 + 2.0*K[2]*molNew_N))*0.5;
	  mol[0] = (molNew_N - mol[1])*0.5;
	  mol[3] = mol[2]*mol[2]/K[0];
	  mol[4] = Ke*mol[2]*mol[1];
	  molNew_sp_O = mol_O - 2*mol[3] - mol[4];
	} 	
}     
