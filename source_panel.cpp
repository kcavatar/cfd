//FIELD : Computational Aerodynamics
//This program computes pressure distribution over an airfoil by invoking source panel method
//Input is done using a text file containing a set of coordinates describing the airfoil
//G R Krishna Chand Avatar, 2016
/* STEPS
   1. Input the text file containing coordinates describing the airfoil 
   2. Find midpoints of the panels
   3. Find the inclination angle of the panels to the x-axis
   4. Apply surface tangency condition of the flow to all the panels
   5. Get source strength distribution for every panel
   6. Find velocity of flow for every panel 
   7. Find coefficient of pressure for every panel
 */
#include<math.h>
#include<stdio.h>
#include<string.h>
#include<iostream>
#include "Eigen/Dense"

#define PI 3.1415
using namespace std;
using namespace Eigen;
double taninv(double x, double y); // Standard tan inverse function tan-1(y/x)
double influence(double beta, double theta, double xi, double yi, double xj, double yj, double plen);  // influence coffecient
double indtanvel(double beta, double theta, double xi, double yi, double xj, double yj, double plen);  // induced tangential velocity
typedef struct {
	double x;
	double y;
	double theta;  //inclination of panel
	double beta;  //ang b/w freestream vel vec and normal to panel
        double len;   //panel length 
        double str;   //panel source strength
        double vel;   //tangential velocity
        double cp;    //coefficient of pressure
} Panel;  // Type Definition of struct block as Panel
Panel panel[100]; // Max no of variables for the Struct Panel is "N-1"
int main()
{
  printf("****  SOURCE PANEL PRESSURE COEFFICIENT CALCULATOR FOR AN AIRFOIL ****\n\n");
  //Initialisation of input variables and file
  double x[100], y[100];
  FILE *p;
  p = fopen("coordinates.txt", "r");
  int N=0;
  double aoa, U;
  printf("Enter freestream velocity (U) in m/s and angle of attack (aoa) in degrees: \t");
  scanf("%lf%lf", &U, &aoa); 
  aoa = aoa*3.1415/180;
                                    /* Read the input coordinates of the end-points of the source panels */
  if(p!=NULL)
  {
   printf("\nThe coordinates of the source panels describing the airfoil as:\n\n");
   printf("\t x\t  y\n");
   while(!feof(p))  //End of file indicator
   {
   	 fscanf(p,"%lf%lf", &x[N], &y[N]);
   	 printf("\t%.6lf\t%.6lf", x[N],y[N]);
   	 N+=1;
   	 printf("\n");
   }
   N-=1;
   printf("\nN = %d\n", N);
   fclose(p);
                                   /*   Finding mid-points of the panels (control points) and their inclination angles    */
   for(int i=0;i<N-1;i++)
   {
  	panel[i].x = (x[i+1] + x[i])/2;
  	panel[i].y = (y[i+1] + y[i])/2;
        panel[i].len = sqrt((x[i+1]-x[i])*(x[i+1]-x[i]) + (y[i+1]-y[i])*(y[i+1]-y[i]));
   }
   
   for(int i=0;i<N;i++)
   {
     double del_x = panel[i+1].x - panel[i].x;
     double del_y = panel[i+1].y - panel[i].y;
     /*if(del_x>=0&&del_y>=0)
           panel[i].theta = atan(del_y/del_x);
     if(del_x<0&&del_y>=0)
           panel[i].theta = PI + atan(del_y/del_x);
     if(del_x<0&&del_y<0)
           panel[i].theta = PI + atan(del_y/del_x);
     if(del_x>=0&&del_y<0)
          panel[i].theta = 2*PI + atan(del_y/del_x);
     */
     panel[i].theta = taninv(del_x, del_y);
     panel[i].beta = (panel[i].theta-aoa) + (PI/2);
   }
                            /* Deriving system of equations for source strenth of all the panels       */ 
                                      /*        COMPUTING INFLUENCE COEFFICIENTS          */
   
   double I[N][N], R[N];
   //MatrixXd A(N,N);
   //VectorXd b(N);
   for(int i=0;i<N;i++)
   {
     int j=0;
     while(j<N)
     {
      if(j==i)
        I[i][j] = 0.5;
      else
        I[i][j] = (0.5/PI)*influence(panel[i].beta, panel[j].theta, panel[i].x, panel[i].y, x[j], y[j], panel[j].len);
      //A (i,j) = I[i][j];
      ++j;
     }
     R[i] = -U*cos(panel[i].beta);
     //b(i) = R[i];
   }
                                   
    //                             Calculation of source strengths of all the panels   

   //VectorXd str = A.fullPivLu().solve(b);
   /* for(unsigned int i=0;i<N;i++)
      panel[i].str = str(i);
   
   */
   
   for(int i=0;i<N;i++)
   {
     panel[i].vel = U*sin(panel[i].beta);
    for(int j=0;j<N;j++)
      if(j!=i)
       panel[i].vel += panel[j].str*0.5/PI + indtanvel(panel[i].beta, panel[j].theta, panel[i].x, panel[i].y, x[j], y[j], panel[j].len);
   }
  
                                           /* CALCULATION OF PRESSURE COEFFICIENTS */
   FILE *p = fopen("press_coeff_dist.txt", "a");
   fprintf("# Panel number \t Pressure coefficient");
   for(int i=0;i<N;i++)
   {
     panel[i].cp = 1.0 - (panel[i].vel/U)*(panel[i].vel/U); 
     fprintf(p,"%d\t%lf", i+1, panel[i].cp);
   }
   fclose(p);
   
  }
  else
  {
   printf("The file couldn't be opened.\n");
  }
  
  printf("\n # No modifications made to the input file.\n");
  
  //getch(); // Keeps the program running for some more time
  return 0;	
}

double influence(double beta, double theta, double xi, double yi, double xj, double yj, double plen)
{
 double A = (xj-xi)*cos(theta) + (yj-yi)*sin(theta);
 double B = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);
 double C = (-cos(beta)*cos(theta) - sin(beta)*sin(theta));
 double D = (xi-xj)*cos(beta) + (yi-yj)*sin(beta);
 double E = (xi-xj)*sin(theta) - (yi-yj)*cos(theta);

 double temp = 0.5*C*log((plen*plen + 2*A*plen + B)/B) + ((D-A*C)/E)*(atan((plen+A)/E) - atan(A/E));
 return temp;
}

double indtanvel(double beta, double theta, double xi, double yi, double xj, double yj, double plen) // induced tangential velocity
{
  double A = (xj-xi)*cos(theta) + (yj-yi)*sin(theta);
  double B = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);
  double C = (-cos(beta)*cos(theta) - sin(beta)*sin(theta));
  double D = (xi-xj)*cos(beta) + (yi-yj)*sin(beta);
  double E = (xi-xj)*sin(theta) - (yi-yj)*cos(theta);
  double temp = ((D-A*C)*0.5/E)*log((plen*plen+2*A*plen+B)/B) - C*(atan((plen+A)/E) - atan(A/E));
  return temp;
}

double taninv(double del_x, double del_y)
{
  if(del_x>=0&&del_y>=0)
     return atan(del_y/del_x);
  if(del_x<0&&del_y>=0)
     return (PI + atan(del_y/del_x));
  if(del_x<0&&del_y<0)
     return (PI + atan(del_y/del_x));
  if(del_x>=0&&del_y<0)
     return (2*PI + atan(del_y/del_x));
}


