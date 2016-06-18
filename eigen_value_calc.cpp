#include<stdio.h>
#include<math.h>

typedef struct{
	float x;
	float y;
} Complex;

int i,j,n; //Global variables

     /*  ----  Declaration of Eigen calculation functions ----- */
void eigen2(int n);
void eigen3(int n);

    /* ------ MAIN FUNCTION ---------*/
int main()
{
	printf("**********   EIGEN-VALUE CALCULATOR  (upto order 3)  ************\n\n");
	//Define variables
	printf("=> Enter the order of square matrix :\t");
    scanf("%d", &n);
    if(n>3) 
	{
	  printf("\n\n    OUT OF BOUNDS!!!!");	
	  return false;
	}
	printf("\n=> Enter the contents of square matrix :\n");
	
	if(n==2)
	 eigen2(n);
	else if(n==3)
	 eigen3(n);
	
	return true;
}

void eigen2(int n)
{
	float a[n][n];
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		scanf("%f", &a[i][j]);
	}
	float D_sq;
	D_sq = pow(a[0][0] + a[1][1], 2) - 4*(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
	float r1, r2; //Roots if D_sq>=0
	printf("\n=> The eigen values of the matrix are :  \n\n");
	if(D_sq>=0)
	{
		r1 = 0.5*((a[0][0] + a[1][1]) + sqrt(D_sq));
		r2 = 0.5*((a[0][0] + a[1][1]) - sqrt(D_sq));
		printf("%f\t%f", r1, r2);
	}
	else
	{
		Complex c1;
		c1.x = 0.5*(a[0][0] + a[1][1]);
		c1.y = 0.5*sqrt(-D_sq);
		printf("\t%.4f + i%.4f\t    %.4f - i%.4f", c1.x, c1.y, c1.x, c1.y);
	}	 
}

void eigen3(int n)
{
	float x[n][n];
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		scanf("%f", &x[i][j]);
	}
	
	printf("\n=> The eigen values of the matrix are :  \n\n");
	//Case : 1  -> All roots are real and equal
	if(a*a/9 == b/3)
	 printf("%.4f\t%.4f\t%.4f", -a/3, -a/3, -a/3);
	//Case : 2  -> All roots are real and different
	 
	
}

