/*

  Adam Majewski
  adammaj1 aaattt o2 dot pl  // o like oxygen not 0 like zero 
  
  
  console program in c programing language 
===============================================================

use potential ( real number) as a measure of exterior





  
  ==============================================
  
  
  Progrm flow: 
  * fill dData ( potential, double)
  * compute rgb color based on the potential (from from dData)  and save the color to rgbData
  * save rgbData to ppm file
  
  ============== Image X ========================
  
  DrawImageOfX -> DrawPointOfX -> ComputeColorOfX 
  
  
  
  

   
  ==========================================

  
  ---------------------------------
  indent d.c 
  default is gnu style 
  -------------------



  c console progam 
  
	export  OMP_DISPLAY_ENV="TRUE"	
  	gcc d.c -lm -Wall -march=native -fopenmp
  	time ./a.out > b.txt


  gcc d.c -lm -Wall -march=native -fopenmp


  time ./a.out

  time ./a.out >a.txt
  
  
  ./g.sh

  ----------------------
  
 real	0m19,809s
user	2m26,763s
sys	0m0,161s


  

*/

#include <stdio.h>
#include <stdlib.h>		// malloc
#include <string.h>		// strcat
#include <math.h>		// M_PI; needs -lm also
#include <complex.h> 		// complex numbers : https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
#include <omp.h>		// OpenMP



/* --------------------------------- global variables and consts ------------------------------------------------------------ */
#define VERSION 20210102

int NumberOfImages = 0;


/*
Representation Function
 http://www.mrob.com/pub/muency/representationfunction.html
*/


// Representation function
typedef enum  {
	Potential = 0,
	Angle = 1,
	Normal = 2
	
		
		} RepresentationFunctionType; 




// transfer function
typedef enum  {
	linear = 0,
	step_linear = 1,
	step_sqrt = 2
		
		} GradientType; 

typedef enum  {
	no = 0,
	average = 1	
		} BlendType; 



// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static unsigned int ixMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int ixMax;	//
static unsigned int iWidth;	// horizontal dimension of array

static unsigned int iyMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int iyMax;	//

static unsigned int iHeight = 10000;	//  




// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
static unsigned int iMax;	// = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array





// see SetCPlane
double radius = 2.5; 
complex double center = -0.75;
double  DisplayAspectRatio  = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)



// c plane = parameter plane
double CxMin ;	//-0.05;
double CxMax ;	//0.75;
double CyMin ;	//-0.1;
double CyMax ;	//0.7;






double PixelWidth;	// =(CxMax-CxMin)/ixMax;
double PixelHeight;	// =(CyMax-CyMin)/iyMax;
double ratio;





// complex numbers of parametr plane 



//  maximal number of iterations
static unsigned long int iterMax = 1000000;	//iHeight*100;
const int iterMax_pot = 400; // potential 

const int iterMax_normal = 2000; //  N in wiki 

/* bail-out value for the bailout test for escaping points
 radius of circle centered ad the origin, exterior of such circle is a target set  */
const double ER_normal = 1000; // big !!!!

double ER = 200.0;		// EscapeRadius for bailout test 
double EscapeRadius=1000000; // = ER big !!!!
 
double ER_POT = 100000.0;  // sqrt(1e24)

double loger; // = log(ER_LSM); // for texture
static double TwoPi=2.0*M_PI; // texture
double MaxFinalRadius;


// potential
double MaxImagePotential;
//
double BoundaryWidth = 3.0; // % of image width  
double distanceMax; //distanceMax = BoundaryWidth*PixelWidth;



/* colors = shades of gray from 0 to 255 */


unsigned char iColorOfExterior = 250;
unsigned char iColorOfInterior = 200;
unsigned char iColorOfInterior1 = 210;
unsigned char iColorOfInterior2 = 180;
unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 30;

// ----------memmory 1D arrays ==================

// The size of array has to be a positive constant integer 
static unsigned int iSize;	// = iWidth*iHeight; 

// array of doubles for better percision
double *dData;


// rgb array = 24bit color = 3 bytes
int iColorSize = 3 ; // RGB = 3*(unsigned char)
unsigned int iSize_rgb; // number of elements in rgb array
unsigned char *rgbData; // for ppm file
//  virtual 2D array of pixels 
// image = file on the disk


 




/* ------------------------------------------ functions -------------------------------------------------------------*/

/**
 * Find maximum between two numbers.
 https://codeforwin.org/2016/02/c-program-to-find-maximum-and-minimum-using-functions.html
 */
double max(double n1, double n2)
{
    return (n1 > n2 ) ? n1 : n2;
}



//---------------------

double min(double n1, double n2)
{
    return (n1 < n2 ) ? n1 : n2;
}


double clip(double d){

	return (d> 1.0) ? 1.0 : d;
}



double frac(double d){

	double fraction = d - ((long)d);
	return fraction;
}




//------------------complex numbers -----------------------------------------------------




double c_arg(complex double z)
{
 double arg;
 arg = carg(z);
 if (arg<0.0) arg+= TwoPi ; 
 return arg; 
}

double c_turn(complex double z)
{
 double arg;
 arg = c_arg(z);
 return arg/TwoPi; 
}


double turn( double x, double y){
	double t = atan2(y,x);
	if ( t<0) t+= TwoPi ;
	return t/TwoPi ;

}



// from screen to world coordinate ; linear mapping
// uses global cons
double GiveCx ( int ix)
{
  return (CxMin + ix * PixelWidth);
}

// uses globaal cons
double GiveCy (int iy) {
  return (CyMax - iy * PixelHeight);
}				// reverse y axis


complex double GiveC( int ix, int iy){
  double Cx = GiveCx(ix);
  double Cy = GiveCy(iy);
	
  return Cx + Cy*I;
	
	


}




int SetCPlane(complex double center, double radius, double a_ratio){

  CxMin = creal(center) - radius*a_ratio;	
  CxMax = creal(center) + radius*a_ratio;	//0.75;
  CyMin = cimag(center) - radius;	// inv
  CyMax = cimag(center) + radius;	//0.7;
  return 0;

}







// ****************** DYNAMICS = trap tests ( target sets) ****************************







/* -----------  array functions = drawing -------------- */

/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i (unsigned int ix, unsigned int iy)
{
  return ix + iy * iWidth;
}


// =============================  tests ============================================


// Check Orientation of c-plane image : mark first quadrant of complex plane 
// it should be in the upper right position
// uses global var :  ...
int CheckCPlaneOrientation(unsigned char A[] )
{
 
	double Cx, Cy; //  C= Cx+Cy*i;
	unsigned i; /* index of 1D array */
	unsigned int ix, iy;		// pixel coordinate 
	
	fprintf(stderr, "compute image CheckOrientation\n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy, i, Cx, Cy) shared(A, ixMax , iyMax) 
	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix){
    			// from screen to world coordinate 
    			Cy = GiveCy(iy);
    			Cx = GiveCx(ix);
	  		i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
	  		if (Cx>0 && Cy>0) A[i]=255-A[i];   // check the orientation of Z-plane by marking first quadrant */
    		}
    	}
   
   
  	return 0;
}














 

// -------------------------- potential========


double ComputePotential(const complex double c){

	double potential = 0.0; // interior
	double s = 0.5;
	complex double z = 0.0;
	double r;
	int iter;
	
	for (iter = 0; iter < iterMax_pot; ++iter){
		
		z = z*z +c; // complex quadratic polynomial
		s *= 0.5;  // 
		r = cabs(z);
		if (r > ER_POT) {break;}
	}
	
	
	
	
	if ( iter == iterMax_pot)
		{ potential = -1.0; } // interior
		else { // exterior and boundary
			potential =  s*log2(r); // log(zn)* 2^(-n)
			potential = fabs(log(potential)); // gives level sets of potential, fabs because log(potential) is < 0
			}
	
	return potential;
	
}


unsigned char ComputePotentialColor(const double potential, const RepresentationFunctionType RepresentationFunction){

	if ( potential >25.0 ){ return 0 ;}// boundary and exterior near boundary = black
     	if ( potential >10.0 ) {return 255;} // 10<potential<25; exterior noisy part 
     	// potential < 10 ; exterior not noisy
     	
     	
     	double p = frac(potential) ; // step function
     	
     	switch(RepresentationFunction){
     	
     		case linear: {p = potential; break;}
     		
     		case step_linear: {p = frac(potential); break;}
     		
     		case step_sqrt: { p = frac(potential); p = sqrt(p); p = 1.0 - p; break;}
     		
     		default: {}
     	
     	}
     	
     	return 255*p; // linear scale  
     							
     							


}



 


 
 
// ****************************************************************************************************
// ****************************** dData **************************************************************
// ***************************************************************************************************


// compute and save  raster point (ix,iy) data
int ComputePoint_dData (double A[], int ix, int iy)
{
  int i;			/* index of 1D array */
  double potential;
  complex double c;


  i = Give_i (ix, iy);		/* compute index of 1D array from indices of 2D array */
  c = GiveC(ix,iy);
  potential = ComputePotential(c);
  
  A[i] = potential;		// 
  
  return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int Fill_dDataArray (double A[])
{
  int ix, iy;		// pixel coordinate 

  	//printf("compute image \n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){ 
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			ComputePoint_dData(A, ix, iy);	//  
  }

  return 0;
}
 

/* 
 roughly speaking the code is

For each pixel:
   x = average change in iteration count of pixels to left and right.
   y = average change in iteration count above and below.
   colourIndex = atan2(x,y)
 
*/ 
unsigned char GiveAngleT( const int i, const double D[]){

	unsigned char t;
	int ix;
	int iy;
	double dx;
	double dy;
	double angle;
	// compute (ix and iy) from i
	// i = ix + iy * iWidth;
	iy = i / iWidth;
	if (iy>iHeight || iy<0) {fprintf(stderr, " bad iy = %d\n", iy);}
	ix = i - iy*iWidth;
	if (ix>iWidth || ix<0) {fprintf(stderr, " bad ix = %d\n", ix);}
	
	
	/*
	compute dx and dy 
	dx = (distance[0][1] - distance[2][1]);
	dy = (distance[1][0] - distance[1][2]);
	*/
	if ( ix-1> -1 && ix+1 < iWidth && iy-1> -1 && iy+1 < iHeight){
		dx = D[Give_i(ix-1, iy)] - D[Give_i(ix+1, iy)];
		dy = D[Give_i(ix, iy-1)] - D[Give_i(ix, iy+1)];
		}
	// compute angle from dx and dy
	//angle = atan2(dy, dx) + M_PI;
	angle = turn(dx, dy);
	angle *= 3.4;
	// scale 
	t = angle *255; 
	return t; 
	
}

double GiveAngleTest( const int i, const double D[]){

	unsigned char t;
	int ix;
	int iy;
	double dx = 0.0;
	double dy = 0.0;
	double angle;
	// compute (ix and iy) from i
	// i = ix + iy * iWidth;
	iy = i / iWidth;
	ix = i - iy*iWidth;
	printf(" i = %d ix = %d iy = %d Give_i(ix,iy) = %d\n", i , ix, iy, Give_i(ix,iy)  );
	
	
	
	/*
	compute dx and dy 
	dx = (distance[0][1] - distance[2][1]);
	dy = (distance[1][0] - distance[1][2]);
	*/
	if ( ix-1> -1 && ix+1 < iWidth && iy-1> -1 && iy+1 < iHeight){
		dx = D[Give_i(ix-1, iy)] - D[Give_i(ix+1, iy)];
		dy = D[Give_i(ix, iy-1)] - D[Give_i(ix, iy+1)];
		}
	/* 
	For clean colouring we need to take colour from nearby pixel at stationaty points.
	if(dx == 0 && dy == 0) {
		dx = (distance[0][0] - distance[2][0]);
		if(dx == 0) {
			dx = (distance[0][2] - distance[2][2]);
			if(dx == 0) {
				dy = (distance[0][0] - distance[0][2]);
				if(dy == 0) dy = (distance[2][0] - distance[2][2]);
			}
		}
	}	
	*/	
	
	
		
		
	// compute angle from dx and dy
	//angle = atan2(dy, dx) + M_PI;
	angle = turn(dx, dy);
	
	printf("angle = %f \n", angle);
	// scale 
	t = angle *255; 
	return t; 
	
}



/* 
 The dot product of two vectors a = [a1, a2, ..., an] and b = [b1, b2, ..., bn] is defined as:[1]
 d = a1b1 + a2b2
*/
double cdot(double complex a, double complex b) {
  return creal(a) * creal(b) + cimag(a) * cimag(b);
}

// 
// output 
// 
double GiveReflection(double complex C, int iMax, double ER) {
  int i = 0; // iteration 

  double complex Z = 0.0; // initial value for iteration Z0
  double complex dC = 0.0; // derivative with respect to c 
  double reflection = -1.0; // inside 

  double h2 = 1.5; // height factor of the incoming light
  double angle = 45.0 / 360.0; // incoming direction of light in turns 
  double complex v = cexp(2.0 * angle * M_PI * I); // = exp(1j*angle*2*pi/360)  // unit 2D vector in this direction
  // incoming light 3D vector = (v.re,v.im,h2)

  double complex u;

  for (i = 0; i < iMax; i++) {
    dC = 2.0 * dC * Z + 1.0;
    Z = Z * Z + C;

    if (cabs(Z) > ER) { // exterior of M set
      u = Z / dC;
      u = u / cabs(u);
      reflection = cdot(u, v) + h2;
      reflection = reflection / (1.0 + h2); // rescale so that t does not get bigger than 1
      if (reflection < 0.0) reflection = 0.0;
      break;

    }
  }

  return reflection;
}




// it do not use data from double array
unsigned char GiveNormalColor(const int i ) {

	double complex c;
	int ix;
	int iy;
	double reflection;
	unsigned char g;
	
	
	// compute (ix and iy) from i
	// i = ix + iy * iWidth;
	iy = i / iWidth;
	if (iy>iHeight || iy<0) {fprintf(stderr, " bad iy = %d\n", iy);}
	ix = i - iy*iWidth;
	if (ix>iWidth || ix<0) {fprintf(stderr, " bad ix = %d\n", ix);}
	// compute c from ix and iy
	c = GiveC(ix,iy);
	// compute color ( shade of gray) from c
	reflection = GiveReflection(c, iterMax_normal, ER_normal);
	g = 255*reflection; // change range
	return g; 
	
	
	


}






 




 
 
// compute color ( shade) of gradient 
unsigned char GiveNonBlendedColor( const int i, const double D[], const double potential, RepresentationFunctionType RepresentationFunction, GradientType Gradient){


	unsigned char g;
	 
	
	switch (RepresentationFunction){
	
		case Potential: { g = ComputePotentialColor(potential, Gradient); break;}
     			
     		case Angle: {g = GiveAngleT(i, D); break;}
     		
     		case Normal: {g = GiveNormalColor(i); break;}
     		
     		default: {}
     		}
     	return g;



}
 

unsigned char GiveBlendedColor( double c1, double c2, BlendType Blend){

	unsigned char t;
	
	switch (Blend){
	
		case average: {t = (c1+c2)/2.0; break;}
		
		
		default: {}
	
	}
	
	return  t;


}


 
unsigned char GiveExteriorColor(const int i, const double D[], const double potential, RepresentationFunctionType RepresentationFunction, GradientType Gradient, BlendType Blend){


	unsigned char t;
	double p;
	double n;
	
	
	if (!Blend )
		{t = GiveNonBlendedColor(i, D, potential,  RepresentationFunction, Gradient); }
		else {
			p = ComputePotentialColor(potential, Gradient); // 
			n = GiveNormalColor(i);
			
			t = GiveBlendedColor(p, n,  Blend);
		
		
		}
			
		
		
			
	return t;

	


} 
 

/*

input :
* int i
* array D of double numbers ( distance)

output : array of rgb colors 

*/
void ComputePointColorAndSave(const int i, const double D[], RepresentationFunctionType RepresentationFunction, GradientType Gradient, BlendType Blend, unsigned char  C[] ){


	
	int iC = i*iColorSize; // compute index of F array
	// color channels from 0 to 255  
	//unsigned char R;
	//unsigned char G;
	//unsigned char B;
	unsigned char t; 
	
	double potential = D[i]; // read 
	
	// compute color
	if (potential<0.0)
		{//t = 255;
		// save color to the rgb array 
		C[iC] 	= 0;
		C[iC+1] = 0;
		C[iC+2] = 127; // blue
		
		
		
		
		} // interior
		else { t = GiveExteriorColor(i, D, potential,RepresentationFunction, Gradient, Blend);
			// save color to the rgb array 
			C[iC] 	= t;
			C[iC+1] = t;
			C[iC+2] = t;
		
		
		
		
		} // exterior
			
		
	
	
	
	
	
	

}







// fill array f using data from d array
// uses global var :  ...
int Fill_rgbData_from_dData (double D[], RepresentationFunctionType  RepresentationFunction, GradientType Gradient, BlendType Blend, unsigned char C[])
{
  int i=0;		// array index 

  	fprintf(stderr, "\nFill_rgbData_from_dData\n");
  	//printf("compute image \n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(i) shared( D, C, iSize)
  	for (i = 0; i < iSize; ++i){
    		//fprintf (stderr, "rgb  %d from %d \r", i, iSize);	//info 
    		ComputePointColorAndSave(i, D, RepresentationFunction, Gradient, Blend, C);	//  
  }
  
  
	
  return 0;
}

 
 
// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************



int Save_PPM( unsigned char A[], double k, char* sName, char* comment )
{
  
  	FILE * fp;
  
  	char name [100]; /* name of file */
  	snprintf(name, sizeof name, "%s", sName); /*  */
	char *filename =strcat(name,".ppm");
	
	char long_comment[200];
  	sprintf (long_comment, "fc(z)=z^2+ c %s", comment);
  
  
  
  	// save image to the pgm file 
  	fp= fopen(filename,"wb"); // create new file,give it a name and open it in binary mode 
  
  	if (!fp ) { fprintf( stderr, "ERROR saving ( cant open) file %s \n", filename); return 1; }
  	// else
  	fprintf(fp,"P6\n%d %d\n255\n",  iWidth, iHeight);  // write header to the file
  	size_t rSize = fwrite(A, sizeof(A[0]), iSize_rgb,  fp);  // write array with image data bytes to the file in one step 
  	fclose(fp); 
  
  	// info 
  	if ( rSize == iSize_rgb) 
  		{
  			printf ("File %s saved ", filename);
  			if (long_comment == NULL || strlen (long_comment) == 0)
   				{printf ("\n"); }
  				else { printf (". Comment = %s \n", long_comment); }
  		}
  		else {printf("wrote %zu elements out of %u requested\n", rSize,  iSize);}
  	
  	
  	// 
	NumberOfImages +=1; // count images using global variable	
    
	return 0;
}







void Test(){

	
	double dx =  PixelWidth/5.0;
	double p;
	complex double c = 0.25;
	int i;
	
	for (i = 0; i < 40; ++i)
		{
		
			p = ComputePotential( c);
			int n = isnormal(p);
			int f = isinf(p);
			printf(" x = %.16f\t p = %f \t isnormal = %d\t isinf = %d\n", creal(c), p,n, f);
			c += dx;
		
		
		}

	
	



}














int PrintInfoAboutProgam()
{
	printf("Number of pgm images = %d \n", NumberOfImages);	
  
  	// display info messages
  	printf ("Numerical approximation of M set for fc(z)= z^2 + c \n");
  	//printf ("iPeriodParent = %d \n", iPeriodParent);
  	//printf ("iPeriodOfChild  = %d \n", iPeriodChild);
  
  	printf ("Image Width = %f in world coordinate\n", CxMax - CxMin);
  	printf ("PixelWidth = %f \n", PixelWidth);
	  
  	printf("for DEM\n");
  	if ( distanceMax<0.0 || distanceMax > ER ) printf("bad distanceMax\n");
	printf("Max distance from exterior to the boundary =  distanceMax = %.16f = %f pixels\n",  distanceMax, BoundaryWidth); 
  	printf("\n");
  
  
  // image corners in world coordinate
  // center and radius
  // center and zoom
  // GradientRepetition
  printf ("Maximal number of iterations = iterMax = %ld \n", iterMax);
  
  
  
  
  printf ("ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  //
  printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
  // OpenMP version is displayed in the console 
  
  
  
  return 0;
}






// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************

int setup ()
{

  fprintf (stderr, "setup start\n");
  
  
  
  
  
	
  /* 2D array ranges */
  
  iWidth = iHeight* DisplayAspectRatio;
  iSize = iWidth * iHeight;	// size = number of points in array 
  iSize_rgb = iSize* iColorSize;
  
  
  // iy
  iyMax = iHeight - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  //ix

  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  
  
  SetCPlane( center, radius,  DisplayAspectRatio );	

  /* Pixel sizes */
  PixelWidth = (CxMax - CxMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
  PixelHeight = (CyMax - CyMin) / iyMax;
  ratio = ((CxMax - CxMin) / (CyMax - CyMin)) / ((double) iWidth / (double) iHeight);	// it should be 1.000 ...
	
  
  //ER2 = ER * ER; // for numerical optimisation in iteration
  
  
  
   	
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  dData = malloc (iSize * sizeof (double));
  rgbData =  malloc (iSize_rgb * sizeof (unsigned char));
  
  	
  if (dData == NULL || rgbData == NULL ){
    fprintf (stderr, " Could not allocate memory");
    return 1;
  }

  
 	
  
  BoundaryWidth = 6.0*iWidth/2000.0  ; //  measured in pixels ( when iWidth = 2000) 
  distanceMax = BoundaryWidth*PixelWidth;
  
  
  
  fprintf (stderr," end of setup \n");
	
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




int end(){


  fprintf (stderr," allways free memory (deallocate )  to avoid memory leaks \n"); // https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free(dData);
  free(rgbData);
 
  PrintInfoAboutProgam();
  return 0;

}

// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main () {
  
  
  
  	setup ();
  
  	  
	Fill_dDataArray(dData);
	//
	Fill_rgbData_from_dData (dData, Potential, step_sqrt, no, rgbData);
	Save_PPM(rgbData, step_linear, "step_sqrt", "potential");
	
	Fill_rgbData_from_dData (dData, Angle, linear, no, rgbData);
	Save_PPM(rgbData, Angle, "angle_linear", "potential angle");
	
  	Fill_rgbData_from_dData (dData, Normal, linear, no, rgbData);
	Save_PPM(rgbData, Angle, "normal", "Normal mapping or slope");
	
	Fill_rgbData_from_dData (dData, Normal, step_sqrt, average, rgbData);
	Save_PPM(rgbData, Angle, "average_sqrt", "average blend = (Potenital_step_sqrt + Normal) / 2 ");
  	
  	
  	end();
  
 	//Test(); 
  	

  	return 0;
}
