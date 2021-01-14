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
	Normal = 1
	
		
		} RepresentationFunctionType; 




// transfer function
typedef enum  {
	linear = 0,
	step_linear = 1,
	step_sqrt = 2
		
		} GradientType; 
		
		
		
		
// https://en.wikipedia.org/wiki/Blend_modes
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
// standard  : ( center = -0.75 and radius = 2.5 )
// main antenna c = -1.75 radius 0.5
// period 3 center = c = -1.754877666246693  radsius +0.038683608637816 
// c = -1.711065402413374     radsius  c = 0.000355134137791


double radius ; //= 0.038683608637816; 
double radius_0 = 2.5;

double zoom ; // = 1/ radius
complex double center ; //= -1.754877666246693;
double  DisplayAspectRatio  = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)



// c plane = parameter plane
double CxMin ;	//-0.05;
double CxMax ;	//0.75;
double CyMin ;	//-0.1;
double CyMax ;	//0.7;



int ExampleNumberMax;

// http://paulbourke.net/fractals/mandelbrot/
// 4 raws of 3 columns
// type arrayName [ x ][ y ];
double examples[][3] = {
	{-0.75, 		0.0, 			2.5}, // standard
	{-1.75, 		0.0, 			0.5},
	{-1.77 ,		0.0, 			0.07 }, // period 3 center 
	{-1.711065402413374, 	0.0, 			0.008}, // good, period 8
	{-1.711161027105541, 	0.0, 			0.0009}, // good, period 8 zoom
	
	{0.365557777904776,	0.613240370349204, 	0.12}, 
	{0.391080345956122, 	0.570677592363374,  	0.01},
	{0.296294860929836,	0.017184282646391,	0.001}, // tune limits
	// period 4 mini M-set
	{-0.170337,		1.06506,  		0.32},
	{-0.170337,		1.06506,  		0.064},
	{-0.170337,		1.06506,  		0.0128}, // mv center
	{-0.170337,		1.06506,  		0.00256},
	{-0.170337,		1.06506,  		0.000512}, // mv center
	// period 6 wake ( 1/6) ? 
	{0.42884,		0.231345, 		0.06}, // mv center
	{0.42884,		0.231345, 		0.01}, // mv center
	{-1.711638937577389,	0.000449229252155, 	0.000001} //    period = 19
	
	
	
	
};



double PixelWidth;	// =(CxMax-CxMin)/ixMax;
double PixelHeight;	// =(CyMax-CyMin)/iyMax;
double ratio;





// complex numbers of parametr plane 



//  maximal number of iterations
static unsigned long int iterMax = 1000000;	//iHeight*100;
const int iterMax_pot = 4000; // potential 

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
double MaxImagePotential = 0.0;
double potential_multiplier;
// limits for potential
double potential_boundary;
double potential_noisy;



//
double BoundaryWidth = 3.0; // % of image width  
double distanceMax; //distanceMax = BoundaryWidth*PixelWidth;



/* colors = shades of gray from 0 to 255 */


unsigned char iColorOfExterior = 250;
unsigned char iColorOfInterior = 127;
unsigned char iColorOfInterior1 = 210;
unsigned char iColorOfInterior2 = 180;
unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 30;
unsigned char iColorOfNoise = 255;

// ----------memmory 1D arrays ==================

// The size of array has to be a positive constant integer 
static unsigned int iSize;	// = iWidth*iHeight; 

// array of doubles for better percision
double *dData1;
double *dData2;


// rgb array = 24bit color = 3 bytes
int iColorSize = 3 ; // RGB = 3*(unsigned char)
unsigned int iSize_rgb; // number of elements in rgb array
unsigned char *rgbData1; // for ppm file
unsigned char *rgbData2; // for ppm file
unsigned char *rgbData3; // for ppm file
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




int SetCPlane(complex double Center, double Radius, double a_ratio){

  // sete up global var
  center = Center;
  radius = Radius;

  CxMin = creal(center) - radius*a_ratio;	
  CxMax = creal(center) + radius*a_ratio;	//0.75;
  CyMin = cimag(center) - radius;	// inv
  CyMax = cimag(center) + radius;	//0.7;
  
  
  return 0;

}


// rows 
#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))


int SetCPlaneFromExamples(const int n, const double a_ratio){

	int nMax = LEN(examples);
	
	
	printf("n = %d \t nMax = %d \n",n,  nMax);
	if (n> nMax)
		{
			SetCPlane(-0.75, 2.5, a_ratio);
			fprintf(stderr, " error n>nMax\n");
			return 1;
		}
		
		
	complex double c = examples[n][0] + I*examples[n][1];
	double r = examples[n][2];
	
	
	SetCPlane(c, r, a_ratio);
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
			
	
	// compute value of global variable
	if (potential >MaxImagePotential ) {MaxImagePotential  = potential;}
	
	return potential;
	
}


unsigned char ComputePotentialColor(const double potential, const GradientType Gradient){

	
	// ranges of potential coputed earlier
	if ( potential > potential_boundary  ){ return iColorOfBoundary ;}// boundary and exterior near boundary = black
     	if ( potential > potential_noisy ) {return iColorOfNoise;} // 10<potential<25; exterior noisy part 
     	// potential < 10 ; exterior not noisy, see below
     	
     	
     	double p ; // local copy of potential
     	
     	switch(Gradient){
     	
     		case linear: {p = potential; break;}
     		
     		case step_linear: {p = frac(potential); break;}
     		
     		case step_sqrt: { p = frac(potential); p = sqrt(p); p = 1.0 - p; break;}
     		
     		default: {}
     	
     	}
     	
     	return 255*p; // change range from [0,1]  to  [0, 255] using linear scale
     							
     							


}



 


 
 





// ****************************************************************************************************
// ****************************** Normal or Slope **************************************************************
// ***************************************************************************************************



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



// ****************************************************************************************************
// ****************************** dData **************************************************************
// ***************************************************************************************************


// compute and save  raster point (ix,iy) data
int ComputePoint_dData (double A[], RepresentationFunctionType  RepresentationFunction, int ix, int iy)
{
  int i;			/* index of 1D array */
  //double potential;
  complex double c;
  double d;
  


  i = Give_i (ix, iy);		/* compute index of 1D array from indices of 2D array */
  c = GiveC(ix,iy);
  
  switch (RepresentationFunction) {
  	case Potential : {d = ComputePotential(c); break;}
  	
  	case Normal : { d = GiveReflection(c, iterMax_normal, ER_normal); break;}
  	
  	//case Angle : {break;}
  	
  	default: {}
  	
  	}
  
  A[i] = d;		// 
  
  return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int Fill_dDataArray (double A[], RepresentationFunctionType  RepresentationFunction)
{
  int ix, iy;		// pixel coordinate 

  	//printf("compute image \n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){ 
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			ComputePoint_dData(A, RepresentationFunction, ix, iy);	//  
  }

  return 0;
}
 




 
// ****************************************************************************************************
// ****************************** RGB *************************************************************
// ***************************************************************************************************



 


 
unsigned char GiveExteriorColor(const int i, const double D[], const double potential, RepresentationFunctionType RepresentationFunction, GradientType Gradient){


	unsigned char g;
	 
	
	switch (RepresentationFunction){
	
		case Potential: { g = ComputePotentialColor(potential, Gradient);  break;}  
     			
     		//case Angle: {g = 255*GiveAngleT(i, D); break;} // !!! needs full double array
     		
     		case Normal: {g = GiveNormalColor(i); break;}
     		
     		default: {}
     		}
     				
	return g;
} 
 

/*

input :
* int i
* array D of double numbers ( distance)

output : array of rgb colors 

*/
void ComputeAndSaveColor(const int i, const double D[], RepresentationFunctionType RepresentationFunction, GradientType Gradient, unsigned char  C[] ){


	
	int iC = i*iColorSize; // compute index of F array
	// color channels from 0 to 255  
	//unsigned char R;
	//unsigned char G;
	//unsigned char B;
	unsigned char t; 
	
	double d = D[i]; 
	// compute color
	if (d<0.0)
		{	// interior = solid blue
			C[iC] 	= 0;
			C[iC+1] = 0;
			C[iC+2] = iColorOfInterior; // blue
		} 
		
		else { 	// exterior = blended gray gradient
			t = GiveExteriorColor(i, D, d, RepresentationFunction, Gradient);
			// save color to the rgb array C
			C[iC] 	= t;
			C[iC+1] = t;
			C[iC+2] = t;
		} 
			
		
}





// fill array f using data from d array
// uses global var :  ...
int Fill_rgbData_from_dData (double D[], RepresentationFunctionType  RepresentationFunction, GradientType Gradient,  unsigned char C[])
{
  int i=0;		// array index 

  	fprintf(stderr, "\nFill_rgbData_from_dData\n");
  	//printf("compute image \n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(i) shared( D, C, iSize)
  	for (i = 0; i < iSize; ++i){
    		//fprintf (stderr, "rgb  %d from %d \r", i, iSize);	//info 
    		ComputeAndSaveColor(i, D, RepresentationFunction, Gradient, C);	//  
  }
  
  
	
  return 0;
}





// *******************************************************************************************
// ********************************** Blend  ****************************
// *********************************************************************************************





 
/* 
Input
* Blend,  see BlendType ( - blend mode)
* 2 colors ( in the same range)

output : color ( in the same range as input colors)
*/
unsigned char GiveBlendedColor(const double c1, const double c2, const BlendType Blend){

	unsigned char t;
	
	switch (Blend){
	
		case average: {t = (c1+c2)/2.0; break;}
		
		
		default: {}
	
	}
	
	return  t;


}






// blend Normal (slope)  and potential(GradientType)
void ComputeAndSaveBlendColor( const unsigned char C1[], const unsigned char C2[], const BlendType Blend, const int i, unsigned char C[]){

	
	
	
	unsigned char t; 
	int iC = i*iColorSize; // compute index of F array
    	double c1 = C1[iC];
    	double c2 = C2[iC];
	//// compute color
	if ( C1[iC+2] == iColorOfInterior && C1[iC]==0) // check for interior : only B  and R from RGB; see ComputeAndSaveColor
		{	// interior = solid blue
			C[iC] 	= 0;
			C[iC+1] = 0;
			C[iC+2] = iColorOfInterior; // blue
		} 
		
		else { 	// exterior = blended gray gradient
			t = GiveBlendedColor( c1 , c2, Blend);
			C[iC] 	= t;
			C[iC+1] = t;
			C[iC+2] = t; 
			
			
		} 
			
}


// 
void MakeBlendImage(const unsigned char C1[], const unsigned char C2[], const BlendType Blend, unsigned char C[]){

	 int i=0;		// array index 

  	fprintf(stderr, "\nFill_rgbData_from_2_dData\n");
  	//printf("compute image \n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(i) shared(  C1, C2, C, iSize)
  	for (i = 0; i < iSize; ++i){
    		
    		ComputeAndSaveBlendColor( C1, C2, Blend, i, C);
  }
  

	
	
}


 
 
// *******************************************************************************************
// ********************************** save A array to ppm file ****************************
// *********************************************************************************************



int Save_PPM( const unsigned char A[], const char* sName, const char* comment, const double radius  )
{
  
  	FILE * fp;
  
  	char name [100]; /* name of file */
  	snprintf(name, sizeof name, "%s_%f", sName, radius); /*  */
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
  
  
  	
  
  //ER2 = ER * ER; // for numerical optimisation in iteration
  
  
  
   	
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  dData1 = malloc (iSize * sizeof (double));
  dData2 = malloc (iSize * sizeof (double));
  rgbData1 =  malloc (iSize_rgb * sizeof (unsigned char));
  rgbData2 =  malloc (iSize_rgb * sizeof (unsigned char));
  rgbData3 =  malloc (iSize_rgb * sizeof (unsigned char));
  
  	
  if (dData1 == NULL || dData2 == NULL  || rgbData1 == NULL || rgbData2 == NULL || rgbData3 == NULL){
    fprintf (stderr, " Could not allocate memory");
    return 1;
  }

  ExampleNumberMax = LEN ( examples);
 	
  
  //BoundaryWidth = 6.0*iWidth/2000.0  ; //  measured in pixels ( when iWidth = 2000) 
  //distanceMax = BoundaryWidth*PixelWidth;
  
  
  
  fprintf (stderr," end of setup \n");
	
  return 0;
  
 } 
  
  
  
  
int local_setup(int example_number)
{

 
  	//SetCPlane( center, radius,  DisplayAspectRatio );	
  	SetCPlaneFromExamples(example_number, DisplayAspectRatio );
  	
  	
  	
  	potential_multiplier = 1+log10(radius_0/radius); //  1+example_number;
  	
  	switch(example_number) {
  	
  	
  		case 4:	; {
  			/*
  			 MaxImagePotential  = 255.0000000000000000 
	 		plane  potential_multiplier = 2777.7777777777778283 
	 		black area : potential > potential_boundary = 25.0* potential_multiplier = 111.0924374808178072  = 0.4356566175718345 * MaxImagePotential 
	 		white area : potential > potential_noisy = 10.0 * potential_multiplier = 44.4369749923271229  = 0.1742626470287338 * MaxImagePotential 
	 		*/
  			iColorOfBoundary = iColorOfNoise; // 
  			potential_boundary = 26.0* potential_multiplier;
  			potential_noisy = 24.0 * potential_multiplier;
  			break;
  			}
  	
  		case 15 : {
  		
  			/*
  			 iterMax_pot = 4000 
	 		ER_POT = 100000.0000000000000000 
	 		MaxImagePotential  = inf 
	 		plane  potential_multiplier = 7.3979400086720375 
	 		black area : potential > potential_boundary =  192.3464402254729748  = 0.0000000000000000 * MaxImagePotential 
	 		white area : potential > potential_noisy  = 177.5505602081288998  = 0.0000000000000000 * MaxImagePotential 
	 		*/
  			iColorOfBoundary = 255;
  			iColorOfNoise = 180 ; // 
  			potential_boundary = 250.0; //35.0* potential_multiplier;
  			potential_noisy = 220.0; //26.0 * potential_multiplier;
  			break;
  		
  		
  		}
  		
  		
  		default : {
  	
  			// automatic limits for potential used for zooming
  			
  			potential_boundary = 25.0* potential_multiplier;
  			potential_noisy = 10.0 * potential_multiplier;
  			}
  			
  		} // switch
	
  /* Pixel sizes */
  PixelWidth = (CxMax - CxMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
  PixelHeight = (CyMax - CyMin) / iyMax;
  ratio = ((CxMax - CxMin) / (CyMax - CyMin)) / ((double) iWidth / (double) iHeight);	// it should be 1.000 ...
 
 return 0; 
  
}  
  

 // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





int PrintInfoAboutProgam(int example_number)
{
	
  
  	// display info messages
  	printf ("Numerical approximation of M set for fc(z)= z^2 + c \n");
  	
  	//printf ("iPeriodParent = %d \n", iPeriodParent);
  	//printf ("iPeriodOfChild  = %d \n", iPeriodChild);
  
  	printf ("Image Width = %f in world coordinate\n", CxMax - CxMin);
  	printf ("PixelWidth = %f  = %.16f * Image Width\n", PixelWidth, PixelWidth/ (CxMax - CxMin));
	  
  	//printf("for DEM\n");
  	//if ( distanceMax<0.0 || distanceMax > ER ) printf("bad distanceMax\n");
	//printf("Max distance from exterior to the boundary =  distanceMax = %.16f = %f pixels\n",  distanceMax, BoundaryWidth); 
  	//printf("\n");
  	
  	
  	
  
  
  	// image corners in world coordinate
  	printf ("example number = %d \n", example_number);
  	printf ("plane center c = ( %.16f ; %.16f ) \n", creal (center), cimag (center));
  	printf ("plane radius = %.16f \n", radius);
  	printf ("plane zoom = 1/radius = %.16f \n", 1.0/radius);
  	
  	
  	printf("\n\n potential \n");
  	printf("\t iterMax_pot = %d \n", iterMax_pot);
  	printf("\t ER_POT = %.16f \n" , ER_POT  ); 
  	
  	printf ("\t MaxImagePotential  = %.16f \n", MaxImagePotential );
  	printf ("\t plane  potential_multiplier = %.16f \n", potential_multiplier );
  	printf("\t black area : potential > potential_boundary =  %.16f  = %.16f * MaxImagePotential \n",potential_boundary, potential_boundary / MaxImagePotential);
  	printf("\t white area : potential > potential_noisy  = %.16f  = %.16f * MaxImagePotential \n", potential_noisy, potential_noisy / MaxImagePotential);
  	printf("\n");
  	
  	
  	// center and radius
  	// center and zoom
  	// GradientRepetition
  	printf ("Maximal number of iterations = iterMax = %ld \n", iterMax);
  
  
  	
  	printf("Number of pgm images = %d \n", NumberOfImages);	
  
  	printf ("ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  	//
  	printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
  	// OpenMP version is displayed in the console 
  
  
  
  return 0;
}




// uses global var
int MakeExampleImages(int example_number){


	 	local_setup(example_number);
  
		// make first input image  	  
		Fill_dDataArray(dData1, Potential);
		//find max potential to update potential limiots
		Fill_rgbData_from_dData (dData1, Potential, step_sqrt, rgbData1);
		//Save_PPM(rgbData1, "potentia_step_sqrt", "potentia_step_sqrt", radius); // bad look
	
	
	
		// make second input image  
		Fill_dDataArray(dData2, Normal);
		Fill_rgbData_from_dData (dData2, Normal, linear, rgbData2);
		Save_PPM(rgbData2, "normal_linear", "normal_linear", radius);
	
	
	
		// make 3-rd image - blend image = mix of previous 2 input images
		MakeBlendImage(rgbData1, rgbData2, average, rgbData3);
		Save_PPM(rgbData3, "average", "average blend = (potential + normal)/2", radius);
	
		PrintInfoAboutProgam(example_number);
		
		return 0;





}












int end(){


  fprintf (stderr," allways free memory (deallocate )  to avoid memory leaks \n"); // https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free(dData1);
  free(dData2);
  free(rgbData1);
  free(rgbData2);
  free(rgbData3);
  
  
  return 0;

}













// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main () {



  	setup ();
  	
  	int example_number = 15;
  	
  	
  	//for (  example_number = 0 ;  example_number < ExampleNumberMax; ++ example_number)
  	{
  
  		MakeExampleImages(example_number);
	
		}
	
	
	
  	
  	
  	end();
  
 	
  	

  	return 0;
}
