/**
*Remember to mount Vitor directory first in order to find qcc
qcc -Wall -O2 source.c -o source -lm

Using libvelt.a:
Now we navigate to the crawler directory and build our executable binary. There are two ways of doing this; we can do
$ gcc -o crawler crawler.c list.c ../util/libtseutil.a
in which we directly specify the path to the the library file; or we can do
$ gcc -o crawler crawler.c list.c -L../util/ -ltseutil


$ qcc -Wall -O2 source.c -o source -lm ./libvelt.a -I/usr/include/python2.7 -lpython2.7
./source
*/

#include "saint-venant.h"
#include "velt.h"



#define LEVEL 8

#define numberOfLayers 5

// Initial film height [m]
double h0 = 1000*7.7e-9;

// Wet area
double wet_area = 4.9e-4;
double yi = -0.014;
double yf = 0.014;
double xi = 0.01;
double xf = 0.045;

// Extended domain
//double domain_lenght = 10*xf; 

// Volumetric flux / Wet area 
// Film height per second [m/s]
double flux = 1000*1.5e-6;

vector topNeumann[];

int main()
{
  
  
  L0 = 2*0.045; // domain_lenght
  size (L0);
  origin (0, - L0/2.);
  G = 9.81;
  N = 1 << LEVEL;
  
  // Number of shear layers
  nl = numberOfLayers;
  // Water dynamic viscosity
  nu = 1.e-3;
  
  dut = topNeumann;
    
  run();
 
  
}


event init (i = 0)
{   
    
  foreach() { 
    topNeumann.x[] = 0.;
    topNeumann.y[] = 0.;
            
    if(x-0.0097+0.47*y-210.15*y*y>0 && x<xf) { 
      h[] = 0.05*h0;
      topNeumann.x[] = get_velocity(y, x);
    } 
  
  }
}


// Ponderate the flux based on each cell area (quadtree)
event source (i++)
{
    
  vector lastLayerVel = ul[numberOfLayers-1];   
    
  //int flag = 1;
  foreach() {    
    if(x-0.0097+0.47*y-210.15*y*y>0 && x<xf) {
      double beta = (Delta*Delta)/wet_area;
      h[] += flux*dt*beta;           
      //if(flag == 1) fprintf (ferr, "%e %e\n", h[], flux*dt*beta); 
      //flag = 0;
    }
    else
      topNeumann.x[] = lastLayerVel.x[]; 
  }  
}


event profiles (t += 10) {
  
  static FILE * fp = fopen ("transient-height-profile.curve", "w");
  fprintf (fp, "#Time %e\n", t);
  
  foreach_leaf() {
#if dimension > 1
    if (fabs(y) < Delta && x < xf)
#endif    
      fprintf (fp, "%e %e \n", x, h[]);
  }
  
  fprintf (fp, "\n");

  fflush (fp);
  
}


event stationary (t = 60) {
  
  FILE * fp = fopen ("stationary-height-profile.curve", "w");
  fprintf (fp, "#Time %e\n", t);
  
  foreach_leaf() {
#if dimension > 1
    if (fabs(y) < Delta)
#endif    
      fprintf (fp, "%e %e \n", x, h[]);
  }
  
  fclose (fp);
  
}


event logfile (i++; t <= 60) {
   fprintf (ferr, "%g %g\n", t, dt); 
}



event vel_bc_distribution (i = 0) {
  output_ppm (topNeumann.x, min = 0., max = 70., file = "vel.ppm", n = 1 << LEVEL);
}


event output (i += 20) {
  static FILE * fp = popen ("ppm2gif > source.gif", "w");
  output_ppm (h, fp, min = 0, max = 4*h0);
}

