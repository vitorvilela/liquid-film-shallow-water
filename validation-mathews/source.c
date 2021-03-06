/*

Using libvelt.a:
Now we navigate to the crawler directory and build our executable binary. 
There are two ways of doing this; we can do
$ gcc -o crawler crawler.c list.c ../util/libtseutil.a
in which we directly specify the path to the the library file; or we can do
$ gcc -o crawler crawler.c list.c -L../util/ -ltseutil


$ qcc -Wall -O2 source.c -o source -lm ./libvelt.a -I/usr/include/python2.7 -lpython2.7
$ qcc -Wall -O2 source.c -o source -lm
./source

*/


#include <stdlib.h>
#include "saint-venant.h"




// Commom Changed Parameters
// Spray Inclination Angle [deg]
#define sprayAngle 45
// Simulation
#define numberOfLayers 10
#define slipFactor 0.03 // Steel-Ice static friction (wiki)
#define gridLevel 6


// Experiment
// Simple Initial Wet Area Model - area of spray collision 
#define wetAreaDiameter 10.e-3
// Wet Area Radius [m]
#define wetAreaRadius 0.5*wetAreaDiameter
// Collision Area [m²]
#define collisionArea M_PI*wetAreaRadius*wetAreaRadius

//Try to find the Wet Area Radius based on nozzle angle and nozzle-wall distance
//#define wetAreaRadius 100.e-3*tan(6*M_PI/180)


// Domain
// Domain Length [m]
#define domainLenght 35.e-3
// Domain Limits [m]
#define yi -0.5*25.e-3
#define yf 0.5*25.e-3
#define xi 0.
#define xf 35.e-3


// Spray 
// Point of Injection [m]
#define pointOfInjection 0.2*domainLenght
// Liquid Density [kg/m³] - Iso-Octane (Trimetilpentano)
#define rho 690
// Liquid Dynamic Viscosity [Pa.s] - Iso-Octane (Trimetilpentano)
#define dynamicViscosity 542e-6
//Liquid Kinematic Viscosity [m²/s] - Iso-Octane (Trimetilpentano)
#define kinematicViscosity dynamicViscosity/rho
// Period of Injection [s]
#define periodOfInjection 8.45e-3
// Injected Mass [kg]
#define experimentInjectedMass 15.6e-6
// Mass Rate [kg/s]
#define massRate experimentInjectedMass/periodOfInjection
// Volume Rate [m³/s]
#define volumeRate massRate/rho


// Simulation
// Experiment resolution: 1.e-6 m
#define filmResolution 1.e-16 
#define timeStep 1.e-3
#define sprayDuration 8.45e-3
#define totalSimulation sprayDuration + 10.e-3


// Plots
#define printProfileAt 1.e-3



// Film Front Displacement
double xFilmFront = 0.;

const scalar slip[] = slipFactor;

vector velocity[];
scalar dropletDiameter[];

// Film Thickness Rate [m/s] - Experiment: 0.034 m/s
double filmThicknessRate = 0.;
double initialWetArea = 0.;

double wetArea = 0.;
double avgDropletVolume = 0.;
double accumulatedThickness = 0.;  
double accumulatedVolume = 0.;
double avgThickness = 0.;
double currentSprayMass = 0.;
double sprayVolume = 0.;
double actualVolume = 0.;



int main()
{
  
  L0 = domainLenght;
  size (L0);
  origin (-pointOfInjection, -L0/2.);
  G = 9.81;
  dry = 1.e-16;
  N = 1 << gridLevel;
  
  // Number of Shear Layers
  nl = numberOfLayers;

  lambda_b = slip;
  
  // Liquid Viscosity
  nu = kinematicViscosity;
    
  dt = timeStep;
      
  run();
   
}


// EXPLICIT
h[right] = dry;
h[left] = dry;
eta[right] = zb[];
eta[left] = zb[];
u.n[right] = u.x[];
u.n[left] = u.x[];
// IMPLICIT
// q.n[right] = q.x[];
// p[right] = dirichlet(0);


event init (i = 0)
{ 
  
  printf("\n*************** INIT ******************");
  
  initialWetArea = 0.;
 
  int counter = 0;
  
  foreach() { 

    velocity.x[] = 0.;
    velocity.y[] = 0.;
    dropletDiameter[] = 0.;
            
    double r = sqrt(x*x + y*y);
    
    // Inside Initial Wet Area
    if(r < wetAreaRadius)   {
      
      initialWetArea += Delta*Delta;    
      
      // Sauter Mean Diameter [m]      
      dropletDiameter[] = -6.75*r*r + 0.0029*r + 2.74e-4;            
         
      double dropletVolume = (M_PI/6)*pow(dropletDiameter[], 3.);
      
      avgDropletVolume += dropletVolume;
      counter++;
            
      // Fourier modes: sin(x * k * pi), where k is the wave number, must cover the range 0 < x < 1 
      double vel = 0.5 * (10423481112.5*pow(r, 4.) + 15379411.8*pow(r, 3.) - 635221.8*pow(r, 2.) - 161.2*r + 20.5); // * fabs(sin((t/sprayDuration) * 10 * M_PI));
      velocity.x[] = vel*sin(sprayAngle*M_PI/180); // + (1-fabs(y)/wetAreaRadius)*vel*cos(sprayAngle*M_PI/180); // vel*sin(sprayAngle*M_PI/180); // + 0.5*vel*cos(sprayAngle*M_PI/180)*(r/wetAreaRadius)*(x/fabs(x));
      velocity.y[] = 0.; //0.5*vel*cos(sprayAngle*M_PI/180)*(y/fabs(y))*(2*(rand()%100)/100);

    }
    
  }  
   
  avgDropletVolume /= counter;
  
  
  filmThicknessRate = volumeRate/initialWetArea;
    
    
  foreach() {  
              
    double r = sqrt(x*x + y*y);
    
    // Inside Initial Wet Area
    if(r < wetAreaRadius) {     
     
      double dropletVolume = (M_PI/6)*pow(dropletDiameter[], 3.); 
      
      double localThicknessSourceRate = (dropletVolume/avgDropletVolume)*((Delta*Delta)/initialWetArea)*filmThicknessRate; 
         
      // Initial film height
      //h[] = localThicknessSourceRate;
      
      // Mass conservation source
      hs[] = localThicknessSourceRate; 
            
      // Momentum conservation source
      foreach_dimension()
        hus.x[] = hs[]*velocity.x[];
           
    }    
  
  }
  
  
  printf("\nmassRate: %e", massRate);
  printf("\nvolumeRate: %e", volumeRate);
  printf("\nfilmThicknessRate: %e", filmThicknessRate);
  printf("\ndt: %f", dt);
  printf("\ninitialWetArea (%e)", initialWetArea);
 
//   getchar();
  
  printf("\n*************** INIT ******************\n\n");
  
  
}


event stopSource (t=sprayDuration)
{   
    
  printf("\n\n*************** STOP SOURCE ******************\n\n");  
  
  foreach() {
    hs[] = 0.;       
    foreach_dimension()
      hus.x[] = 0;      
  }
    
}



event motion_stats (i++) {
    
    printf("\n*************** STATS ******************");  
    
    double averageFilmVelocity = 0.;
    int counter = 0;
    double min = 1.;
    double max = 0.;
    
    foreach() {        
              
        if (h[] > filmResolution) {
	    counter++;
	    double vel = sqrt(u.x[]*u.x[] + u.y[]*u.y[]);
	    averageFilmVelocity += vel;           
        
	    if (h[] < min) min = h[];
	    if (h[] > max) max = h[];
	    avgThickness += h[];
	                
	    if (x > xFilmFront) xFilmFront = x;  
	}
	
	if (h[] > dry) {
	  actualVolume += (Delta*Delta)*h[]; 
	}
                               
    }
        
    if (counter != 0) {
      averageFilmVelocity /= counter;
      avgThickness /= counter;
      
      printf("\nmin: %e, max: %e, avgThickness: %e", min, max, avgThickness);
      printf("\nAverage Film Velocity: %e m/s", averageFilmVelocity);
      printf("\nX Film Front at: %e", xFilmFront);
    }
    
    printf("\n*************** STATS ******************\n\n");
    
}


event vol_conservation (i=end) {
    
    printf("\n*************** CONSERVATION ******************");
    
    printf("\n\nSimulation Volume: %e", actualVolume);
    printf("\n\nTheoretical Volume: %e", volumeRate*sprayDuration);
    
    printf("\n*************** CONSERVATION ******************\n\n");
    
}

// event initial_profile (i = 0) {
//       
//   static FILE * fp = fopen ("initial-velocity-profile.curve", "w");
//   fprintf (fp, "#Time %e\n", t);
//   
//   foreach_leaf() {
// #if dimension > 1
//     if (fabs(y) < Delta)
// #endif    
//       fprintf (fp, "%e %e \n", x, topNeumann.x[]);
//   }
//   
//   fprintf (fp, "\n");
// 
//   fflush (fp);
//   
// }


event transient (t += printProfileAt; t <= totalSimulation) {
      
  static FILE * fp = fopen ("transient-thickness-profile.curve", "w");
  fprintf (fp, "#Time %e\n", t);
  
  foreach_leaf() {
#if dimension > 1
    if (fabs(y) < Delta)
#endif    
      fprintf (fp, "%e %e \n", x, h[]);
  }
  
  fprintf (fp, "\n");

  fflush (fp);
  
}


event stationary (t = totalSimulation) {
  
  FILE * fp = fopen ("stationary-thickness-profile.curve", "w");
  fprintf (fp, "#Time %e\n", t);
  
  foreach_leaf() {
#if dimension > 1
    if (fabs(y) < Delta)
#endif    
      fprintf (fp, "%e %e \n", x, h[]);
  }
  
  fclose (fp);
  
}



event logfile (i++; t <= totalSimulation) {
   fprintf (ferr, "%g %g\n", t, dt); 
}


event start_distributions (i = 0) { 
  output_ppm (dropletDiameter, file = "d.ppm", n = 1 << 10); 
  output_ppm (h, file = "h.ppm", n = 1 << 10);
  output_ppm (u.x, file = "vx.ppm", n = 1 << 10);
  output_ppm (u.y, file = "vy.ppm", n = 1 << 10);
}


// event distributions (t+=1.e-3; t<=totalSimulation) {   
//   char name[25];
//   sprintf(name, "h%.3f.ppm", t);
//   output_ppm (h, file = name, n = 1 << 10);
//   sprintf(name, "vx%.3f.ppm", t);
//   output_ppm (u.x, file = name, n = 1 << 10);
//   sprintf(name, "vy%.3f.ppm", t);
//   output_ppm (u.y, file = name, n = 1 << 10);
// }

event spray_distributions (t=sprayDuration) {   
  char name[25];
  sprintf(name, "h%f.ppm", t);
  output_ppm (h, file = name, n = 1 << 10);
  sprintf(name, "vx%f.ppm", t);
  output_ppm (u.x, file = name, n = 1 << 10);
  sprintf(name, "vy%f.ppm", t);
  output_ppm (u.y, file = name, n = 1 << 10);
}

event final_distributions (t=totalSimulation) {   
  char name[25];
  sprintf(name, "h%f.ppm", t);
  output_ppm (h, file = name, n = 1 << 10);
  sprintf(name, "vx%f.ppm", t);
  output_ppm (u.x, file = name, n = 1 << 10);
  sprintf(name, "vy%f.ppm", t);
  output_ppm (u.y, file = name, n = 1 << 10);
}



event output_gifs (i += 10; t<=totalSimulation) {
  static FILE * fp1 = popen ("ppm2gif > thickness.gif", "w");  
  output_ppm (h, fp1, n = 1 << 10);
  
  static FILE * fp2 = popen ("ppm2gif > velocity-x.gif", "w");  
  output_ppm (u.x, fp2, n = 1 << 10);
  
  static FILE * fp3 = popen ("ppm2gif > velocity-y.gif", "w"); 
  output_ppm (u.y, fp3, n = 1 << 10);
  
//   scalar l[];
//   foreach()
//     l[] = level;
//   static FILE * fp4 = popen ("ppm2gif > level.gif", "w");  
//   output_ppm (l, fp4);  

}



/*

event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-5}, level+1);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}*/

