#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef pi
#define pi 3.141593
#endif

#ifndef EPSILON
#define EPSILON 0.00000001
#endif


/*
 * Hammersley point sets. Deterministic and look random.
 * Base p = 2, which is especially fast for computation.
 * https://www.cse.cuhk.edu.hk/~ttwong/papers/udpoint/udpoints.html
 */
void SphereHammersley(float *xyz_sphere, int n)
{
  float p, t, st, phi, phirad;
  int k, kk, pos;

  for (k=0, pos=0 ; k<n ; k++)
  {
    t = 0;
    for (p=0.5, kk=k ; kk ; p*=0.5, kk>>=1)
      if (kk & 1)                           // kk mod 2 == 1
	t += p;
    t = 2.0 * t  - 1.0;                     // map from [0,1] to [-1,1]

    phi = (k + 0.5) / n;                    // a slight shift
    phirad =  phi * 2.0 * pi;             // map to [0, 2 pi)

    st = sqrt(1.0-t*t);
    xyz_sphere[pos] = st * cos(phirad);
    xyz_sphere[pos+n] = st * sin(phirad);
    xyz_sphere[(pos++)+2*n] = t;
  }
}

/*
 * Fibbonnacci points
 * https://people.sc.fsu.edu/~jburkardt/c_src/sphere_fibonacci_grid/sphere_fibonacci_grid.html
 */

void SphereFibbonnacci(float *xyz_sphere, int n, int it)
{
  double cphi, i_r8, ng_r8, r8_phi, sphi, theta;
  const double r8_pi = 3.141592653589793;
  int sampling_scale = 2;
  int j, pos;

  if(it==1)//means that the image is a fisheye, we need to change the sampling to get half a sphere
  {
    sampling_scale=1;
  }

  r8_phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;
  ng_r8 = ( double ) ( n );

  for ( j = 0, pos=0; j < n; j++ )
  {
    i_r8 = (double)(-n+1+sampling_scale*j);
    theta = 2.0 * r8_pi * i_r8 / r8_phi;
    sphi = i_r8 / ng_r8;
    cphi = sqrt ( ( ng_r8 + i_r8 ) * ( ng_r8 - i_r8 ) ) / ng_r8;

    xyz_sphere[pos] = cphi * sin ( theta );
    xyz_sphere[pos+n] = cphi * cos ( theta );
    xyz_sphere[(pos++)+2*n] = sphi;
  }
}


void xyz_3Dsphere_to_xy_2Dplanar(float *xyz_sphere, int n, int h, int w, float* x, float *y) {

    float x3, y3, z3, phic, thetac;

    for (int i=0; i<n; i++) {
        x3     = xyz_sphere[i];
        y3     = xyz_sphere[i+n];
        z3     = xyz_sphere[i+2*n];
        phic   = acos(z3);
        thetac = atan2(y3,x3);
   	    y[i]   = ((phic*h/(pi)))-1;
        x[i]   = ((thetac)*w/(2*pi))-1;

        if (x[i]<0)
            x[i] = w+x[i];

    }
}

void seeds_sp_sampling_hammersley(int n, int h, int w, float* x, float *y) {

    float *xyz_sphere = (float *)malloc(n*3*sizeof(float));

    SphereHammersley(xyz_sphere, n);
    xyz_3Dsphere_to_xy_2Dplanar(xyz_sphere, n, h, w, x, y);

    free(xyz_sphere);
}

void seeds_sp_sampling_fibbonnacci(int n, int h, int w, int it, float* x, float *y) {

    float *xyz_sphere = (float *)malloc(n*3*sizeof(float));

    SphereFibbonnacci(xyz_sphere, n, it);
    xyz_3Dsphere_to_xy_2Dplanar(xyz_sphere, n, h, w, x, y);

    free(xyz_sphere);
}
