/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 * @(#)multiply.cc
 * @author Karl Ljungkvist <karl.ljungkvist@it.uu.se>
 *
 */

#include <sys/time.h>
#include <cstdlib>
#include <cstdio>

typedef double number;


double get_time()
{
  struct timeval timeval_time;
  gettimeofday(&timeval_time,NULL);
  return (double)timeval_time.tv_sec + (double)timeval_time.tv_usec*1e-6;
}

template <int dim>
void mult_elem_centric(number *dst, const number *src, const number *a,
                       const number *A, const int intervals)
{
  const int intervals_in_color = intervals >> 1;
  const int N = intervals + 1;
  const int ndofs = 1<<dim;

  number u[ndofs];
  number v[ndofs];

  // loop over all colors
  for (int color1=0; color1 <2; color1++) {
    for (int color2=0; color2 <2; color2++) {
      for (int color3=0; color3 <2; color3++) {


        // loop over all elements
        for (int k1=0; k1 < intervals_in_color; k1++) {
          for (int k2=0; k2 < intervals_in_color; k2++) {
            for (int k3=0; k3 < intervals_in_color; k3++) {

              const int K = (color3 + k3*2) + intervals*((color2 + k2*2) + intervals*(color1 + k1*2));

              //---------------------------------------------------------------
              // element-local operation
              //---------------------------------------------------------------

              // read local dof values
              for (int i1=0; i1 <2; i1++) {
                for (int i2=0; i2 <2; i2++) {
                  for (int i3=0; i3 <2; i3++) {
                    const int global_idx = (i3 + color3 + k3*2) + N*((i2 + color2 + k2*2) + N*(i1 + color1 + k1*2));
                    u[i3+2*(i2+2*i1)] = src[global_idx];
                  }
                }
              }


              for(int i=0; i<ndofs; i++) {
                number tmp = 0;

                for(int j=0; j<ndofs; j++) {
                  tmp += a[i*ndofs+j]*u[j];
                }

                v[i] = tmp*A[K];
              }



              // write back local dof values
              for (int i1=0; i1 <2; i1++) {
                for (int i2=0; i2 <2; i2++) {
                  for (int i3=0; i3 <2; i3++) {
                    const int global_idx = (i3 + color3 + k3*2) + N*((i2 + color2 + k2*2) + N*(i1 + color1 + k1*2));
                    dst[global_idx] += v[i3+2*(i2+2*i1)];
                  }
                }
              }

            }
          }
        }


      }
    }
  }

}

template <int dim>
void mult_node_centric(number *dst, const number *src, const number *a, const number *A, const int intervals)
{
  const int interior_nodes = intervals-1;
  const int N = intervals+1;

  for (int n1=0; n1 < interior_nodes; n1++) {
    for (int n2=0; n2 < interior_nodes; n2++) {
      for (int n3=0; n3 < interior_nodes; n3++) {

        const int my_idx = (n3+1) + N*((n2+1)  + N*(n1+1));


        number tmp = 0;
        // go through all neighbors
        for (int i1=0; i1 <3; i1++) {
          for (int i2=0; i2 <3; i2++) {
            for (int i3=0; i3 <3; i3++) {
              const int K = 0; // FIXME!
              const number a = 1.0; // FIXME!

              const int global_idx = (n3+i3) + N*((n2+i2)  + N*(n1+i1));
              tmp += src[global_idx]*a*A[K];
            }
          }
        }


        dst[my_idx] = tmp;

      }
    }
  }

}

int main(int argc, char *argv[])
{
  const int dim = 3;
  const int nlocdofs = 1<<dim;

  const int Ninterv = 1<<8;
  const int N = 1+Ninterv;

  const int Nelems = Ninterv * (dim>1?Ninterv:1) * (dim>2?Ninterv:1);
  const int Ndofs = N * (dim>1?N:1) * (dim>2?N:1);

  number *u = new number[Ndofs];
  number *v = new number[Ndofs];

  number *A = new number[Nelems];

  number *a = new number[nlocdofs*nlocdofs];

  double t = get_time();
  mult_elem_centric<dim>(u,v,a,A,Ninterv);
  printf("elapsed time: %g s\n",get_time()-t);

  t = get_time();
  mult_node_centric<dim>(u,v,a,A,Ninterv);
  printf("elapsed time: %g s\n",get_time()-t);


  delete[] u;
  delete[] v;
  delete[] A;
  delete[] a;

  return 0;
}
