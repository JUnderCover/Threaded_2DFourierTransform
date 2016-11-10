// Threaded two-dimensional Discrete FFT transform
// Jai Chauhan
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

#define SWAP(a,b) tempr = (a); (a) = (b); (b) = tempr

Complex* data;
int w;
unsigned long w_ulong; 
int h;
int N; // number of points in transform
int numthreads = 16;
int rowseach;
int numelements;
unsigned long myrowlength;
int count;
int start;
int end = 0;
pthread_mutex_t startMutex;
pthread_mutex_t  endMutex;
pthread_cond_t exitCond;
pthread_mutex_t countMutex;
pthread_mutex_t exitMutex;
// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;


void fastTransform(int s, unsigned long len) 
{
  double* myrow = new double[len * 2];
  int counter = 0;
  for(int i = 0; i < len*2; i+=2)
    {
      myrow[i] = data[s + counter].real;
      myrow[i+1] = data[s + counter].imag;
      counter++; 
    }
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  n = len<<1;
  j = 1;
  for(i=1;i<n;i+=2)
    {
      if(j>i) 
	{
	  SWAP(myrow[j-1],myrow[i-1]);
	  SWAP(myrow[j],myrow[i]);
	}
      m = n>>1;
      while((m>=2)&&(j>m))
	{
	  j = j - m;
	  m = m >> 1;
	}
      j = j + m;
    }
  mmax = 2;
  while(n>mmax)
    {
      istep = mmax<<1;
      theta = (-2*M_PI/mmax);
      wtemp = sin(0.5 * theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for(m=1;m<mmax;m+=2)
	{
	  for(i=m; i<=n; i+=istep)
	    {
	      j = i + mmax;
	      tempr = wr*myrow[j-1] - wi*myrow[j];
	      tempi = wr*myrow[j] + wi*myrow[j-1];
	      myrow[j-1] = myrow[i-1] - tempr;
	      myrow[j]  = myrow[i] - tempi;
	      myrow[i-1] = myrow[i-1] + tempr;
	      myrow[i] = myrow[i] + tempi;
	    }
	  wtemp = wr;
	  wr = wr*wpr-wi*wpi+wr;
	  wi = wi*wpr+wtemp*wpi+wi;
	}
      mmax = istep;
    }
  double real; 
  double imag;
  counter = 0;
  for(int i = 0; i < len * 2; i+=2)
    {
      real = myrow[i];
      imag = myrow[i+1];
      data[s + counter] = Complex(real,imag);
      counter++;
    }
}
// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.


//I have no clue what this does, see my own implementation 
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}


// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main
void MyBarrier_Init()// you will likely need some parameters)
{
}

// Each thread calls MyBarrier after completing the row-wise DFT
void MyBarrier() // Again likely need parameters
{
}
                    
void Transform1D(Complex* h, int N)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
}

void* Transform2DThread(void* v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  unsigned long  threadnum = (unsigned long)v;
  pthread_mutex_lock(&startMutex);
  for(int i = 0; i < rowseach; i++)
    {
      start = threadnum * numelements + i * w; 
      fastTransform(start,w_ulong);
    }
  pthread_mutex_unlock(&startMutex);
  pthread_mutex_lock(&countMutex); 
  count--;
  if (count == 0)
    {
      count = numthreads;
      pthread_mutex_unlock(&countMutex);
      pthread_mutex_lock(&exitMutex);
      pthread_cond_signal(&exitCond);
      pthread_mutex_unlock(&exitMutex);
    }
  else
    {
      pthread_mutex_unlock(&countMutex);
    }
  return 0;
}

void Transform2D(const char* inputFN)
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Create the global pointer to the image array data
  // Create 16 threads
  // Wait for all threads complete
  // Write the transformed data
  //initialize globals
  data = image.GetImageData();
  w = image.GetWidth();
  w_ulong = (unsigned long) w;
  h = image.GetHeight();
  N = w * h;
  rowseach = h / numthreads;
  numelements = rowseach * w;
  myrowlength = (unsigned long) numelements;
  pthread_cond_init(&exitCond,0);
  pthread_mutex_init(&exitMutex,0);
  pthread_mutex_init(&startMutex,0);
  pthread_mutex_init(&endMutex,0);
  pthread_mutex_init(&countMutex,0);
  count = numthreads;
  pthread_mutex_lock(&exitMutex);
  for (int i = 0; i <numthreads; i++)
    {
      pthread_t pt;
      printf("Creating thread %d\n",i);
      pthread_create(&pt,0,Transform2DThread,(void*)i);
    }
    pthread_cond_wait(&exitCond,&exitMutex);
    cout<<"1D transform complete" << endl;
    image.SaveImageData("MyAfter1D.txt",data,w,h);
    //transpose data array
    int n = sqrt(w*h);
    Complex temp;
    for(int i = 0; i < n; i++)
      {
	for(int j = i + 1; j < n; j++)
	  {
	    temp = data[i*n + j];
	    data[i*n+j] = data[j*n+i];
	    data[j*n+i] = temp;
	  }
      }
    //reinitialize
    pthread_cond_init(&exitCond,0);
    pthread_mutex_init(&exitMutex,0); 
    pthread_mutex_lock(&exitMutex);
    //recreate threads
    for (int i = 0; i <numthreads; i++)
      {
	pthread_t pt;
	printf("Creating thread %d\n",i);
	pthread_create(&pt,0,Transform2DThread,(void*)i);
      }
    pthread_cond_wait(&exitCond,&exitMutex);
    cout << "Done with 2d" << endl;
    //retranspose you foolx
    for(int i = 0; i < n; i++)
      {
	for(int j = i + 1; j < n; j++)
	  {
	    temp = data[i*n + j];
	    data[i*n+j] = data[j*n+i];
	    data[j*n+i] = temp;
	  }
      }
    image.SaveImageData("MyAfter2D.txt",data,w,h);
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
