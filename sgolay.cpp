/*
Implements the sgolay filter design. 
B = SGOLAY(K,F) designs a Savitzky-Golay (polynomial) FIR smoothing
filter B.  The polynomial order, K, must be less than the frame size,
F, and F must be odd.  

Note that if the polynomial order K equals F-1, no smoothing will occur.

Author: Olalekan Ogunmolu

*/
#include <cmath>
#include <vector>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

using namespace std;

/*Function Prototypes*/
void sgolay(int k, int F);
void vander(const emxArray_real_T *v, emxArray_real_T *A);


/*
void sgolay(k, F)
{
	//compute projection matrix B
	void vander(const double v[5], double A[25]);
	
}
*/
void vander(const emxArray_real_T *v, emxArray_real_T *A)
{
  int jtilecol;
  int ibtile;
  int k;

  for (jtilecol = 0; jtilecol < 5; jtilecol++) 
  {
    ibtile = jtilecol * 5;
    for (k = 0; k < 5; k++) 
    {
      A[ibtile + k] = v[k];
    }
  }

  for (jtilecol = 0; jtilecol < 5; jtilecol++) 
  {
    A[20 + jtilecol] = 1.0;
  }

  for (jtilecol = 0; jtilecol < 5; jtilecol++) 
  {
    for (k = 0; k < 4; k++) 
    {
      A[(jtilecol - (k + 1) * 5) + 20] *= A[(jtilecol - k * 5) + 20];
    }
  }
  cout << "A: " << A << endl;
}

int main(int argc, char** argv)
{
  double v[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
  double A[25];
	cout << "v: " << v << endl;
	vander(v, A);
}
/*
typedef struct emxArray_real_T
{
    double *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
} emxArray_real_T;
*/