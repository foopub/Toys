#include <stdio.h>

/*Start by reading some initial conditions from a csv file.
 *Define a simple five point derivative, where possible.
 * 
 * */

enum {size=5,slots=3};		//More slots can be added 
//typedef struct { 	//Use an explicit vector type vs n-array?
//	int x; 		//This might be more efficient, but it also
//	int y;		//might just compile to the same binary?
//} vector_t;	

void derivative (double vectors[slots][size], int n) {
/* n-point derivative for arbitrarily spaced points. Generally:
 *
 * First find the n-degree polynomial using Lagrange interpolation,
 * (S_m and P_m mean the summation and product from i=0 to n
 * and j=0 to n but j!=i ):
 *
 * p(x)=S_i P_j {y_i * (x - x_j)/(x_i - x_j)}
 *
 * which can be differentiated with respect to x (k=0 to n, k!=j):
 *
 * p'(x)= S_i P_j S_k {y_i * ((x - x_j) / (x_i - x_j)) / (x - x_k))}
 * 
 * where the extra k term ALWAYS CANCELS one of the j terms; if x = x_k,
 * there's no problem. This gives a total of n*(n-1) additive terms for 
 * the n point derivative. */

	if (n>size) {		//Check if there's enough points 
		printf("n exceeds the number of points!");
		return;
	}

	double *x = &vectors[0][0];
	double *y = &vectors[1][0];
	double factor, deriv, start, end;	//placeholders
	int mid = (n-1)/2;

	for (int p=0; p<size; p++) {		//the point number
	/* First determine the relative start and end points.
	 * For the end points an asymmetric derivative is needed...
	 * doesn't work for even n yet :( */
		if (p<mid) {			//forward
			start = 0-p;
			end = n-p;
		} else if (size-p-1<mid) {	//backward
			start = size-n-p;
			end = size-p;
		} else { 			//regular
			start = -mid;
			end = mid+1;
		}
		vectors[2][p] = 0, deriv = 0;	//clear old values
//printf("\n");
		for (int i=start; i<end; i++) { 	
//printf("P %d, i %d,\t", p, i);
			if (i!=0) { 
				factor = 1;
			} else {				
				factor = 0;
			}
			for (int j=start; j<end; j++) {	//the product 
				if (j==i) continue;
				if (i==0) {
//printf("\nJ: %d x[p]: %f x[p+j]: %f ", j, x[p], x[p+j]);
				factor += 1/(x[p]-x[p+j]);
				continue;
				}
				factor /= (x[p+i]-x[p+j]);
				if (j==0) continue;
				factor *= (x[p]-x[p+j]);
			}
//printf("Factor is %f\n", factor);
			deriv += y[p+i]*factor;
		}
		vectors[2][p] = deriv;
	}
	return;
}

void readCSV (double vectors[slots][size]) {
	/* add some function to read CSV file here */

}

int main () {
	double vectors[slots][size];
	readCSV(vectors);
	derivative(vectors, 5);
	for (int i=0; i<size; i++) {
		printf("%f ", vectors[2][i]);
		printf("\n");
	}
	return 0;
}
