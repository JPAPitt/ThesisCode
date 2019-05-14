// C file that does all the hard work of the Serre solver


//----------Libraries--------------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//------------Constants-----------
const double i24 = 1.0/24.0;
const double i12 = 1.0/12.0;
const double i3 = 1.0/3.0;
const double i8 = 1.0/8.0;
const double i48 = 1.0/48.0;


//------------Reconstruction Functions-----------


double minmod(double a, double b, double c)
{
	//Minmod function (Defintion in Thesis (3.1) p33)

    if((a > 0) && (b>0) && (c>0))
    {
        return fmin(a,fmin(b,c));
    }
    else if((a < 0) && (b<0) && (c<0))
    {
        return fmax(a,fmax(b,c));
    }
        return 0.0;
}


//------------Matrix Solve Functions-----------
void TDMA(double *a, double *b, double *c, double *d, int n, double *x)
{
	
	// Solve matrix equation Ax = d, where A is a triadiagonal matrix with a, b and c as the lower,middle and upper diagonal respectively. n is the length of vector and size of matrix (n x n)

	//Thomas Algorithm to solve tri-diagonal matrix (https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm))
	//(https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm)
	
	double *alpha = malloc(n*sizeof(double));
	double *beta = malloc(n*sizeof(double));

    alpha[0] = c[0] / b[0];
    beta[0] = d[0] / b[0];


    int i;
    double m;
    for(i=1; i<n-1; i++)
    {
        m = 1.0 / (b[i] - a[i-1]*alpha[i-1]);
        alpha[i] = c[i]*m;
        beta[i] = (d[i] - a[i-1]*beta[i-1])*m;   
        
    }

    m = 1.0 / (b[n-1] - a[n-2]*alpha[n-2]);
    beta[n-1] = (d[n-1] - a[n-2]*beta[n-2])*m;

    x[n-1] = beta[n-1];

    for (i=n-2; i > -1; i--)
    {
        x[i] = beta[i] - alpha[i]*x[i+1];  
   
    }

    free(alpha);
    free(beta);

}

//------------Miscellaneous Functions-----------


void conc(double *a , double *b, double *c,int n,int m ,int k, double *d)
{
	// concatenates list a,b and c into list d. a has size n, b has size m, c has size k and d has size n+m+k
    memcpy(d,a,n*sizeof(double));
    memcpy(d+n,b,m*sizeof(double));
    memcpy(d+n+m,c,k*sizeof(double));
}



//------------Solution Algorithm Functions-----------

void addBC(double *h, double *h0, double *h1, int nBC, int n, int nBCs, double *nh)
{ 
    //Function that adds ghost cells to vector of quantity h (size n). Ghost cells given by h0 (left) and h1 (right) both have length nBCs (used for both velocity and height)
	// nh is new quantity vector with size n + 2*nBC with ghost cells, as only nBC points are needed to evolve system (nBC <= nBCs)
    int j,k;


    double *hb = malloc(nBC*sizeof(double));
    double *he = malloc(nBC*sizeof(double));

    //front end
    //i keeps track of big array
    //j keeps track of small bc arrays
    //k keeps track of small new bc arrays
    j = nBCs-1;

    for(k = nBC -1; k > -1 ; k--)
    {
        hb[k] = h0[j];
        j--;
    }

    //back end
    j = 0;
    for(k = 0; k < nBC ; k++)
    {
        he[k] = h1[j];
        j++;

    }

    //concatenate h and its boundary conditions into nh

    conc(hb,h,he,nBC,n,nBC,nh);

    free(hb);
    free(he);  

}


void evolveh(double *h, double *u, double *pu, double g, double dx, double dt, int nBC, int n,double *nh)
{
    //Solve the conservation of mass equation
	// h,u is the vector of h and u with ghost cells at boundaries (size n + 2 nBC)
	// pu is the vector of u at the previous time step with ghost cells over boundaries
	// g is acceleration due to gravity
	// dx is spatial resolution
	// dt is temporal resolution
	// n is size of vectors without ghost cells
	// nBC is size of ghost cells to left and right

	// nh is h at next time step without ghost cells

    double idx = 1.0 / dx;
    int i = nBC - 1;
    double hiph,himh, uiph,uimh;

    for (i = nBC ; i < n +nBC;i++)
    {
        hiph = 0.5*(h[i+1] + h[i]) - 0.5*(dt*idx)*(u[i+1]*h[i+1] - u[i]*h[i]);
        himh = 0.5*(h[i] + h[i-1]) - 0.5*(dt*idx)*(u[i]*h[i] - u[i-1]*h[i-1]);

		uiph = 0.5*(0.5*(u[i+1] + pu[i+1]) + 0.5*(u[i] + pu[i]));
		uimh = 0.5*(0.5*(u[i] + pu[i]) + 0.5*(u[i-1] + pu[i-1]));

        nh[i -nBC] = h[i] - dt*idx*(hiph*uiph -  himh*uimh); 
 
    }    


   
}
void evolveu(double *h, double *u, double *pu, double g, double dx, double dt, int nBC, int n,double *fu)
{

    //Solve the conservation of momentum equation
	// h,u is the vector of h and u with ghost cells at boundaries (size n + 2 nBC)
	// pu is the vector of u at the previous time step with ghost cells over boundaries
	// g is acceleration due to gravity
	// dx is spatial resolution
	// dt is temporal resolution
	// n is size of vectors without ghost cells
	// nBC is size of ghost cells to left and right

	// fu is u at next time step without ghost cells

    double idx = 1.0 / dx;
    double *a = malloc((n-1)*sizeof(double));
    double *b = malloc(n*sizeof(double));
    double *c = malloc((n-1)*sizeof(double));
    double *f = malloc(n*sizeof(double));

    int i,j;
    double nu,nh,nux,nuxx,nuxxx,nhx,pux,puxx,F,S;


	// Middle of the domain
    for (i =1;i < n-1 ; i++)
    {
        //i for diagonals, j for the BC u,h vectors 
        j = i + nBC;
        nhx = 0.5*idx*(h[j+1] - h[j-1]);
        
        a[i-1] = 0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];
        b[i] = h[j] + 2*i3*idx*idx*h[j]*h[j]*h[j];
        c[i] = -0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];

        nu = u[j];
        nh = h[j];
        nux = 0.5*idx*(u[j+1] - u[j-1]);
        nuxx = idx*idx*(u[j+1] - 2*u[j] + u[j-1]);
        nuxxx = 0.5*idx*idx*idx*(u[j+2] - 2*u[j+1]  + 2*u[j-1] - u[j-2]);
        pux = 0.5*idx*(pu[j+1] - pu[j-1]);
        puxx = idx*idx*(pu[j+1] - 2*pu[j] + pu[j-1]);
        F = nu*nh*nux + g*nh*nhx + nh*nh*nhx*nux*nux + i3*nh*nh*nh*nux*nuxx 
            - nu*nh*nh*nhx*nuxx -i3*nh*nh*nh*nu*nuxxx;

        S = 2*dt*F - pu[j]*nh + nh*nh*nhx*pux + i3*nh*nh*nh*puxx;
        f[i] = -S;

    }

    //Boundaries

	// left
    i = 0;
    j = i + nBC;

    nhx = 0.5*idx*(h[j+1] - h[j-1]);


    b[i] = h[j] + 2*i3*idx*idx*h[j]*h[j]*h[j];
    c[i] = -0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];

    nu = u[j];
    nh = h[j];
    nux = 0.5*idx*(u[j+1] - u[j-1]);
    nuxx = idx*idx*(u[j+1] - 2*u[j] + u[j-1]);
    nuxxx = 0.5*idx*idx*idx*(u[j+2] - 2*u[j+1]  + 2*u[j-1] - u[j-2]);
    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    pux = 0.5*idx*(pu[j+1] - pu[j-1]);
    puxx = idx*idx*(pu[j+1] - 2*pu[j] + pu[j-1]);
    F = nu*nh*nux + g*nh*nhx + nh*nh*nhx*nux*nux + i3*nh*nh*nh*nux*nuxx 
        - nu*nh*nh*nhx*nuxx -i3*nh*nh*nh*nu*nuxxx;

    S = 2*dt*F - pu[j]*nh + nh*nh*nhx*pux + i3*nh*nh*nh*puxx;
    f[i] = -S - u[j-1]*(0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j]);


	// right
    i = n-1;
    j = i + nBC;

    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    a[i-1] = 0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j];
    b[i] = h[j] + 2*i3*idx*idx*h[j]*h[j]*h[j];

    nu = u[j];
    nh = h[j];
    nux = 0.5*idx*(u[j+1] - u[j-1]);
    nuxx = idx*idx*(u[j+1] - 2*u[j] + u[j-1]);
    nuxxx = 0.5*idx*idx*idx*(u[j+2] - 2*u[j+1]  + 2*u[j-1] - u[j-2]);
    nhx = 0.5*idx*(h[j+1] - h[j-1]);
    pux = 0.5*idx*(pu[j+1] - pu[j-1]);
    puxx = idx*idx*(pu[j+1] - 2*pu[j] + pu[j-1]);
    F = nu*nh*nux + g*nh*nhx + nh*nh*nhx*nux*nux + i3*nh*nh*nh*nux*nuxx 
        - nu*nh*nh*nhx*nuxx -i3*nh*nh*nh*nu*nuxxx;

    S = 2*dt*F - pu[j]*h[j] + nh*nh*nhx*pux + i3*nh*nh*nh*puxx;
    f[i] = -S - u[j+1]*(-0.5*idx*h[j]*h[j]*nhx - i3*idx*idx*h[j]*h[j]*h[j]);


	// Solve the generated tridiagonal matrix to get u at next time
    TDMA(a,b,c,f,n,fu);

	// free temporary variables
    free(a);
    free(b);
    free(c); 
    free(f);  
   
}



void evolvewrap(double *u, double *h, double *pubc, double *h0, double *h1, double *u0, double *u1, double g, double dx, double dt, int nBC, int n, int nBCs)
{
	// Wrapper function to solve the Serre equations over a single time step

	// h,u is the vector of h and u without ghost cells at boundaries (size n)
	// pubc is the vector of u at the previous time step with ghost cells over boundaries
	// h0,h1 are the boundary condition of h to the left and right respectively, size (nBCs). 
	// hu,u1 are the boundary condition of u to the left and right respectively, size (nBCs). 
	// g is acceleration due to gravity
	// dx is spatial resolution
	// dt is temporal resolution
	// n is size of vectors without ghost cells
	// nBC is size of ghost cells that will be added to h and u to left and right (nBC <= nBCs)

    double *ubc = malloc((n + 2*nBC)*sizeof(double));
    double *hbc = malloc((n + 2*nBC)*sizeof(double));

    addBC(h, h0, h1,nBC, n,nBCs,hbc);
    addBC(u, u0, u1,nBC, n,nBCs,ubc);
    evolveu(hbc, ubc,pubc,g,dx,dt,nBC,n,u); 
    memcpy(pubc,ubc,(n + 2*nBC)*sizeof(double));


	addBC(u, u0, u1,nBC, n,nBCs,ubc);

	evolveh(hbc, ubc,pubc,g,dx, dt,nBC, n,h); 
    free(ubc);
    free(hbc);

}


//------------Python Interface Functions-----------

double *mallocPy(int n)
{
	//Allow python to allocate C list of doubles with length n
    double *x = malloc(n*sizeof(double));
    return x;
}

void writetomem(double *x, int i , double f)
{
	//Allow python to write to C list x at location i with new double f.
    x[i] = f;

}

double readfrommem(double*x,int i)
{
	//Allow python to read from C list x at location i.
    return x[i];
}

void deallocPy(double *x)
{
	//Allow python to deallocate memory x.
    free(x);
}

