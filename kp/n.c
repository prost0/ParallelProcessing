#include  <stdio.h>
#include <stdlib.h>
#include  <mpi.h>
/***/
//#define  MX 20
//#define  MY 20
#define  NITER  10000
#define  STEPITER  100

#define  min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define  max(X, Y) (((X) > (Y)) ? (X) : (Y))
					 
int  main(int  argc, char  **argv)
{
    
    int i, j, k, n, m, mx,
       	lprev, rprev, lnext, rnext;
    double t1, t2;
    int MX, MY;
    FILE *fp;
    FILE *txt_fp;
    
    if ( argc != 2 )
    {
	fprintf( stderr, "Usage: %s <nrows> <ncolumns>\n", argv[0] );
	    return (-1);
    }
    
    FILE *finput;
    printf("Before scanf\n");
    printf("File name = %s\n", argv[1]);
    finput = fopen(argv[1], "r");
    fscanf(finput, "%d%d", &MX, &MY);
    
    printf("MX = %d , MY = %d\n", MX, MY);
    MPI_Init(&argc, &argv);
    mx = MX;
 
    int *l = (int *)malloc(mx * sizeof(int));
    int *r = (int *)malloc(mx * sizeof(int));    
    
    for(i = 0; i < mx; i++){
	fscanf(finput, "%d%d", &l[i], &r[i]);
    }    
    
    fclose(finput);

    double **f = (double**)malloc(mx * sizeof(double*));
    double **df = (double**)malloc(mx * sizeof(double*));
    for(i = 0; i < mx; i++) {
        f[i] = (double*)malloc(((r[i] - l[i]) + 1) * sizeof(double));
	f[i] = &f[i][0] - l[i];
        df[i] = (double*)malloc(((r[i] - l[i]) + 1) * sizeof(double));
	df[i] = &df[i][0] - l[i];
    }
			
    printf("solving task on %d bу %d grid\n",  mx,  MY  );	
    fflush(stdout);
   /*  Initial  conditions:  */	

    t1 = MPI_Wtime();
	    
    for (i = 0; i < mx; i++)
    {
        for (j = l[i]; j <= r[i] ; j++)
        {
            f[i][j] = df[i][j] = 0;
	    
            if ((i == 0)  || (j == l[i])) f[i][j] = 1.0;
            else if ((i == (mx - 1)) ||
			    (j == r[i])) f[i][j] = 1.0;
         }
	
	 if(i > 0 && i < (mx - 1))
	 {
	     for(k = l[i]; k <= l[i - 1]; ++k)
	         f[i][k] = 1;
	     for(k = r[i - 1]; k <= r[i]; ++k)
	       	f[i][k] = 1;
	     for(k = l[i]; k <= l[i + 1]; ++k)
	 	f[i][k] = 1;
	     for(k = r[i + 1]; k <= r[i]; ++k)
	  	f[i][k] = 1;	
	 }
    }

    
    /*Iteration  loop:	*/
 
    for (n = 0; n < NITER; n++)
    {
        if (!(n % STEPITER))
            printf("Iteration %d\n", n);
            /*EXCHANGE ТНЕ FIRST AND LAST ROWS HERE! ! ! !*/
	    /*Step of calculation starts here : */
	for (i = 1; i < (mx - 1); i++)
	{
	    for (j = max(max(l[i], l[i - 1]), l[i + 1]) + 1; j < min(min(r[i - 1], r[i]), r[i + 1]); j++)
	    {	
		if(f[i][j] != 1)
	            df[i][j] = (f[i][j + 1] + f[i][j - 1] +
			        f[i - 1][j] + f[i + 1][j]) * 0.25 - f[i][j];
	    }
	}
	for (i = 1; i < (mx - 1); i++)
	{
	   for (j = max(max(l[i], l[i - 1]), l[i + 1]) + 1; j < min(min(r[i - 1], r[i]), r[i + 1]); j++)
	   {
		if(f[i][j] != 1)
	            f[i][j] += df[i][j];
	   }
	}
    }
    
    
    t2 = MPI_Wtime();
 /*  Calculation is done, F  array is а  result: */
    fp = fopen("fn.dat", "w");
    txt_fp = fopen("fn.txt", "w");

    for(i = 0; i < mx; i++)
    {
	printf("Write %d\n", i);
	fwrite(f[i] + l[i], r[i] - l[i] + 1, sizeof(double), fp);
	
	for(k = 0; k < l[i]; ++k) 
	    fprintf(txt_fp, "     ");
	
	for(k = l[i]; k <= r[i]; k++)
	    fprintf(txt_fp, "%.2f ", f[i][k]);
	
    	fprintf(txt_fp, "\n");	
    }
	    
    fclose(fp);
    fclose(txt_fp);
    
    printf("Time = %f\n", t2 - t1);
    MPI_Finalize();
    return  0;
}
