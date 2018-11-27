#include  <stdio.h>
#include <stdlib.h>
/***/
//#define  MX 30
//#define  MY 15
#define  NITER  10000
#define  STEPITER  100

//static  float	f[MX][MY];
//static  float  df[MX][MY];
/***/


int  main(int  argc, char  **argv)
{
    int i, j, n, m, MX, MY;
    FILE *fp;
    FILE *txt_fp;
    
    if ( argc != 3 )
    {
	fprintf(stderr, "Usage: %s <nrows> <ncolumns>\n", argv[0]);
	return (-1);
    }
    
    MX = (int)atol( argv[1] );
    MY = (int)atol( argv[2] );
    

    double **f = (double**)malloc(MX * sizeof(double*));
    double **df = (double**)malloc(MX * sizeof(double*));
    for(i = 0; i < MX; i++) {
	f[i] = (double*)malloc(MY * sizeof(double));
	df[i] = (double*)malloc(MY * sizeof(double));
    }
    
    //float f[MX][MY];
    //float df[MX][MY];

    printf("Solving heat conduction task on %d bу %d grid\n", MX,  MY  );	
    fflush(stdout);
   /*  Initial  conditions:  */	
    for (i = 0; i < MX; i++)
    {
        for (j = 0; j < MY; j++)
        {
            f[i][j] = df[i][j] = 0.0;
            if ((i == 0) || (j == 0)) f[i][j] = 1.0;
            else if ( (i == (MX - 1)) || (j == (MY - 1)) ) f[i][j] = 0.5;
         }
    }
    /*Iteration  loop:	*/
 
    for (n = 0; n < NITER; n++)
    {
        if (!(n % STEPITER))
            printf("Iteration %d\n", n);
            /*EXCHANGE ТНЕ FIRST AND LAST ROWS HERE! ! ! !*/
        /*Step of calculation starts here : */
	for (i = 1; i < (MX - 1); i++)
	{
	   for (j = 1; j < (MY - 1); j++)
	   {
	      df[i][j] = ( f[i][j+1] + f[i][j-1] + f[i-1][j] +
			       f[i+1][j] ) * 0.25 - f[i][j];;
	   }
	}
	for (i = 1; i < (MX - 1); i++)
	{
	   for (j = 1; j < (MY - 1); j++)
	   {
	       f[i][j] += df[i][j];
	   }
	}
    }
 /*  Calculation is done, F  array is а  result: */
    fp = fopen("fn.dat", "w");
    txt_fp = fopen("fn.txt", "w");
    for(i = 1; i < (MX - 1); i++)
    {
        fwrite(f[i] + 1, MY - 2, sizeof(f[0][0]), fp);
    	for(j = 1; j < (MY - 1); j++)
		fprintf(txt_fp, "%.2f ", f[i][j]);
        
	fprintf(txt_fp, "\n");
    }
    
    fclose(fp);
    fclose(txt_fp);
    return  0;
}
