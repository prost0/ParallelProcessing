#include  <stdio.h>
#include <stdlib.h>
#include  <mpi.h>
/***/
//#define  MX 20
//#define  MY 20
#define  NITER  10000
#define  STEPITER  100
#define MIN_ERR 0.01

//static  float	f[MX][MY];
//static  float  df[MX][MY];
/***/


int  main(int  argc, char  **argv)
{
    
    MPI_Status status;
   /***/
    int i, j, k, n, m, mx, size, rank;
    double t1, t2;
    int MX, MY;
    FILE *fp;
    FILE *txt_fp;
    MPI_Init(&argc, &argv);
    if ( argc != 3 )
    {
	fprintf( stderr, "Usage: %s <nrows> <ncolumns>\n", argv[0] );
	return (-1);
    }
    MX = (int)atol( argv[1]);
    MY = (int)atol( argv[2]);

    
    //float f[MX][MY];
    //float df[MX][MY];
    t1 = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mx = ((MX - 2) + size - 1) / size;
    if (rank == (size - 1))
    {
	//mx = MX - (mx * (size - 1));
        mx = (MX - 2) - (size - 1) * mx;
    }
    mx += 2;
    if(mx < 0) {
	mx = 2;
    }
    double (*f)[MY] = (double(*)[MY])malloc(mx * MY * sizeof(double*));
    double (*df)[MY] = (double(*)[MY])malloc(mx * MY *  sizeof(double*));
    
    int *r_lengthes = (int *)malloc(size * sizeof(int));
    int *r_displacements = (int *)malloc(size * sizeof(int));
    int *s_lengthes = (int *)malloc(size * sizeof(int));
    int *s_displacements = (int *)malloc(size * sizeof(int));
			
    printf("%d of %d is solving heat conduction task on %d bу %d grid\n", rank,  size,  mx,  MY  );	
    fflush(stdout);
   /*  Initial  conditions:  */	
    for (i = 0; i < mx; i++)
    {
        for (j = 0; j < MY; j++)
        {
            f[i][j] = df[i][j] = 0;
            if (((i == 0) && (rank == 0)) || (j == 0)) f[i][j] = 1.0;
            else if (((i == (mx - 1)) && (rank == (size - 1))) || (j == (MY - 1))) f[i][j] = 0.5;
         }
    }
    
    /*Prepare the lengthes and displacements: */
    for (i = 0; i < size; i++)
    {
	if(i==(rank-1))
	{
	    r_lengthes[i] = s_lengthes[i] = MY - 2;
	    r_displacements[i] = 1;
	    s_displacements[i] = MY + 1;
	} else if (i == (rank + 1))
	{
	    r_lengthes[i] = s_lengthes[i] = MY - 2;
	    r_displacements[i] = MY*(mx-1)+1;
	    s_displacements[i] = MY*(mx-2)+1;
	}
	else
            r_lengthes[i] = s_lengthes[i] = r_displacements[i] = 
		    s_displacements[i] = 0;
    }
    
    double *recv_buf = (double *)malloc((size) * sizeof(double));
    /*Iteration  loop:	*/
    double dfs = 99999999.;
    double sumdfs = 999999999.;
    for (n = 0; sumdfs > MIN_ERR; n++)
    {
	dfs = 0;
        if (!(n % STEPITER))
            printf("Iteration %d\n", n);
        /*Do transfers*/
	MPI_Alltoallv(&f[0][0],s_lengthes,
			s_displacements,
			MPI_DOUBLE,
			&f[0][0], r_lengthes,
			r_displacements,
			MPI_DOUBLE,
			MPI_COMM_WORLD);
    	/*Step of calculation starts here : */
	
	for (i = 1; i < (mx - 1); i++)
	{
	    for (j = 1; j < (MY - 1); j++)
	    {
	        df[i][j] = (f[i][j + 1] + f[i][j - 1] + f[i - 1][j] +
			       	f[i + 1][j]) * 0.25 - f[i][j];
		dfs += df[i][j];
	    }
	}
	for (i = 1; i < (mx - 1); i++)
	{
	   for (j = 1; j < (MY - 1); j++)
	   {
	       f[i][j] += df[i][j];
	   }
	}
	MPI_Allgather(&dfs, 1, MPI_DOUBLE, 
			recv_buf, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	sumdfs = 0;
	int pr;
	for(pr = 0; pr < size; ++pr)
	{
	    sumdfs += recv_buf[pr];
	}
    }

    t2 = MPI_Wtime();
 /*  Calculation is done, F  array is а  result: */
    
    if (rank == 0)
    {
        fp = fopen("fp.dat", "w");
	txt_fp = fopen("fp.txt", "w");
	fclose(fp);
	fclose(txt_fp);
    }
   
    
    for(j = 0; j < size; j++)
    {
	MPI_Barrier(MPI_COMM_WORLD);

	if(j == rank)
	{
	    fp = fopen("fp.dat", "a");
	    txt_fp = fopen("fp.txt", "a"); 
	   
	    // fprintf(txt_fp, "rank = %d, size = %d\n", rank, size);
	    for(i = 1; i < (mx - 1); i++)
	    {
		fwrite(f[i] + 1, MY - 2, sizeof(f[0][0]), fp);
		for(k = 1; k < (MY - 1); k++)
		    fprintf(txt_fp, "%.2f ", f[i][k]);
		
	    	fprintf(txt_fp, "\n");	
	    }

	    fclose(fp);
	    fclose(txt_fp);
        }
    }
    
    printf("Time = %f\n", t2 - t1);
    MPI_Finalize();
    return  0;
}
