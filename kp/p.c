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
    
    MPI_Status status;
    int i, j, k, n, m, mx, size, rank, number_of_elements = 0, message_length,
       	lprev, rprev, lnext, rnext, mxt;
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
    
    int *ltemp = (int *)malloc( MX * sizeof(int));
    int *rtemp = (int *)malloc( MX * sizeof(int));
    for(i = 0; i < MX; ++i) 
	fscanf(finput, "%d%d", &ltemp[i], &rtemp[i]);
    
    printf("MX = %d , MY = %d\n", MX, MY);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mx = ((MX - 2) + size - 1) / size;
    mxt = mx;
    if (rank == (size - 1))
    {
        mx = (MX - 2) - (size - 1) * mx;
    }
    mx += 2;

    
    if(mx < 0) {
	mx = 2;
    }
 
    int *l = (int *)malloc(mx * sizeof(int));
    int *r = (int *)malloc(mx * sizeof(int));
    
    int buf_size = (1000 * 2 * (MY-2)*sizeof(double)) + MPI_BSEND_OVERHEAD;
    double* buf = (double *)malloc(buf_size); 
    MPI_Buffer_attach(buf, buf_size);

    if(buf == NULL)
	fprintf(stderr, "Cannot allocate memory\n");
    
    if(rank > 0){
	l[0] = ltemp[rank * mxt - 1 + 1];
	r[0] = rtemp[rank * mxt - 1 + 1];
    } else {
	l[0] = ltemp[0];
	r[0] = rtemp[0];
    }
    if(rank < size - 1) {
	l[mx - 1] = ltemp[(rank + 1) * (mxt) + 1];
	r[mx - 1] = rtemp[(rank + 1) * (mxt) + 1];
    } else {
	l[mx - 1] = ltemp[MX - 1];
	r[mx - 1] = rtemp[MX - 1];
    }
    
    for(i = 1; i < mx - 1; ++i){
	l[i] = ltemp[rank * mxt + i - 1 + 1]; 
	r[i] = rtemp[rank * mxt + i - 1 + 1];
    }
   

    for( i = 0; i < mx; ++i){
	printf("rank = %d, i = %d, l[i] = %d, r[i] = %d\n", rank, i, l[i], r[i]); 
    }
    /*
    for(i = 0; i < mx; ++i)
    {
	fscanf(finput, "%d%d", &l[i], &r[i]);
	number_of_elements += (r[i] - l[i]) + 1;
    }
    printf("l[0] = %d, r[0] = %d \n", l[0], r[0]);
    
    for(j = 0; j < size; j++)
    {
	MPI_Barrier(MPI_COMM_WORLD);
	if(j == rank)
	{
	    if(rank == 0)
		fscanf(finput, "%d%d", &l[0], &r[0]);
	    for(i = 1; i < (mx - 1); i++){
		fscanf(finput, "%d%d", &l[i], &r[i]);
	    }
	    if(rank == size - 1)
		fscanf(finput, "%d%d", &l[mx - 1], &r[mx - 1]);
        }
    }
    fclose(finput);

 	if(rank != 0)
	{
            MPI_Bsend(&l[1], 1, MPI_INT,
		       rank - 1, 100, MPI_COMM_WORLD);
	    MPI_Recv( &l[0], 1, MPI_INT, rank-1,
	     	      MPI_ANY_TAG, MPI_COMM_WORLD,
	              &status );
	}
	if(rank != size - 1)
	{    
	    MPI_Bsend(&l[mx - 2], 1, MPI_INT, rank+1, 100,MPI_COMM_WORLD);
	    MPI_Recv( &l[mx - 1], 1, MPI_INT,
	               rank+1, MPI_ANY_TAG,
		       MPI_COMM_WORLD, &status ); 
        }
	if(rank != 0)
	{
            MPI_Bsend(&r[1], 1, MPI_INT,
		       rank - 1, 100, MPI_COMM_WORLD);
	    MPI_Recv( &r[0], 1, MPI_INT, rank-1,
	     	      MPI_ANY_TAG, MPI_COMM_WORLD,
	              &status );
	}
	if(rank != size - 1)
	{    
	    MPI_Bsend(&r[mx - 2], 1, MPI_INT, rank+1, 100,MPI_COMM_WORLD);
	    MPI_Recv( &r[mx - 1], 1, MPI_INT,
	               rank+1, MPI_ANY_TAG,
		       MPI_COMM_WORLD, &status ); 
        }
	   
    for(i = 0; i < mx; ++i)
    {
	printf("rank=%d, i = %d, l[i]=%d, r[i]=%d\n", rank,i, l[i], r[i]);
    }
    */
    
    printf("After bsend l[0] = %d, r[0] = %d \n", l[0], r[0]);
    double **f = (double**)malloc(mx * sizeof(double*));
    double **df = (double**)malloc(mx * sizeof(double*));
    for(i = 0; i < mx; i++) {
        f[i] = (double*)malloc(((r[i] - l[i]) + 1) * sizeof(double));
	f[i] = &f[i][0] - l[i];
        df[i] = (double*)malloc(((r[i] - l[i]) + 1) * sizeof(double));
	df[i] = &df[i][0] - l[i];
    }
			
    printf("%d of %d is solvingtask on %d bу %d grid\n", rank,  size,  mx,  MY  );	
    fflush(stdout);
   /*  Initial  conditions:  */	

    t1 = MPI_Wtime();
    for (i = 0; i < mx; i++)
    {
        for (j = l[i]; j <= r[i] ; j++)
        {
            f[i][j] = df[i][j] = 0;
	    
            if (((i == 0) && (rank == 0)) || (j == l[i])) f[i][j] = 1.0;
            else if (((i == (mx - 1)) && (rank == (size - 1))) ||
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
	if(rank != 0)
	{
	    message_length = min(r[0] - l[0], r[1] - l[1]) - 1;
            MPI_Bsend(&f[1][max(l[1], l[0]) + 1],
		       message_length, MPI_DOUBLE,
		       rank - 1, 100, MPI_COMM_WORLD);
	    MPI_Recv( &f[0][max(l[0], l[1]) + 1],
		       message_length, MPI_DOUBLE, rank-1,
	     	      MPI_ANY_TAG, MPI_COMM_WORLD,
	              &status );
	}
	if(rank != size - 1)
	{
            message_length = min(r[mx - 1] - l[mx - 1], r[mx - 2] - l[mx - 2]) - 1;
	    MPI_Bsend(&f[mx - 2][max(l[mx - 2], l[mx - 1]) + 1],
			    message_length, MPI_DOUBLE,
			    rank+1, 100,MPI_COMM_WORLD);
	    MPI_Recv( &f[mx - 1][max(l[mx - 1], l[mx - 2]) + 1],
		       message_length, MPI_DOUBLE,
	               rank+1, MPI_ANY_TAG,
		       MPI_COMM_WORLD, &status ); 
        }
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
	    if(j == 0) {
		    
	        fwrite(&f[0][0] + l[0], r[0] - l[0] + 1, sizeof(double), fp);
		for(k = 0; k < l[0]; ++k) 
		    fprintf(txt_fp, "     ");
		for(k = l[0]; k <= r[0]; k++)
		    fprintf(txt_fp, "%.2f ", f[0][k]);
		
	    	fprintf(txt_fp, "\n");
	    }
	    for(i = 1; i < (mx - 1); i++)
	    {
		fwrite(&f[i][0] + l[i], r[i] - l[i] + 1, sizeof(double), fp);
		for(k = 0; k < l[i]; ++k) 
		    fprintf(txt_fp, "     ");
		for(k = l[i]; k <= r[i]; k++)
		    fprintf(txt_fp, "%.2f ", f[i][k]);
		
	    	fprintf(txt_fp, "\n");	
	    }
	    
            if(j == size - 1) {

		fwrite(&f[mx - 1][0] + l[mx - 1], r[mx - 1] - l[mx - 1] + 1, sizeof(double), fp);
		for(k = 0; k < l[mx - 1]; ++k) 
		    fprintf(txt_fp, "     ");
		for(k = l[mx - 1]; k <= r[mx - 1]; k++)
		    fprintf(txt_fp, "%.2f ", f[mx - 1][k]);
		
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
