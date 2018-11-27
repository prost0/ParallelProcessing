#include  <stdio.h>
#include <stdlib.h>
/***/
//#define  MX 20
//#define  MY 20
#define  R 512

#define  min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define  max(X, Y) (((X) > (Y)) ? (X) : (Y))
					 
int  main(int  argc, char  **argv)
{
    FILE *f;
    f = fopen("inputbest.txt", "w");
    int i , j;
    fprintf(f, "%d %d\n", R, R);
    for(int i = 0; i < R / 2; ++i)
	  fprintf(f, "%d %d\n", R / 2 - i - 1, R / 2 + i);
    
    for(int i = R / 2 - 1; i >= 0; --i)
	  fprintf(f, "%d %d\n", R / 2 - i - 1, R / 2 + i);
   
    fclose(f);
    return  0;
}
