// Thomas Wilson
// gcc canny.c -o canny -lm
// canny.exe garb34.pgm output1.pgm output2.pgm 1 0
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PICSIZE 256
#define MAXMASK 100

int pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
int edgeflag[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double ival[PICSIZE][PICSIZE];
unsigned char cand[PICSIZE][PICSIZE];

int main(int argc, char **argv) {
    int i,j,p,q,s,t,mr,centx,centy;
    double  maskvalx,maskvaly,sumx,sumy,sig,maxival,minival,maxval,ZEROTOL;
    FILE *fo1, *fo2,*fp1;
    char *foobar;

    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

    argc--; argv++;
    foobar = *argv;
    fo1=fopen(foobar,"wb");

    argc--; argv++;
    foobar = *argv;
    fo2=fopen(foobar,"wb");

    argc--; argv++;
    foobar = *argv;
    sig = atof(foobar);

    argc--; argv++;
    foobar = *argv;
    ZEROTOL = atof(foobar);

    mr = (int)(sig * 3);
    centx = (MAXMASK / 2);
    centy = (MAXMASK / 2);

    char header[256];
    fgets(header, 256, fp1); 
    fgets(header, 256, fp1);
    fgets(header, 256, fp1);

    for (i=0;i<256;i++) {
        for (j=0;j<256;j++) {
            pic[i][j] = getc (fp1);
        }
    }

    for (p = -mr; p <= mr; p++) {
        for (q = -mr; q <= mr; q++) {
            maskvalx = (-(double)q / (sig*sig)) * exp(-(double)(p*p + q*q) / (2.0 * sig * sig));
            maskvaly = (-(double)p / (sig*sig)) * exp(-(double)(p*p + q*q) / (2.0 * sig * sig));

            (maskx[p+centy][q+centx]) = maskvalx;
            (masky[p+centy][q+centx]) = maskvaly;
        }
    }

    for (i = mr; i <= 255 - mr; i++) {
        for (j = mr; j <= 255 - mr; j++) {
            sumx = 0.0;
            sumy = 0.0;
            for (p=-mr;p<=mr;p++) {
                for (q=-mr;q<=mr;q++) {
                    sumx += pic[i+p][j+q] * maskx[p+centy][q+centx];
                    sumy += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
            }
            outpicx[i][j] = sumx;
            outpicy[i][j] = sumy;
        }
    }

    maxival = 0;
    for (i = mr; i < 256 - mr; i++) {
        for (j = mr; j < 256 - mr; j++) {
            ival[i][j]=sqrt((double)((outpicx[i][j] * outpicx[i][j]) + (outpicy[i][j] * outpicy[i][j])));
            if (ival[i][j] > maxival) {
                maxival = ival[i][j];
            }
        }
    }


    for (i = 0; i < PICSIZE; i++) {
        for (j = 0; j < PICSIZE; j++) {
            cand[i][j] = 0;
        }
    }
        
    for (i = mr; i < 256 - mr; i++) {
        for (j = mr; j < 256 - mr; j++) {
            
            double gx = outpicx[i][j];

            if (gx == 0.0) {
                gx = 0.00001;
            }

            double slope = outpicy[i][j]/gx;

            if (slope <= 0.4142 && slope > -0.4142) {
                if (ival[i][j] > ival[i][j-1] && ival[i][j] > ival[i][j+1]) {
                    cand[i][j] = 255;
                }
            } else if (slope <= 2.4142 && slope > .4142) {
                if (ival[i][j] > ival[i-1][j-1] && ival[i][j] > ival[i+1][j+1]) {
                    cand[i][j] = 255;
                }
            } else if (slope <= -0.4142 && slope > -2.4142) {
                if (ival[i][j] > ival[i+1][j-1] && ival[i][j] > ival[i-1][j+1]) {
                    cand[i][j] = 255;
                }
            } else {
                if (ival[i][j] > ival[i-1][j] && ival[i][j] > ival[i+1][j]) {
                    cand[i][j] = 255;
                }
            }
        }
    }
    
    fprintf(fo1, "P5\n");
    fprintf(fo1, "%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo1, "255\n");

    fprintf(fo2, "P5\n");
    fprintf(fo2, "%d %d\n", PICSIZE, PICSIZE);
    fprintf(fo2, "255\n");

    for (i = 0; i < 256; i++) {
        for (j = 0; j < 256; j++) {
            ival[i][j] = (ival[i][j] / maxival) * 255;            
            fprintf(fo1,"%c",(char)((int)(ival[i][j])));
            fprintf(fo2,"%c",(char)((int)(cand[i][j])));
        }
    }

    fclose(fp1);
    fclose(fo1);
    fclose(fo2);
    return 0;
}
