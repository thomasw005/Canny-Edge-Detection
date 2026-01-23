// Thomas Wilson
// gcc canny.c -o canny -lm
// canny.exe garb34.pgm output1.pgm output2.pgm output3.pgm 1 0.05
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Constants
#define PICSIZE 256
#define MAXMASK 100

// Globals
int pic[PICSIZE][PICSIZE];
double outpicx[PICSIZE][PICSIZE];
double outpicy[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double ival[PICSIZE][PICSIZE];
unsigned char cand[PICSIZE][PICSIZE];
unsigned char finaledge[PICSIZE][PICSIZE];
unsigned char peaks[PICSIZE][PICSIZE];

int main(int num_args, char **args) {
    FILE *fo1, *fo2, *fo3, *fp1;
    char *temp;

    // Read in arguments.
    num_args--; args++;
    temp = *args;
    fp1=fopen(temp,"rb");

    num_args--; args++;
    temp = *args;
    fo1=fopen(temp,"wb");

    num_args--; args++;
    temp = *args;
    fo2=fopen(temp,"wb");

    num_args--; args++;
    temp = *args;
    fo3=fopen(temp,"wb");

    num_args--; args++;
    temp = *args;
    double sig = atof(temp);

    num_args--; args++;
    temp = *args;
    double percentage = atof(temp);

    int hist[256];
    int mask_radius = (int)(sig * 3);
    int centx = (MAXMASK / 2);
    int centy = (MAXMASK / 2);
    
    // Get rid of header.
    char header[256];
    fgets(header, 256, fp1); 
    fgets(header, 256, fp1);
    fgets(header, 256, fp1);
    
    // Initialize arrays.
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            pic[i][j] = getc(fp1);

            outpicx[i][j] = 0.0;
            outpicy[i][j] = 0.0;
            ival[i][j]    = 0.0;
            cand[i][j]    = 0;
            peaks[i][j]   = 0;
            finaledge[i][j]=0;
            hist[i] = 0;
        }
    }

    // Compute gaussian.
    double maskvalx, maskvaly;
    for (int p = -mask_radius; p <= mask_radius; p++) {
        for (int q = -mask_radius; q <= mask_radius; q++) {
            maskvalx = (-(double)q / (sig*sig)) * exp(-(double)(p*p + q*q) / (2.0 * sig * sig));
            maskvaly = (-(double)p / (sig*sig)) * exp(-(double)(p*p + q*q) / (2.0 * sig * sig));

            (maskx[p+centy][q+centx]) = maskvalx;
            (masky[p+centy][q+centx]) = maskvaly;
        }
    }

    // Convolve with gaussian.
    double sumx, sumy;
    for (int i = mask_radius; i <= 255 - mask_radius; i++) {
        for (int j = mask_radius; j <= 255 - mask_radius; j++) {
            sumx = 0.0;
            sumy = 0.0;
            for (int p = -mask_radius; p <= mask_radius; p++) {
                for (int q = -mask_radius; q <= mask_radius; q++) {
                    sumx += pic[i+p][j+q] * maskx[p+centy][q+centx];
                    sumy += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
            }
            outpicx[i][j] = sumx;
            outpicy[i][j] = sumy;
        }
    }

    // Compute magnitude of gradient.
    double maxival = 0;
    for (int i = mask_radius; i < 256 - mask_radius; i++) {
        for (int j = mask_radius; j < 256 - mask_radius; j++) {
            ival[i][j]=sqrt((double)((outpicx[i][j] * outpicx[i][j]) + (outpicy[i][j] * outpicy[i][j])));
            if (ival[i][j] > maxival) {
                maxival = ival[i][j];
            }
        }
    }

    // Scale magnitudes to 0-255
    for (int i = mask_radius; i < 256 - mask_radius; i++) {
      for (int j = mask_radius; j < 256 - mask_radius; j++) {
        ival[i][j] = (ival[i][j] / maxival) * 255.0;
      }
   }
    
   // Find peaks.
    for (int i = mask_radius; i < 256 - mask_radius; i++) {
        for (int j = mask_radius; j < 256 - mask_radius; j++) {
            
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

    // Auto Threshold
    int HI = 0;
    double LO = 0;

    // Build histogram from magnitudes.
    for (int i = mask_radius; i < 256 - mask_radius; i++) {
        for (int j = mask_radius; j < 256 - mask_radius; j++) {
            int mag = (int)(ival[i][j] + 0.5);
            if (mag < 0) mag = 0;
            if (mag > 255) mag = 255;
            hist[mag]++;
        }
    }

    // Compute cutoff from percentage*rows*cols.
    int numPixels = (256 - 2*mask_radius) * (256 - 2*mask_radius);
    int cutoff = (int)(percentage * numPixels);  
    int areaOfTops = 0;

    // Loop downwards until we reach cutoff.
    int hival = 255;
    for (hival = 255; hival >= 1; hival--) {
        areaOfTops += hist[hival];
        if (areaOfTops > cutoff) break;
    }
    HI = hival;
    LO = 0.35 * HI;

    // Save state of peaks for later writing to file.
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            peaks[i][j] = cand[i][j];
        }
    }
      
    // Apply Double Threshold
    for (int i = mask_radius; i < 256 - mask_radius; i++) {
        for (int j = mask_radius; j < 256 - mask_radius; j++) {
            if (cand[i][j] == 255) {
                if (ival[i][j] > HI) {
                    cand[i][j] = 0;
                    finaledge[i][j] = 255;
                }
                else if (ival[i][j] < LO) {
                    cand[i][j] = 0;
                    finaledge[i][j] = 0;
                }
            }
        }
    }

    int moretodo = 1;
    while (moretodo) {
        moretodo = 0;
        for (int i = mask_radius; i < 256 - mask_radius; i++) {
            for (int j = mask_radius; j < 256 - mask_radius; j++) {
                if (cand[i][j] == 255) {
                    for (int p = -1; p <= 1; p++) {
                        for (int q = -1; q <= 1; q++) {
                            if (finaledge[i+p][j+q] == 255) {
                                cand[i][j] = 0;
                                finaledge[i][j] = 255;
                                moretodo = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Write headers to file outputs.
    fprintf(fo1, "P5\n%d %d\n255\n", PICSIZE, PICSIZE);
    fprintf(fo2, "P5\n%d %d\n255\n", PICSIZE, PICSIZE);
    fprintf(fo3, "P5\n%d %d\n255\n", PICSIZE, PICSIZE);

    // Write values to files.
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            fprintf(fo1,"%c",(char)((int)(ival[i][j])));
            fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
            fprintf(fo3, "%c", (unsigned char)finaledge[i][j]);
        }
    }

    // Close files.
    fclose(fp1);
    fclose(fo1);
    fclose(fo2);
    fclose(fo3);
    return 0;
}
