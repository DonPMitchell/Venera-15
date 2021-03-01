//
//  Venera15.cpp : Defines the entry point for the console application.
//  D.P. Michell  2016/02/12.
//
#include "stdafx.h"
#include "Venera15.h"
#define     PLOTDAY 19
#define     PLOTMONTH 12
#define     PLOTSIZE 1024

Complex     rgcIQ[256];
char        rgnBadSAR[3200];
char        rgnBadBand[3200];
RadarLook   rgrlSurvey[3200];       // full radar swath, interleaved from both tape recorders
int         rgnRadiometry15[400], rgnRadiometry16[400];
int         rgnAltimetry15[400], rgnAltimetry16[400];
unsigned char    rgnPolar[PLOTSIZE][PLOTSIZE];
unsigned char    rgnSolar[PLOTSIZE][PLOTSIZE];
char        szRecord[80];

//
//  Radar looks are alternately written to two tape records, and they are transmitted as
//  two files that must be interleaved.
//
static void
LoadSurvey()
{
    FILE *pf1, *pf2;
    char szName[512];
    int i, j;
    size_t n1, n2;
    static RadarLook rl, rlZero;

    sprintf(szName, "%s1610_50.dat", DIR);          // Test swath, October 16, 1983
    pf1 = fopen(szName, "rb");
    if (pf1 == 0)
        printf("ERROR: Failed to open %s\n", szName);
    sprintf(szName, "%s1610_51.dat", DIR);
    pf2 = fopen(szName, "rb");
    if (pf2 == 0)
        printf("ERROR: Failed to open %s\n", szName);
    n1 = n2 = 0;
    for (j = 0; j < 1600; j++) {
        rl = rlZero;
        n1 += fread(&rl, sizeof(rl), 1, pf1);
        rgrlSurvey[2*j+0] = rl;
        rl = rlZero;
        n2 += fread(&rl, sizeof(rl), 1, pf2);
        rgrlSurvey[2*j+1] = rl;
    }
    printf("LoadSurvey: %d and %d records\n", n1, n2);
    //printf("  First record loaded, n = %d\n", rgrlSurvey[0].rhHeader.nLine);
    fclose(pf1);
    fclose(pf2);
    //
    //  Mark bad records
    //
    for (j = 0; j < 3200; j++) {
        rgnBadSAR[j] =  rgrlSurvey[j].rgiqSAR[9][0].n == 0 && rgrlSurvey[j].rgiqSAR[9][1].n == 0 && 
                        rgrlSurvey[j].rgiqSAR[9][2].n == 0 && rgrlSurvey[j].rgiqSAR[9][0].n == 0;
        rgnBadBand[j] = rgrlSurvey[j].rgs1.rgn[0][11] == 0 && rgrlSurvey[j].rgs1.rgn[1][11] == 0;
        //if ((j%2 == 0 && j >= n1) || (j%2 == 1 && j >= n2))
        //    rgnBadSAR[j] = rgnBadBand[j] = 1;
    }
    //for (i = 0; i < 33; i++)
    //    rgrlSurvey[0].rhHeader.rgn[i] = 255;
    //rgrlSurvey[0].rgiqSAR[0][0].n = 0;
    pf1 = fopen("survey.raw", "wb");
    fwrite(rgrlSurvey, sizeof(rgrlSurvey), 1, pf1);
    fclose(pf1);
}

//
//  Decode the 8-bit Venera encoding of complex signal sample (4-bits real, 4-bits imaginary)
//
static void
BuildTableIQ()
{
    IQ iq;
    int i, n, nRev;
    Complex c;

    iq.n = 7 + 128;
    printf("%d = %d + %di, size = %d\n",iq.n, iq.re,iq.im, sizeof(iq));
    for (i = 0; i < 16; i++) {
        iq.re = i;
        printf("%f ", (float(iq.re) - 8.0f)/8.0f);
    }
    printf(" IQ decoding\n");
    for (i = 0; i < 256; i++) {
        iq.n = i;
        n = iq.re;
        nRev = n;
        c.re = (float(nRev) - 7.0f)/8.0f;   // -8.0 vs -7.5 seems to make no differenc
        n = iq.im;
        nRev = n;
        c.im = (float(nRev) - 7.0f)/8.0f;
        rgcIQ[i] = c;
    }
}
//
//  Plot altimeter and radiometry data
//
void
Plot(float fLat, float fLong, unsigned char n)
{
    float r, x, y;

    r = float(PLOTSIZE/2) - float(PLOTSIZE/2)*fLat/90.0;
    x = r*cos(fLong*D_DTR);
    y = r*sin(fLong*D_DTR);
    rgnPolar[int(x + PLOTSIZE/2)][int(y + PLOTSIZE/2)] = n;
    rgnPolar[int(x + PLOTSIZE/2)+1][int(y + PLOTSIZE/2)] = n;
    rgnPolar[int(x + PLOTSIZE/2)][int(y + PLOTSIZE/2)+1] = n;
    rgnPolar[int(x + PLOTSIZE/2)+1][int(y + PLOTSIZE/2)+1] = n;
}

void
Plot2(float fLat, float fLong, unsigned char n)
{
    float r, x, y;

    r = float(PLOTSIZE/2) - float(PLOTSIZE/2)*fLat/90.0;
    x = r*cos(fLong*D_DTR);
    y = r*sin(fLong*D_DTR);
    rgnPolar[int(x + PLOTSIZE/2)][int(y + PLOTSIZE/2)] = n;
}

void
Plot3(float xData, float yData, unsigned char n)
{
    float r, x, y;

    x = float(PLOTSIZE/2) + float(PLOTSIZE/2)*xData/1.49597870660e8;
    y = float(PLOTSIZE/2) + float(PLOTSIZE/2)*yData/1.49597870660e8;
    rgnSolar[int(x)][int(y)] = n;
    rgnSolar[int(x)+1][int(y)] = n;
    rgnSolar[int(x)][int(y)+1] = n;
    rgnSolar[int(x)+1][int(y)+1] = n;
}

void
PlotGrid()
{
    float fLat, fLong;
    int i;

    for (i = 0; i < 79; i++)
        szRecord[i] = ' ';
    for (fLat = 0; fLat <= 90.0; fLat += 30.0) {
        for (fLong = 0.0; fLong <= 360.0; fLong += 1.0) {
            Plot2(fLat, fLong, 64);
        }
    }
    for (fLat = 0; fLat <= 90.0; fLat += 1.0) {
        for (fLong = 0.0; fLong < 360.0; fLong += 45.0) {
            Plot2(fLat, fLong, 64);
        }
    }
}
//
//  Load Omega-V radiometry data
//
void
LoadRadiometry()
{
    FILE *pf;
    int nMonth, nDay, nSpacecraft, nX, nY, nKV, nKH;
    float fLat, fLong, fIncidence;
    int n15, n16, i, nTotal15, nTotal16, nBoth;
    double f15, f16;
    static int rgn15[400], rgn16[400];

    n15 = n16 = 0;
    f15 = f16 = 0.0;
    pf = fopen("mpi-radiometry.dat", "rb");
    if (pf == 0)
        printf("Failed to open radiometry data\n");
    while (fscanf(pf, "%d%d%d%d%d%f%f%d%d%f", &nMonth, &nDay, &nSpacecraft, &nX, &nY,
        &fLong, &fLat, &nKV, &nKH, &fIncidence) == 10) {
        i = nDay + 31*nMonth;
        if (nSpacecraft == 15) {
            n15++;
            rgn15[i]++;
            f15 += fIncidence;
            if (nMonth == PLOTMONTH && nDay == PLOTDAY)
                Plot(fLat, fLong + 180.0, 128);
        }
        if (nSpacecraft == 16) {
            n16++;
            rgn16[i]++;
            f16 += fIncidence;
            if (nMonth == PLOTMONTH && nDay == PLOTDAY)
                Plot(fLat, fLong + 180.0, 64);
        }
    }
    nTotal15 = nTotal16 = nBoth = 0;
    for (i = 0; i < 400; i++) {
        rgnRadiometry15[i] = rgn15[i];
        rgnRadiometry16[i] = rgn16[i];
        if (rgn15[i])
            nTotal15++;
        if (rgn16[i])
            nTotal16++;
        if (rgn15[i] && rgn16[i])
            nBoth++;
    }
    printf("Radiometry Data:\n");
    printf("  Venera-15 %10d points, %f angle\n", n15, f15/double(n15));
    printf("  Venera-16 %10d points, %f angle\n", n16, f16/double(n16));
    printf("  Surveys on Venera-15: %4d, Venera-16: %4d,  Both: %4d\n", nTotal15, nTotal16, nBoth);
    fclose(pf);
}
//
//  Load radar altimetry data
//
void
LoadAltimetry()
{
    FILE *pf;
    int i, nOrbit, nData, nSpacecraft, nYear, nMonth, nDay, nSurveys, n15, n16, n;
    float fLat, fLong, fRadius;
    int nTotal15, nTotal16, nBoth;
    static int rgn15[400], rgn16[400];

    pf = fopen("mpi-venus-alt.dat", "rb");
    if (pf == 0)
        printf("Failed to open radio altimetry data\n");
    nSurveys = n15 = n16 = 0;
    while (fscanf(pf, "%d%d VENERA %d%d%d%d", &nOrbit, &nData, &nSpacecraft, &nYear, &nMonth, &nDay) == 6) {
        nSurveys++;
        n = nDay + 31*nMonth;
        for (i = 0; i < nData; i++) {
            fscanf(pf, "%f%f%f", &fLat, &fLong, &fRadius);
            if (nSpacecraft == 15) {
                if (nMonth == PLOTMONTH && nDay == PLOTDAY)
                    Plot(fLat, fLong, 255);
            } else {
                if (nMonth == PLOTMONTH && nDay == PLOTDAY)
                    Plot(fLat, fLong, 190);
            }
        }
        if (nSpacecraft == 15) {
            n15++;
            rgn15[n]++;
        } else {
            n16++;
            rgn16[n]++;
        }
    }
        nTotal15 = nTotal16 = nBoth = 0;
    for (i = 0; i < 400; i++) {
        rgnAltimetry15[i] = rgn15[i];
        rgnAltimetry16[i] = rgn16[i];
        if (rgn15[i])
            nTotal15++;
        if (rgn16[i])
            nTotal16++;
        if (rgn15[i] && rgn16[i])
            nBoth++;
    }
    printf("Radio Altimetry:\n");
    printf("  %4d Surveys, Venera-15: %4d, Venera-16: %4d\n", nSurveys, n15, n16);
    printf("  Surveys on Venera-15: %4d, Venera-16: %4d,  Both: %4d\n", nTotal15, nTotal16, nBoth);
}
//
//  Load KS-18 cosmic ray data
//
void
LoadCosmicRays()
{
    FILE *pf;
    char sz[256];
    int nUTC, n, nYear, nMonth, nDay, nHour;
    float fP30, fP60, fP65, fP11, fP8, fA24, fP70, fP9, fP10, fA22, x, y, z;

    pf = fopen("vg_ks-18-6v.dat", "rb");
    if (pf == 0)
        printf("Failed to open cosmic-ray data\n");
    fgets(sz, 256, pf); // skip header lines
    fgets(sz, 256, pf);
    n = 0;
    while (fscanf(pf, "%d%f%f%f%f%f%f%f%f%f%f%f%f%f", &nUTC, &fP30, &fP60, &fP65, &fP11, &fP8,
        &fA24, &fP70, &fP9, &fP10, &fA22, &x, &y, &z) == 14) {
        n++;
        nHour  = nUTC % 100; nUTC /= 100;
        nDay   = nUTC % 100; nUTC /= 100;
        nMonth = nUTC % 100; nUTC /= 100;
        nYear  = nUTC % 100;
        Plot3(x, y, 255);
    }
    printf("Cosmic Rays: %d telemetry readings\n", n);
    fclose(pf);
}

void
CommonRadiometryAltimetry()
{
    int i;

    printf("Surveys with both radiometry and altimetry:\n");
    printf("  Venera-15: ");
    for (i = 0; i < 400; i++) {
        if (rgnAltimetry15[i] && rgnRadiometry15[i])
            printf("%d/%d, ", i/31, i%31);
    }
    printf("\n  Venera-16: ");
    for (i = 0; i < 400; i++) {
        if (rgnAltimetry16[i] && rgnRadiometry16[i])
            printf("%d/%d, ", i/31, i%31);
    }
        printf("\n");
}
int _tmain(int argc, _TCHAR* argv[])
{
    RadarLook rl;
    FILE *pf;

    PlotGrid();
    LoadRadiometry();
    LoadAltimetry();
    LoadCosmicRays();
    CommonRadiometryAltimetry();
    pf = fopen("AltRad.raw", "wb");
    fwrite(rgnPolar, sizeof(rgnPolar), 1, pf);
    fclose(pf);
    pf = fopen("CosmicRay.raw", "wb");
    fwrite(rgnSolar, sizeof(rgnSolar), 1, pf);
    fclose(pf);
    return 1;

    printf("Record size = %d (should be 4096)\n", sizeof(rl));
    printf("   Header size = %d, offset = %d\n", sizeof(rl.rhHeader), (char *)&(rl.rhHeader) - (char *)&rl);
    printf("   SAR data size = %d, offset = %d\n", sizeof(rl.rgiqSAR), (char *)&rl.rgiqSAR - (char *)&rl);
    LoadSurvey();
    BuildTableIQ();
    ProcessStripeImages();
    InitializeOrbitalData();
    ProcessFrameImages();
	return 1;
}

