//
//  Process synthetic-aperture radar data in the primary frame survey of Venera-15
//  D.P. Mitchell  2016/02/12.
//
#pragma implicit sqrt(), cos(), sin()
#include "stdafx.h"
#include "Venera15.h"
#define D_2PI           6.28318530717958647692528676655900576839433879875021
#define DOPPLERSHIFT    (D_2PI/double(20))

static Complex          rgcHologram[3200][20][127];
static Complex          rgcDoppler[3200][20][127];
static Complex          rgcDoppler2[3200][20][127];
static Complex          rgcDoppler3[3200][40][127];
static Complex          rgcShift[3200][20][127];
static Complex          rgcPulse[3200][20][127];
static Complex          rgcPulsePhase[3200][20][127];
static Complex          rgcPhase[3200][20][127];
static Complex          rgcPhase2[3200][20][127];
static Complex          rgcPhase3[3200][40][127];
static Complex          rgcCorrelated[3200][20][127];

static double           rgfSynth2[3220][127];
static double           rgfSynth3[6440][127];
static double           rgfWeight[6500][127];
static double           rgfReal3[3200][127];

static unsigned char    rgnHologram[3200][3200];
static unsigned char    rgnDoppler[3200][3200];
static unsigned char    rgnDoppler2[3200][3200];
static unsigned char    rgnDoppler3[3200][6000];    // 40 doppler channels
static unsigned char    rgnShift[3200][3200];
static unsigned char    rgnPulse[3200][3200];
static unsigned char    rgnPhase[3200][3200];
static unsigned char    rgnPhase2[3200][3200];
static unsigned char    rgnPhase3[3200][6000];
static unsigned char    rgnCorrelated[3200][3200];
static unsigned char    rgnSynth2[3220][127];
static unsigned char    rgnSynth3[6440][127];
static unsigned char    rgnReal3[3200][127];
static unsigned char    rgnFigure_Frame[1800][3200];

static int              nCode127 = 127;
static double           rgfMcode[127];

static int
Shift127()
{
    int n1, n4, n7;

    n4 = (nCode127 >> 3) & 1;       // Venera-15 taps are 3 (could be 0, 2, 3, 5)
    n7 = (nCode127 >> 6) & 1;
    n1 = n4 ^ n7;
    nCode127 = (nCode127 << 1) | n1;
    return n1;
}

static void
InitializeCode()
{
    int i;
    double fSum;

    fSum = 0.0;
    for (i = 0; i < 127; i++) {
        rgfMcode[i] = double(2*Shift127()) - 1.0;
        fSum += rgfMcode[i];
    }
    printf("Sum of M-code %f\n", fSum);
}
//
//  circular correlation with the M-code pattern, to compress CW signal into pulse reflection
//
void
PulseCompression(Complex rgcSrc[127], Complex rgcDst[127])
{
    int i, j;
    Complex cSum;

    for (i = 0; i < 127; i++) {
        cSum = 0.0;
        for (j = 0; j < 127; j++)
            cSum = cSum + rgcSrc[j]*rgfMcode[(i + j) % 127];
        rgcDst[i] = cSum;
    }
}

void
PulseCompression2()
{
    int i, j, k, ii;
    double f;
    Complex c;
    FILE *pf;
    static Complex rgcSrc[22*127], rgcDst[20*127];

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        ii = 0;
        for (i = 0; i < 127; i++)
            rgcSrc[ii++] = rgcHologram[j][0][i];
        for (k = 0; k < 20; k++)
            for (i = 0; i < 127; i++)
                rgcSrc[ii++] = rgcHologram[j][k][i];
        for (i = 0; i < 127; i++)
            rgcSrc[ii++] = rgcHologram[j][19][i];

        for (ii = 0; ii < 20*127; ii++) {
            c = 0.0;
            for (i = 0; i < 127; i++) {
                c = c + rgcSrc[ii + i + 63]*rgfMcode[i];
            }
            rgcDst[ii] = c;
        }
        ii = 0;
        for (k = 0; k < 20; k++) {
            for (i = 0; i < 127; i++) {
                rgcCorrelated[j][k][i] = c = rgcDst[ii++];
                f = 0.0625*c.abs() - 0.1;
                if (f > 1.0)
                    rgnCorrelated[j][i+127*k] = 255;
                if (f < 0.0)
                    rgnCorrelated[j][i+127*k] = 0;
                else
                    rgnCorrelated[j][i+127*k] = int(255.0*f);
            }
        }
    }
    pf = fopen("correlated.raw", "wb");
    fwrite(rgnCorrelated, sizeof(rgnCorrelated), 1, pf);
    fclose(pf);
}

Complex
cis(double f)
{
    return Complex(cos(f), sin(f));
}

void
SlowTransform(Complex x[20][127], Complex S[20][127])
{
    int tau, nu, p, G_rbo;
    Complex c, rgc[20][127];

    G_rbo = 0;
    for (tau = 0; tau < 127; tau++) {
        for (nu = 0; nu < 20; nu++) {
            c = 0;
            for (p = 0; p < 20; p++) {
                c = c + x[p][tau]*cis(-D_2PI*(nu - 9 + G_rbo)*(tau + 127*p)/2540.0);
            }
            rgc[nu][tau] = c;
        }
    }
    for (tau = 0; tau < 127; tau++)
        for (nu = 0; nu < 20; nu++)
            S[nu][tau] = rgc[nu][tau];
}

void
DopplerAnalyze2()
{
    int i, j, k;
    Complex c;
    double f;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        SlowTransform(rgcCorrelated[j], rgcDoppler2[j]);
        for (i = 0; i < 127; i++) {
            for (k = 0; k < 20; k++) {
                c = rgcDoppler2[j][k][i];
                f = 0.0125*c.abs();
                if (f > 1.0)
                    rgnDoppler2[j][i + 127*k] = 255;
                else
                    rgnDoppler2[j][i + 127*k] = int(255.0*f);
            }
        }

    }
    pf = fopen("doppler2.raw", "wb");
    fwrite(rgnDoppler2, sizeof(rgnDoppler2), 1, pf);
    fclose(pf);
}

void
PhaseShift2()
{
    int i, j, k, iPhase;
    double x, y, f;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = A + B*x + C*x*x;
        iPhase = int(127.0*y/63.0 + 0.5) + 98 - 70;
        iPhase = int(2.0*y + 0.5) + 98 - 70;             // This works better.  Ratio of pulse duration.
        if (j == 7 || j == 207)
            printf("x = %f, y = %f, iPhase_SAR = %d, ", x, y, iPhase);
        for (k = 0; k < 20; k++) {
            for (i = 0; i < 127; i++) {
                rgcPhase2[j][k][(i + iPhase)%127] = rgcDoppler2[j][k][i];
            }
            for (i = 0; i < 127; i++) {
                f = 0.0125*rgcPhase2[j][k][i].abs();
                if (f > 1.0)
                    rgnPhase2[j][i + 127*k] = 255;
                else
                    rgnPhase2[j][i + 127*k] = int(255.0*f);
            }
        }
    }
    pf = fopen("phase2.raw", "wb");
    fwrite(rgnPhase2, sizeof(rgnPhase2), 1, pf);
    fclose(pf);
}

void
ApertureSynthesis2()
{
    int i, j, k, jDoppler;
    double fPower, fBright;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        for (k = 6; k <= 14; k++) {
            jDoppler = -int(floor(double(k - 10)*rgod[j].fFramePerDoppler + 0.5));
            if (j == 40 && k == 6)
                printf("jDoppler = %d\n", jDoppler);
            for (i = 0; i < 127; i++) {
                rgfSynth2[j+10+jDoppler][i] += rgcPhase2[j][k][i].power();
                rgfWeight[j+10+jDoppler][i] += 1.0;
            }
        }
    }
    for (j = 0; j < 3220; j++) {
        for (i = 0; i < 127; i++) {
            if (rgfWeight[j][i])
                fPower = rgfSynth2[j][i]/rgfWeight[j][i];
            else
                fPower = 0.0;
            fBright = 0.00625*sqrt(fPower) - 0.1;
            if (fBright < 0.0)
                rgnSynth2[j][i] = 0;
            else if (fBright < 1.0)
                rgnSynth2[j][i] = int(255.0*fBright);
            else
                rgnSynth2[j][i] = 255;
        }
    }
    pf = fopen("synth2.raw", "wb");
    fwrite(rgnSynth2, sizeof(rgnSynth2), 1, pf);
    fclose(pf);
}
//
//  Extract more azumuth resolution
//
void
SlowTransform3(Complex x[20][127], Complex S[40][127])
{
    int tau, nu, p, G_rbo, i, j;
    Complex c, rgc[40][127], xx[40][127];

    G_rbo = 0;                              // Soviet's added a heterodyne correction
    for (j = 0; j < 127; j++) {
        for (i = 0; i < 20; i++)
            xx[i][j] = x[i][j];
        for (; i < 40; i++)
            xx[i][j] = 0.0;
    }
    for (tau = 0; tau < 127; tau++) {
        for (nu = 0; nu < 40; nu++) {
            c = 0;
            for (p = 0; p < 40; p++) {
                c = c + xx[p][tau]*cis(-D_2PI*(nu - 19 + G_rbo)*(tau + 127*p)/5080.0);
            }
            rgc[nu][tau] = c;
        }
    }
    for (tau = 0; tau < 127; tau++)
        for (nu = 0; nu < 40; nu++)
            S[nu][tau] = rgc[nu][tau];
}

void
DopplerAnalyze3()
{
    int i, j, k;
    Complex c;
    double f;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        SlowTransform3(rgcCorrelated[j], rgcDoppler3[j]);
        for (i = 0; i < 127; i++) {
            for (k = 0; k < 40; k++) {
                c = rgcDoppler3[j][k][i];
                f = 0.0125*c.abs() - 0.1;
                if (f > 1.0)
                    rgnDoppler3[j][i + 127*k] = 255;
                else if (f < 0.0)
                    rgnDoppler3[j][i + 127*k] = 0;
                else
                    rgnDoppler3[j][i + 127*k] = int(255.0*f);
            }
        }

    }
    pf = fopen("doppler3.raw", "wb");
    fwrite(rgnDoppler3, sizeof(rgnDoppler3), 1, pf);
    fclose(pf);
}

void
PhaseShift3()
{
    int i, j, k, iPhase;
    double x, y, f;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = A + B*x + C*x*x;
        iPhase = int(2.0*y + 0.5) + 98 - 70;             // This works better.  Ratio of pulse duration.
        if (j == 7 || j == 207)
            printf("x = %f, y = %f, iPhase_SAR = %d, ", x, y, iPhase);
        for (k = 0; k < 40; k++) {
            for (i = 0; i < 127; i++) {
                rgcPhase3[j][k][(i + iPhase)%127] = rgcDoppler3[j][k][i];
                rgcPulsePhase[j][k][(i + iPhase)%127] = rgcCorrelated[j][k][i];
            }
            for (i = 0; i < 127; i++) {
                f = 0.0125*rgcPhase3[j][k][i].abs() - 0.1;
                if (f > 1.0)
                    rgnPhase3[j][i + 127*k] = 255;
                else if (f < 0.0)
                    rgnPhase3[j][i + 127*k] = 0;
                else
                    rgnPhase3[j][i + 127*k] = int(255.0*f);
            }
        }
    }
    pf = fopen("phase3.raw", "wb");
    fwrite(rgnPhase3, sizeof(rgnPhase3), 1, pf);
    fclose(pf);
}

void
ApertureSynthesis3()
{
    int i, j, jj, k, jDoppler;
    double fPower, fBright, fDoppler;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        jj = 2*j;
        for (k = 12; k <= 28; k++) {
            fDoppler = double(k - 20)*rgod[j].fFramePerDoppler;
            jDoppler = -int(floor(fDoppler + 0.5));
            for (i = 0; i < 127; i++) {
                rgfSynth3[jj+20+jDoppler][i] += rgcPhase3[j][k][i].power();
                rgfWeight[jj+20+jDoppler][i] += 1.0;
            }
        }
    }
    for (j = 0; j < 6440; j++) {
        for (i = 0; i < 127; i++) {
            if (rgfWeight[j][i])
                fPower = rgfSynth3[j][i]/rgfWeight[j][i];
            else
                fPower = 0.0;
            fBright = 0.00625*sqrt(fPower) - 0.1;
            if (fBright < 0.0)
                rgnSynth3[j][i] = 0;
            else if (fBright < 1.0)
                rgnSynth3[j][i] = int(255.0*fBright);
            else
                rgnSynth3[j][i] = 255;
        }
    }
    pf = fopen("synth3.raw", "wb");
    fwrite(rgnSynth3, sizeof(rgnSynth3), 1, pf);
    fclose(pf);
}

void
RealAperture3()
{
    int i, j, k;
    double fPower, fBright;
    FILE *pf;

    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        for (i = 0; i < 127; i++) {
            fPower = 0.0;
            for (k = 0; k < 20; k++) {
                fPower += rgcPulsePhase[j][k][i].power();
            }
            fBright = 0.00625*sqrt(fPower) - 0.1;
            if (fBright < 0.0)
                rgnReal3[j][i] = 0;
            else if (fBright < 1.0)
                rgnReal3[j][i] = int(255.0*fBright);
            else
                rgnReal3[j][i] = 255;
        }
    }
    pf = fopen("real3.raw", "wb");
    fwrite(rgnReal3, sizeof(rgnReal3), 1, pf);
    fclose(pf);
}
//
//  Build figure for book.  Full Synth3 swath is 6267 pixels long.  6300 = 1575 x 4.
//
void
FigureFrameSurvey()
{
    int i, j, k;
    double f;
    FILE *pf;

    for (j = 0; j < 1800; j++)
        for (i = 0; i < 3200; i++)
            rgnFigure_Frame[j][i] = 180;     // light-gray background
    for (j = 0; j < 1575; j++) {
        for (i = 0; i < 127; i++) {
            rgnFigure_Frame[j+10][i +  10    ] =    rgrlSurvey[j + 600].rgiqSAR[10][i].re << 4;
            rgnFigure_Frame[j+10][i +  30+127] = rgnCorrelated[j + 600][i + 10*127];
            rgnFigure_Frame[j+10][i +  50+254] =     rgnPhase3[j + 600][i + 10*127];
            // rgnFigure_Frame[j+10][i +  70+381] =      rgnReal3[j + 0][i];
            rgnFigure_Frame[j+10][i +  70+381] =      rgnSynth3[j + 1200 + 20][i];
            //rgnFigure_Frame[j+10][i + 110+635] =      rgnSynth3[j + 20+1575][i];
            //rgnFigure_Frame[j+10][i + 130+762] =      rgnSynth3[j + 20+2*1575][i];
            //rgnFigure_Frame[j+10][i + 150+889] =      rgnSynth3[j + 20+3*1575][i];
        }
    }
    for (k = 0; k < 40; k++) {
        for (i = 0; i < 127; i++) {
            f = 0.0125*rgcPhase3[600+240][k][i].abs() - 0.1;
            if (f < 0.0)
                rgnFigure_Frame[10+k+240-20][i + 90+508] = 0;
            else if (f < 1.0)
                rgnFigure_Frame[10+k+240-20][i + 90+508] = int(255.0*f);
            else
                rgnFigure_Frame[10+k+240-20][i + 90+508] = 255;
        }
    }
    pf = fopen("Figure_Frame.raw", "wb");
    fwrite(rgnFigure_Frame, sizeof(rgnFigure_Frame), 1, pf);
    fclose(pf);
}

void
ProcessFrameImages()
{
    int i, j, k, iPhase;
    double f, x, y;
    Complex c, cShift, rgc[20];
    FILE *pf;

    printf("Processing frame-survey images\n");
    InitializeCode();
    //
    //  Decode 8-bit coherent samples into complex numbers
    //
    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        for (k = 0; k < 20; k++) {
            for (i = 0; i < 127; i++) {
                rgcHologram[j][k][i] = c = rgcIQ[rgrlSurvey[j].rgiqSAR[k][i].n];
                f = c.abs();
                if (f > 1.0)
                    rgnHologram[j][i + 127*k] = 255;
                else
                    rgnHologram[j][i + 127*k] = int(255.0*f);
                rgnHologram[j][i + 127*k] = rgrlSurvey[j].rgiqSAR[k][i].re << 4;
            }
        }
    }
    pf = fopen("holo.raw", "wb");
    fwrite(rgnHologram, sizeof(rgnHologram), 1, pf);
    fclose(pf);
    //
    //  Test using 40 doppler channels
    //
    PulseCompression2();
    DopplerAnalyze3();
    PhaseShift3();
    ApertureSynthesis3();
    RealAperture3();
    FigureFrameSurvey();
    return;
    //
    //  Test doing pulse compression first, using a running correlation
    //
    PulseCompression2();
    DopplerAnalyze2();
    PhaseShift2();
    ApertureSynthesis2();
    return;
    //
    //  Fourier transform into 20 Doppler-shift channels
    //
    /*
    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        for (i = 0; i < 127; i++) {
            for (k = 0; k < 20; k++)
                rgc[k] = rgcHologram[j][k][i];
            fft(rgc, 20);
            for (k = 0; k < 20; k++) {
                rgcDoppler[j][k][i] = c = rgc[k];
                f = 0.125*c.abs();
                if (f > 1.0)
                    rgnDoppler[j][i + 127*k] = 255;
                else
                    rgnDoppler[j][i + 127*k] = int(255.0*f);
            }
        }
    }
    */
    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        SlowTransform(rgcHologram[j], rgcDoppler[j]);
        for (i = 0; i < 127; i++) {
            for (k = 0; k < 20; k++) {
                c = rgcDoppler[j][k][i];
                f = 0.125*c.abs();
                if (f > 1.0)
                    rgnDoppler[j][i + 127*k] = 255;
                else
                    rgnDoppler[j][i + 127*k] = int(255.0*f);
            }
        }

    }
    pf = fopen("doppler.raw", "wb");
    fwrite(rgnDoppler, sizeof(rgnDoppler), 1, pf);
    fclose(pf);
    //
    //  Frequency shift the Doppler channels before pulse compression
    //
    /*
    cShift = cis(3.0*DOPPLERSHIFT);
    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        for (k = 0; k < 20; k++)
            for (i = 0; i < 127; i++) {
                rgcDoppler[j][k][i] = rgcDoppler[j][k][i] * cShift;
            }
    }
    */
    //
    //  Pulse compression
    //
    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        for (k = 0; k < 20; k++) {
            PulseCompression(rgcDoppler[j][k], rgcPulse[j][k]);
            for (i = 0; i < 127; i++) {
                c = rgcPulse[j][k][i];
                f = 0.0125*c.abs();
                if (f > 1.0)
                    rgnPulse[j][i + 127*k] = 255;
                else
                    rgnPulse[j][i + 127*k] = int(255.0*f);
            }
        }
    }
    pf = fopen("pulse.raw", "wb");
    fwrite(rgnPulse, sizeof(rgnPulse), 1, pf);
    fclose(pf);
    //
    //  Phase correction, based on quartic curve derived from the band-survey images
    //
    for (j = 0; j < 3200; j++) {
        if (rgnBadSAR[j])
            continue;
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = A + B*x + C*x*x;
        iPhase = int(127.0*y/63.0 + 0.5) + 98;
        iPhase = int(2.0*y + 0.5) + 98;             // This works better.  Ratio of pulse duration.
        if (j == 7 || j == 207)
            printf("x = %f, y = %f, iPhase_SAR = %d, ", x, y, iPhase);
        for (k = 0; k < 20; k++) {
            for (i = 0; i < 127; i++) {
                rgcPhase[j][k][i] = rgcPulse[j][k][(i + iPhase)%127];
            }
            for (i = 0; i < 127; i++) {
                f = 0.0125*rgcPhase[j][k][i].abs();
                if (f > 1.0)
                    rgnPhase[j][i + 127*k] = 255;
                else
                    rgnPhase[j][i + 127*k] = int(255.0*f);
            }
        }
    }
    pf = fopen("phase.raw", "wb");
    fwrite(rgnPhase, sizeof(rgnPhase), 1, pf);
    fclose(pf);
}