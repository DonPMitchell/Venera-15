//
//  Process the on-board band images from the OKB MEI radar system
//  D.P. Mitchell  2016/02/12.
//
#include "stdafx.h"
#include "Venera15.h"

static double           rgfStrip[3200][63*3];       // Band image, which has three Doppler-shift channels
static double           rgfPhase[3200][63*3];       // Phase-corrected image
static double           rgfDoppler[3200][63];
static double           rgfGain[3200];
static double           rgfAntennaGain[63];
static unsigned char    rgnStrip[3200][200*20];
static int              rgnStripMin[3200];
static int              rgnStripMinUnmod[3200];
static unsigned char    rgnCurveFit[3200][3200];
static unsigned char    rgnFigure_Band[1500][1500];

double                  A, B, C, fPericenter = 1734.8;

static double
Det(double a00, double a01, double a02,
    double a10, double a11, double a12,
    double a20, double a21, double a22)
{
    return a00*a11*a22 + a01*a12*a20 + a02*a10*a21 - a02*a11*a20 - a00*a12*a21 - a01*a10*a22;
}

static void
CurveFit()
{
    int i, j, iMod, nCheck;
    double a, b, c, d, x, y, xMin;
    double sy, sx, sx2, sxy, sx3, sx4, sx2y, sn;
    FILE *pf;
    
    printf("  Curve Fitting to find phase shifts\n");
    iMod = 0;
    nCheck = 2;
    for (j = 0; j < 3100; j++) {
        if (rgnBadBand[j])
            continue;
        if (nCheck < 0 && rgnStripMin[j] - rgnStripMin[j-1] > 40) {
            iMod -= 63;
            nCheck = 20;
        } else if (nCheck < 0 && rgnStripMin[j] - rgnStripMin[j-1] < -40) {
            iMod += 63;
            nCheck = 20;
        }
        --nCheck;
        rgnStripMinUnmod[j] = rgnStripMin[j] + iMod;
        rgnCurveFit[j][rgnStripMinUnmod[j]] = 128;
    }
    //
    //  LSE fit phase offset to curve, y = a + bx + cx**2
    //  where x = (j - fPericenter)**2
    //
    sy = sx = sx2 = sxy = sx3 = sx4 = sx2y = sn = 0.0;
    for (j = 0; j < 3100; j++) {
        if (rgnBadBand[j])
            continue;
        x = (double(j) - fPericenter)*(double(j) - fPericenter);  
        y = rgnStripMinUnmod[j];
        sy += y;
        sx += x;
        sx2 += x*x;
        sxy += x*y;
        sx3 += x*x*x;
        sx4 += x*x*x*x;
        sx2y += x*x*y;
        sn += 1.0;
    }
    a = Det(sy, sx, sx2, sxy, sx2, sx3, sx2y, sx3, sx4);
    b = Det(sn, sy, sx2, sx, sxy, sx3, sx2, sx2y, sx4);
    c = Det(sn, sx, sy, sx, sx2, sxy, sx2, sx3, sx2y);
    d = Det(sn, sx, sx2, sx, sx2, sx3, sx2, sx3, sx4);
    a /= d;
    b /= d;
    c /= d;
    A = a;
    B = b;
    C = c;
    //
    //  Use the curve to fix mod 63 errors in StripMinUnmod[], and then compute better fit
    //
    for (j = 0; j < 3100; j++) {
        if (rgnBadBand[j])
            continue;
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = A + B*x + C*x*x;
        i = int(y + 0.5);
        if (rgnStripMinUnmod[j] < i - 32)
            rgnStripMinUnmod[j] += 63;
        if (rgnStripMinUnmod[j] > i + 32)
            rgnStripMinUnmod[j] -= 63;
        rgnCurveFit[j][200 + rgnStripMinUnmod[j]] = 255;
    }
    sy = sx = sx2 = sxy = sx3 = sx4 = sx2y = sn = 0.0;
    for (j = 0; j < 3100; j++) {
        if (rgnBadBand[j])
            continue;
        x = (double(j) - fPericenter)*(double(j) - fPericenter);  
        y = rgnStripMinUnmod[j];
        sy += y;
        sx += x;
        sx2 += x*x;
        sxy += x*y;
        sx3 += x*x*x;
        sx4 += x*x*x*x;
        sx2y += x*x*y;
        sn += 1.0;
    }
    a = Det(sy, sx, sx2, sxy, sx2, sx3, sx2y, sx3, sx4);
    b = Det(sn, sy, sx2, sx, sxy, sx3, sx2, sx2y, sx4);
    c = Det(sn, sx, sy, sx, sx2, sxy, sx2, sx3, sx2y);
    d = Det(sn, sx, sx2, sx, sx2, sx3, sx2, sx3, sx4);
    a /= d;
    b /= d;
    c /= d;
    A = a;
    B = b;
    C = c;
    for (j = 0; j < 3100; j++) {
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = A + B*x + C*x*x;
        i = int(y + 0.5);
        rgnCurveFit[j][i] = 255;
    }
    //
    //  Calculate a good pericenter by fitting parabola to minimum
    //
    sy = sx = sx2 = sxy = sx3 = sx4 = sx2y = sn = 0.0;
    for (j = 0; j < 3100; j++) {
        if (rgnBadBand[j])
            continue;
        x = double(j);  
        y = rgnStripMinUnmod[j];
        sy += y;
        sx += x;
        sx2 += x*x;
        sxy += x*y;
        sx3 += x*x*x;
        sx4 += x*x*x*x;
        sx2y += x*x*y;
        sn += 1.0;
    }
    a = Det(sy, sx, sx2, sxy, sx2, sx3, sx2y, sx3, sx4);
    b = Det(sn, sy, sx2, sx, sxy, sx3, sx2, sx2y, sx4);
    c = Det(sn, sx, sy, sx, sx2, sxy, sx2, sx3, sx2y);
    d = Det(sn, sx, sx2, sx, sx2, sx3, sx2, sx3, sx4);
    a /= d;
    b /= d;
    c /= d;
    xMin = -b/(2.0*c);
    printf("a = %g, b = %g, c = %g\n", A, B, C);
    printf("jPericenter delta = %f\n", xMin);           // 1734.747, inserted in fPericenter above
   
    pf = fopen("CurveFit.raw", "wb");
    fwrite(rgnCurveFit, sizeof(rgnCurveFit), 1, pf);
    fclose(pf);
}
//
//  Make figures for book
//
void
FigureBandSurvey()
{
    int i, j;
    FILE *pf;

    for (j = 0; j < 1500; j++)
        for (i = 0; i < 1500; i++)
            rgnFigure_Band[j][i] = 180;     // light-gray background
    for (j = 0; j < 1040; j++) {
        for (i = 0; i < 189; i++) {
            rgnFigure_Band[j+10][i +  10    ] = rgnStrip[j + 800][i];
            rgnFigure_Band[j+10][i +  30+189] = rgnStrip[j + 800][i + 200];
            rgnFigure_Band[j+10][i +  50+378] = rgnStrip[j + 800][i + 400];
        }
        for (i = 0; i < 63; i++) {
            rgnFigure_Band[j+10][i +  70+567] = rgnStrip[j +  800][i + 600];
            rgnFigure_Band[j+10][i +  90+630] = rgnStrip[j +    0][i + 700];
            rgnFigure_Band[j+10][i +  92+693] = rgnStrip[j + 1040][i + 700];
            rgnFigure_Band[j+10][i +  94+756] = rgnStrip[j + 2080][i + 700];
            rgnFigure_Band[j+10][i + 104+819] = 0;
        }
    }
    pf = fopen("Figure_Band.raw", "wb");
    fwrite(rgnFigure_Band, sizeof(rgnFigure_Band), 1, pf);
    fclose(pf);
}
//
//  Calculate orbital parameters
//
static void
Orbital()
{
    double a, e, p, fPeriod, fEnergy, fCosine;
    double rMin, rMax, rPerigeeRange, rPericenter, rVenus, rRange, x, y;
    double aa, bb, cc;
    int j;

    printf("Orbital Parameters\n");
    a = V_SEMIMAJOR*1000.0;                         // Defined in Venera-15 Oct 16 telemetry file
    e = V_ECCENTRIC;
    rPericenter = a*(1.0 - e);
    p = a*(1 - e*e);                                // Semi-latus rectum
    fPeriod = D_2PI*sqrt(a*a*a/C_VENUS_GM);
    printf("  Period = %f seconds\n", fPeriod);
    fEnergy = -C_VENUS_GM/(2.0*a);
    fCosine = cos(D_2PI*V_RADARANGLE/360.0);
    rVenus = 1000.0 * V_VENUSRADIUS;
    //
    //  use law of cosines to find the range to target at pericenter. 
    //  rVenus*rVenus = rRange*rRange + rPericenter*rPericenter - 2.0*rRange*rPericenter*fCosine
    //
    aa = 1.0;
    bb = -2.0*rPericenter*fCosine;
    cc = rPericenter*rPericenter - rVenus*rVenus;
    rPerigeeRange = (-bb - sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);
    printf("  Radius = %f, pericenter = %f, range = %f\n", rVenus, rPericenter, rPerigeeRange);
    printf("  Perigee = %f\n", rPericenter - rVenus);
    //
    //  Use phase shift to calculate change of range to target spot on surface.
    //  We know the pericenter from orbital elements, but later we calculate height from radio-range
    //
    rMin = 1000000000000.0;
    rMax = 0.0;
    for (j = 0; j < 3200; j++) {
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = B*x + C*x*x; 
        rgod[j].rPhaseRange = rRange = rPerigeeRange - 0.5*C_SPEEDOFLIGHT * y * V_BANDPULSE;    // 0.5, twice distance
    }
}
//
//  Process the OKB MEI band-survey data
//
void
ProcessStripeImages()
{
    int i, j, k, iMin, iPhase, ii, n;
    double f, fSum, fMin, fMiddle, rgf[63];
    double x, y;
    FILE *pf;

    printf("Processing Band-Image data\n");
    //
    //  Correct for automatic gain by dividing by the average
    //
    for (j = 0; j < 3200; j++) {
        fSum = 0.0;
        for (i = 0; i < 3*63; i++) {
            rgfStrip[j][i] = f = float(rgrlSurvey[j].rgs1.rgn[i/63][i%63]);
            fSum += f;
        }
        if (fSum)
            rgfGain[j] = (0.5f/1.25f)*30000.0f/fSum;
        if (j == 7)
            printf("Typical gain: %f\n", rgfGain[j]);
        for (i = 0; i < 3*63; i++) {
            rgfStrip[j][i] = f = float(rgrlSurvey[j].rgs1.rgn[i/63][i%63]) * rgfGain[j];
            rgnStrip[j][i] = unsigned char(0.5f*float(rgrlSurvey[j].rgs1.rgn[i/63][i%63]));
            if (f <= 255.0f)
                rgnStrip[j][i+200] = int(f);
            else
                rgnStrip[j][i+200] = 255;
        }
    }
    //
    //  Phase correction, distace of orbit to surface fits a quartic equations very closely
    //
    for (j = 0; j < 3200; j++) {
        if (rgnBadBand[j]) {
            rgnStrip[j][192] = rgnStrip[j][193] = 255;
            continue;
        }
        fMin = 1000000.0;
        iMin = 0;
        for (i = 63; i < 2*63; i++) {
            f = 0.0;
            for (k = -6; k <= 6; k++)
                f += rgfStrip[j][i+k];
            if (f < fMin) {
                fMin = f;
                iMin = i;
            }
        }
        rgnStripMin[j] = iMin;
        if (j && rgnBadBand[j-1])
            rgnStripMin[j-1] = iMin;
        rgnCurveFit[j][iMin] = 255;
    }
    CurveFit();
    for (j = 0; j < 3200; j++) {
        x = (double(j) - fPericenter)*(double(j) - fPericenter);
        y = A + B*x + C*x*x;
        iPhase = int(y + 0.5);
        for (i = 0; i < 3*63; i++) {
            k = i/63;
            ii = (i + iPhase)%63;
            f = rgfStrip[j][k*63 + ii];
            rgfPhase[j][i] = f;
            if (f <= 255.0f)
                rgnStrip[j][i+400] = int(f);
            else
                rgnStrip[j][i+400] = 255;
        }
    }
    //
    //  Average the Doppler channels together for less speckle noise
    //
    for (j = 0; j < 3200-4; j++) {
        n = 0;
        for (i = 0; i < 63; i++)
            rgf[i] = 0.0;
        for (k = 0; k < 3; k++) {
            if (rgnBadBand[j+2*k])
                continue;
            n++;
            for (i = 0; i < 63; i++)
                rgf[i] += rgfPhase[j+2*k][i+63*k];
        }
        if (n)
            for (i = 0; i < 63; i++) {
                f = rgf[i] / float(n);
                rgfDoppler[j][i] = f;
                if (f <= 255.0f)
                    rgnStrip[j][i+600] = int(f);
                else
                    rgnStrip[j][i+600] = 255;
            }
    }
    //
    //  Correct for antenna gain pattern (Bockstein's method)
    //
    for (j = JGAIN; j < 3200-JGAIN; j++) {
        for (i = 0; i < 63; i++)
            rgfAntennaGain[i] = 0.0;
        n = 0;
        for (k = -JGAIN; k <= JGAIN; k++) {
            if (rgnBadBand[j+k])
                continue;
            n++;
            for (i = 0; i < 63; i++)
                rgfAntennaGain[i] += rgfDoppler[j+k][i];
        }
        if (n) {
            for (i = 0; i < 63; i++)
                rgfAntennaGain[i] /= float(n);
        } else
            continue;
        fMiddle = (rgfAntennaGain[30] + rgfAntennaGain[31] + rgfAntennaGain[32])/3.0f;
        if (j == 400)
            printf("fMiddle = %f (use instead of 128.0 Bockstein constant\n", fMiddle);
        if (j == JGAIN) {   
            for (k = 0; k < JGAIN; k++)
                for (i = 0; i < 63; i++) {
                    f = 90.0f*rgfDoppler[k][i]/rgfAntennaGain[i];
                    if (f <= 255.0f)
                        rgnStrip[k][i+700] = int(f);
                    else
                        rgnStrip[k][i+700] = 255;
                }
        }
        for (i = 0; i < 63; i++) {
            f = 90.0f*rgfDoppler[j][i]/rgfAntennaGain[i];
            if (f <= 255.0f)
                rgnStrip[j][i+700] = int(f);
            else
                rgnStrip[j][i+700] = 255;
        }
    }
    //
    //  Correct for antenna gain (my method)
    //
    for (j = 0; j < 1600; j++)
        for (i = 0; i < 63; i++)
            rgfAntennaGain[i] += rgfDoppler[j][i];
    f = 0.0;
    for (i = 31-2; i <= 31+2; i++)
        f += rgfAntennaGain[i]/5.0f;
    for (i = 0; i < 63; i++)
        rgfAntennaGain[i] /= f;
    for (j = 0; j < 3200; j++) {
        for (i = 0; i < 63; i++) {
            f = rgfDoppler[j][i] / rgfAntennaGain[i];
            if (f <= 255.0f)
                rgnStrip[j][i+800] = int(f);
            else
                rgnStrip[j][i+800] = 255;
        }
    }
    //
    //  Use phase shift to calculate distance from illuminated stop on surface
    //
    Orbital();
    //
    //  Write out the final image
    //
    pf = fopen("Stripe.raw", "wb");
    fwrite(rgnStrip, sizeof(rgnStrip), 1, pf);
    fclose(pf);
    FigureBandSurvey();
}

