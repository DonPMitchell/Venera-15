#include "stdafx.h"
#include "Venera15.h"
#include "Vector.h"
#define SOLARDAY    86400.0                 // seconds

struct KeplerianElements {
    double  fEpoch;             // time of parameter definitions
    double  fInclination;       // i, inclinarion of orbital plane
    double  fAscendingNode;     // OMEGA, longitude of ascending node
    double  fEccentricity;      // e, eccentricity of elliptical orbit
    double  fPerigee;           // omega, argument of perigee (not longitude)
    double  fMeanAnomaly;       // angle or revolution at epoch
    double  fMeanMotion;        // degrees per solar day (should be revolutions?)
};

static int rgnDays[] = { 31,28,31,30,31,30,31,31,30,31,30,31 };

OrbitalData     rgod[3200];
Vector          vOrbitalAxis;       // vector perpendicular to orbital plane
Vector          vAscendingNode;     // vector pointing at ascending node in orbital plane
Vector          vBinormal;          // vector perpendicular to axis and node

double
JulianDay(int nYear, int nMonth, int nDay, int nHour, int nMinute, double fSecond)
{
    int jd;

    if (nMonth < 1 || nMonth > 12 || nDay < 1 || nDay > rgnDays[nMonth - 1])
        printf("BAD JULIAN DAY ARGUMENTS %d %d %d\n", nYear, nMonth, nDay);
    jd = ( 1461 * ( nYear + 4800 + ( nMonth - 14 ) / 12 ) ) / 4 +
         ( 367 * ( nMonth - 2 - 12 * ( ( nMonth - 14 ) / 12 ) ) ) / 12 -
         ( 3 * ( ( nYear + 4900 + ( nMonth - 14 ) / 12 ) / 100 ) ) / 4 +
         nDay - 32075;
    return jd + (double(nHour-12) + (double(nMinute) + fSecond/60.0)/60.0)/24.0;
}

double
CubeRoot(double x)
{
    double t, t2;
    unsigned *pnT;

    if (x == 0.0)
        return 0.0;
    pnT = (unsigned *)&t;
    t = fabs(x);
    pnT[1] = pnT[1]/3 + 715094163;  // plus 2/3 of one as integer
    if (x < 0.0)
        t = -t;
    t2 = x/(t*t);
    t = t*(t + t2 + t2)/(t + t + t2);               // Halley-method step
    t = (t + t + x/(t*t))/3.0;                      // Newton-method steps
    t = (t + t + x/(t*t))/3.0;
    return t;
}

Vector
LocateKeplerian(double fJulianDay, KeplerianElements &ke, double fGM, int j = -1)
{
    double fMeanAnomaly, fEccentricAnomaly, fE, t;
    double x, y, r, fTrueAnomaly, fA, fCosN, fSinN, fCosVW, fSinVW, fCosI, fSinI;
    Vector v;

    //
    //  Calculate mean, eccentric and true anomalies
    //
    t = fJulianDay - ke.fEpoch;
    fMeanAnomaly = fmod((ke.fMeanAnomaly + ke.fMeanMotion*t), 360.0);
    if (fMeanAnomaly < 0.0)
        fMeanAnomaly += 360.0;
    fMeanAnomaly *= D_DTR;
    fEccentricAnomaly = fMeanAnomaly + ke.fEccentricity*sin(fMeanAnomaly)
                        * (1.0 + ke.fEccentricity*cos(fMeanAnomaly));
    do {
        fE = fEccentricAnomaly;
        fEccentricAnomaly = fE - (fE - ke.fEccentricity*sin(fE) - fMeanAnomaly)
                            / (1.0 - ke.fEccentricity*cos(fE));
    } while (fabs(fE - fEccentricAnomaly) > 1.0e-8);
    fE = fEccentricAnomaly;

    fA = CubeRoot(fGM/(ke.fMeanMotion*ke.fMeanMotion*D_DTR*D_DTR/(SOLARDAY*SOLARDAY)));
    x = fA * (cos(fEccentricAnomaly) - ke.fEccentricity);
    y = fA * (sqrt(1.0 - ke.fEccentricity*ke.fEccentricity)*sin(fEccentricAnomaly));
    fTrueAnomaly = atan2(y, x);
    if (j >= 0)
        rgod[j].fTrueAnomaly = fTrueAnomaly/D_DTR;        // A side effect, but its most accurate to save it now
    r = sqrt(x*x + y*y);
    fCosN = cos(ke.fAscendingNode*D_DTR);
    fSinN = sin(ke.fAscendingNode*D_DTR);
    fCosVW = cos(ke.fPerigee*D_DTR + fTrueAnomaly);
    fSinVW = sin(ke.fPerigee*D_DTR + fTrueAnomaly);
    fCosI = cos(ke.fInclination*D_DTR);
    fSinI = sin(ke.fInclination*D_DTR);
    v.x = r*(fCosN*fCosVW - fSinN*fSinVW*fCosI);
    v.y = r*(fSinN*fCosVW + fCosN*fSinVW*fCosI);
    v.z = r*(fSinVW*fSinI);
    return v;
}

Vector
VelocityKeplerian(double fJulianDay, KeplerianElements &ke, double fGM)
{
    Vector vP1; 
    Vector vP2;
    Vector vP3;
    Vector vP4;
    double tDelta;
    
    tDelta = 0.125/ke.fMeanMotion;           // time to go 0.25 degrees
    vP1 = LocateKeplerian(fJulianDay - tDelta, ke, fGM);
    vP2 = LocateKeplerian(fJulianDay + tDelta, ke, fGM);
    vP3 = LocateKeplerian(fJulianDay - 2.0*tDelta, ke, fGM);
    vP4 = LocateKeplerian(fJulianDay + 2.0*tDelta, ke, fGM);
    //return (vP2 - vP1)*(0.5/(tDelta*SOLARDAY));
    return (vP2*8.0 - vP1*8.0 + vP3 - vP4)*(1.0/(12.0*tDelta*SOLARDAY));
}

void
OrbitalBasis(KeplerianElements &ke)
{
    Vector u, v, w;
    double c, s;

    //
    //  Ascending node is on the equator (z = 0 plane)
    //
    vAscendingNode = u = Vector(cos(ke.fAscendingNode*D_DTR), sin(ke.fAscendingNode*D_DTR), 0.0);
    //
    //  Rotate z axis about node by angle of inclination
    //
    c = cos(ke.fInclination*D_DTR);
    s = sin(ke.fInclination*D_DTR);
    vOrbitalAxis = v = Vector(u.x*u.z*(1.0 - c) + u.y*s, u.y*u.z*(1.0 - c) - u.x*s, c + u.z*u.z*(1.0 - c));
    vBinormal = CrossProduct(u, v);
    //
    //  Test
    //
    v = LocateKeplerian(ke.fEpoch, ke, C_VENUS_GM);
    printf("Basis: dot with axis = %f\n", DotProduct(v, vOrbitalAxis));
    v = LocateKeplerian(ke.fEpoch + 0.33, ke, C_VENUS_GM);
    printf("       dot with axis = %f\n", DotProduct(v, vOrbitalAxis));
}

/*
struct KeplerianElements {
    double  fEpoch;             // time of parameter definitions
    double  fInclination;       // i, inclinarion of orbital plane
    double  fAscendingNode;     // OMEGA, longitude of ascending node
    double  fEccentricity;      // e, eccentricity of elliptical orbit
    double  fPerigee;           // omega, argument of perigee (not longitude)
    double  fMeanAnomaly;       // angle or revolution at epoch
    double  fMeanMotion;        // degrees per solar day (should be revolutions?)
};
*/

KeplerianElements keVenera15 = {
    2445622.5,      // Oct 15, 1983
    109.9200,       // Orbital inclination
    351.8140,       // Longitude of ascending node
    0.82080360,     // Eccentricity
    113.6510,       // Argument of pericenter
    320.149516,     // Mean anomaly at epoch (9343.9 seconds before pericenter)
    368.484444      // Mean Motion (degrees per day)
};

//
//  Build orbit data associated with radar records
//
void
InitializeOrbitalData()
{
    int i, j;
    double a, p, e, t, fJulian, r, fVel, aa, bb, cc, rRadius, rRange, fCosine;
    double fSinPhi, fLat, fCosPhi, fHorz, fVert, fDeltai, fCosi, fSini, fRange, fVelRange;
    double fDoppler, fLastDoppler, rgfDoppler[3], fDegPerFrame, fDopplerPerFrame, fFramePerChannel;
    Vector vSpacecraft, vNorm;
    Vector vVel, vVelHorz;

    OrbitalBasis(keVenera15);
    a = 1000.0*V_SEMIMAJOR;
    e = V_ECCENTRIC;
    p = a*(1.0 - e*e);
    rRadius = C_VENUSRADIUS;                 // JPL's modern value for radius of Venus
    fCosine = cos(D_DTR*V_RADARANGLE);        // Cosine of look angle of spacecraft (10 degrees)
    for (j = 0; j < 3200; j++) {
        t = 0.300*(double(j) - fPericenter);                            // seconds from pericenter
        rgod[j].tPericenter = t;
        fJulian = JulianDay(1983, 10, 15, 0, 0, 9343.9 + t);
        vSpacecraft = LocateKeplerian(fJulian, keVenera15, C_VENUS_GM, j);           // coordinates of spacecraft
        r = vSpacecraft.abs();
        rgod[j].rAltitude = r - rRadius;
        rgod[j].x = vSpacecraft.x;
        rgod[j].y = vSpacecraft.y;
        rgod[j].z = vSpacecraft.z;
        fVel = C_VENUS_GM*(2.0/r - 1.0/a);
        vVel = VelocityKeplerian(fJulian, keVenera15, C_VENUS_GM);      // velocity of spaceraft
        rgod[j].vx = vVel.x;
        rgod[j].vy = vVel.y;
        rgod[j].vz = vVel.z;
        //
        //  Law of cosines: rRadius*rRadius = rRange*rRange + r*r - 2.0*rRange*r*fCosine
        //
        aa = 1.0;
        bb = -2.0*r*fCosine;
        cc = r*r - rRadius*rRadius;
        rRange = (-bb - sqrt(bb*bb - 4.0*aa*cc))/(2.0*aa);
        rgod[j].rRange = rRange;
        fSinPhi = rRange*sin(D_2PI*V_RADARANGLE/360.0)/rRadius;     // Law of sines
        rgod[j].fPhi = asin(fSinPhi)/D_DTR;
        //
        //  Calculate Doppler-shift outsets
        //
        fCosPhi = cos(rgod[j].fPhi*D_DTR);
        vNorm = vSpacecraft *(1.0/vSpacecraft.abs());   // unit vector from planet center to spacecraft
        fVert = -DotProduct(vVel, vNorm);
        vVelHorz = vVel - vNorm*fVert;
        fHorz = vVelHorz.abs();
        fLastDoppler = 0.0;
        for (i = -1; i <= 1; i++) {
            fDeltai = double(i)/50.0;
            fCosi = cos(fDeltai*D_DTR);
            fSini = sin(fDeltai*D_DTR);
            fRange = sqrt(rRadius*rRadius + r*r - 2.0*rRadius*r*fCosi*fCosPhi);
            fVelRange = (fHorz*rRadius*fSini*fCosPhi - fVert*(r - rRadius*fCosi*fCosPhi))/fRange;
            fDoppler = -2.0*fVelRange/V_WAVELENGTH;
            rgfDoppler[i+1] = fDoppler;
            if (j == 0 || j == 1735)
                printf("  Di = %f, VR = %f, Dop = %f, Del = %f\n", fDeltai, fVelRange, fDoppler, fDoppler-fLastDoppler);
            fLastDoppler = fDoppler;
        }
        rgod[j].fDopplerPerDegree = (rgfDoppler[0] - rgfDoppler[2])/(2.0/50.0);
        //
        //  Print some states at beginning, pericenter, end of survey
        //
        if (j == 0 || j == 1735 || j == 3130) {
            printf("OrbitalData [%d]\n", j);
            printf("  h = %f, orbital velocity %f %f\n", rgod[j].rAltitude,
                sqrt(fVel), sqrt(vVel.x*vVel.x + vVel.y*vVel.y + vVel.z*vVel.z));
            printf("  fVert = %f, fHorz = %f\n", fVert, fHorz);
            printf("  v = %f, %f, %f\n", rgod[j].x, rgod[j].y, rgod[j].z);
            printf("  True Anomaly = %f, polar angle Phi = %f\n", rgod[j].fTrueAnomaly, rgod[j].fPhi);
            fLat = acos(vSpacecraft.z/r);
            printf("   Latitude = %f\n", 90.0 - fLat/D_DTR);
        }
    }
    //
    //  process true anomaly, time and doppler shift
    //
    for (j = 0; j < 3200-1; j++) {
        fDegPerFrame = rgod[j].fTrueAnomaly - rgod[j+1].fTrueAnomaly;
        fDopplerPerFrame = rgod[j].fDopplerPerDegree * fDegPerFrame;
        fFramePerChannel = 0.5*511.811 / fDopplerPerFrame;
        if (j == 0 || j == 1735 || j == 3130) {
            printf("DopplerChannels: DPF = %f, FPC = %f\n", fDopplerPerFrame, fFramePerChannel);
        }
        rgod[j].fFramePerDoppler = fFramePerChannel;
    }
}
//
//  Transform the Phi, i image space to time-delay/Doppler-shift coordinates
//
DelayDoppler
SurfaceToSignal(OrbitalData &od, double fTrueAnomaly, double fLatitude)
{
    DelayDoppler dd;
    Vector vSpacecraft, vSurface;
    double fCos;

    vSpacecraft = Vector(od.x, od.y, od.z);
    vSurface = vSpacecraft *(C_VENUSRADIUS/vSpacecraft.abs());
    fCos = cos(fLatitude*D_DTR);
    vSurface = Vector(cos(fTrueAnomaly*D_DTR)*fCos, sin(fTrueAnomaly*D_DTR)*fCos, sin(fLatitude*D_DTR))*C_VENUSRADIUS;
    return dd;
}