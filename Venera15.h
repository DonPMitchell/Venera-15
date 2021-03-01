//
//  Venera15 - analyze Venera-15 telemetry
//  D.P. Mitchell  2016/02/12.
//
#pragma once
#define DIR     "E:\\Pictures\\Soviet Image Data\\Venera 15 and 16\\"
#define JGAIN   31
#define D_2PI   6.28318530717958647692528676655900576839433879875021
#define D_DTR   (D_2PI/360.0)
//
//  Modern parameters
//
#define C_SPEEDOFLIGHT  2.99792458e+8
#define C_JPL_GM2       0.724345248616270270e-9         // GM Venus, newer JPL value
#define C_SOV_GM2       0.7243452486e-9                 // Soviet value (EPM2004)
#define C_JPL_AU2M      448485855946515.83566006340476058e+9    // au**3/day**2 to m**3/sec**2
#define C_VENUS_GM      (C_JPL_AU2M*C_JPL_GM2)   // DE-405 && 403, m**3/sec**2
#define C_VENUSORBIT    1.08208930e+11  // Meters
#define C_VENUSRADIUS   6.05184e+6      // Meters
#define C_VENUSDAY     -243.018484      // Solar Days (retrograde)
#define C_VENUSYEAR     0.61519726      // Sidereal Years
//
//  Parameters defined by Soviets
//
#define V_PULSEFREQUENCY    650000.0
#define V_MASTERFREQUENCY   (V_PULSEFREQUENCY * 80.0)       // 52 MHz crystal master oscillator
#define V_CARRIERFREQUENCY  (V_MASTERFREQUENCY * 72.0)      // 3744 MHz
#define V_WAVELENGTH        (C_SPEEDOFLIGHT/V_CARRIERFREQUENCY)
#define V_SARPULSE      (1.0/V_PULSEFREQUENCY)              // about 1.54 microseconds
#define V_SARCYCLE      (127.0/V_PULSEFREQUENCY)            // M-code cycle time, about 195 microseconds
#define V_BANDPULSE     (2.0/V_PULSEFREQUENCY)
#define V_SPEEDOFLIGHT  299792.46           // used by Soviets
#define V_VENUSRADIUS   6051.00             // used by Soviets
#define V_VENUS_GM2     3.2485961e14        // Russian, Venera-15
#define V_RADARANGLE    -10.0               // look angle of radar on Oct 16 survey
#define V_SEMIMAJOR     38848.68            // semimajor axis
#define V_ECCENTRIC     0.82080360          // eccentricity of orbit
#define V_PERICENTER    (V_SEMIMAJOR*(1.0 - V_ECCENTRIC))       // 6961.544 km
#define V_PERIGEE       (V_PERICENTER - V_VENUSRADIUS)          //  910.544 km

//
//  8-bit representation of coherent signal samples
//
struct IQ {
    union {
        struct {
            unsigned char   im : 4;
            unsigned char   re : 4;
        };
        unsigned char n;
    };
};
//
//  Venera-15/16 telemetry record
//
struct RawHeader {
    //unsigned short  nLine;
    unsigned char   rgn[33];
};

struct RawStripe {
    unsigned short  rgn[3][63];             // doppler channels, range data.
};

struct RadarLook {
    RawHeader       rhHeader;               // automatic gain levels, M-code, etc.
    IQ              rgiqSAR[20][127];       // radio hologram (inphase, quadrature samples).
    IQ              rgiqALT[14][31];        // altimeter data. First sample is zero. REVIEW: two zero bytes
    RawStripe       rgs1;                   // band-survey image stripe
    RawStripe       rgs2;                   // second copy of some kind of band-survey data.
    unsigned char   rgnFill[332];           // zeros
};
//
//  Orbital data record
//
struct OrbitalData {
    double      rAltitude;          // distance from surface
    double      rRange;             // range to radar target
    double      rPhaseRange;        // range calculated from radar signal phase
    double      fTrueAnomaly;       // angle from pericenter, spacecraft longitude
    double      tPericenter;        // time from pericenter, seconds
    double      fPhi;               // angle from orbital plane, latitude
    double      x, y, z;            // coordinates of spacecraft
    double      vx, vy, vz;         // velocity  of spacecraft
    double      fDopplerPerDegree;
    double      fFramePerDoppler;   // delta of true anomaly for each Doppler channel
};
//
//  Delay-Doppler coordinates in units of the 127 x 20 radio hologram
//
struct DelayDoppler {
    double tDelay;                  // Radar signal delay, in units of 1/650000 Hz, the pulse time
    double fDoppler;                // Radar Doppler frequency shift, in units
};
//
//  Simple complex arithmetic class
//
struct Complex {
    double          re;
    double          im;

                    Complex(double f1, double f2) : re(f1), im(f2) {};
                    Complex(double f) : re(f), im(0.0) {};
                    Complex() : re(0.0f), im(0.0f) {};
    Complex         operator +(Complex c) { return Complex(re+c.re, im+c.im); };
    Complex         operator -(Complex c) { return Complex(re-c.re, im-c.im); };
    Complex         operator *(double f)  { return Complex(re*f, im*f); };
    Complex         operator *(Complex c) { return Complex(re*c.re - im*c.im, re*c.im + im*c.re); }
    double          abs()   { return sqrt(re*re + im*im); };
    double          power() { return     (re*re + im*im); };
};

RadarLook       rgrlSurvey[];               // full radar swath, interleaved from both tape recorders
Complex         rgcIQ[];
char            rgnBadSAR[];
char            rgnBadBand[];
extern double   A, B, C, fPericenter;       // curve fit to equal-phase curve in band-survey images
OrbitalData     rgod[];

void ProcessStripeImages();
void ProcessFrameImages();
void fft(Complex [], int);
void InitializeOrbits();
void InitializeOrbitalData();