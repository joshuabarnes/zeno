/*
 * BODYTAGS.H: definitions for N-body snapshot files.
 */

#ifndef _bodytags_h
#define _bodytags_h

#define BodyTag  "Body"

/* 
 * Item tags for SnapShot components.
 */

#define SnapShotTag		"SnapShot"

#define   ParametersTag		"Parameters"
#define     NBodyTag            "NBody"
#define     NGasTag		"NGas"
#define     TimeTag             "Time"

#define   ParticlesTag		"Particles"
#define     PosTag              "Position"
#define     VelTag              "Velocity"
#define     MassTag             "Mass"
#define	    PhiTag	        "Potential"
#define	    AccTag	        "Acceleration"
#define     SmoothTag           "SmoothLength"
#define     RhoTag              "Density"
#define     EntFuncTag		"EntropyFunc"
#define     UinternTag		"Uinternal"
#define     UdotIntTag		"UdotInternal"
#define     UdotRadTag		"UdotRadiation"
#define     UdotVisTag          "UdotViscosity"
#define     TauTag              "OpticalDepth"
#define     TypeTag             "BodyType"
#define     BirthTag		"BirthDate"
#define     DeathTag		"DeathDate"
#define     KeyTag              "Key"
#define     AuxTag              "Aux"
#define     AuxVecTag           "AuxVec"

#define   DiagnosticsTag	"Diagnostics"
#define     EnergyTag           "Energy"
#define     KETensorTag         "KETensor"
#define     PETensorTag		"PETensor"
#define     AMVectorTag		"AMVector"
#define     CMPhaseTag	        "CMPhaseSpace"

/*
 * Tags from older implementations.
 */

#define     NobjTag             "Nobj"
#define     PhaseTag            "PhaseSpace"

#endif /* ! _bodytags_h */
