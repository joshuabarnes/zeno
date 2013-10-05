#include "stdinc.h"
#include "mathfns.h"
#include "vectdefs.h"
#define SafeSelect TRUE
#include "phatbody.h"

#define pos      Pos(_p)
#define x        PosX(_p)
#define y        PosY(_p)
#define z        PosZ(_p)
#define vel      Vel(_p)
#define vx       VelX(_p)
#define vy       VelY(_p)
#define vz       VelZ(_p)
#define m        Mass(_p)
#define phi      Phi(_p)
#define acc      Acc(_p)
#define ax       AccX(_p)
#define ay       AccY(_p)
#define az       AccZ(_p)
#define smooth   Smooth(_p)
#define rho      Rho(_p)
#define entf     EntFunc(_p)
#define uint     Uintern(_p)
#define udot     UdotInt(_p)
#define udotrad  UdotRad(_p)
#define udotvis  UdotVis(_p)
#define tau      Tau(_p)
#define type     Type(_p)
#define birth    Birth(_p)
#define death    Death(_p)
#define key      Key(_p)
#define aux      Aux(_p)
#define auxv     AuxVec(_p)
#define auxvx    AuxVecX(_p)
#define auxvy    AuxVecY(_p)
#define auxvz    AuxVecZ(_p)

#define r        absv(pos)
#define R        rsqrt(x*x + y*y)
#define v        absv(vel)
#define vr       (dotvp(pos,vel)/r)
#define vt       rsqrt(dotvp(vel,vel) - rsqr(vr))
#define etot     (phi + 0.5*dotvp(vel,vel))
#define jx       (vy*z - vz*y)
#define jy       (vz*x - vx*z)
#define jz       (vx*y - vy*x)
#define jtot     rsqrt(jx*jx + jy*jy + jz*jz)

#define Value(b)  SelectFloat(b, phatbody[NewBodyFields+0].offset)

void extendbody(void)
{
    new_field(&phatbody[NewBodyFields+0], "f", "Value");
    new_field(&phatbody[NewBodyFields+1], NULL, NULL);
}

void computemap(bodyptr _q, bodyptr _p, real t, int i, int n)
{
    Value(_q) = (x);
}
