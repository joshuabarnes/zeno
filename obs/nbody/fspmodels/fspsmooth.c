/*
 * FSPSMOOTH.C: smooth mass distribution of finite spherical profile.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "fsp.h"

string defv[] = {		";Smooth mass distribution of FSP",
    "in=???",			";Input file with FSP",
    "out=???",			";Output file for FSP",
    "eps=0.025",		";Plummer smoothing length",
    "kappa=1.5",		";Interpolation parameter",
    "VERSION=1.1",		";Josh Barnes  26 September 1994",
    NULL,
};

void fspsmooth(fsprof *, real, real);

int main(int argc, string argv[])
{
    stream istr, ostr;
    fsprof *fsp;

    initparam(argv, defv);
    istr = stropen(getparam("in"), "r");
    get_history(istr);
    fsp = get_fsprof(istr);
    fspsmooth(fsp, getdparam("eps"), getdparam("kappa"));
    ostr = stropen(getparam("out"), "w");
    put_history(ostr);
    put_fsprof(ostr, fsp);
    strclose(ostr);
    return (0);
}

void fspsmooth(fsprof *fsp, real eps, real kappa)
{
    int i;
    real c0, r, f, d;

    eprintf("[fspsmooth:  rho_eps(0) = %f  r[0]/eps = %f]\n",
	    fsp->density[0] * rpow(eps / fsp->radius[0], fsp->alpha),
	    fsp->radius[0] / eps);
    c0 = rpow(2.0/3.0, kappa/fsp->alpha);
    for (i = 0; i < fsp->npoint; i++) {
	r = fsp->radius[i];
	f = rpow(1 + c0 * rpow(eps/r, kappa), fsp->alpha/kappa);
	d = - fsp->alpha * c0 * (rpow(eps/r, kappa) / r) *
	      rpow(1 + c0 * rpow(eps/r, kappa), fsp->alpha/kappa - 1);
	fsp->density[i] = f * fsp->density[i] +
			    d * fsp->mass[i] / (4 * PI * rsqr(r));
	fsp->mass[i] = f *fsp->mass[i];
    }
    fsp->alpha = 0.0;
}
