/*
 * SNAPSET.C: set body variables to given expressions.
 */

#include "stdinc.h"
#include "getparam.h"
#include "vectdefs.h"
#include "bodytags.h"

string defv[] = {               ";Set body variables in snapshot",
    "in=???",                   ";Input snapshot file name",
    "out=???",                  ";Output snapshot file name",
    "times=all",                ";Range of times to process",
    "x=",			";Expression for X position",
    "y=",			";Expression for Y position",
    "z=",			";Expression for Z position",
    "vx=",			";Expression for X velocity",
    "vy=",			";Expression for Y velocity",
    "vz=",			";Expression for Z velocity",
    "m=",			";Expression for mass",
    "phi=",                     ";Expression for potential",
    "ax=",			";Expression for X acceleration",
    "ay=",			";Expression for Y acceleration",
    "az=",			";Expression for Z acceleration",
    "smooth=",			";Expression for smoothing length",
    "rho=",			";Expression for density",
    "entf=",			";Expression for entropy function",
    "uint=",			";Expression for internal energy",
    "udot=",			";Expression for derivative of uint",
    "udotrad=",			";Expression for radiative losses",
    "udotvis=",			";Expression for viscous heating",
    "tau=",			";Expression for optical depth",
    "type=",			";Expression for body type",
    "birth=",			";Expression for birth date",
    "death=",			";Expression for death date",
    "key=",			";Expression for key value",
    "aux=",                     ";Expression for aux value",
    "auxvx=",                   ";Expression for X aux vector",
    "auxvy=",                   ";Expression for Y aux vector",
    "auxvz=",                   ";Expression for Z aux vector",
    "require=",                 ";List of input items necessary",
    "produce=",                 ";List of output items produced",
    "passall=true",		";If true, pass on all input data",
    "seed=",			";Seed for random number generator",
    "VERSION=2.7",              ";Josh Barnes  21 May 2012",
    NULL,
};

/*
 * Define mapping from expressions for body fields to macro names.
 * Note: this table is derived from ~/nbody/library/mapdefs.h, and
 * should be updated whenever new standard fields are defined.
 */

local string mapdefs[][2] = {
    { "x",       "PosX"    },
    { "y",       "PosY"    },
    { "z",       "PosZ"    },
    { "vx",      "VelX"    },
    { "vy",      "VelY"    },
    { "vz",      "VelZ"    },
    { "m",       "Mass"    },
    { "phi",     "Phi"     },
    { "ax",      "AccX"    },
    { "ay",      "AccY"    },
    { "az",      "AccZ"    },
    { "smooth",  "Smooth"  },
    { "rho",     "Rho"     },
    { "entf",    "EntFunc" },
    { "uint",    "Uintern" },
    { "udot",    "UdotInt" },
    { "udotrad", "UdotRad" },
    { "udotvis", "UdotVis" },
    { "tau",     "Tau"     },
    { "type",    "Type"    },
    { "birth",   "Birth"   },
    { "death",   "Death"   },
    { "key",     "Key"     },
    { "aux",     "Aux"     },
    { "auxvx",   "AuxVecX" },
    { "auxvy",   "AuxVecY" },
    { "auxvz",   "AuxVecZ" },
    { NULL,      NULL      }
};

#define MAXMAP  (sizeof(mapdefs) / (2 * sizeof(string)))

void buildmap(string, string *, string *, string *, string, int);

void execmap(string);

int main(int argc, string argv[])
{
    string prog, names[MAXMAP], exprs[MAXMAP];
    int i, j;

    initparam(argv, defv);
    prog = tempnam("/tmp", "sm");
    for (i = j = 0; mapdefs[i][0] != NULL; i++) {
	if (! strnull(getparam(mapdefs[i][0]))) {
	    exprs[j] = getparam(mapdefs[i][0]);
	    names[j] = mapdefs[i][1];
	    j++;
	}
    }
    exprs[j] = names[j] = NULL;
    eprintf("[%s: defined %d variable%s]\n", getargv0(), j, j!=1 ? "s" : "");
    buildmap(prog, names, exprs, NULL, Precision, NDIM);
    execmap(prog);
    if (unlink(prog) != 0)
        error("%s: can't unlink %s\n", getargv0(), prog);
    return (0);
}

#include <unistd.h>
#include <sys/wait.h>

void execmap(string prog)
{
    int mappid, mapstat;
    char histbuf[512];

    mappid = fork();
    if (mappid == 0) {                          /* if this is child process */
	sprintf(histbuf, "HISTORY=%s", *ask_history());
        execl(prog, getargv0(), getparam("in"), getparam("out"),
	      getparam("times"), getparam("require"), getparam("produce"),
	      getparam("passall"), getparam("seed"), histbuf, NULL);
        error("%s: execl %s failed\n", getargv0(), prog);
    }
    while (wait(&mapstat) != mappid)
        eprintf("[%s: waiting on subprocess %d]\n", getargv0(), mappid);
}
