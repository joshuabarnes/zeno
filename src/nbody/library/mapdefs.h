/*
 * MAPDEFS.H: Mapping between identifiers used in expressions and
 *            corresponding macro names.
 */

local string mapdefs[][2] = {
    { "pos",     "Pos"     },
    { "x",       "PosX"    },
    { "y",       "PosY"    },
    { "z",       "PosZ"    },
    { "vel",     "Vel"     },
    { "vx",      "VelX"    },
    { "vy",      "VelY"    },
    { "vz",      "VelZ"    },
    { "m",       "Mass"    },
    { "phi",     "Phi"     },
    { "acc",     "Acc"     },
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
    { "auxv",    "AuxVec"  },
    { "auxvx",   "AuxVecX" },
    { "auxvy",   "AuxVecY" },
    { "auxvz",   "AuxVecZ" },
    { NULL,      NULL      }
};

/*
 * Derived quantities.
 */

local string expdefs[][2] = {
    { "r",       "absv(pos)"                        },
    { "R",       "rsqrt(x*x + y*y)"                 },
    { "v",       "absv(vel)"                        },
    { "vr",      "(dotvp(pos,vel)/r)"               },
    { "vt",      "rsqrt(dotvp(vel,vel) - rsqr(vr))" },
    { "etot",    "(phi + 0.5*dotvp(vel,vel))"       },
    { "jx",      "(vy*z - vz*y)"                    },
    { "jy",      "(vz*x - vx*z)"                    },
    { "jz",      "(vx*y - vy*x)"                    },
    { "jtot",    "rsqrt(jx*jx + jy*jy + jz*jz)"     },
    { "taux",    "(ay*z - az*y)"                    },
    { "tauy",    "(az*x - ax*z)"                    },
    { "tauz",    "(ax*y - ay*x)"                    },
    { "tautot",  "rsqrt(ax*ax + ay*ay + az*az)"     },
    { NULL,       NULL                              }
};
