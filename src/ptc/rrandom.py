"""R's default random number generator, so that results match the R package.

R uses the Mersenne-Twister generator, seeded by ``set.seed()`` with its own
scrambling, and draws normal numbers by inversion (``qnorm`` of a uniform
with extra precision). This module reproduces ``set.seed(seed)`` followed by
``rnorm(n)`` exactly, using numpy's Mersenne-Twister for the 32-bit draws.
"""

import math

import numpy as np

_BIG = 134217728  # 2^27, as in R's snorm.c (INVERSION)
_I2_32M1 = 2.328306437080797e-10  # 1 / (2^32 - 1)


class RRandom:
    """Random numbers as produced by R after ``set.seed(seed)``."""

    def __init__(self, seed):
        self.set_seed(seed)

    def set_seed(self, seed):
        # Same scrambling as R's set.seed (RNG.c: Randomize / RNG_Init)
        s = int(seed) & 0xFFFFFFFF
        for _ in range(50):
            s = (69069 * s + 1) & 0xFFFFFFFF
        key = np.empty(625, dtype=np.uint32)
        for j in range(625):
            s = (69069 * s + 1) & 0xFFFFFFFF
            key[j] = s
        # key[0] is R's position counter (mti); R sets it to 624, so the
        # first draw regenerates the state, exactly like numpy with pos 624
        self._bitgen = np.random.MT19937(0)
        self._bitgen.state = {"bit_generator": "MT19937",
                              "state": {"key": key[1:], "pos": 624}}

    def unif_rand(self, n):
        """n uniform numbers in (0, 1), as R's unif_rand()."""
        u = self._bitgen.random_raw(n).astype(np.float64) * 2.3283064365386963e-10
        # R's fixup(): keep the values strictly inside (0, 1)
        u[u <= 0.0] = 0.5 * _I2_32M1
        u[1.0 - u <= 0.0] = 1.0 - 0.5 * _I2_32M1
        return u

    def rnorm(self, n):
        """n standard normal numbers, as R's rnorm(n)."""
        u = self.unif_rand(2 * n)
        # R draws two uniforms per normal: u = (int)(BIG*u1) + u2
        p = (np.floor(_BIG * u[0::2]) + u[1::2]) / _BIG
        return qnorm(p)


def qnorm(p):
    """Standard normal quantiles with R's algorithm (AS 241, qnorm.c)."""
    p = np.asarray(p, dtype=np.float64)
    q = p - 0.5
    val = np.empty_like(p)

    central = np.abs(q) <= 0.425
    if central.any():
        qc = q[central]
        r = 0.180625 - qc * qc
        val[central] = (qc * (((((((r * 2509.0809287301226727 +
                                    33430.575583588128105) * r + 67265.770927008700853) * r +
                                  45921.953931549871457) * r + 13731.693765509461125) * r +
                                1971.5909503065514427) * r + 133.14166789178437745) * r +
                              3.387132872796366608)
                        / (((((((r * 5226.495278852854561 +
                                 28729.085735721942674) * r + 39307.89580009271061) * r +
                               21213.794301586595867) * r + 5394.1960214247511077) * r +
                             687.1870074920579083) * r + 42.313330701600911252) * r + 1.0))

    tail = ~central
    if tail.any():
        pt = p[tail]
        qt = q[tail]
        # r = min(p, 1 - p); R computes 1 - p as (0.5 - p + 0.5)
        r = np.where(qt > 0, 0.5 - pt + 0.5, pt)
        # math.log is the C library log() that R uses; numpy's log can differ
        # in the last bit
        r = np.sqrt(-np.array([math.log(x) for x in r]))
        v = np.empty_like(r)
        near = r <= 5.0
        if near.any():
            rr = r[near] - 1.6
            v[near] = ((((((((rr * 7.7454501427834140764e-4 +
                              .0227238449892691845833) * rr + .24178072517745061177) *
                            rr + 1.27045825245236838258) * rr +
                           3.64784832476320460504) * rr + 5.7694972214606914055) *
                         rr + 4.6303378461565452959) * rr +
                        1.42343711074968357734)
                       / (((((((rr *
                                1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                               rr + .0151986665636164571966) * rr +
                              .14810397642748007459) * rr + .68976733498510000455) *
                            rr + 1.6763848301838038494) * rr +
                           2.05319162663775882187) * rr + 1.0))
        far = ~near
        if far.any():
            rr = r[far] - 5.0
            v[far] = ((((((((rr * 2.01033439929228813265e-7 +
                             2.71155556874348757815e-5) * rr +
                            .0012426609473880784386) * rr + .026532189526576123093) *
                          rr + .29656057182850489123) * rr +
                         1.7848265399172913358) * rr + 5.4637849111641143699) *
                       rr + 6.6579046435011037772)
                      / (((((((rr *
                               2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                              rr + 1.8463183175100546818e-5) * rr +
                             7.868691311456132591e-4) * rr + .0148753612908506148525)
                           * rr + .13692988092273580531) * rr +
                          .59983220655588793769) * rr + 1.0))
        v = np.where(qt < 0.0, -v, v)
        val[tail] = v
    return val
