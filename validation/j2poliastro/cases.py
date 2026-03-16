from hapsira.core.perturbations import J2_perturbation
from hapsira.core.propagation import func_twobody
from astropy import units as u
from astropy.time import Time
from hapsira.bodies import Earth
from hapsira.twobody import Orbit
from hapsira.twobody.propagation import CowellPropagator
import numpy as np

print("Python constants")
print(f"mu={Earth.k.value}")
print(f"J2={Earth.J2.value}")
print(f"R={Earth.R.to(u.km).value}")


def f(t0, state, k):
    du_kep = func_twobody(t0, state, k)
    ax, ay, az = J2_perturbation(
        t0, state, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value
    )
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_ad


orbits = []
orbits.append(
    (90, Orbit.from_classical(
        Earth,
        a=6778.137 * u.km,
        ecc=0.001 * u.one,
        inc=51.6 * u.deg,
        raan=0 * u.deg,
        argp=0 * u.deg,
        nu=0 * u.deg,
    ))
)

orbits.append(
    (1440, Orbit.from_classical(
        Earth,
        a=6778.137 * u.km,
        ecc=0.001 * u.one,
        inc=51.6 * u.deg,
        raan=0 * u.deg,
        argp=0 * u.deg,
        nu=0 * u.deg,
    ))
)

orbits.append(
    (1440, Orbit.from_classical(
        Earth,
        a=126778.137 * u.km,
        ecc=0.5 * u.one,
        inc=51.6 * u.deg,
        raan=0 * u.deg,
        argp=0 * u.deg,
        nu=0 * u.deg,
    ))
)

for i, orbit in enumerate(orbits):
    print(f"Orbit {i + 1}, {orbit[0]}min")
    print("r0 =", orbit[1].r)
    print("v0 =", orbit[1].v)
    orb_j2 = orbit[1].propagate(orbit[0] * u.min, method=CowellPropagator(f=f))
    print("r1 =", orb_j2.r)
    print("v1 =", orb_j2.v)
