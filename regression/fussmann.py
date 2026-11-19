"""
Analytic prompt redeposition (Fussmann's model) and the self-sputtering
chain, for checking Flan's sputtering runs.

Prompt redeposition
-------------------
A sputtered atom leaves a wall with speed v and ionizes after a path L.
The ion it becomes gyrates about the field with gyroradius rho = v / Omega
(Omega = q B / m for the singly charged ion). If its first gyro-orbit
reaches the wall, it lands there within one gyration: a prompt return.
For a field parallel to the wall, a cosine distribution of directions and
a fixed L, the fraction that comes back promptly is

    r(x) = 1 / (1 + x^2),     x = L / rho.

Both L and rho go as v, so with L the mean ionization length,
lambda = v / (n_e S), the ratio is independent of the speed:

    p = lambda / rho = Omega / (n_e S) = q B / (m n_e S).

The ionization length is exponentially distributed about lambda, so the
fraction Flan should see (it samples ionization as a Poisson process) is
r averaged over that distribution:

    <r>(p) = (1 / p) int_0^inf exp(-x / p) / (1 + x^2) dx.

1 / (1 + p^2) is only the fixed-length result. The two differ a lot at
p ~ 1 (0.50 against 0.62). The model leaves out an electric field (no
sheath), motion along B (none for B parallel to the wall), and further
ionization during the first gyration. Reimold (2011), section 2.1.1.3,
gives the model and the exponential average, and attributes the
derivation to R. Dux; the usual citation is Fussmann et al. (1995).
geometry_check() reproduces r(x) by Monte Carlo from the geometry, as a
check on the formula.

Self-sputtering chain
---------------------
If a fraction <r> of the sputtered atoms comes back promptly and each one
sputters Y_self atoms (and Y_D atoms are sputtered per incident ion), the
net yield, counting every generation of the chain, is (Reimold 2011, eq.
2.18, for one projectile species)

    Y_eff = (1 - <r>) Y_D / (1 - <r> Y_self),

which runs away when <r> Y_self >= 1.

References
----------
- G. Fussmann et al., "High-Z elements as target materials in fusion
  devices", in Plasma Physics and Controlled Nuclear Fusion Research 1994,
  IAEA, Vol. 2, p. 143 (1995).
- F. Reimold, "Numerical simulations of tungsten impurity transport in
  ASDEX-Upgrade with EIRENE", IPP Report 17/27 (2011).

Usage: python fussmann.py prints a table of <r>(p) and the geometry check.
"""
import numpy as np
from scipy import integrate

E_CHARGE = 1.602176634e-19	# C
AMU = 1.66053906660e-27		# kg


def p_param(B, mass_amu, ne, S, charge=1):
	"""
	Ratio of the mean ionization length to the gyroradius.

	Inputs
	 - B (float)        : Field strength (T).
	 - mass_amu (float) : Mass of the sputtered atom (amu).
	 - ne (float)       : Electron density (m^-3).
	 - S (float)        : Ionization rate coefficient of the atom (m^3/s).
	 - charge (int)     : Charge of the ion it becomes (1).

	Returns p = q B / (m n_e S).
	"""
	return charge * E_CHARGE * B / (mass_amu * AMU * ne * S)


def r_fixed(x):
	"""
	Fraction promptly redeposited for a fixed ionization length x rho:
	1 / (1 + x^2).
	"""
	x = np.asarray(x, dtype=float)
	return 1.0 / (1.0 + x * x)


def r_exp(p):
	"""
	Fraction promptly redeposited for an exponentially distributed
	ionization length with mean p rho: r_fixed averaged over it.

	Inputs
	 - p (float) : Mean ionization length over the gyroradius, > 0.
	"""
	val, _ = integrate.quad(lambda x: np.exp(-x / p) / (1.0 + x * x), 0.0,
		np.inf)
	return val / p


def geometry_check(x, n=1_000_000, seed=1):
	"""
	Fraction promptly redeposited for a fixed ionization length, by Monte
	Carlo from the geometry, to check r_fixed.

	The wall is z = 0 with its normal along +z and B along y. An atom
	leaves along (sin t cos f, sin t sin f, cos t), with cos^2 t uniform
	(a cosine distribution), and ionizes at height x rho cos t. The ion
	gyrates in the x-z plane with radius rho_perp = rho sqrt(1 - (sin t
	sin f)^2) about a center offset by rho sin t cos f in z (the sign of
	the offset depends on the charge and the field direction; f is
	uniform, so either covers both). It comes back if the orbit reaches
	z = 0.

	Inputs
	 - x (float)  : Ionization length over the gyroradius at full speed.
	 - n (int)    : Number of atoms.
	 - seed (int) : Random seed.
	"""
	rng = np.random.default_rng(seed)
	cos_t = np.sqrt(rng.random(n))
	sin_t = np.sqrt(1.0 - cos_t ** 2)
	f = 2.0 * np.pi * rng.random(n)
	vx = sin_t * np.cos(f)
	vperp = np.sqrt(vx ** 2 + cos_t ** 2)
	return np.mean(x * cos_t + vx < vperp)


def y_eff(r, y_d, y_self):
	"""
	Net yield with prompt redeposition and self-sputtering, summed over
	every generation: (1 - r) y_d / (1 - r y_self). Returns inf when
	r y_self >= 1 (runaway).
	"""
	if r * y_self >= 1.0:
		return np.inf
	return (1.0 - r) * y_d / (1.0 - r * y_self)


if __name__ == "__main__":
	print("   p    <r> (exponential)   1/(1+p^2) (fixed length)")
	for p in [0.1, 0.3, 1.0, 2.55, 3.0, 10.0]:
		print("{:6.2f}  {:10.4f}          {:10.4f}".format(p, r_exp(p),
			float(r_fixed(p))))
	print("\nGeometry check of r_fixed:")
	print("   x    Monte Carlo   1/(1+x^2)")
	for x in [0.3, 1.0, 3.0]:
		print("{:5.1f}   {:8.4f}    {:8.4f}".format(x, geometry_check(x),
			float(r_fixed(x))))
