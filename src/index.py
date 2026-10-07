from pathlib import Path

import numpy as np
from scipy.optimize import brentq

MU_EARTH = 398600.0  # km^3/s^2
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "outputs" / "text"
TAU = 2.0 * np.pi


def heading(problem: str) -> None:
    """Print the heading used to separate each problem's output."""
    print(f"\n------\n {problem}\n------")


def write_output(problem: str, lines: list[str]) -> None:
    """Print and save a problem's text output for inclusion in Typst."""
    output = "\n".join(lines) + "\n"
    print(output, end="")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / f"s{problem}.txt").write_text(output, encoding="utf-8")


def angle(first: np.ndarray, second: np.ndarray, normal: np.ndarray) -> float:
    """Return the oriented angle about a unit normal in [0, 2 pi)."""
    return np.arctan2(np.dot(normal, np.cross(first, second)), first @ second) % TAU


def geometry(r: np.ndarray, v: np.ndarray) -> dict:
    """Compute two-body invariants from Cartesian position and velocity."""
    radius = np.linalg.norm(r)
    momentum = np.cross(r, v)
    eccentricity = np.cross(v, momentum) / MU_EARTH - r / radius
    energy = v @ v / 2.0 - MU_EARTH / radius
    return {
        "r": radius,
        "vr": r @ v / radius,
        "H": momentum,
        "N": np.cross([0.0, 0.0, 1.0], momentum),
        "ev": eccentricity,
        "e": np.linalg.norm(eccentricity),
        "a": -MU_EARTH / (2.0 * energy),
        "energy": energy,
    }


def classical(g: dict, r: np.ndarray) -> dict:
    """Classical elements for noncircular, inclined ellipses and hyperbolas."""
    momentum, node, eccentricity = g["H"], g["N"], g["ev"]
    if np.linalg.norm(node) < 1e-12 or g["e"] < 1e-12:
        raise ValueError("Classical node/perigee angles are undefined.")
    normal = momentum / np.linalg.norm(momentum)
    theta = angle(eccentricity, r, normal)
    elements = {
        "e": g["e"],
        "a": g["a"],
        "i": np.arctan2(np.linalg.norm(momentum[:2]), momentum[2]),
        "Omega": np.arctan2(node[1], node[0]) % TAU,
        "omega": angle(node, eccentricity, normal),
        "theta": theta,
    }
    if g["e"] < 1.0:
        anomaly = (
            np.arctan2(
                np.sqrt(1.0 - g["e"] ** 2) * np.sin(theta), g["e"] + np.cos(theta)
            )
            % TAU
        )
        elements.update(E=anomaly, M=(anomaly - g["e"] * np.sin(anomaly)) % TAU)
    elif g["e"] > 1.0:
        signed_theta = np.arctan2(np.sin(theta), np.cos(theta))
        anomaly = hyperbolic_anomaly(g["e"], signed_theta)
        elements.update(
            theta_signed=signed_theta,
            F=anomaly,
            M_h=g["e"] * np.sinh(anomaly) - anomaly,
            varpi=(elements["Omega"] + elements["omega"]) % TAU,
        )
    else:
        raise ValueError("Parabolic orbits require a different element set.")
    return elements


def hyperbolic_anomaly(e: float, theta: float) -> float:
    """Signed hyperbolic anomaly F; inbound values are negative, never wrapped."""
    denominator = 1.0 + e * np.cos(theta)
    if e <= 1.0 or denominator <= 0.0:
        raise ValueError(
            "A true anomaly on the physical hyperbolic branch is required."
        )
    return np.arcsinh(np.sqrt(e * e - 1.0) * np.sin(theta) / denominator)


def equinoctial_basis(p: float, q: float) -> tuple[np.ndarray, np.ndarray]:
    """Orthonormal in-plane axes, nonsingular at i = 0."""
    scale = 1.0 + p * p + q * q
    f = np.array([1.0 - p * p + q * q, 2.0 * p * q, -2.0 * p]) / scale
    g = np.array([2.0 * p * q, 1.0 + p * p - q * q, 2.0 * q]) / scale
    return f, g


def equinoctial(g: dict, r: np.ndarray) -> dict:
    """Direct conversion for ellipses, with a stated hyperbolic extension.

    h = e sin(Omega + omega), k = e cos(Omega + omega),
    p = tan(i/2) sin(Omega), q = tan(i/2) cos(Omega).
    L and F are true and eccentric longitude; lambda is mean longitude.
    beta = 1 / (1 + sqrt(1 - h^2 - k^2)).
    For hyperbolas, lambda = varpi + M_h, with varpi = atan2(h, k)
    in [0, 2 pi), M_h = e sinh(F) - F, and no wrapping of lambda.
    This extends the mean-longitude definition, not the elliptic F/beta formulas.
    The posigrade convention is singular only at i = pi.
    """
    momentum = g["H"]
    denominator = np.linalg.norm(momentum) + momentum[2]
    if denominator <= 1e-12 or g["e"] == 1.0:
        raise ValueError(
            "A nonparabolic orbit with inclination below 180 deg is required."
        )
    p, q = momentum[0] / denominator, -momentum[1] / denominator
    f, b = equinoctial_basis(p, q)
    h, k = g["ev"] @ b, g["ev"] @ f
    longitude = np.arctan2(r @ b, r @ f) % TAU
    if g["e"] < 1.0:
        beta = 1.0 / (1.0 + np.sqrt(1.0 - h * h - k * k))
        eccentric_longitude = longitude + 2.0 * np.arctan2(
            h * np.cos(longitude) - k * np.sin(longitude),
            1.0 / beta + k * np.cos(longitude) + h * np.sin(longitude),
        )
        mean_longitude = (
            eccentric_longitude
            + h * np.cos(eccentric_longitude)
            - k * np.sin(eccentric_longitude)
        ) % TAU
    else:
        varpi = np.arctan2(h, k) % TAU
        anomaly = hyperbolic_anomaly(g["e"], longitude - varpi)
        mean_longitude = varpi + g["e"] * np.sinh(anomaly) - anomaly
    return {
        "a": g["a"],
        "h": h,
        "k": k,
        "p": p,
        "q": q,
        "lambda": mean_longitude,
        "L": longitude,
    }


def reconstruct(eq: dict) -> tuple[np.ndarray, np.ndarray]:
    """Invert elliptic elements or the unwrapped hyperbolic extension."""
    a, h, k, p, q, longitude = (eq[key] for key in ("a", "h", "k", "p", "q", "lambda"))
    eccentricity = np.hypot(h, k)
    if a < 0.0 and eccentricity > 1.0:
        varpi = np.arctan2(h, k) % TAU
        mean_anomaly = longitude - varpi
        # e sinh(F) - F is odd and strictly increasing for e > 1.
        bound = max(1.0, np.arcsinh(abs(mean_anomaly) / eccentricity) + 1.0)
        while eccentricity * np.sinh(bound) - bound < abs(mean_anomaly):
            bound *= 2.0
        anomaly = brentq(
            lambda f: eccentricity * np.sinh(f) - f - mean_anomaly,
            -bound,
            bound,
            xtol=1e-14,
        )
        scale = -a
        beta = np.sqrt(eccentricity**2 - 1.0)
        position = scale * np.array(
            [eccentricity - np.cosh(anomaly), beta * np.sinh(anomaly)]
        )
        rate = np.sqrt(MU_EARTH / scale**3) / (eccentricity * np.cosh(anomaly) - 1.0)
        velocity = scale * rate * np.array([-np.sinh(anomaly), beta * np.cosh(anomaly)])
        f, b = equinoctial_basis(p, q)
        # Perigee and transverse directions, expressed in the equinoctial basis.
        basis = np.column_stack(
            ((k * f + h * b) / eccentricity, (-h * f + k * b) / eccentricity)
        )
        return basis @ position, basis @ velocity
    if a <= 0.0 or eccentricity >= 1.0:
        raise ValueError(
            "Semimajor axis and eccentricity must describe a nonparabolic conic."
        )
    anomaly = brentq(
        lambda f: f + h * np.cos(f) - k * np.sin(f) - longitude,
        longitude - eccentricity - 1e-12,
        longitude + eccentricity + 1e-12,
        xtol=1e-14,
    )
    beta = 1.0 / (1.0 + np.sqrt(1.0 - h * h - k * k))
    transform = np.array(
        [
            [1.0 - beta * h * h, beta * h * k],
            [beta * h * k, 1.0 - beta * k * k],
        ]
    )
    position = a * (transform @ [np.cos(anomaly), np.sin(anomaly)] - [k, h])
    rate = np.sqrt(MU_EARTH / a**3) / (1.0 - k * np.cos(anomaly) - h * np.sin(anomaly))
    velocity = a * rate * (transform @ [-np.sin(anomaly), np.cos(anomaly)])
    basis = np.column_stack(equinoctial_basis(p, q))
    return basis @ position, basis @ velocity


def verify(
    r: np.ndarray, v: np.ndarray, g: dict, eq: dict, coe: dict | None
) -> list[str]:
    """Check state reconstruction and independent classical identities."""
    recovered_r, recovered_v = reconstruct(eq)
    np.testing.assert_allclose(recovered_r, r, rtol=0.0, atol=1e-7)
    np.testing.assert_allclose(recovered_v, v, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(
        np.dot(g["H"], g["H"]) / MU_EARTH,
        g["a"] * (1.0 - g["e"] ** 2),
        rtol=1e-12,
    )
    if coe is not None:
        varpi = coe["Omega"] + coe["omega"]
        expected = [
            coe["e"] * np.sin(varpi),
            coe["e"] * np.cos(varpi),
            np.tan(coe["i"] / 2.0) * np.sin(coe["Omega"]),
            np.tan(coe["i"] / 2.0) * np.cos(coe["Omega"]),
        ]
        np.testing.assert_allclose(
            [eq[key] for key in ("h", "k", "p", "q")], expected, atol=1e-12
        )
        if g["e"] < 1.0:
            delta = eq["lambda"] - varpi - coe["M"]
            np.testing.assert_allclose(
                [np.cos(delta), np.sin(delta)], [1.0, 0.0], atol=1e-12
            )
        else:
            np.testing.assert_allclose(
                eq["lambda"], varpi % TAU + coe["M_h"], atol=1e-12
            )
            # Independent state-vector relation for the signed hyperbolic anomaly.
            np.testing.assert_allclose(
                np.sinh(coe["F"]),
                (r @ v) / (g["e"] * np.sqrt(MU_EARTH * -g["a"])),
                atol=1e-12,
            )
    return [
        f"position reconstruction error = {np.linalg.norm(recovered_r - r):.3e} km",
        f"velocity reconstruction error = {np.linalg.norm(recovered_v - v):.3e} km/s",
    ]


def element_lines(elements: dict) -> list[str]:
    """Format numerical elements; the report reads these key/value lines."""
    angles = {
        "i",
        "Omega",
        "omega",
        "theta",
        "theta_signed",
        "E",
        "M",
        "lambda",
        "L",
        "varpi",
    }
    lines = []
    for key, value in elements.items():
        # if key == "lambda" and elements["a"] < 0.0:
        #     lines.append(f"lambda = {value:.10f} rad")
        # elif key in angles:
        if key in angles:
            lines.append(f"{key} = {np.degrees(value):.6f} deg")
        elif key == "a":
            lines.append(f"a = {value:.6f} km")
        else:
            lines.append(f"{key} = {value:.10f}")
    return lines


def solve_state(number: str, r: list[float], v: list[float], both: bool = True) -> None:
    """Save classical/equinoctial results and Cartesian diagnostics."""
    r, v = np.array(r, dtype=float), np.array(v, dtype=float)
    g = geometry(r, v)
    eq = equinoctial(g, r)
    coe = classical(g, r) if both else None
    diagnostics = [
        f"r = {g['r']:.6f} km",
        f"vr = {g['vr']:.6f} km/s",
        f"energy = {g['energy']:.9f} km^2/s^2",
        f"H vector = {np.array2string(g['H'], precision=6)} km^2/s",
        f"N vector = {np.array2string(g['N'], precision=6)} km^2/s",
        f"e vector = {np.array2string(g['ev'], precision=10)}",
    ]
    checks = verify(r, v, g, eq, coe)
    if both:
        assert coe is not None
        write_output(number + "a", element_lines(coe) + diagnostics)
        write_output(number + "b", element_lines(eq) + checks)
    else:
        write_output(
            number, element_lines(eq) + [f"e = {g['e']:.10f}"] + diagnostics + checks
        )


def problem_04() -> None:
    """Curtis 4.6: choose the inbound branch of true anomaly."""
    heading("p04")
    r = np.array([-6000.0, -1000.0, -5000.0])
    ev = np.array([-0.4, -0.5, -0.6])
    cosine = r @ ev / (np.linalg.norm(r) * np.linalg.norm(ev))
    principal = np.arccos(np.clip(cosine, -1.0, 1.0))
    theta = TAU - principal
    assert np.sin(theta) < 0.0
    write_output(
        "04",
        [
            f"r = {np.linalg.norm(r):.6f} km",
            f"e = {np.linalg.norm(ev):.10f}",
            f"e dot r = {ev @ r:.6f} km",
            f"cos theta = {cosine:.10f}",
            f"principal = {np.degrees(principal):.6f} deg",
            f"theta = {np.degrees(theta):.6f} deg",
        ],
    )


def problem_05() -> None:
    """Curtis 4.7: inward radial motion fixes the orbit normal's sign."""
    heading("p05")
    r = np.array([-6600.0, -1300.0, -5200.0])
    ev = np.array([-0.4, -0.5, -0.6])
    # e x r = e r sin(theta) H_hat; inbound means sin(theta) < 0.
    cross = np.cross(r, ev)
    normal = cross / np.linalg.norm(cross)
    inclination = np.arctan2(np.linalg.norm(normal[:2]), normal[2])
    # Reconstruct a compatible velocity to verify the selected direction.
    semilatus = np.linalg.norm(r) + ev @ r
    momentum = np.sqrt(MU_EARTH * semilatus)
    velocity = MU_EARTH / momentum * np.cross(normal, ev + r / np.linalg.norm(r))
    assert r @ velocity < 0.0
    np.testing.assert_allclose(np.cross(r, velocity) / momentum, normal, atol=1e-12)
    write_output(
        "05",
        [
            f"r cross e = {np.array2string(cross, precision=6)} km",
            f"cross magnitude = {np.linalg.norm(cross):.6f} km",
            f"cos i = {normal[2]:.10f}",
            f"i = {np.degrees(inclination):.6f} deg",
            f"compatible radial velocity = {r @ velocity / np.linalg.norm(r):.6f} km/s",
        ],
    )


def main() -> None:
    #########
    # p01 #
    #########
    heading("p01")
    solve_state("01", [2500, 16000, 4000], [-3, -1, 5])

    #########
    # p02 #
    #########
    heading("p02")
    solve_state("02", [0, 0, -13000], [4, 5, 6])

    #########
    # p03 #
    #########
    heading("p03")
    solve_state("03", [6500, -7500, -2500], [4, 3, -3])

    #########
    # p04 #
    #########
    problem_04()

    #########
    # p05 #
    #########
    problem_05()

    #########
    # p06 #
    #########
    heading("p06")
    solve_state("06a", [3698.3455, -34232.4473, 0], [2.1482, -1.4944, 0], both=False)
    solve_state("06b", [-21082.0848, 36515.242, 0], [-2.6627, -1.5373, 0], both=False)
    solve_state(
        "06c",
        [-15843.4562, -2247.7766, 21200.4627],
        [1.4544, -3.5188, 0.7138],
        both=False,
    )


if __name__ == "__main__":
    main()
