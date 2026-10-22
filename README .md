# Circuit Transient Response Calculator

A Python command-line program that calculates analytical transient-response equations for **RC, RL, and series RLC circuits** with a constant DC, sinusoidal, or exponential voltage source and specified initial conditions.

This project combines circuit analysis with object-oriented programming. Each circuit class inherits from a common `Circuit` class and implements its own `calculate_response()` method.

## Features

- Calculate the time constant and capacitor-voltage equation for an RC circuit.
- Calculate the time constant and inductor-current equation for an RL circuit.
- Classify a series RLC response as overdamped, critically damped, or underdamped.
- Calculate the series RLC capacitor-voltage equation using the initial voltage and current.
- Choose a DC, sinusoidal, or exponential voltage source for each circuit.
- Include initial conditions and matching exponential decay rates in the analytical response.
- Enter component values through an interactive menu and perform multiple calculations.

## Requirements

- Python 3.10 or newer
- NumPy

## Installation and Usage

Open a terminal in the repository folder and install NumPy:

```bash
python -m pip install numpy
```

Run the program:

```bash
python main.py
```

Choose a circuit from the menu:

```text
1. RC Circuit
2. RL Circuit
3. Series RLC Circuit
4. Exit
```

After entering the component values, select a source:

```text
1. Constant DC: Vs(t) = Vs
2. Sinusoidal:  Vs(t) = A sin(2π f t + φ)
3. Exponential: Vs(t) = A exp(k t)
```

For a sinusoid, enter amplitude in volts, frequency in Hz, and phase in degrees.
Frequency must be nonnegative. For an exponential, enter amplitude in volts and
rate `k` in inverse seconds: negative means decay, positive means growth, and zero
means DC. Sources apply for `t >= 0`; initial voltage/current are specified at `t = 0`.

Enter values in SI units:

| Input | Meaning | Unit |
| --- | --- | --- |
| `R` | Resistance | ohm |
| `L` | Inductance | henry |
| `C` | Capacitance | farad |
| `Vs` | Constant source voltage for t ≥ 0 | volt |
| `V0` | Initial capacitor voltage | volt |
| `I0` | Initial inductor current | ampere |

For example, enter `0.0001` or `1e-4` for a 100 µF capacitor. Use positive values for resistance, inductance, and capacitance.

## Time-Varying Sources

The program solves these equations with the entered initial conditions:

```text
RC:         RC Vc′ + Vc = Vs(t)
RL:         L iL′ + R iL = Vs(t)
Series RLC: LC Vc″ + RC Vc′ + Vc = Vs(t)
```

The printed result uses `Re[ ... ]`, meaning take the real part, and `j` denotes
the imaginary unit. Each term has the form `coefficient * t^n * exp(exponent*t)`.
Complex exponentials express sinusoidal responses compactly. The result includes
both the forced response and the natural response needed to satisfy initial
conditions. If an exponential source rate matches a natural root, the forced
term includes `t` (or `t^2` for a repeated root). Root comparisons use a relative
tolerance of `1e-12`.

Example: RC with `R=1`, `C=1`, source `2 exp(-t)`, and `V0=0` gives
`Vc(t) = 2 t exp(-t)` V. Sinusoidal responses include the startup transient,
not just the steady-state sinusoid.

For programmatic evaluation, import the circuit classes and `VoltageSource`:

```python
from main import RCCircuit, VoltageSource

source = VoltageSource.sinusoidal(amplitude=10, frequency=50, phase_degrees=0)
response = RCCircuit(R=1000, C=1e-4, Vs=source, V0=0).response()
print(response.equation())
print(response.value(0.01))  # capacitor voltage at 10 ms
```

Passing a number as `Vs` continues to mean a constant DC voltage.

## Circuit Equations for DC Sources

### RC Circuit

The capacitor voltage approaches the source voltage with time constant τ = RC:

```text
Vc(t) = Vs + (V0 − Vs)e^(−t/τ)
```

### RL Circuit

The inductor current approaches Vs/R with time constant τ = L/R:

```text
iL(t) = Vs/R + (I0 − Vs/R)e^(−t/τ)
```

### Series RLC Circuit

The damping factor and undamped natural angular frequency are:

```text
α = R/(2L)
ω₀ = 1/√(LC)
```

The capacitor-voltage response depends on their relationship:

| Condition | Response | Capacitor voltage |
| --- | --- | --- |
| α > ω₀ | Overdamped | Vs + A·e^(s₁t) + B·e^(s₂t) |
| α = ω₀ | Critically damped | Vs + (A + Bt)e^(−αt) |
| α < ω₀ | Underdamped | Vs + e^(−αt)[A·cos(ωdt) + B·sin(ωdt)] |

Here, s₁ and s₂ = −α ± √(α² − ω₀²), and ωd = √(ω₀² − α²). The program calculates A and B from Vc(0) = V0 and dVc/dt at t = 0 equal to I0/C.

## Example: RC Charging

Choose circuit option `1`, enter R and C, then choose source option `1` (DC) and enter Vs and V0:

```text
R  = 1000 ohm
C  = 0.0001 F
Vs = 10 V
V0 = 0 V
```

The calculated result is:

```text
Time Constant = 0.1 s
Vc(t) = Re[ (10)*exp((0)*t) + (-10)*exp((-10)*t) ] V
```

The capacitor starts at 0 V and approaches 10 V. At one time constant, its voltage is approximately 6.32 V. The program prints the equation; the value at a specific time can be evaluated separately.

## Object-Oriented Design

| Class | Role |
| --- | --- |
| `Circuit` | Abstract base class (`ABC`); requires subclasses to implement `response()` |
| `RCCircuit` | RC capacitor-voltage response |
| `RLCircuit` | RL inductor-current response |
| `SeriesRLC` | Series RLC capacitor-voltage response and damping classification |

The circuit subclasses implement the abstract `response()` method. The shared
`calculate_response()` method prints the result; subclasses extend it with time
constants or damping information using `super()`.

`VoltageSource` is the parent of `DCSource`, `SinusoidalSource`, and
`ExponentialSource`. Circuits accept any of these source objects, so the menu
can select a source without changing the circuit equations. `Response` is a
regular class with an explicit constructor, methods for evaluation, and equation
formatting.

NumPy supplies exponentials, complex square roots, phase conversion, approximate
comparisons, and array evaluation. For example, `response.value(np.linspace(0, 1, 100))`
returns 100 response values (after `import numpy as np`).

Invalid inputs raise `ValueError`. The menu catches these exceptions and lets the
user retry. Numerical range errors are handled separately; NumPy overflow during
response evaluation raises `FloatingPointError`. EOF or Ctrl+C exits cleanly.

## Assumptions and Current Limitations

- Components are ideal. Sources use the selected DC, sinusoidal, or exponential form for t ≥ 0.
- The RC and RL models use a series resistor with the energy-storage element; the RLC model is a series circuit.
- Initial conditions are specified at t = 0. Series RLC current is positive in the direction that charges the capacitor's positive terminal.
- The program prints analytical equations; it does not yet plot responses or generate time-series data.
- CLI inputs must be finite numbers; R, L, and C must be positive. Invalid input returns to the circuit menu.
- Critical damping and source/root matching use a relative tolerance of `1e-12`; extremely close roots can be numerically sensitive.
- Very large exponential growth rates or evaluation times can exceed floating-point range; evaluation raises a numerical exception.

## Planned Improvements

- Plot voltage and current against time.
- Calculate response values at user-selected times.
- Export calculated response data to CSV.

## Tests

Run from the repository root:

```bash
python -m unittest discover -s tests -v
```

Tests check DC regressions, initial conditions, the governing differential equations,
matching exponential rates (including repeated roots), invalid components, and all
nine circuit/source combinations through the CLI.
