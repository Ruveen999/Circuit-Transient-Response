# Circuit Transient Response Calculator

A Python command-line program that calculates analytical transient-response equations for **RC, RL, and series RLC circuits** with a constant DC source and specified initial conditions.

This project combines circuit analysis with object-oriented programming. Each circuit class inherits from a common `Circuit` class and implements its own `calculate_response()` method.

## Features

- Calculate the time constant and capacitor-voltage equation for an RC circuit.
- Calculate the time constant and inductor-current equation for an RL circuit.
- Classify a series RLC response as overdamped, critically damped, or underdamped.
- Calculate the series RLC capacitor-voltage equation using the initial voltage and current.
- Enter component values through an interactive menu and perform multiple calculations.

## Requirements

- Python 3
- NumPy

## Installation and Usage

Download the repository and open a terminal in its folder. Save the program as `main.py`, then install NumPy:

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

## Circuit Equations

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

Choose option `1` and enter:

```text
R  = 1000 ohm
C  = 0.0001 F
Vs = 10 V
V0 = 0 V
```

The calculated result is:

```text
Time Constant = 0.1 s
Vc(t) = 10 − 10e^(−t/0.1) V
```

The capacitor starts at 0 V and approaches 10 V. At one time constant, its voltage is approximately 6.32 V. The program prints the equation; the value at a specific time can be evaluated separately.

## Object-Oriented Design

| Class | Role |
| --- | --- |
| `Circuit` | Common parent class with a placeholder response method |
| `RCCircuit` | RC capacitor-voltage response |
| `RLCircuit` | RL inductor-current response |
| `SeriesRLC` | Series RLC capacitor-voltage response and damping classification |

The subclasses override the same method, demonstrating inheritance and method overriding.

## Assumptions and Current Limitations

- Components are ideal, and the source voltage remains constant for t ≥ 0.
- The RC and RL models use a series resistor with the energy-storage element; the RLC model is a series circuit.
- Initial conditions are specified at t = 0. Series RLC current is positive in the direction that charges the capacitor's positive terminal.
- The program prints analytical equations; it does not yet plot responses or generate time-series data.
- Inputs are not currently validated. Non-numeric entries or zero component values can cause errors.
- Critical damping currently uses exact floating-point equality. A tolerance-based comparison would improve classification near the critical condition.

## Planned Improvements

- Validate component values and handle invalid input.
- Use a tolerance-based critical-damping comparison.
- Plot voltage and current against time.
- Calculate response values at user-selected times.
- Export calculated response data to CSV.
