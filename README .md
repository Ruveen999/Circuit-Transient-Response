# Circuit Transient Response Calculator

A simple Python program for series RC, RL, and RLC circuits. It prints the
capacitor-voltage or inductor-current equation, including initial conditions.

## Run

Requires Python 3 and NumPy:

```bash
python -m pip install numpy
python main.py
```

Choose a circuit, enter component values, select a source, and enter the initial
voltage/current. Choose `4` to exit.

## Sources

| Source | Voltage for t >= 0 | Inputs |
| --- | --- | --- |
| DC | Vs(t) = V | Voltage in volts |
| Sinusoidal | Vs(t) = A sin(2πft + phase) | Amplitude in volts, frequency in Hz, phase in degrees |
| Exponential | Vs(t) = A exp(kt) | Amplitude in volts, rate k in inverse seconds |

Negative `k` means exponential decay; positive `k` means growth. Zero means DC.
Use positive R, L, and C values in ohms, henries, and farads. Frequency must be
nonnegative. Initial capacitor voltage is in volts and initial inductor current
is in amperes.

Enter finite numbers at number prompts. There is no exception handling: entering
text instead of a number stops the program. Simple `if` conditions check menu
choices and component signs. Very large values can overflow NumPy calculations.

## Code Structure

- `VoltageSource` stores the source type and its values using a plain constructor.
- `Circuit` is a regular parent class with a method to print the answer.
- `RCCircuit`, `RLCircuit`, and `SeriesRLC` inherit from `Circuit`. Each has a
  constructor and its own `calculate_response()` method, like the original code.
- `source_response()` is an ordinary function that calculates the source part of
  the answer. It returns the equation, its starting value, and its starting slope.
- `read_source()` reads source details. `main()` runs the menu.

There are no abstract classes, decorators, exception handlers, or complex numbers
in the application. NumPy provides sine, cosine, square root, and phase conversion.

Start reading at `main()`, then read `RCCircuit.calculate_response()`. Read the
source function next, followed by RL and RLC. The circuit methods follow three steps:

1. Calculate the part caused by the source.
2. Calculate constants so the complete answer has the entered initial conditions.
3. Print the equation.

## Circuit Equations

```text
RC:         RC*Vc' + Vc = Vs(t)
RL:         L*iL' + R*iL = Vs(t)
Series RLC: LC*Vc'' + RC*Vc' + Vc = Vs(t)
```

The source function uses `a*y'' + b*y' + c*y = Vs(t)` so these three circuits
can share the same source calculation. For RC, for example, `a=0`, `b=RC`, `c=1`.

For a sinusoidal source, the source part is `A*sin(omega*t + phase) +
B*cos(omega*t + phase)`. For an exponential source, it is normally
`A*exp(k*t)`. When `k` matches a natural decay rate, it needs a factor of `t`;
a repeated match needs `t**2`. Comparisons use a relative tolerance of `1e-12`.

RC and RL have one natural decay term. Series RLC has overdamped, critically
damped, and underdamped cases, calculated with separate `if`/`elif`/`else` branches.
The RLC initial voltage slope is `I0/C`.

The printed equations use `exp(x)` for e to the power x and `**2` for squaring.
All equations use real numbers. Components are ideal and the output applies for
`t >= 0`. Near-matching roots can be numerically sensitive. Plotting and CSV
export are not included.

## Example: RC Charging

Choose circuit `1` and enter `R=1000`, `C=0.0001`. Choose DC source `1`, enter
`Vs=10`, and set `V0=0`.

The program prints a time constant of `0.1 s` and an equation equivalent to:

```text
Vc(t) = 10 - 10*exp(-t/0.1) V
```

For an exponential example, use `R=1`, `C=1`, amplitude `2`, rate `-1`, and
`V0=0`. The result is `Vc(t) = 2*t*exp(-t)` V.

## Tests

```bash
python -m unittest discover -s tests -v
```

The five tests check DC examples, initial conditions, the circuit differential
equations, matching exponential rates, menu conditions, and all nine circuit/source
combinations. The test file uses extra testing tools; these are separate from the
basic application code in `main.py`.
