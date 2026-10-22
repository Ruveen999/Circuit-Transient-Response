"""Analytical responses of series RC, RL and RLC circuits for t >= 0."""
from abc import ABC, abstractmethod
import math
import numpy as np


class VoltageSource:
    """Represent a voltage as Re[coefficient * exp(exponent * t)]."""

    def __init__(self, coefficient, exponent=0):
        if not np.isfinite(coefficient) or not np.isfinite(exponent):
            raise ValueError('Source parameters must be finite.')
        self.coefficient = coefficient
        self.exponent = exponent

    def value(self, t):
        with np.errstate(over='raise', invalid='raise'):
            return np.real(self.coefficient * np.exp(self.exponent * np.asarray(t)))

    @classmethod
    def sinusoidal(cls, amplitude, frequency, phase_degrees=0):
        return SinusoidalSource(amplitude, frequency, phase_degrees)

    @classmethod
    def exponential(cls, amplitude, rate):
        return ExponentialSource(amplitude, rate)


class DCSource(VoltageSource):
    def __init__(self, voltage):
        super().__init__(voltage)


class SinusoidalSource(VoltageSource):
    def __init__(self, amplitude, frequency, phase_degrees=0):
        if not all(np.isfinite(v) for v in (amplitude, frequency, phase_degrees)):
            raise ValueError('Source parameters must be finite.')
        if frequency < 0:
            raise ValueError('Frequency must be nonnegative.')
        self.amplitude = amplitude
        self.frequency = frequency
        self.phase_degrees = phase_degrees
        coefficient = amplitude * np.exp(1j * np.deg2rad(phase_degrees)) / 1j
        super().__init__(coefficient, 2j * np.pi * frequency)


class ExponentialSource(VoltageSource):
    def __init__(self, amplitude, rate):
        self.amplitude = amplitude
        self.rate = rate
        super().__init__(amplitude, rate)


class Response:
    """Store analytical terms and evaluate them at one time or a NumPy array."""

    def __init__(self, terms, label, unit):
        # Each term is (coefficient, power of t, exponential rate).
        self.terms = terms
        self.label = label
        self.unit = unit

    def value(self, t, derivative=0):
        if not isinstance(derivative, int) or derivative < 0:
            raise ValueError('Derivative order must be a nonnegative integer.')
        times = np.asarray(t, dtype=float)
        if not np.all(np.isfinite(times)) or np.any(times < 0):
            raise ValueError('Evaluation times must be finite and nonnegative.')
        total = np.zeros_like(times, dtype=complex)
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            for coefficient, power, exponent in self.terms:
                for j in range(min(power, derivative) + 1):
                    total += (coefficient * math.comb(derivative, j)
                              * math.factorial(power) / math.factorial(power - j)
                              * times ** (power - j) * exponent ** (derivative - j)
                              * np.exp(exponent * times))
        return total.real.item() if total.ndim == 0 else total.real

    def equation(self):
        def number(z):
            z = complex(z)
            return f'{z.real:.12g}' if z.imag == 0 else f'{z:.12g}'
        parts = []
        for coefficient, power, exponent in self.terms:
            polynomial = '' if power == 0 else ('*t' if power == 1 else f'*t^{power}')
            parts.append(f'({number(coefficient)}){polynomial}*exp(({number(exponent)})*t)')
        return f'{self.label}(t) = Re[ ' + ' + '.join(parts) + f' ] {self.unit}'



class Circuit(ABC):
    """Shared validation and output; subclasses provide the circuit response."""

    @abstractmethod
    def response(self):
        """Return the complete analytical response, including initial conditions."""
        pass

    def calculate_response(self):
        response = self.response()
        print(response.equation())
        return response

    @staticmethod
    def _particular(source, a2, a1, a0, gain, label, unit):
        """Solve P(D)y = gain*source, including simple and repeated root matching."""
        s = source.exponent
        p = a2 * s * s + a1 * s + a0
        scale = max(abs(a2 * s * s), abs(a1 * s), abs(a0))
        if abs(p) > 1e-12 * scale:
            term = (gain * source.coefficient / p, 0, s)
        else:
            dp = 2 * a2 * s + a1
            if abs(dp) > 1e-12 * max(abs(2 * a2 * s), abs(a1)):
                term = (gain * source.coefficient / dp, 1, s)
            else:
                term = (gain * source.coefficient / (2 * a2), 2, s)
        return Response([term], label, unit)

    @staticmethod
    def positive(*values):
        if any(not np.isfinite(v) or v <= 0 for v in values):
            raise ValueError('Component values must be finite and positive.')

    @staticmethod
    def source(value):
        return value if isinstance(value, VoltageSource) else DCSource(value)


class RCCircuit(Circuit):
    def __init__(self, R, C, Vs, V0):
        self.positive(R, C)
        self.R, self.C, self.Vs, self.V0 = R, C, self.source(Vs), V0

    def response(self):
        tau = self.R * self.C
        result = self._particular(self.Vs, 0, tau, 1, 1, 'Vc', 'V')
        result.terms.append((self.V0 - result.value(0), 0, -1 / tau))
        return result

    def calculate_response(self):
        print('Time Constant =', self.R * self.C, 's')
        return super().calculate_response()


class RLCircuit(Circuit):
    """Series RL circuit (name retained for compatibility)."""
    def __init__(self, R, L, Vs, I0):
        self.positive(R, L)
        self.R, self.L, self.Vs, self.I0 = R, L, self.source(Vs), I0

    def response(self):
        result = self._particular(self.Vs, 0, self.L, self.R, 1, 'iL', 'A')
        result.terms.append((self.I0 - result.value(0), 0, -self.R / self.L))
        return result

    def calculate_response(self):
        print('Time Constant =', self.L / self.R, 's')
        return super().calculate_response()


class SeriesRLC(Circuit):
    def __init__(self, R, L, C, Vs, V0, I0):
        self.positive(R, L, C)
        self.R, self.L, self.C = R, L, C
        self.Vs, self.V0, self.I0 = self.source(Vs), V0, I0

    def damping(self):
        alpha = self.R / (2 * self.L)
        omega0 = 1 / np.sqrt(self.L * self.C)
        if np.isclose(alpha, omega0, rtol=1e-12, atol=0):
            return 'Critically Damped', alpha, omega0
        return ('Overdamped' if alpha > omega0 else 'Underdamped'), alpha, omega0

    def response(self):
        result = self._particular(self.Vs, self.L * self.C, self.R * self.C,
                            1, 1, 'Vc', 'V')
        y0 = self.V0 - result.value(0)
        dy0 = self.I0 / self.C - result.value(0, 1)
        kind, alpha, omega0 = self.damping()
        if kind == 'Critically Damped':
            result.terms.extend([(y0, 0, -alpha), (dy0 + alpha * y0, 1, -alpha)])
        else:
            delta = np.emath.sqrt(alpha * alpha - omega0 * omega0)
            s1, s2 = -alpha + delta, -alpha - delta
            A = (dy0 - s2 * y0) / (s1 - s2)
            result.terms.extend([(A, 0, s1), (y0 - A, 0, s2)])
        return result

    def calculate_response(self):
        kind, alpha, omega0 = self.damping()
        print('alpha =', alpha)
        print('omega0 =', omega0)
        print('Response Type =', kind)
        return super().calculate_response()


def read_number(prompt):
    value = float(input(prompt))
    if not np.isfinite(value):
        raise ValueError('Inputs must be finite numbers.')
    return value


def read_source():
    print('Source: 1. Constant DC  2. Sinusoidal  3. Exponential')
    choice = input('Choose source: ').strip()
    if choice == '1':
        return DCSource(read_number('Enter Vs (V): '))
    if choice == '2':
        return SinusoidalSource(read_number('Enter amplitude (V): '),
                                        read_number('Enter frequency (Hz): '),
                                        read_number('Enter phase (degrees): '))
    if choice == '3':
        return ExponentialSource(read_number('Enter amplitude (V): '),
                                         read_number('Enter rate k (1/s), negative for decay: '))
    raise ValueError('Invalid source choice.')


def main():
    while True:
        print('\n1. RC Circuit\n2. RL Circuit\n3. Series RLC Circuit\n4. Exit')
        try:
            choice = input('Enter your choice: ').strip()
            if choice == '4':
                print('Exiting program.')
                break
            if choice not in ('1', '2', '3'):
                print('Invalid choice.')
                continue
            R = read_number('Enter R (ohm): ')
            if choice == '1':
                C = read_number('Enter C (F): ')
                source = read_source()
                circuit = RCCircuit(R, C, source, read_number('Enter V0 (V): '))
            elif choice == '2':
                L = read_number('Enter L (H): ')
                source = read_source()
                circuit = RLCircuit(R, L, source, read_number('Enter I0 (A): '))
            else:
                L = read_number('Enter L (H): ')
                C = read_number('Enter C (F): ')
                source = read_source()
                circuit = SeriesRLC(R, L, C, source, read_number('Enter V0 (V): '),
                                    read_number('Enter I0 (A): '))
            circuit.calculate_response()
        except ValueError as error:
            print('Invalid input:', error)
        except (FloatingPointError, OverflowError, ZeroDivisionError) as error:
            print('Calculation failed: values exceed the supported numerical range.', error)
        except (EOFError, KeyboardInterrupt):
            print('Exiting program.')
            break


if __name__ == '__main__':
    main()
