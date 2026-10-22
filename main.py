"""Analytical responses of series RC, RL and RLC circuits for t >= 0."""
import numpy as np


class VoltageSource:
    """Represent a voltage as Re[coefficient * exp(exponent * t)]."""

    def __init__(self, coefficient, exponent=0):
        self.coefficient = coefficient
        self.exponent = exponent

    def value(self, t):
        return np.real(self.coefficient * np.exp(self.exponent * np.asarray(t)))



class DCSource(VoltageSource):
    def __init__(self, voltage):
        super().__init__(voltage)


class SinusoidalSource(VoltageSource):
    def __init__(self, amplitude, frequency, phase_degrees=0):
        self.amplitude = amplitude
        self.frequency = frequency
        self.phase_degrees = phase_degrees
        # sin(angle) = the real part of -j * exp(j * angle).
        # Python writes the imaginary unit as 1j.
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
        times = np.asarray(t, dtype=float)
        total = np.zeros_like(times, dtype=complex)
        for coefficient, power, rate in self.terms:
            # The polynomial can only be 1, t, or t squared.
            if power == 0:
                polynomial = np.ones_like(times)
                first_derivative = 0
                second_derivative = 0
            elif power == 1:
                polynomial = times
                first_derivative = 1
                second_derivative = 0
            else:
                polynomial = times ** 2
                first_derivative = 2 * times
                second_derivative = 2

            # Product rule for polynomial * exp(rate * t).
            if derivative == 1:
                polynomial = first_derivative + rate * polynomial
            elif derivative == 2:
                polynomial = (second_derivative + 2 * rate * first_derivative
                              + rate ** 2 * polynomial)

            total += coefficient * polynomial * np.exp(rate * times)

        # Complex numbers are used to represent sine waves. The physical
        # voltage/current is always the real part of the answer.
        if total.ndim == 0:
            return total.real.item()
        return total.real

    def equation(self):
        parts = []
        for coefficient, power, rate in self.terms:
            term = f'({coefficient:.6g})'
            if power == 1:
                term += '*t'
            elif power == 2:
                term += '*t^2'
            term += f'*exp(({rate:.6g})*t)'
            parts.append(term)

        expression = ' + '.join(parts)
        return f'{self.label}(t) = Re[ {expression} ] {self.unit}'



def make_source(value):
    # A number means DC; a source object can also be passed directly.
    if isinstance(value, VoltageSource):
        return value
    return DCSource(value)


def forced_response(source, second_order, first_order, constant, label, unit):
    """Find the source-driven part of the response.

    The circuit equation is:
    second_order*y'' + first_order*y' + constant*y = source voltage.
    Try y = coefficient * exp(rate*t), then substitute into that equation.
    """
    rate = source.exponent
    amplitude = source.coefficient
    denominator = second_order * rate ** 2 + first_order * rate + constant
    scale = max(abs(second_order * rate ** 2), abs(first_order * rate), abs(constant))

    if abs(denominator) > 1e-12 * scale:
        coefficient = amplitude / denominator
        power = 0
    else:
        # Matching a natural root requires an extra factor of t.
        slope = 2 * second_order * rate + first_order
        slope_scale = max(abs(2 * second_order * rate), abs(first_order))
        if abs(slope) > 1e-12 * slope_scale:
            coefficient = amplitude / slope
            power = 1
        else:
            # A repeated natural root requires t squared instead.
            coefficient = amplitude / (2 * second_order)
            power = 2

    return Response([(coefficient, power, rate)], label, unit)


class Circuit:
    """Parent class with the output method shared by the three circuits."""

    def calculate_response(self):
        response = self.response()
        print(response.equation())
        return response



class RCCircuit(Circuit):
    def __init__(self, R, C, Vs, V0):
        self.R = R
        self.C = C
        self.Vs = make_source(Vs)
        self.V0 = V0

    def response(self):
        tau = self.R * self.C
        result = forced_response(self.Vs, 0, tau, 1, 'Vc', 'V')
        # Add a decaying term so the capacitor starts at V0.
        initial_difference = self.V0 - result.value(0)
        result.terms.append((initial_difference, 0, -1 / tau))
        return result

    def calculate_response(self):
        print('Time Constant =', self.R * self.C, 's')
        return super().calculate_response()


class RLCircuit(Circuit):
    """Series RL circuit (name retained for compatibility)."""
    def __init__(self, R, L, Vs, I0):
        self.R = R
        self.L = L
        self.Vs = make_source(Vs)
        self.I0 = I0

    def response(self):
        result = forced_response(self.Vs, 0, self.L, self.R, 'iL', 'A')
        # Add a decaying term so the inductor current starts at I0.
        initial_difference = self.I0 - result.value(0)
        result.terms.append((initial_difference, 0, -self.R / self.L))
        return result

    def calculate_response(self):
        print('Time Constant =', self.L / self.R, 's')
        return super().calculate_response()


class SeriesRLC(Circuit):
    def __init__(self, R, L, C, Vs, V0, I0):
        self.R = R
        self.L = L
        self.C = C
        self.Vs = make_source(Vs)
        self.V0 = V0
        self.I0 = I0

    def damping(self):
        alpha = self.R / (2 * self.L)
        omega0 = 1 / np.sqrt(self.L * self.C)
        if np.isclose(alpha, omega0, rtol=1e-12, atol=0):
            return 'Critically Damped', alpha, omega0
        if alpha > omega0:
            return 'Overdamped', alpha, omega0
        return 'Underdamped', alpha, omega0

    def response(self):
        result = forced_response(self.Vs, self.L * self.C, self.R * self.C, 1, 'Vc', 'V')
        # Subtract the forced response from the required initial conditions.
        voltage_difference = self.V0 - result.value(0)
        slope_difference = self.I0 / self.C - result.value(0, 1)
        kind, alpha, omega0 = self.damping()
        if kind == 'Critically Damped':
            result.terms.extend([(voltage_difference, 0, -alpha),
                                 (slope_difference + alpha * voltage_difference, 1, -alpha)])
        else:
            delta = np.emath.sqrt(alpha * alpha - omega0 * omega0)
            s1, s2 = -alpha + delta, -alpha - delta
            A = (slope_difference - s2 * voltage_difference) / (s1 - s2)
            B = voltage_difference - A
            result.terms.append((A, 0, s1))
            result.terms.append((B, 0, s2))
        return result

    def calculate_response(self):
        kind, alpha, omega0 = self.damping()
        print('alpha =', alpha)
        print('omega0 =', omega0)
        print('Response Type =', kind)
        return super().calculate_response()


def read_number(prompt):
    return float(input(prompt))


def read_source():
    print('Source: 1. Constant DC  2. Sinusoidal  3. Exponential')
    choice = input('Choose source: ').strip()
    if choice == '1':
        return DCSource(read_number('Enter Vs (V): '))
    if choice == '2':
        amplitude = read_number('Enter amplitude (V): ')
        frequency = read_number('Enter frequency (Hz): ')
        phase = read_number('Enter phase (degrees): ')
        if frequency < 0:
            print('Frequency must be nonnegative.')
            return None
        return SinusoidalSource(amplitude, frequency, phase)
    if choice == '3':
        amplitude = read_number('Enter amplitude (V): ')
        rate = read_number('Enter rate k (1/s), negative for decay: ')
        return ExponentialSource(amplitude, rate)
    print('Invalid source choice.')
    return None


def main():
    while True:
        print('\n1. RC Circuit\n2. RL Circuit\n3. Series RLC Circuit\n4. Exit')
        choice = input('Enter your choice: ').strip()
        if choice == '4':
            print('Exiting program.')
            break
        if choice not in ('1', '2', '3'):
            print('Invalid choice.')
            continue
        R = read_number('Enter R (ohm): ')
        if R <= 0:
            print('R must be positive.')
            continue
        if choice == '1':
            C = read_number('Enter C (F): ')
            if C <= 0:
                print('C must be positive.')
                continue
            source = read_source()
            if source is None:
                continue
            V0 = read_number('Enter V0 (V): ')
            circuit = RCCircuit(R, C, source, V0)
        elif choice == '2':
            L = read_number('Enter L (H): ')
            if L <= 0:
                print('L must be positive.')
                continue
            source = read_source()
            if source is None:
                continue
            I0 = read_number('Enter I0 (A): ')
            circuit = RLCircuit(R, L, source, I0)
        else:
            L = read_number('Enter L (H): ')
            if L <= 0:
                print('L must be positive.')
                continue
            C = read_number('Enter C (F): ')
            if C <= 0:
                print('C must be positive.')
                continue
            source = read_source()
            if source is None:
                continue
            V0 = read_number('Enter V0 (V): ')
            I0 = read_number('Enter I0 (A): ')
            circuit = SeriesRLC(R, L, C, source, V0, I0)
        circuit.calculate_response()


if __name__ == '__main__':
    main()
