import numpy as np


class VoltageSource:
    def __init__(self, source_type, amplitude, frequency=0, phase=0, rate=0):
        self.source_type = source_type
        self.amplitude = amplitude
        self.frequency = frequency      # Hz
        self.phase = phase              # degrees
        self.rate = rate                # exponential rate in 1/s


def source_response(source, a, b, c):
    # Find the source part of the answer for a*y'' + b*y' + c*y = Vs(t).
    # Return its equation, value at t=0, and slope at t=0.
    amplitude = source.amplitude

    if source.source_type == 'DC':
        voltage = amplitude / c
        return str(voltage), voltage, 0

    elif source.source_type == 'Sinusoidal':
        omega = 2 * np.pi * source.frequency
        phase = np.deg2rad(source.phase)

        # The answer has the form A*sin(omega*t + phase) + B*cos(omega*t + phase).
        sine_factor = c - a * omega ** 2
        cosine_factor = b * omega
        denominator = sine_factor ** 2 + cosine_factor ** 2
        A = amplitude * sine_factor / denominator
        B = -amplitude * cosine_factor / denominator

        equation = f'({A})*sin({omega}*t + {phase}) + ({B})*cos({omega}*t + {phase})'
        initial_value = A * np.sin(phase) + B * np.cos(phase)
        initial_slope = omega * (A * np.cos(phase) - B * np.sin(phase))
        return equation, initial_value, initial_slope

    else:  # Exponential source: amplitude * exp(rate*t)
        rate = source.rate
        denominator = a * rate ** 2 + b * rate + c
        scale = max(abs(a * rate ** 2), abs(b * rate), abs(c))

        if abs(denominator) > 1e-12 * scale:
            A = amplitude / denominator
            equation = f'({A})*exp({rate}*t)'
            return equation, A, rate * A

        # When the source rate matches a natural decay rate, include t.
        slope = 2 * a * rate + b
        slope_scale = max(abs(2 * a * rate), abs(b))
        if abs(slope) > 1e-12 * slope_scale:
            A = amplitude / slope
            equation = f'({A})*t*exp({rate}*t)'
            return equation, 0, A

        # A repeated natural decay rate needs t squared.
        A = amplitude / (2 * a)
        equation = f'({A})*t**2*exp({rate}*t)'
        return equation, 0, 0


class Circuit:
    def show_response(self, name, equation, unit):
        print(name + '(t) = ' + equation + ' ' + unit)


class RCCircuit(Circuit):
    def __init__(self, R, C, source, V0):
        self.R = R
        self.C = C
        self.source = source
        self.V0 = V0

    def calculate_response(self):
        tau = self.R * self.C

        # RC*Vc' + Vc = Vs(t)
        source_equation, source_start, source_slope = source_response(self.source, 0, tau, 1)

        # Choose A so the complete answer starts at V0.
        A = self.V0 - source_start
        equation = source_equation + f' + ({A})*exp(-t/{tau})'

        print('Time Constant =', tau, 's')
        self.show_response('Vc', equation, 'V')
        return equation


class RLCircuit(Circuit):
    # This class represents a series RL circuit.
    def __init__(self, R, L, source, I0):
        self.R = R
        self.L = L
        self.source = source
        self.I0 = I0

    def calculate_response(self):
        tau = self.L / self.R

        # L*iL' + R*iL = Vs(t)
        source_equation, source_start, source_slope = source_response(self.source, 0, self.L, self.R)

        A = self.I0 - source_start
        equation = source_equation + f' + ({A})*exp(-t/{tau})'

        print('Time Constant =', tau, 's')
        self.show_response('iL', equation, 'A')
        return equation


class SeriesRLC(Circuit):
    def __init__(self, R, L, C, source, V0, I0):
        self.R = R
        self.L = L
        self.C = C
        self.source = source
        self.V0 = V0
        self.I0 = I0

    def calculate_response(self):
        alpha = self.R / (2 * self.L)
        omega0 = 1 / np.sqrt(self.L * self.C)

        # LC*Vc'' + RC*Vc' + Vc = Vs(t)
        source_equation, source_start, source_slope = source_response(
            self.source, self.L * self.C, self.R * self.C, 1)

        # The full response must start at V0 with slope I0/C.
        voltage_start = self.V0 - source_start
        slope_start = self.I0 / self.C - source_slope

        print('alpha =', alpha)
        print('omega0 =', omega0)

        if np.isclose(alpha, omega0, rtol=1e-12, atol=0):
            print('Response Type = Critically Damped')
            A = voltage_start
            B = slope_start + alpha * A
            natural_equation = f'(({A}) + ({B})*t)*exp(-{alpha}*t)'

        elif alpha > omega0:
            print('Response Type = Overdamped')
            s1 = -alpha + np.sqrt(alpha ** 2 - omega0 ** 2)
            s2 = -alpha - np.sqrt(alpha ** 2 - omega0 ** 2)
            A = (slope_start - s2 * voltage_start) / (s1 - s2)
            B = voltage_start - A
            natural_equation = f'({A})*exp({s1}*t) + ({B})*exp({s2}*t)'

        else:
            print('Response Type = Underdamped')
            omega_d = np.sqrt(omega0 ** 2 - alpha ** 2)
            A = voltage_start
            B = (slope_start + alpha * A) / omega_d
            natural_equation = f'exp(-{alpha}*t)*(({A})*cos({omega_d}*t) + ({B})*sin({omega_d}*t))'

        equation = source_equation + ' + ' + natural_equation
        self.show_response('Vc', equation, 'V')
        return equation


def read_source():
    print('Source: 1. Constant DC  2. Sinusoidal  3. Exponential')
    choice = input('Choose source: ')

    if choice == '1':
        voltage = float(input('Enter Vs (V): '))
        return VoltageSource('DC', voltage)

    elif choice == '2':
        amplitude = float(input('Enter amplitude (V): '))
        frequency = float(input('Enter frequency (Hz): '))
        phase = float(input('Enter phase (degrees): '))
        if frequency < 0:
            print('Frequency must be nonnegative.')
            return None
        return VoltageSource('Sinusoidal', amplitude, frequency, phase)

    elif choice == '3':
        amplitude = float(input('Enter amplitude (V): '))
        rate = float(input('Enter rate k (1/s), negative for decay: '))
        return VoltageSource('Exponential', amplitude, rate=rate)

    else:
        print('Invalid source choice.')
        return None


def main():
    while True:
        print('\n1. RC Circuit')
        print('2. RL Circuit')
        print('3. Series RLC Circuit')
        print('4. Exit')
        choice = input('Enter your choice: ')

        if choice == '4':
            print('Exiting program.')
            break
        if choice not in ('1', '2', '3'):
            print('Invalid choice.')
            continue

        R = float(input('Enter R (ohm): '))
        if R <= 0:
            print('R must be positive.')
            continue

        if choice == '1':
            C = float(input('Enter C (F): '))
            if C <= 0:
                print('C must be positive.')
                continue
            source = read_source()
            if source is None:
                continue
            V0 = float(input('Enter V0 (V): '))
            circuit = RCCircuit(R, C, source, V0)

        elif choice == '2':
            L = float(input('Enter L (H): '))
            if L <= 0:
                print('L must be positive.')
                continue
            source = read_source()
            if source is None:
                continue
            I0 = float(input('Enter I0 (A): '))
            circuit = RLCircuit(R, L, source, I0)

        else:
            L = float(input('Enter L (H): '))
            C = float(input('Enter C (F): '))
            if L <= 0 or C <= 0:
                print('L and C must be positive.')
                continue
            source = read_source()
            if source is None:
                continue
            V0 = float(input('Enter V0 (V): '))
            I0 = float(input('Enter I0 (A): '))
            circuit = SeriesRLC(R, L, C, source, V0, I0)

        circuit.calculate_response()


if __name__ == '__main__':
    main()
