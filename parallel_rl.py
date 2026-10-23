# Transient response of a parallel RL circuit (R and L in parallel).
# The source can be a voltage source or a current source.

import numpy as np


def rnd(x):
    # This fun rounds up every numerical value upto 4decimal points throughout the whole program.
    return float(round(x, 4))


class parallel_RL():
    # This class stores the circuit values and calculates the response.
    #   R, L        : resistor (ohm) and inductor (H) connected in parallel
    #   Rs          : series resistance of a voltage source (ohm)
    #                 Rs = 0 means an ideal voltage source connected directly across R and L
    #   I0          : inductor current at t = 0 (A)
    #   source_type : "voltage" or "current"
    #   wave        : "dc"  -> K
    #                 "exp" -> K*e^(-a*t)
    #                 "sin" -> K*sin(w*t + phase)
    #   K           : source amplitude (V for a voltage source, A for a current source)
    #   a           : decay rate of an exponential source (1/s), negative means growth
    #   w           : angular frequency of a sinusoidal source (rad/s)
    #   phase_deg   : phase of a sinusoidal source (degrees)

    def __init__(self, R, L, Rs, I0, source_type, wave, K, a, w, phase_deg):
    # here we are saving the values into the variables we need in the program and also checking if the values are correct or not.
        if R <= 0 or L <= 0:
            raise ValueError("R and L must be positive")
        if Rs < 0:
            raise ValueError("Rs cannot be negative")
        if wave == "sin" and w <= 0:
            raise ValueError("w must be greater than 0 for a sinusoidal source")

        self.R = R
        self.L = L
        self.Rs = Rs if source_type == "voltage" else 0   # Rs is only used by a voltage source
        self.I0 = I0
        self.source_type = source_type
        self.wave = wave
        self.K = K
        self.a = a
        self.w = w
        self.phase = np.radians(phase_deg)   

    def calculate_response(self):
        # Chooses the correct formulas for the selected source and prints
        # iL(t) (inductor current), v(t) (voltage across R and L) and iR(t) (resistor current).
        # An ideal voltage source (Rs = 0) has no exponential transient,
        # so it is handled in its own function below.
        if self.source_type == "voltage" and self.Rs == 0:
            self.ideal_voltage_response()
            return

        R, L, I0 = self.R, self.L, self.I0
        K, a, w, ph = self.K, self.a, self.w, self.phase

        # Step 1: turn the source into an equivalent current source Keq
        # feeding Rth in parallel with L, so every case uses the same equation
        if self.source_type == "current":
            Rth = R                              # the inductor only sees R
            Keq = K                              # the source current itself
        else:
            Rth = R * self.Rs / (R + self.Rs)    # Rs in parallel with R
            Keq = K / self.Rs                    # voltage source converted to a current source

        alpha = Rth / L          # decay rate of the circuit (alpha = Rth/L = 1/tau)
        tau = 1 / alpha          # time constant that is L/R

        print("Source =", self.source_type, "source, wave type =", self.wave)
        print("Rth =", rnd(Rth), "ohm")
        print("alpha =", rnd(alpha), " tau =", rnd(tau), "s")

        # Step 2: choose the formulas depending on the waveform of the source
        if self.wave == "dc":
            print("Response Type = Constant DC source (natural + forced)")
            A = I0 - Keq                         # iL(0) = Keq + A
            Av = Rth * (Keq - I0)                # v(t) = L*diL/dt
            print("iL(t) =", rnd(Keq), "+ (", rnd(A), ")*e^(", rnd(-alpha), "t) A")
            print("v(t)  = (", rnd(Av), ")*e^(", rnd(-alpha), "t) V")
            print("iR(t) = (", rnd(Av / R), ")*e^(", rnd(-alpha), "t) A")

        elif self.wave == "exp":
            if np.isclose(a, alpha):
                # source decays at the same rate as the circuit, so a t*e^(-alpha*t) term appears
                # (np.isclose is used instead of == because these are decimal numbers) (yo chai aba ai ley prefer gareko)
                print("Response Type = Exponential source, decay rate equal to alpha")
                Av = Rth * (Keq - I0)
                Bv = -Rth * alpha * Keq 
                print("iL(t) = (", rnd(I0), "+ (", rnd(alpha * Keq), ")t )*e^(", rnd(-alpha), "t) A")
                print("v(t)  = (", rnd(Av), "+ (", rnd(Bv), ")t )*e^(", rnd(-alpha), "t) V")
                print("iR(t) = (", rnd(Av / R), "+ (", rnd(Bv / R), ")t )*e^(", rnd(-alpha), "t) A")
            else:
                print("Response Type = Exponential source (natural + forced)")
                C = alpha * Keq / (alpha - a)    # forced part of iL
                print("iL(t) = (", rnd(C), ")*e^(", rnd(-a), "t) + (", rnd(I0 - C), ")*e^(", rnd(-alpha), "t) A")
                print("v(t)  = (", rnd(Rth * (Keq - C)), ")*e^(", rnd(-a), "t) + (", rnd(-Rth * (I0 - C)), ")*e^(", rnd(-alpha), "t) V")
                print("iR(t) = (", rnd(Rth * (Keq - C) / R), ")*e^(", rnd(-a), "t) + (", rnd(-Rth * (I0 - C) / R), ")*e^(", rnd(-alpha), "t) A")

        else:
            print("Response Type = Sinusoidal source (transient + steady state)")
            M = alpha * Keq / np.hypot(alpha, w)   # amplitude of the steady state iL, hypot(x, y) = sqrt(x*x + y*y)
            theta = ph - np.arctan2(w, alpha)      # phase of the steady state iL
            B = I0 - M * np.sin(theta)             # coefficient of the decaying transient
            deg = np.degrees(theta)                # numpy converts radians back to degrees for printing
            print("iL(t) = (", rnd(M), ")sin(", rnd(w), "t + (", rnd(deg), "deg)) + (", rnd(B), ")*e^(", rnd(-alpha), "t) A")
            print("v(t)  = (", rnd(L * M * w), ")cos(", rnd(w), "t + (", rnd(deg), "deg)) + (", rnd(-Rth * B), ")*e^(", rnd(-alpha), "t) V")
            print("iR(t) = (", rnd(L * M * w / R), ")cos(", rnd(w), "t + (", rnd(deg), "deg)) + (", rnd(-Rth * B / R), ")*e^(", rnd(-alpha), "t) A")

    def ideal_voltage_response(self):
        # Used only when Rs = 0 (ideal voltage source directly across R and L)
        # The voltage across L is forced by the source: v(t) = source
        # So iL(t) = I0 + (1/L) * integral of v(t), and iR(t) = v(t)/R
        R, L, I0 = self.R, self.L, self.I0
        K, a, w, ph = self.K, self.a, self.w, self.phase

        print("Source = ideal voltage source (Rs = 0), wave type =", self.wave)
        print("Response Type = Forced response only (no exponential transient)")

        if self.wave == "dc" or (self.wave == "exp" and np.isclose(a, 0)):
            print("iL(t) =", rnd(I0), "+ (", rnd(K / L), ")t A")
            print("v(t)  =", rnd(K), "V")
            print("iR(t) =", rnd(K / R), "A")
            print("Note: iL keeps growing, in a real circuit the wire/source resistance would limit it")

        elif self.wave == "exp":
            D = K / (L * a)                      # integral of K*e^(-a*t) divided by L
            print("iL(t) = (", rnd(I0 + D), ") + (", rnd(-D), ")*e^(", rnd(-a), "t) A")
            print("v(t)  = (", rnd(K), ")*e^(", rnd(-a), "t) V")
            print("iR(t) = (", rnd(K / R), ")*e^(", rnd(-a), "t) A")

        else:
            deg = np.degrees(ph)                 # phase in degrees for printing
            print("iL(t) = (", rnd(I0 + K * np.cos(ph) / (L * w)), ") + (", rnd(-K / (L * w)), ")cos(", rnd(w), "t + (", rnd(deg), "deg)) A")
            print("v(t)  = (", rnd(K), ")sin(", rnd(w), "t + (", rnd(deg), "deg)) V")
            print("iR(t) = (", rnd(K / R), ")sin(", rnd(w), "t + (", rnd(deg), "deg)) A")


def ask_choice(prompt, options):
    # to show a menu
    print(prompt)
    for key, label in options.items():
        print("  " + key + ".", label)
    while True:
        choice = input("Enter choice: ").strip()
        if choice in options:
            return choice
        print("Invalid choice, try again.")


def ask_number(prompt, default=None):
    # Asks for a number and keeps asking until the user types a valid one
    # iff a default is given pressing Enter uses it
    if default is not None:
        prompt = prompt + " [" + str(default) + "]"
    while True:
        text = input(prompt + ": ").strip()
        if text == "" and default is not None:
            return default
        try:
            return float(text)
        except ValueError:
            print("Please enter a number.")


def main():
    # Asks the user for the source and circuit values, then prints the response
    print(" Parallel RL circuit transient response \n")

    source_type = {"1": "voltage", "2": "current"}[
        ask_choice("Source type:", {"1": "Voltage source", "2": "Current source"})]

    wave = {"1": "dc", "2": "exp", "3": "sin"}[
        ask_choice("\nSource waveform:", {
            "1": "Constant DC       K",
            "2": "Exponential       K*e^(-a*t)",
            "3": "Sinusoidal        K*sin(w*t + phase)"})]

    unit = "V" if source_type == "voltage" else "A"
    print()
    K = ask_number("Amplitude K (" + unit + ")")
    a, w, phase_deg = 0.0, 0.0, 0.0
    if wave == "exp":
        a = ask_number("Decay rate a (1/s, negative for growth)")
    if wave == "sin":
        w = ask_number("Angular frequency w (rad/s)")
        phase_deg = ask_number("Phase (degrees)", 0)

    R = ask_number("Resistance R (ohm)")
    L = ask_number("Inductance L (H)")
    Rs = 0
    if source_type == "voltage":
        Rs = ask_number("Series resistance Rs of the source (ohm, 0 = ideal source across R and L)", 0)
    I0 = ask_number("Initial inductor current I0 (A)", 0)

    try:
        circuit = parallel_RL(R, L, Rs, I0, source_type, wave, K, a, w, phase_deg)
    except ValueError as err:
        print("\nInvalid input:", err)
        return

    print()
    circuit.calculate_response()


main()