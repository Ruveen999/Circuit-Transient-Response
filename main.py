import numpy as np

class Circuit:
    def calculate_response(self):
        pass
#I made three classes for RC, RL and RLC circuit.

class RCCircuit(Circuit):
    def __init__(self, R, C, Vs, V0):
        self.R = R
        self.C = C
        self.Vs = Vs
        self.V0 = V0

    def calculate_response(self):
        tau = self.R * self.C
        Vf = self.Vs
        A = self.V0 - Vf
        print("Time Constant =", tau, "s")
        print("Vc(t) =", Vf, "+ (", A, ") * e^(-t/", tau, ") V")


class RLCircuit(Circuit):
    def __init__(self, R, L, Vs, I0):
        self.R = R
        self.L = L
        self.Vs = Vs
        self.I0 = I0

    def calculate_response(self):
        tau = self.L / self.R
        If = self.Vs / self.R
        A = self.I0 - If
        print("Time Constant =", tau, "s")
        print("iL(t) =", If, "+ (", A, ") * e^(-t/", tau, ") A")


class SeriesRLC(Circuit):
    def __init__(self, R, L, C, Vs, V0, I0):
        self.R = R
        self.L = L
        self.C = C
        self.Vs = Vs
        self.V0 = V0
        self.I0 = I0

    def calculate_response(self):
        alpha = self.R / (2 * self.L)   # formula for alpha is = R/2L
        omega0 = 1 / np.sqrt(self.L * self.C) #formula for omega(angular frequency) is 1/(root(L*C)). I used a numpy funciton to calculate the square root
        Vf = self.Vs
        dV0 = self.I0 / self.C

        print("alpha =", alpha)
        print("omega0 =", omega0)

        if alpha > omega0:
            print("Response Type = Overdamped")
            s1 = -alpha + np.sqrt(alpha * alpha - omega0 * omega0)
            s2 = -alpha - np.sqrt(alpha * alpha - omega0 * omega0)
            A = (dV0 - s2 * (self.V0 - Vf)) / (s1 - s2)
            B = (self.V0 - Vf) - A
            print("s1 =", s1, " s2 =", s2)
            print("Vc(t) =", Vf, "+ (", A, ")*e^(", s1, "t) + (", B, ")*e^(", s2, "t) V")

        elif alpha == omega0:
            print("Response Type = Critically Damped")
            A = self.V0 - Vf
            B = dV0 + alpha * A
            print("Vc(t) =", Vf, "+ (", A, "+ (", B, ")t )*e^(-", alpha, "t) V")

        else:
            print("Response Type = Underdamped")
            omega_d = np.sqrt(omega0 * omega0 - alpha * alpha) #omegad is angular frequency in underdamped condition
            A = self.V0 - Vf
            B = (dV0 + alpha * A) / omega_d
            print("omega_d =", omega_d)
            print("Vc(t) =", Vf, "+ e^(-", alpha, "t)[(", A, ")cos(", omega_d, "t) + (", B, ")sin(", omega_d, "t)] V")


def main():
    while True:
        print("\n1. RC Circuit")
        print("2. RL Circuit")
        print("3. Series RLC Circuit")
        print("4. Exit")

        choice = int(input("Enter your choice: "))

        if choice == 1:
            R = float(input("Enter R (ohm): "))
            C = float(input("Enter C (F): "))
            Vs = float(input("Enter Vs (V): "))
            V0 = float(input("Enter V0 (V): "))
            circuit = RCCircuit(R, C, Vs, V0)
            circuit.calculate_response()

        elif choice == 2:
            R = float(input("Enter R (ohm): "))
            L = float(input("Enter L (H): "))
            Vs = float(input("Enter Vs (V): "))
            I0 = float(input("Enter I0 (A): "))
            circuit = RLCircuit(R, L, Vs, I0)
            circuit.calculate_response()

        elif choice == 3:
            R = float(input("Enter R (ohm): "))
            L = float(input("Enter L (H): "))
            C = float(input("Enter C (F): "))
            Vs = float(input("Enter Vs (V): "))
            V0 = float(input("Enter V0 (V): "))
            I0 = float(input("Enter I0 (A): "))
            circuit = SeriesRLC(R, L, C, Vs, V0, I0)
            circuit.calculate_response()

        elif choice == 4:
            print("Exiting program.")
            break

        else:
            print("Invalid choice.")


main()