import cmath
import math
import subprocess
import sys
import unittest

import numpy as np

from main import (Circuit, DCSource, SinusoidalSource, ExponentialSource,
                  RCCircuit, RLCircuit, SeriesRLC, VoltageSource)


class ResponseTests(unittest.TestCase):
    def test_abstract_circuit(self):
        with self.assertRaises(TypeError):
            Circuit()

        class IncompleteCircuit(Circuit):
            pass

        with self.assertRaises(TypeError):
            IncompleteCircuit()

    def test_numpy_arrays_and_source_classes(self):
        times = np.linspace(0, 1, 11)
        sources = [DCSource(4), SinusoidalSource(4, 2, 30), ExponentialSource(4, -2)]
        expected = [np.full_like(times, 4), 4 * np.sin(4 * np.pi * times + np.pi / 6),
                    4 * np.exp(-2 * times)]
        for source, values in zip(sources, expected):
            np.testing.assert_allclose(source.value(times), values, atol=1e-14)
            response = RCCircuit(1, 1, source, 0).response()
            np.testing.assert_allclose(response.value(times), [response.value(t) for t in times])
            np.testing.assert_allclose(response.value(times, 1) + response.value(times), values,
                                       atol=1e-12)

    def test_exceptions_and_cli_recovery(self):
        for source in (lambda: DCSource(float('nan')),
                       lambda: ExponentialSource(1, float('inf'))):
            with self.assertRaises(ValueError):
                source()
        response = RCCircuit(1, 1, 1, 0).response()
        for t in (-1, float('nan')):
            with self.assertRaises(ValueError):
                response.value(t)
        with self.assertRaises(FloatingPointError):
            ExponentialSource(1, 1000).value(1)
        p = subprocess.run([sys.executable, 'main.py'], input='1\noops\n4\n',
                           text=True, capture_output=True, timeout=10, check=True)
        self.assertIn('Invalid input:', p.stdout)
        self.assertIn('Exiting program.', p.stdout)

    def test_initial_conditions_and_differential_equations(self):
        sources = [VoltageSource(10), VoltageSource.sinusoidal(4, 0.7, 35),
                   VoltageSource.exponential(3, -0.8),
                   VoltageSource.exponential(2, 0.4),
                   VoltageSource.sinusoidal(4, 0, 90)]
        count = 0
        for source in sources:
            circuits = [RCCircuit(3, 0.4, source, 2), RLCircuit(3, 0.4, source, 2)]
            circuits += [SeriesRLC(r, 1, 1, source, 2, 0.3) for r in (1, 2, 3)]
            for circuit in circuits:
                with self.subTest(source=source, circuit=type(circuit), R=circuit.R):
                    response = circuit.response()
                    self.assertAlmostEqual(response.value(0), 2)
                    if isinstance(circuit, SeriesRLC):
                        self.assertAlmostEqual(response.value(0, 1), 0.3)
                    for t in (0, 0.1, 0.8, 2):
                        y, dy = response.value(t), response.value(t, 1)
                        if isinstance(circuit, RCCircuit):
                            lhs = circuit.R * circuit.C * dy + y
                        elif isinstance(circuit, RLCircuit):
                            lhs = circuit.L * dy + circuit.R * y
                        else:
                            lhs = circuit.L * circuit.C * response.value(t, 2) + circuit.R * circuit.C * dy + y
                        rhs = (source.coefficient * cmath.exp(source.exponent * t)).real
                        self.assertAlmostEqual(lhs, rhs, places=9)
                    count += 1
        self.assertEqual(count, 25)

    def test_matching_exponential_rates(self):
        # First-order, distinct second-order roots, and a repeated second-order root.
        for circuit in [RCCircuit(1, 1, VoltageSource.exponential(2, -1), 0),
                        RLCircuit(1, 1, VoltageSource.exponential(2, -1), 0),
                        SeriesRLC(3, 1, 0.5, VoltageSource.exponential(2, -1), 0, 0)]:
            response = circuit.response()
            self.assertEqual(response.terms[0][1], 1)
            self.assertTrue(math.isfinite(response.value(1)))
            for t in (0, 0.5, 1):
                y, dy = response.value(t), response.value(t, 1)
                if isinstance(circuit, SeriesRLC):
                    lhs = circuit.L * circuit.C * response.value(t, 2) + circuit.R * circuit.C * dy + y
                else:
                    lhs = dy + y
                self.assertAlmostEqual(lhs, 2 * math.exp(-t))
        for cls in (RCCircuit, RLCircuit):
            response = cls(1, 1, VoltageSource.exponential(2, -1), 0).response()
            self.assertAlmostEqual(response.value(1), 2 / math.e)
        response = SeriesRLC(2, 1, 1, VoltageSource.exponential(2, -1), 0, 0).response()
        self.assertEqual(response.terms[0][1], 2)
        self.assertAlmostEqual(response.value(1), 1 / math.e)
        self.assertAlmostEqual(response.value(0), 0)
        self.assertAlmostEqual(response.value(0, 1), 0)
        for t in (0, 0.5, 1):
            self.assertAlmostEqual(response.value(t, 2) + 2 * response.value(t, 1) + response.value(t), 2 * math.exp(-t))

    def test_dc_regression(self):
        self.assertAlmostEqual(RCCircuit(1000, 0.0001, 10, 0).response().value(0.1), 10 * (1 - math.exp(-1)))
        self.assertAlmostEqual(RLCircuit(10, 2, 20, 0).response().value(0.2), 2 * (1 - math.exp(-1)))

    def test_invalid_components(self):
        for bad in (0, -1, float('inf'), float('nan')):
            with self.assertRaises(ValueError):
                RCCircuit(bad, 1, 10, 0)
        with self.assertRaises(ValueError):
            VoltageSource.sinusoidal(1, -1)

    def test_cli_all_sources_and_circuits(self):
        for circuit, components, initial in [('1', ['1', '1'], ['0']),
                                              ('2', ['1', '1'], ['0']),
                                              ('3', ['2', '1', '1'], ['0', '0'])]:
            for source in (['1', '10'], ['2', '10', '1', '30'], ['3', '10', '-1']):
                data = '\n'.join([circuit] + components + source + initial + ['4', ''])
                p = subprocess.run([sys.executable, 'main.py'], input=data, text=True,
                                   capture_output=True, timeout=10, check=True)
                self.assertIn('(t) = Re[', p.stdout)
                self.assertIn('Exiting program.', p.stdout)
                self.assertNotIn('Invalid input:', p.stdout)
                self.assertNotIn('nan', p.stdout)


if __name__ == '__main__':
    unittest.main()
