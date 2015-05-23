using System;
using System.Linq;

namespace GaussianQF
{
    public static class Poly
    {
        public static double Evaluate(double[] polynomial, double x)
        {
            double result = 0;
            for (int i = 0; i < polynomial.Length; i++)
            {
                result += Math.Pow(x, i) * polynomial[i];
            }
            return result;
        }
        public static double[] Derivative(double[] polynomial)
        {
            var derivative = new double[polynomial.Length - 1];
            for(int i=1; i<polynomial.Length ;i++)
            {
                derivative[i - 1] = polynomial[i] * i;
            }
            return derivative;
        }
        public static double[] RootsFinding(double[] polynomial, double a, double b, int maxIteration = 1000, double tolerance = 1E-15)
        {
            //if (polynomial.Length == 4)
            //{
            //    var f = polynomial.Reverse().ToArray();
            //    double Q = (f[1] * f[1] - 3 * f[2]) / 9;
            //    double R = (2 * f[1] * f[1] * f[1] - 9 * f[1] * f[2] + 27 * f[3]) / 54;
            //    double A = Math.Acos(R / Math.Pow(Q * Q * Q, 0.5)) / 3;

            //    return new double[] {
            //    -2 * Math.Pow(Q,0.5) * Math.Cos(A) - f[1]/3,
            //    -2 * Math.Pow(Q,0.5) * Math.Cos(A+2.0/3*Math.PI) - f[1] / 3,
            //    -2 * Math.Pow(Q,0.5) * Math.Cos(A-2.0/3*Math.PI) - f[1] / 3};
            //}
            var roots = new double[polynomial.Length - 1];

            //localize roots
            int n = 50;
            while (n < maxIteration)
            {
                var h = (b - a) / (n - 1);
                int nextRoot = 0;
                for (int i = 0; i < n; i++)
                {
                    var left = a + i * h;
                    var right = a + (i + 1) * h;

                    if (Evaluate(polynomial, left) * Evaluate(polynomial, right) < 0)
                    {
                        roots[nextRoot] = left;
                        nextRoot++;
                        if (nextRoot == (polynomial.Length - 1)) return FindRoot(polynomial, roots, h, tolerance, maxIteration);
                    }
                }
                n *= 20;
            }
            throw new ApplicationException("Не удалось локализовать корни");
        }
        public static double[] FindRoot(double[] polynomial, double[] rootsLeft, double h, double tolerance, int maxIteration)
        {
            for (int i = 0; i < rootsLeft.Length; i++)
            {
                var root = NewtonRaphsonMethod(polynomial, rootsLeft[i], rootsLeft[i] + h, tolerance, maxIteration);
                if (root < rootsLeft[i] || root > rootsLeft[i] + h || double.IsNaN(root))
                    root = BisectionMethod(polynomial, rootsLeft[i], rootsLeft[i] + h, tolerance);
                if (root < rootsLeft[i] || root > rootsLeft[i] + h)
                    throw new ApplicationException("Не удалось уточнить корень");
                rootsLeft[i] = root;
            }

            return rootsLeft;
        }
        public static double BisectionMethod(double[] polynomial, double a, double b, double epsilon)
        {
            double x1 = a;
            double x2 = b;
            double fb = Poly.Evaluate(polynomial, b);
            while (Math.Abs(x2 - x1) > epsilon)
            {
                double midpt = 0.5 * (x1 + x2);
                if (fb * Poly.Evaluate(polynomial, midpt) > 0)
                    x2 = midpt;
                else
                    x1 = midpt;
            }
            return x2 - (x2 - x1) * Poly.Evaluate(polynomial, x2) / (Poly.Evaluate(polynomial, x2) - Poly.Evaluate(polynomial, x1));
        }
        public static double NewtonRaphsonMethod(double[] polynomial, double a, double b, double epsilon, int maxIteration)
        {
            int i = 0;
            var derivative = Derivative(polynomial);
            double f0 = Poly.Evaluate(polynomial, a);
            double x = a;
            while (Math.Abs(Poly.Evaluate(polynomial, x)) > epsilon)
            {
                x -= f0 / Poly.Evaluate(derivative, x);
                f0 = Poly.Evaluate(polynomial, x);
                i++;
                if (i >= maxIteration) return double.NaN;
            }
            return x;
        }

    }
}
