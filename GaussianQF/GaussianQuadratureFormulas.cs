using System;

namespace GaussianQF
{
    public class GaussianQuadratureFormulas : IQuadratureFormula
    {
        private double left;
        private double right;
        private Function function;
        private WeightFunction weights;
        private SystemSolver systemSolver;
        private PolynomialSolver polynomialSolver;
        private double[] poly;
        private double[] nodes;
        private double[] coefficients;
        private QFType type;
        private int parts;
        public delegate double WeightFunction(int k, double a, double b);
        public delegate double[] SystemSolver(double[,] matrix);
        public delegate double[] PolynomialSolver(double[] poly);
        public delegate double Function(double x);
        public int Parts { get { return parts; } set {parts = value; } }
        public double H { get; set; }
        public int N { get; set; }
        public GaussianQuadratureFormulas(int n, double a, double b, QFType type, Function function, WeightFunction weights, SystemSolver systemSolver, PolynomialSolver polynomialSolver, int parts = 10)
        {
            this.left = a;
            this.right = b;
            this.weights = weights;
            this.N = n;
            this.systemSolver = systemSolver;
            this.polynomialSolver = polynomialSolver;
            this.function = function;
            this.type = type;
            this.parts = parts;
        }
        public double CalculateIntegral()
        {
            if (type == QFType.Simple)
            {
                GenerateAndSolveSystem();
                FindPolynomialRoots();
                FindCoefficients();
                return Summ();
            }
            else if (type == QFType.Complex)
            {
                var start = left;
                var end = right;
                double result = 0;

                H = (right - left) / parts;
                for (int i = 0; i < parts; i++)
                {
                    left = start + i * H;
                    right = start + (i + 1) * H;

                    GenerateAndSolveSystem();
                    FindPolynomialRoots();
                    FindCoefficients();
                    result += Summ();
                    //Console.WriteLine("From {0} to {1} = {2}", left, right, Summ());
                }

                left = start;
                right = end;
                return result;
            }
            else throw new NotImplementedException();

        }
        private double Summ()
        {
            double summ = 0;
            for (int i = 0; i < N; i++)
            {
                summ += coefficients[i] * function(nodes[i]);
            }
            return summ;
        }
        private void GenerateAndSolveSystem()
        {
            var matrix = new double[N, N + 1];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    matrix[i, j] = weights(i + j, left, right);
                }
                matrix[i, N] = -weights(i + N, left, right);
            }
            poly = systemSolver(matrix);
            Array.Resize<double>(ref poly, N + 1);
            poly[N] = 1;
        }
        private void FindPolynomialRoots()
        {
            nodes = polynomialSolver(poly);
        }
        private void FindCoefficients()
        {
            var matrix = new double[N, N + 1];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    matrix[i, j] = Math.Pow(nodes[j], i);
                }
                matrix[i, N] = weights(i, left, right);
            }
            coefficients = systemSolver(matrix);
        }
    }
}
