using System;

namespace GaussianQF
{
    public class QuadratureFormulas : IQuadratureFormula
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
        private QFMethod method;
        private int parts;
        public delegate double WeightFunction(int k, double a, double b);
        public delegate double[] SystemSolver(double[,] matrix);
        public delegate double[] PolynomialSolver(double[] poly, double a, double b);
        public delegate double Function(double x);
        public int Parts { get { return parts; } set { parts = value; } }
        public double H { get; set; }
        public int N { get; set; }
        public QuadratureFormulas(int n, double a, double b, Function function, WeightFunction weights, QFMethod method, QFType type, int parts = 10)
        {
            this.left = a;
            this.right = b;
            this.weights = weights;
            this.N = n;
            this.function = function;
            this.type = type;
            this.method = method;
            this.parts = parts;
            this.polynomialSolver = (x, left, right) => Poly.RootsFinding(x, left, right, 1000000, 1E-15);
            this.systemSolver = Infrastructure.GaussSolve;
        }
        public void SetPolynomialSolver(PolynomialSolver polynomialSolver)
        {
            this.polynomialSolver = polynomialSolver;
        }
        public void SetSystemSolver(SystemSolver systemSolver)
        {
            this.systemSolver = systemSolver;
        }
        public double CalculateIntegral()
        {
            if (type == QFType.Simple)
            {
                GenerateNodes();
                FindCoefficients();
                return Sum();
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

                    GenerateNodes();
                    FindCoefficients();
                    result += Sum();
                }
                left = start;
                right = end;
                return result;
            }
            else throw new NotImplementedException();

        }
        private double Sum()
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
            nodes = polynomialSolver(poly, left, right);
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
            foreach (var c in coefficients) if (c < 0) throw new ApplicationException("Появился отрицательный коэффициент! Погрешность начинает расти!");
        }
        private void GenerateEquidistantNodes()
        {
            nodes = new double[N];
            var dx = (right - left) / (N - 1);
            var node = left;
            for (int i = 0; i < N; i++)
            {
                nodes[i] = node;
                node += dx;
            }
        }
        private void GenerateNodes()
        {
            if (method == QFMethod.Gauss)
            {
                GenerateAndSolveSystem();
                FindPolynomialRoots();
            }
            else if (method == QFMethod.NewtonCotes)
            {
                GenerateEquidistantNodes();
            }
            else
            {
                throw new NotImplementedException();
            }
        }
    }
}
