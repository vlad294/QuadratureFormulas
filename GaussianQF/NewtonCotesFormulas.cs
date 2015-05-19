using System;

namespace GaussianQF
{
    public class NewtonCotesFormulas : IQuadratureFormula
    {
        private double left;
        private double right;
        private Function function;
        private WeightFunction weights;
        private SystemSolver systemSolver;
        private double[] nodes;
        private double[] coefficients;
        private int parts;
        private QFType type;

        public delegate double WeightFunction(int k, double a, double b);
        public delegate double[] SystemSolver(double[,] matrix);
        public delegate double Function(double x);
        public int Parts { get { return parts; } set { parts = value; } }
        public double H { get; set; }
        public int N { get; set; }
        public NewtonCotesFormulas(int n, double a, double b, QFType type, Function function, WeightFunction weights, SystemSolver systemSolver, int parts = 10)
        {
            this.left = a;
            this.right = b;
            this.weights = weights;
            this.N = n;
            this.systemSolver = systemSolver;
            this.function = function;
            this.parts = parts;
            this.type = type;
        }
        public void GenerateEquidistantNodes()
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
        private void GenerateAndSolveSystem()
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
        private double Summ()
        {
            double summ = 0;
            for (int i = 0; i < N; i++)
            {
                var c = function(nodes[i]);
                summ += coefficients[i] * function(nodes[i]);
            }
            return summ;
        }
        public double CalculateIntegral()
        {
            if (type == QFType.Simple)
            {
                GenerateEquidistantNodes();
                GenerateAndSolveSystem();
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

                    GenerateEquidistantNodes();
                    GenerateAndSolveSystem();
                    result += Summ();
                   // Console.WriteLine("от {0} до {1} = {2}", left, right, Summ());
                }

                left = start;
                right = end;
                return result;
            }
            else throw new NotImplementedException();
        }
    }

}
