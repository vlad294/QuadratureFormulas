using System;

namespace GaussianQF
{
    /// <summary>
    /// Вычисление интегралов "с особенностью"
    /// Квадратурные формулы Ньютона-Котеса и Гаусса
    /// </summary>
    public class QuadratureFormulas : IQuadratureFormula
    {
        /// <summary>
        /// Левая граница интегрирования
        /// </summary>
        private double left;

        /// <summary>
        /// Правая граница интегрирования
        /// </summary>
        private double right;

        /// <summary>
        /// Функция (без весовой функции)
        /// </summary>
        private Function function;

        /// <summary>
        /// Моменты весовой функции
        /// </summary>
        private WeightFunction weights;

        /// <summary>
        /// Метод решения СЛАУ
        /// </summary>
        private SystemSolver systemSolver;

        /// <summary>
        /// Метод нахождения корней полинома
        /// </summary>
        private PolynomialSolver polynomialSolver;

        /// <summary>
        /// Узловой полином КФ Гаусса
        /// </summary>
        private double[] poly;

        /// <summary>
        /// Узлы квадратурной формулы
        /// </summary>
        private double[] nodes;

        /// <summary>
        /// Коэффициенты квадратурной формулы
        /// </summary>
        private double[] coefficients;

        /// <summary>
        /// Тип: простая или составная
        /// </summary>
        private QFType type;

        /// <summary>
        /// Способ выбора узлов
        /// Равноотстоящие (Ньютона-Котеса)
        /// Гаусс
        /// </summary>
        private QFMethod method;

        /// <summary>
        /// Количество участков интегрирования
        /// </summary>
        private int parts;

        /// <summary>
        /// Количество участков интегрирования
        /// </summary>
        public int Parts { get { return parts; } set { parts = value; } }

        /// <summary>
        /// Шаг сетки составной квадратурной формулы
        /// </summary>
        public double H { get; set; }

        /// <summary>
        /// Количество точек квадратурной формулы
        /// </summary>
        public int N { get; set; }

        /// <summary>
        /// Создание квадратурной формулы
        /// </summary>
        /// <param name="n">Количество точек</param>
        /// <param name="a">Левая граница</param>
        /// <param name="b">Правая граница</param>
        /// <param name="function">Функция</param>
        /// <param name="weights">Функция, возвращающая k-ый момент весовой функции</param>
        /// <param name="method">Ньютон-Котес (равноотстоящие узлы) или Гаусс</param>
        /// <param name="type">Простая либо составная квадратурная формула</param>
        /// <param name="parts">Количество участков разбиения составной квадратурной формулы</param>
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
            SetPolynomialSolver((x, left, right) => Poly.RootsFinding(x, left, right, 1000, 1E-15));
            SetSystemSolver(Infrastructure.GaussSolve);
        }

        /// <summary>
        /// Изменить метод поиска корней полинома
        /// </summary>
        /// <param name="polynomialSolver">Метод поиска корней полинома на промежутке</param>
        public void SetPolynomialSolver(PolynomialSolver polynomialSolver)
        {
            this.polynomialSolver = polynomialSolver;
        }
        
        /// <summary>
        /// Изменить метод решения СЛАУ
        /// </summary>
        /// <param name="systemSolver">Метод решения СЛАУ</param>
        public void SetSystemSolver(SystemSolver systemSolver)
        {
            this.systemSolver = systemSolver;
        }
        
        /// <summary>
        /// Вычисление значения интеграла
        /// </summary>
        /// <returns></returns>
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
        
        /// <summary>
        /// Сумма коэффициентов квадратурной формулы на значение функции в узлах
        /// </summary>
        /// <returns></returns>
        private double Sum()
        {
            double summ = 0;
            for (int i = 0; i < N; i++)
            {
                summ += coefficients[i] * function(nodes[i]);
            }
            return summ;
        }
        
        /// <summary>
        /// Составление системы для узлового полинома квадратурной формулы Гаусса
        /// </summary>
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

        /// <summary>
        /// Вычисление корней полинома
        /// </summary>
        private void FindPolynomialRoots()
        {
            nodes = polynomialSolver(poly, left, right);
        }
        
        /// <summary>
        /// Вычисление коэффициентов квадратурной формулы
        /// </summary>
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
        
        /// <summary>
        /// Генерация равноотстоящих узлов (квадратурная формула Ньютона-Котеса)
        /// </summary>
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
        
        /// <summary>
        /// Вычисление узлов квадратурной формулы
        /// </summary>
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
