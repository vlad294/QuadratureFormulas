using System;

namespace GaussianQF
{
    static class Infrastructure
    {
        public delegate double Function(double x);
        public static void Richardson(IQuadratureFormula qf, int L, int q, double tolerance, int maxr)
        {
            var initialPartsCount = qf.Parts;
            Console.WriteLine("Считаем интеграл с наперед заданной точностью {0}", tolerance);
            Console.Write("Определяем параметры m и количество интервалов разбиения. Строим последовательность mi:");
            double mi = 0;
            double miOld = 0;
            while (Math.Abs(mi - miOld) > 0.1 || miOld == 0)
            {
                var s0 = qf.CalculateIntegral(); //
                qf.Parts *= L;
                var s1 = qf.CalculateIntegral();
                qf.Parts *= L;
                var s2 = qf.CalculateIntegral();

                miOld = mi;
                qf.Parts = initialPartsCount * L;

                mi = Math.Log((Math.Abs(s0 - s1) / Math.Abs(s1 - s2))) / Math.Log(q);
                Console.Write("{0} ", mi);
            }
            int parts = qf.Parts;
            int m = (int)Math.Ceiling(mi);
            Console.WriteLine("\n\nПараметр m={0}, количество интервалов разбиения={1}\n", m, parts);

            var integrals = new double[maxr + 2];
            var H = new double[maxr + 2];
            var prevSystemSize = 0;
            for (int r = 0; r < maxr; r++)
            {
                Console.WriteLine("Составляем систему для r={0}", r);
                //нужно посчитать 2+r интегралов
                for (int i = 0; i < r + 2; i++)
                {
                    if (i <= prevSystemSize && (i + r) != 0) continue; //уже посчитанные интегралы не считаем
                    qf.Parts = parts * (int)Math.Pow(q, i);
                    integrals[i] = qf.CalculateIntegral();
                    H[i] = qf.H;
                    Console.WriteLine("{2}-ый интеграл по {0} частям равен {1} с шагом {3}", qf.Parts, integrals[i], i + 1, H[i]);
                    prevSystemSize = i;
                }
                var matrix = new double[r + 2, r + 3];
                for (int i = 0; i < 2 + r; i++)
                {
                    for (int j = 0; j <= r; j++)
                    {
                        matrix[i, j] = Math.Pow(H[i], j + m);
                    }
                    matrix[i, r + 1] = -1;
                    matrix[i, r + 2] = -integrals[i];
                }
                var answer = Infrastructure.GaussSolve(matrix);
                var epsilon = Math.Abs(answer[r + 1] - integrals[r + 1]);
                Console.WriteLine("R={0}, Sj={1}, epsilon={2}\n", answer[r + 1], integrals[r + 1], epsilon);
                if (epsilon <= tolerance)
                {
                    Console.WriteLine("Заданная точность достигнута");
                    Console.WriteLine("Значение интеграла: {0} при разбиении на {1} частей", integrals[r + 1], qf.Parts);
                    return;
                }
            }
        }
        public static double CalculateByDefinition(Function f, double start, double end, int n)
        {
            var H = (end - start) / (n - 1);
            double summ = 0;
            for (int i = 0; i < n; i++)
            {
                var point = start + i * H + H / 2;
                summ += f(point) * H;
            }
            return summ;

        }
        public static double[] GaussSolve(double[,] M)
        {
            // размерность прямоугольной матрицы M
            int n = M.GetUpperBound(0) + 1;

            for (int col = 0; col + 1 < n; col++) if (M[col, col] == 0)
                // проверка на нулевые коэффициенты
                {
                    // поиск ненулевого коэффициента
                    int swapRow = col + 1;
                    for (; swapRow < n; swapRow++) if (M[swapRow, col] != 0) break;

                    if (M[swapRow, col] != 0) // найден ненулевой коэффициент
                    {
                        // смена значений местами
                        double[] tmp = new double[n + 1];
                        for (int i = 0; i < n + 1; i++)
                        {
                            tmp[i] = M[swapRow, i];
                            M[swapRow, i] = M[col, i];
                            M[col, i] = tmp[i];
                        }
                    }
                    else return null;
                }
            // исключение
            for (int sourceRow = 0; sourceRow + 1 < n; sourceRow++)
            {
                for (int destRow = sourceRow + 1; destRow < n; destRow++)
                {
                    double df = M[sourceRow, sourceRow];
                    double sf = M[destRow, sourceRow];
                    for (int i = 0; i < n + 1; i++)
                        M[destRow, i] = M[destRow, i] * df - M[sourceRow, i] * sf;
                }
            }
            // обратный ход гаусса
            for (int row = n - 1; row >= 0; row--)
            {
                double f = M[row, row];
                if (f == 0) return null;

                for (int i = 0; i < n + 1; i++) M[row, i] /= f;
                for (int destRow = 0; destRow < row; destRow++)
                { M[destRow, n] -= M[destRow, row] * M[row, n]; M[destRow, row] = 0; }
            }
            // формирование ответа
            var x = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = M[i, n];
            }
            return x;
        }
        public static double[] SeidelSolve(double[,] matrix, double[] right,
                              double relaxation, int iterations,
                              double[] lo, double[] hi)
        {
            // Validation omitted
            var x = right;
            double delta;

            // Gauss-Seidel with Successive OverRelaxation Solver
            for (int k = 0; k < iterations; ++k)
            {
                for (int i = 0; i < right.Length; ++i)
                {
                    delta = 0.0f;

                    for (int j = 0; j < i; ++j)
                        delta += matrix[i, j] * x[j];
                    for (int j = i + 1; j < right.Length; ++j)
                        delta += matrix[i, j] * x[j];

                    delta = (right[i] - delta) / matrix[i, i];
                    x[i] += relaxation * (delta - x[i]);
                    // Project the solution within the lower and higher limits
                    if (x[i] < lo[i])
                        x[i] = lo[i];
                    if (x[i] > hi[i])
                        x[i] = hi[i];
                }
            }
            return x;
        }
    }
    public enum QFType
    {
        Simple,
        Complex
    }
}
