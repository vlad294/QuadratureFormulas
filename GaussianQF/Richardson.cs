using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GaussianQF
{
    class Richardson
    {
        public static double Aitken(IQuadratureFormula qf, int q, double tolerance, out int parts)
        {
            var initialPartsCount = qf.Parts;
            Console.WriteLine("Считаем интеграл с наперед заданной точностью {0}", tolerance);
            Console.Write("Определяем параметры m и количество интервалов разбиения. \nСтроим последовательность mi:");
            double mi = 0, miOld = 0;
            while (Math.Abs(mi - miOld) > 0.1 || miOld == 0)
            {
                var s0 = qf.CalculateIntegral();
                qf.Parts *= q;
                var s1 = qf.CalculateIntegral();
                qf.Parts *= q;
                var s2 = qf.CalculateIntegral();
                miOld = mi;
                initialPartsCount *= q;
                qf.Parts = initialPartsCount;
                mi = Math.Log((Math.Abs(s0 - s1) / Math.Abs(s1 - s2))) / Math.Log(q);
                Console.Write("{0} ", mi);
            }
            parts = qf.Parts;
            Console.WriteLine("\n\nСтепень при главном члене погрешности m={0}, количество интервалов разбиения={1}\n", mi, qf.Parts);
            return mi;
        }
        public static double Calculate(IQuadratureFormula qf, int L, int q, double tolerance, int maxr)
        {
            int parts;
            double m = Aitken(qf, q, tolerance, out parts);
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
                    qf.Parts = parts * (int)Math.Pow(L, i);
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
                    return integrals[r + 1];
                }
            }
            throw new ApplicationException("Try increase maxr parameter");
        }
    }
}
