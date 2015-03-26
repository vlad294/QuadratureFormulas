using System;

namespace GaussianQF
{
    class Program
    {

        static void Main(string[] args)
        {
            var summ = Infrastructure.CalculateByDefinition(x => F(x) * Math.Pow(x - 2.5, -2.0 / 3.0), 2.5, 3.3, 1000);
            Console.WriteLine("Интеграл по определению с n={0} равен {1}\n", 1000, summ);

            SimpleR();

            var qf = new GaussianQuadratureFormulas(9, 0, 0.8, QFType.Simple, FInNewVars, WeightsVariant9, Infrastructure.GaussSolve, x => Poly.RootsFinding(x, 0, 0.8, 10000, 1E-15));
            Console.WriteLine("Интеграл по КФ Гаусса = {0}\n", qf.CalculateIntegral());

            var qf2 = new NewtonCotesFormulas(3, 0, 0.8, QFType.Simple, FInNewVars, WeightsVariant9, Infrastructure.GaussSolve);
            Console.WriteLine("Интеграл по КФ Ньютона Котеса = {0}\n", qf2.CalculateIntegral());

            Console.WriteLine("\nПроцесс Ричардсона для КФ Ньютона-Котеса");
            IQuadratureFormula qf4 = new NewtonCotesFormulas(3, 0, 0.8, QFType.Complex, FInNewVars, WeightsVariant9, Infrastructure.GaussSolve, 1);
            Infrastructure.Richardson(qf4, 2, 2, 1E-10, 6);

            Console.WriteLine("\nПроцесс Ричардсона для КФ Гаусса");

            IQuadratureFormula qf3 = new GaussianQuadratureFormulas(3, 0, 0.8, QFType.Complex, FInNewVars, WeightsVariant9, Infrastructure.GaussSolve, x => Poly.RootsFinding(x, 0, 0.8, 10000, 1E-15), 1);
            Infrastructure.Richardson(qf3, 2, 2, 1E-14, 10);
        }

        public static double WeightsVariant9(int k, double a, double b)
        {
            Func<double, double> weight = t => (3.0 / (3.0 * k + 1.0)) * Math.Pow(t, ((3.0 * k + 1.0) / 3.0));
            return weight(b) - weight(a);
        }

        public static double F(double x)
        {
            return 3.0 * Math.Cos(1.5 * x) * Math.Exp(x / 4.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 4 * x;
        }
        public static double FInNewVars(double t)
        {
            return F(t + 2.5);
        }
        public static void SimpleR()
        {
            // d^3(3.0 * Cos(1.5 * x) * e^(x / 4.0) + 4.0 * Sin(3.5 * x) * e^(-3 * x) + 4 * x)/dx^3 from 2.5 to 3.3 
            // |M3| = 30

            // d^6(3.0 * Cos(1.5 * x) * e^(x / 4.0) + 4.0 * Sin(3.5 * x) * e^(-3 * x) + 4 * x)/dx^6 from 2.5 to 3.3
            // |M6| = 70
            double m3 = 30; // делить на 3! - факториал
            //integral of  |x(x-0.4)(x-0.8)/((x)^(2/3))| dx from 0 to 0.8
            double int1 = 0.033;
            double r1 = m3 / 6 * int1;
            Console.WriteLine("Оценка для одиночной КФ Ньютона-Котеса: {0}\n", r1);
            //w(x) = x(x-0.4)(x-0.8)
            double m6 = 70; //делить на 6!
            // integral of  |(x(x-0.4)(x-0.8))^2/((x)^(2/3))| dx from 0 to 0.8
            double int2 = 0.0006;
            double r2 = m6 / 144 * int2;
            Console.WriteLine("Оценка для одиночной КФ Гаусса: {0}\n", r2);
        }
    }


}
