using System;

namespace GaussianQF
{
    class Program
    {

        static void Main(string[] args)
        {
            //var summ = Infrastructure.CalculateByDefinition(x => F19(x) * Math.Pow(x - 1.5, -1.0 / 5.0), 1.5, 2.3, 10000);
            //Console.WriteLine("Интеграл по определению с n={0} равен {1}\n", 10000, summ);

            ////SimpleR();

            var qf = new GaussianQuadratureFormulas(3, 0, 0.5, QFType.Simple, F8InNewVars, WeightsVariant8, Infrastructure.GaussSolve, (x, left, right) => Poly.RootsFinding(x, left, right, 1000000, 1E-15));
            Console.WriteLine("Интеграл по КФ Гаусса = {0}\n", qf.CalculateIntegral());

            var qf2 = new NewtonCotesFormulas(3, 0, 0.5, QFType.Simple, F8InNewVars, WeightsVariant8, Infrastructure.GaussSolve);
            Console.WriteLine("Интеграл по КФ Ньютона Котеса = {0}\n", qf2.CalculateIntegral());


            Console.WriteLine("\nПроцесс Ричардсона для КФ Ньютона-Котеса");
            IQuadratureFormula qf4 = new NewtonCotesFormulas(3, 0, 0.5, QFType.Complex, F8InNewVars, WeightsVariant8, Infrastructure.GaussSolve, 2);
            Infrastructure.Richardson(qf4, 2, 2, 1E-6, 6);

            Console.WriteLine("\nПроцесс Ричардсона для КФ Гаусса");

            IQuadratureFormula qf3 = new GaussianQuadratureFormulas(3, 0, 0.5, QFType.Complex, F8InNewVars, WeightsVariant8, Infrastructure.GaussSolve, (x, left, right) => Poly.RootsFinding(x, left, right, 1000000, 1E-15), 2);
            Infrastructure.Richardson(qf3, 2, 2, 1E-10, 10);
            //Console.ReadLine();
        }


        public static double WeightsVariant9(int k, double a, double b)
        {
            Func<double, double> weight = t => (3.0 / (3.0 * k + 1.0)) * Math.Pow(t, ((3.0 * k + 1.0) / 3.0));

            return weight(b) - weight(a);
        }
        public static double WeightsVariant1(int k, double a, double b)
        {
            Func<double, double> weight = t => (3.0 / (3.0 * k + 1.0)) * Math.Pow(t, ((3.0 * k + 1.0) / 3.0));
            return weight(b) - weight(a);
        }
        public static double WeightsVariant3(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 4.0 / 5.0)) * Math.Pow(t, k + 4.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double WeightsVariant13(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 4.0 / 5.0)) * Math.Pow(t, k + 4.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double WeightsVariant19(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 3.0 / 5.0)) * Math.Pow(t, k + 3.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double WeightsVariant8(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 2.0 / 5.0)) * Math.Pow(t, k + 2.0 / 5.0);
            return weight(b) - weight(a);
        }

        public static double F(double x)
        {
            return 3.0 * Math.Cos(1.5 * x) * Math.Exp(x / 4.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 4 * x;
        }
        public static double F2(double x)
        {
            return 2.0 * Math.Cos(2.5 * x) * Math.Exp(x / 3.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + x;
        }
        public static double F3(double x)
        {
            return 2.5 * Math.Cos(2 * x) * Math.Exp(2.0 * x / 3.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 3 * x;
        }
        public static double F13(double x)
        {
            return 2.0 * Math.Cos(3.5 * x) * Math.Exp(5.0 * x / 3.0) + 3.0 * Math.Sin(1.5 * x) * Math.Exp(-4 * x) + 3;
        }
        public static double F19(double x)
        {
            return 0.5 * Math.Cos(3 * x) * Math.Exp(2 * x / 5) + 4 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 3 * x;
        }
        public static double F5(double x)
        {
            return Math.Cos(1.5 * x) * Math.Exp(2 * x / 3) + 3 * Math.Sin(5.5 * x) * Math.Exp(-2 * x) + 3 * x;
        }
        public static double F8(double x)
        {
            return 3.7 * Math.Cos(1.5 * x) * Math.Exp(-4.0 * x / 3.0) + 2.4 * Math.Sin(4.5 * x) * Math.Exp(2.0 * x / 3.0) + 4;
        }
        public static double F8InNewVars(double t)
        {
            return F8(-t + 2.3);
        }
        public static double F13InNewVars(double t)
        {
            return F13(t + 1.5);
        }
        public static double F3InNewVars(double t)
        {
            return F3(t + 0.1);
        }
        public static double F2InNewVars(double t)
        {
            return F2(t + 1.5);
        }
        public static double FInNewVars(double t)
        {
            return F(t + 2.5);
        }
        public static double F19InNewVars(double t)
        {
            return F19(t + 1.1);
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
