using System;

namespace GaussianQF
{
    class Program
    {
        static void Main(string[] args)
        {
            var summ = Infrastructure.CalculateByDefinition(x => F15(x) * Math.Pow(x - 1.1, -4.0 / 5.0), 1.1, 2.3, 1000);
            Console.WriteLine("Интеграл по определению с n={0} равен {1}\n", 10000, summ);

            IQuadratureFormula qf = new QuadratureFormulas(9, 0, 1.2, F15InNewVars, WeightsVariant15, QFMethod.Gauss, QFType.Simple);
            Console.WriteLine("Интеграл по КФ Гаусса = {0}\n", qf.CalculateIntegral());

            IQuadratureFormula qf2 = new QuadratureFormulas(3, 0, 1.2, F15InNewVars, WeightsVariant15, QFMethod.NewtonCotes, QFType.Simple);
            Console.WriteLine("Интеграл по КФ Ньютона Котеса = {0}\n", qf2.CalculateIntegral());

            Console.WriteLine("\nПроцесс Ричардсона для КФ Ньютона-Котеса");
            IQuadratureFormula qf4 = new QuadratureFormulas(3, 0, 1.2, F15InNewVars, WeightsVariant15, QFMethod.NewtonCotes, QFType.Complex, 2);
            Richardson.Calculate(qf4, 2, 2, 1E-6, 6);

            Console.WriteLine("\nПроцесс Ричардсона для КФ Гаусса");

            IQuadratureFormula qf3 = new QuadratureFormulas(3, 0, 1.2, F15InNewVars, WeightsVariant15, QFMethod.Gauss, QFType.Complex, 2);
            Richardson.Calculate(qf3, 2, 2, 1E-6, 10);
        }
        #region Вариант 1
        public static double WeightsVariant1(int k, double a, double b)
        {
            Func<double, double> weight = t => (3.0 / (3.0 * k + 1.0)) * Math.Pow(t, ((3.0 * k + 1.0) / 3.0));
            return weight(b) - weight(a);
        }
        public static double F(double x)
        {
            return 2.0 * Math.Cos(2.5 * x) * Math.Exp(x / 3.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + x;
        }
        public static double FInNewVars(double t)
        {
            return F(t + 1.5);
        }
        #endregion
        #region Вариант 3
        public static double WeightsVariant3(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 4.0 / 5.0)) * Math.Pow(t, k + 4.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double F3(double x)
        {
            return 2.5 * Math.Cos(2 * x) * Math.Exp(2.0 * x / 3.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 3 * x;
        }
        public static double F3InNewVars(double t)
        {
            return F3(t + 0.1);
        }
        #endregion
        #region Вариант 5
        public static double WeightsVariant5(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 5.0 / 7.0)) * Math.Pow(t, k + 5.0 / 7.0);
            return weight(b) - weight(a);
        }
        public static double F5(double x)
        {
            return Math.Cos(1.5 * x) * Math.Exp(2 * x / 3) + 3 * Math.Sin(5.5 * x) * Math.Exp(-2 * x) + 2;
        }
        public static double F5InNewVars(double t)
        {
            return F5(t + 2.5);
        }
        #endregion
        #region Вариант 8
        public static double WeightsVariant8(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 2.0 / 5.0)) * Math.Pow(t, k + 2.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double F8(double x)
        {
            return 3.7 * Math.Cos(1.5 * x) * Math.Exp(-4.0 * x / 3.0) + 2.4 * Math.Sin(4.5 * x) * Math.Exp(2.0 * x / 3.0) + 4;
        }
        public static double F8InNewVars(double t)
        {
            return F8(-t + 2.3);
        }
        #endregion
        #region Вариант 9
        public static double WeightsVariant9(int k, double a, double b)
        {
            Func<double, double> weight = t => (3.0 / (3.0 * k + 1.0)) * Math.Pow(t, ((3.0 * k + 1.0) / 3.0));

            return weight(b) - weight(a);
        }
        public static double F9(double x)
        {
            return 3.0 * Math.Cos(1.5 * x) * Math.Exp(x / 4.0) + 4.0 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 4 * x;
        }
        public static double F9InNewVars(double t)
        {
            return F9(t + 2.5);
        }
        #endregion
        #region Вариант 13
        public static double WeightsVariant13(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 4.0 / 5.0)) * Math.Pow(t, k + 4.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double F13(double x)
        {
            return 2.0 * Math.Cos(3.5 * x) * Math.Exp(5.0 * x / 3.0) + 3.0 * Math.Sin(1.5 * x) * Math.Exp(-4 * x) + 3;
        }
        public static double F13InNewVars(double t)
        {
            return F13(t + 1.5);
        }
        #endregion
        #region Вариант 15
        public static double F15(double x)
        {
            return 3.5 * Math.Cos(0.7 * x) * Math.Exp(-5.0 * x / 3.0) + 2.4 * Math.Sin(5.5 * x) * Math.Exp(-3.0 * x / 4.0) + 5;
        }
        public static double F15InNewVars(double t)
        {
            return F15(t + 1.1);
        }
        public static double WeightsVariant15(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 1.0 / 5.0)) * Math.Pow(t, k + 1.0 / 5.0);
            return weight(b) - weight(a);
        }
        #endregion
        #region Вариант 19
        public static double WeightsVariant19(int k, double a, double b)
        {
            Func<double, double> weight = t => (1.0 / (k + 3.0 / 5.0)) * Math.Pow(t, k + 3.0 / 5.0);
            return weight(b) - weight(a);
        }
        public static double F19(double x)
        {
            return 0.5 * Math.Cos(3 * x) * Math.Exp(2 * x / 5) + 4 * Math.Sin(3.5 * x) * Math.Exp(-3 * x) + 3 * x;
        }
        public static double F19InNewVars(double t)
        {
            return F19(t + 1.1);
        }
        #endregion
    }


}
