
namespace GaussianQF
{
    public delegate double WeightFunction(int k, double a, double b);
    public delegate double[] SystemSolver(double[,] matrix);
    public delegate double[] PolynomialSolver(double[] poly, double a, double b);
    public delegate double Function(double x);
    interface IQuadratureFormula
    {
        double CalculateIntegral();
        int Parts { get; set; }
        double H { get; set; }
        int N { get; set; }
        void SetPolynomialSolver(PolynomialSolver polynomialSolver);
        void SetSystemSolver(SystemSolver systemSolver);
    }
    public enum QFType
    {
        Simple,
        Complex
    }
    public enum QFMethod
    {
        NewtonCotes,
        Gauss
    }

}
