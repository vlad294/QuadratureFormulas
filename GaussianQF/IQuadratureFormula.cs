
namespace GaussianQF
{
    interface IQuadratureFormula
    {
        double CalculateIntegral();
        int Parts { get; set; }
        double H { get; set; }
        int N { get; set; }
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
