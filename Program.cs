using CMAIBOVSII_Lab2;

var n = 50;
var A = SomeMatrix.CreateSparseMatrix(n);

var b = new double[n];
b[0] = 1;
b[n - 1] = -1;

Console.WriteLine("метод сопряженных градиентов");
var solution = SomeMatrix.ConjugateGradientMethod(A, b);
foreach (var item in solution)
{
    Console.WriteLine(item);
}

Console.WriteLine();
Console.WriteLine("предобусловленный метод сопряженных градиентов");
solution = SomeMatrix.PreconditionedConjugateGradientMethod(A, b);
foreach (var item in solution)
{
    Console.WriteLine(item);
}