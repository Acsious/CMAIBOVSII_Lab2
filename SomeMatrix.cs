namespace CMAIBOVSII_Lab2
{
    internal class SomeMatrix
    {
        public static double[][] CreateSparseMatrix(int n)
        {
            double[][] A = new double[n][];
            for (int i = 0; i < n; i++)
            {
                A[i] = new double[n];
                for (int j = 0; j < n; j++)
                {
                    A[i][j] = 0;
                }
            }
            for (int i = 1; i < n - 1; i++)
            {
                A[i][i] = 2;
                A[i][i + 1] = -1;
                A[i][i - 1] = -1;
            }
            A[0][0] = 2;
            A[0][1] = -1;
            A[n - 1][n - 1] = 2;
            A[n - 1][n - 2] = -1;
            return A;
        }

        public static double[] ConjugateGradientMethod(double[][] A, double[] b, double tolerance = 1e-6, int maxIterations = 500)
        {
            var n = A.GetLength(0);
            var x = new double[n];
            var r = VectorSubtract(b, MatrixVectorMultiply(A, x));
            var p = r;

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                var alpha = VectorScalarMultiply(r, r) / VectorScalarMultiply(p, MatrixVectorMultiply(A, p));
                x = VectorAdd(x, ScalarAndVectorMultiply(alpha, p));
                var r_new = VectorSubtract(r, ScalarAndVectorMultiply(alpha, MatrixVectorMultiply(A, p)));

                if (Math.Sqrt(VectorScalarMultiply(r_new, r_new)) < tolerance)
                {
                    Console.WriteLine($"Кол-во итераций: {iteration + 1}");
                    return x;
                }

                var beta = VectorScalarMultiply(r_new, r_new) / VectorScalarMultiply(r, r);
                p = VectorAdd(r_new, ScalarAndVectorMultiply(beta, p));
                r = r_new;
            }

            Console.WriteLine("максимум итераций");
            return x;
        }

        public static double[] PreconditionedConjugateGradientMethod(double[][] A, double[] b, double tolerance = 1e-6, int maxIterations = 500)
        {
            int n = A.GetLength(0);
            var x = new double[n];
            var r = VectorSubtract(b, MatrixVectorMultiply(A, x)); // невязка
            var z = SolvePreconditioner(A, r);
            var p = z; //направление

            for (int iteration = 0; iteration < maxIterations; iteration++)
            {
                var alpha = VectorScalarMultiply(r, z) / VectorScalarMultiply(p, MatrixVectorMultiply(A, p));
                x = VectorAdd(x, ScalarAndVectorMultiply(alpha, p));
                var r_new = VectorSubtract(r, ScalarAndVectorMultiply(alpha, MatrixVectorMultiply(A, p)));

                if (Math.Sqrt(VectorScalarMultiply(r_new, r_new)) < tolerance)
                {
                    Console.WriteLine($"Кол-во итераций: {iteration + 1}");
                    return x;
                }

                var z_new = SolvePreconditioner(A, r_new); 
                var beta = VectorScalarMultiply(r_new, z_new) / VectorScalarMultiply(r, z);
                p = VectorAdd(z_new, ScalarAndVectorMultiply(beta, p)); 

                r = r_new;
                z = z_new;
            }

            Console.WriteLine("максимум итераций");
            return x;
        }

        public static double[] MatrixVectorMultiply(double[][] A, double[] x)
        {
            int rows = A.GetLength(0);
            int cols = A[0].GetLength(0);
            double[] result = new double[rows];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    result[i] += A[i][j] * x[j];
                }
            }
            return result;
        }

        public static double VectorScalarMultiply(double[] x, double[] y)
        {
            if (x.Length != y.Length)
            {
                throw new ArgumentException("векторы разной длины");
            }
            double result = 0;
            for (int i = 0; i < x.Length; i++)
            {
                result += x[i] * y[i];
            }
            return result;
        }

        public static double[] VectorSubtract(double[] x, double[] y)
        {
            if (x.Length != y.Length)
            {
                throw new ArgumentException("векторы разной длины");
            }
            double[] result = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                result[i] = x[i] - y[i];
            }
            return result;
        }

        public static double[] VectorAdd(double[] x, double[] y)
        {
            if (x.Length != y.Length)
            {
                throw new ArgumentException("векторы разной длины");
            }
            double[] result = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                result[i] = x[i] + y[i];
            }
            return result;
        }

        public static double[] ScalarAndVectorMultiply(double scalar, double[] x)
        {
            double[] result = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                result[i] = scalar * x[i];
            }
            return result;
        }

        // предобуславливатель якоби
        public static double[] SolvePreconditioner(double[][] A, double[] r)
        {
            int n = A.GetLength(0);
            double[] z = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (A[i][i] != 0)  
                {
                    z[i] = r[i] / A[i][i]; // диагональ
                }
                else
                {
                    z[i] = 0;
                }
            }
            return z;
        }
    }
}
