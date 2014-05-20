/*
 * Copyright (C) 2014. Daniel Asarnow
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package mdscale;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: da
 * Date: 10/8/11
 * Time: 2:27 AM
 */
public class Data {

    public static void doubleCenter(double[][] matrix) {
        int n = matrix[0].length;
        int k = matrix.length;

        for (int j = 0; j < k; j++) {
            double avg = 0.0D;
            for (int i = 0; i < n; i++) avg += matrix[j][i];
            avg /= n;
            for (int i = 0; i < n; i++) matrix[j][i] -= avg;
        }

        for (int i = 0; i < n; i++) {
            double avg = 0.0D;
            for (int j = 0; j < k; j++) avg += matrix[j][i];
            avg /= matrix.length;
            for (int j = 0; j < k; j++) matrix[j][i] -= avg;
        }
    }

    public static void multiply(double[][] matrix, double factor) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] *= factor;
            }
        }
    }

    public static void squareEntries(double[][] matrix) {
        int n = matrix[0].length;
        int k = matrix.length;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = Math.pow(matrix[i][j], 2.0D);
            }
        }
    }

    public static void normalize(double[][] x) {
        for (int i = 0; i < x.length; i++) normalize(x[i]);
    }

    public static double normalize(double[] x) {
        double norm = Math.sqrt(prod(x, x));
        for (int i = 0; i < x.length; i++) x[i] /= norm;
        return norm;
    }

    public static double prod(double[] x, double[] y) {
        double result = 0.0D;
        int length = Math.min(x.length, y.length);
        for (int i = 0; i < length; i++) result += x[i] * y[i];
        return result;
    }

/*    public static void eigen(double[][] matrix, double[][] evecs, double[] evals) {
        double eps = 1.0E-06D;
        int maxiter = 100;
        int d = evals.length;
        int n = matrix.length;
        for (int m = 0; m < d; m++) {
            if (m > 0)
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        matrix[i][j] -= evals[(m - 1)] * evecs[(m - 1)][i] * evecs[(m - 1)][j];
            for (int i = 0; i < n; i++)
                evecs[m][i] = Math.random();
            Data.normalize(evecs[m]);

            double r = 0.0D;

            for (int iter = 0; (Math.abs(1.0D - r) > 1.0E-06D) && (iter < 100); iter++) {
                double[] q = new double[n];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++)
                        q[i] += matrix[i][j] * evecs[m][j];
                }
                evals[m] = Data.prod(evecs[m], q);
                Data.normalize(q);
                r = Math.abs(Data.prod(evecs[m], q));
                evecs[m] = q;
            }
        }
    }*/

    public static void eigen(double[][] matrix, double[][] evecs, double[] evals) {
        int d = evals.length;
        int k = matrix.length;
        double d1 = 1.E-05D;
        double r = 0.0D;
        for (int m = 0; m < d; m++) {
            double eps;
            evals[m] = Data.normalize(evecs[m]);
        }
        int iterations = 0;
        while (r < 0.9999900000000001D) {
            double[][] tempOld = new double[d][k];

            for (int m = 0; m < d; m++) {
                for (int i = 0; i < k; i++) {
                    tempOld[m][i] = evecs[m][i];
                    evecs[m][i] = 0.0D;
                }
            }

            for (int m = 0; m < d; m++) {
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < k; j++)
                        evecs[m][j] += matrix[i][j] * tempOld[m][i];
            }
            for (int m = 0; m < d; m++) {
                for (int p = 0; p < m; p++) {
                    double fac = Data.prod(evecs[p], evecs[m]) / Data.prod(evecs[p], evecs[p]);
                    for (int i = 0; i < k; i++) evecs[m][i] -= fac * evecs[p][i];
                }
            }

            for (int m = 0; m < d; m++) evals[m] = Data.normalize(evecs[m]);
            r = 1.0D;
            for (int m = 0; m < d; m++)
            {
                r = Math.min(Math.abs(Data.prod(evecs[m], tempOld[m])), r);
            }

            iterations++;
        }
    }


    public static void randomize(double[][] matrix) {
        Random random = new Random();
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] = random.nextDouble();
            }
        }
    }

    public static int[] landmarkIndices(double[][] matrix) {
        int k = matrix.length;
        int n = matrix[0].length;
        int[] result = new int[k];
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n; j++) {
                if (matrix[i][j] == 0.0D) {
                    result[i] = j;
                }
            }
        }
        return result;
    }

    public static double[][] copyMatrix(double[][] matrix) {
        double[][] copy = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                copy[i][j] = matrix[i][j];
            }
        }
        return copy;
    }
}
