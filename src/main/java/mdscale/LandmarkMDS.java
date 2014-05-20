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

/**
 * Created by IntelliJ IDEA.
 * User: da
 * Date: 10/8/11
 * Time: 2:43 AM
 */
public class LandmarkMDS {
    private double[][] x;
    private double[][] d;
    private double[][] w;

    public LandmarkMDS(double[][] d, double[][] x, double[][] w)
    {
        this.x = x;
        this.d = d;
        this.w = w;
    }

    public LandmarkMDS(double[][] d, double[][] x)
    {
        this.x = x;
        this.d = d;
        this.w = weightMatrix(d, 0.0D);
    }

    public double[][] getDissimilarities()
    {
        return this.d;
    }

    public double[][] getWeights()
    {
        return this.w;
    }

    public double[][] getPositions()
    {
        return this.x;
    }

    public void setDissimilarities(double[][] d)
    {
        this.d = d;
    }

    public void setWeights(double[][] w)
    {
        this.w = w;
    }

    public void setPositions(double[][] x)
    {
        this.x = x;
    }

    public String iterate()
    {
        return iterate(1);
    }

    public String iterate(int n)
    {
        return majorize(this.x, this.d, this.w, n, 0, 0);
    }

    public String iterate(int iter, int threshold, int timeout)
    {
        return majorize(this.x, this.d, this.w, iter, threshold, timeout);
    }

    public double getStress()
    {
        return stress(this.d, this.w, this.x);
    }

    public double getNormalizedStress()
    {
        return normalizedStress(this.d, this.w, this.x);
    }

    public static double[][] weightMatrix(double[][] D, double exponent)
    {
        int n = D[0].length;
        int k = D.length;
        double[][] result = new double[k][n];
        for (int i = 0; i < k; i++)
            for (int j = 0; j < n; j++)
                if (D[i][j] > 0.0D)
                    result[i][j] = Math.pow(D[i][j], exponent);
        return result;
    }

    public static String majorize(double[][] x, double[][] d, double[][] w, int iter, int threshold, int timeout)
    {
        String report = "";
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;
        int[] index = Data.landmarkIndices(d);
        double[] wSum = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                wSum[i] += w[j][i];
            }
        }
        double eps = Math.pow(10.0D, -threshold);
        long time = System.nanoTime();
        if (iter == 0)
            iter = 10000000;
        for (int c = 0; c < iter; c++) {
            double change = 0.0D;
            double magnitude = 0.0D;
            for (int i = 0; i < n; i++) {
                double[] xnew = new double[dim];
                for (int j = 0; j < k; j++) {
                    double inv = 0.0D;
                    for (int m = 0; m < dim; m++) {
                        inv += Math.pow(x[m][i] - x[m][index[j]], 2.0D);
                    }
                    if (inv != 0.0D) inv = Math.pow(inv, -0.5D);
                    for (int m = 0; m < dim; m++) {
                        xnew[m] += w[j][i] * (x[m][index[j]] +
                                d[j][i] * (x[m][i] - x[m][index[j]]) * inv);
                    }
                }
                if (wSum[i] != 0.0D) {
                    for (int m = 0; m < dim; m++) {
                        change += Math.pow(xnew[m] / wSum[i] - x[m][i], 2.0D);
                        magnitude += Math.pow(x[m][i], 2.0D);
                        x[m][i] = (xnew[m] / wSum[i]);
                    }
                }
            }
            change = Math.sqrt(change / magnitude);
            long timediff = (System.nanoTime() - time) / 1000000L;

            if ((timeout > 0) && (timediff > timeout)) {
                return c + 1 + " iterations, " +
                        timediff + " milliseconds, " + change + " relative change";
            }
            if ((threshold > 0) && (change < eps)) {
                return c + 1 + " iterations, " +
                        timediff + " milliseconds, " + change + " relative change";
            }
            if ((iter > 0) && (c >= iter - 1)) {
                report = c + 1 + " iterations, " +
                        timediff + " milliseconds, " + change + " relative change";
            }
        }
        return report;
    }

    public static void majorize(double[][] x, double[][] d, double[][] w, int iter)
    {
        int n = x[0].length;
        int dim = x.length;
        double[] wSum = new double[n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) wSum[i] += w[i][j];
        for (int c = 0; c < iter; c++)
            for (int i = 0; i < n; i++) {
                double[] xnew = new double[dim];
                for (int j = 0; j < n; j++) {
                    double inv = 0.0D;
                    for (int k = 0; k < dim; k++) inv += Math.pow(x[k][i] - x[k][j], 2.0D);
                    if (inv != 0.0D) inv = Math.pow(inv, -0.5D);
                    for (int k = 0; k < dim; k++) xnew[k] += w[i][j] * (x[k][j] + d[i][j] * (x[k][i] - x[k][j]) * inv);
                }
                if (wSum[i] != 0.0D)
                    for (int k = 0; k < dim; k++) x[k][i] = (xnew[k] / wSum[i]);
            }
    }

    public static double stress(double[][] d, double[][] w, double[][] x) {
        double result = 0.0D;
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;
        int[] index = Data.landmarkIndices(d);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n; j++) {
                double dist = 0.0D;
                for (int m = 0; m < dim; m++)
                    dist += Math.pow(x[m][index[i]] - x[m][j], 2.0D);
                result += w[i][j] * Math.pow(d[i][j] - Math.sqrt(dist), 2.0D);
            }
        }
        return result;
    }

    public static double normalizedStress(double[][] d, double[][] w, double[][] x) {
        double result = 0.0D;
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;
        int[] index = Data.landmarkIndices(d);
        double sum = 0.0D;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n; j++) {
                double dist = 0.0D;
                for (int m = 0; m < dim; m++)
                    dist += Math.pow(x[m][index[i]] - x[m][j], 2.0D);
                result += w[i][j] * Math.pow(d[i][j] - Math.sqrt(dist), 2.0D);
                sum += w[i][j] * Math.pow(d[i][j], 2.0D);
            }
        }
        return result / sum;
    }
}
