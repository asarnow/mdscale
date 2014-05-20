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
 * A class implementing Stress Minimization by Majorizing a Complex Function (SMACOF).
 *
 * Created by IntelliJ IDEA.
 * User: da
 * Date: 10/8/11
 * Time: 2:43 AM
 */

public class SMACOF {
    private double[][] x;
    private double[][] d;
    private double[][] w;

    /**
     * Construct a new SMACOF instance.
     * @param d distance matrix
     * @param x initial coordinate matrix
     * @param w weights matrix
     */
    public SMACOF(double[][] d, double[][] x, double[][] w)
    {
        this.x = x;
        this.d = d;
        this.w = w;
    }

    /**
     * Construct a SMACOF instance without weights.
     * @param d distance matrix
     * @param x initial coordinate matrix
     */
    public SMACOF(double[][] d, double[][] x)
    {
        this.x = x;
        this.d = d;
//        this.w = weightMatrix(d, 0.0D);
        this.w = null;
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

    /**
     * Perform 1 majorization iteration using this SMACOF instance.
     * @return report
     */
    public String iterate()
    {
        return iterate(1);
    }

    /**
     * Perform n majorization iterations using this SMACOF instance.
     * @param n number of iterations
     * @return report
     */
    public String iterate(int n)
    {
        if (this.w != null)
            return majorize(this.x, this.d, this.w, n, 0, 0);
        return majorize(this.x, this.d, n, 0, 0);
    }

    /**
     * Perform majorization iterations until the maximum number of iterations
     * is reached, the maximum runtime has elapsed, or the change in
     * normalized stress falls below the threshold, whichever comes first.
     * @param iter maximum number of iterations
     * @param threshold threshold for change in normalized stress
     * @param timeout maximum runtime in milliseconds
     * @return report
     */
    public String iterate(int iter, int threshold, int timeout)
    {
        if (this.w !=null)
            return majorize(this.x, this.d, this.w, iter, threshold, timeout);
        return majorize(this.x, this.d, iter, threshold, timeout);
    }

    /**
     * Compute the absolute stress for this SMACOF instance.
     * @return stress
     */
    public double getStress()
    {
        if (this.w != null)
            return stress(this.d, this.w, this.x);
        return stress(this.d,this.x);
    }

    /**
     * Compute the normalized stress for this SMACOF instance.
     * @return normalized stress
     */
    public double getNormalizedStress()
    {
        if (this.w !=null)
            return normalizedStress(this.d, this.w, this.x);
        return normalizedStress(this.d,this.x);
    }

    /**
     * Element-wise matrix exponentiation for self-weighting of distances.
     * @param D distance matrix or initial weights
     * @param exponent power to raise each element the matrix
     * @return exponentiated weights
     */
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

    /**
     * SMACOF algorithm (weighted).
     * @param x coordinates matrix
     * @param d distance matrix
     * @param w weights matrix
     * @param iter maximum iterations
     * @param threshold halting threshold for change in normalized stress
     * @param timeout maximum runtime in milliseconds
     * @return report
     */
    public static String majorize(double[][] x, double[][] d, double[][] w, int iter, int threshold, int timeout)
    {
        String report = "";
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;

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
                        inv += Math.pow(x[m][i] - x[m][j], 2.0D);
                    }
                    if (inv != 0.0D) inv = Math.pow(inv, -0.5D);
                    for (int m = 0; m < dim; m++) {
                        xnew[m] += w[j][i] * (x[m][j] +
                                d[j][i] * (x[m][i] - x[m][j]) * inv);
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

    /**
     * SMACOF algorithm (unweighted).
     * @param x coordinates matrix
     * @param d distance matrix
     * @param iter maximum iterations
     * @param threshold halting threshold for change in normalized stress
     * @param timeout maximum runtime in milliseconds
     * @return report
     */
    public static String majorize(double[][] x, double[][] d, int iter, int threshold, int timeout)
        {
            String report = "";
            int n = x[0].length;
            int k = d.length;
            int dim = x.length;

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
                            inv += Math.pow(x[m][i] - x[m][j], 2.0D);
                        }
                        if (inv != 0.0D) inv = Math.pow(inv, -0.5D);
                        for (int m = 0; m < dim; m++) {
                            xnew[m] += (x[m][j] +
                                    d[j][i] * (x[m][i] - x[m][j]) * inv);
                        }
                    }
                    for (int m = 0; m < dim; m++) {
                        change += Math.pow(xnew[m] / n - x[m][i], 2.0D);
                        magnitude += Math.pow(x[m][i], 2.0D);
                        x[m][i] = (xnew[m] / n);
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

    /**
     * Bare SMACOF algorithm. Convenient for reading the algorithm.
     * @param x coordinates matrix
     * @param d distance matrix
     * @param w weights matrix
     * @param iter number of iterations
     */
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

    /**
     * Compute the absolute stress between a weighted distance matrix and a configuration of coordinates.
     * @param d distance matrix
     * @param w weights matrix
     * @param x coordinates matrix
     * @return stress
     */
    public static double stress(double[][] d, double[][] w, double[][] x) {
        double result = 0.0D;
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;

        for (int i = 0; i < k; i++) {
            for (int j = i+1; j < n; j++) {
                double dist = 0.0D;
                for (int m = 0; m < dim; m++)
                    dist += Math.pow(x[m][i] - x[m][j], 2.0D);
                result += w[i][j] * Math.pow(d[i][j] - Math.sqrt(dist), 2.0D);
            }
        }
        return result;
    }

    /**
     * Compute the absolute stress between a distance matrix and a configuration of coordinates.
     * @param d distance matrix
     * @param x coordinates matrix
     * @return stress
     */
    public static double stress(double[][] d, double[][] x) {
        double result = 0.0D;
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;

        for (int i = 0; i < k; i++) {
            for (int j = i+1; j < n; j++) {
                double dist = 0.0D;
                for (int m = 0; m < dim; m++)
                    dist += Math.pow(x[m][i] - x[m][j], 2.0D);
                result += Math.pow(d[i][j] - Math.sqrt(dist), 2.0D);
            }
        }
        return result;
    }

    /**
     * Compute the normlized stress between a weighted distance matrix and a configuration of coordinates.
     * @param d distance matrix
     * @param w weights matrix
     * @param x coordinate matrix
     * @return normalized stress
     */
    public static double normalizedStress(double[][] d, double[][] w, double[][] x) {
        double result = 0.0D;
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;

        double sum = 0.0D;
        for (int i = 0; i < k; i++) {
            for (int j = i+1; j < n; j++) {
                double dist = 0.0D;
                for (int m = 0; m < dim; m++)
                    dist += Math.pow(x[m][i] - x[m][j], 2.0D);
                result += w[i][j] * Math.pow(d[i][j] - Math.sqrt(dist), 2.0D);
                sum += w[i][j] * Math.pow(d[i][j], 2.0D);
            }
        }
        return result / sum;
    }

    /**
     * Return the normlized stress between a distance matrix and a configuration of coordinates.
     * @param d distance matrix
     * @param x coordinate matrix
     * @return normalized stress
     */
    public static double normalizedStress(double[][] d, double[][] x) {
        double result = 0.0D;
        int n = x[0].length;
        int k = d.length;
        int dim = x.length;

        double sum = 0.0D;
        for (int i = 0; i < k; i++) {
            for (int j = i+1; j < n; j++) {
                double dist = 0.0D;
                for (int m = 0; m < dim; m++)
                    dist += Math.pow(x[m][i] - x[m][j], 2.0D);
                result += Math.pow(d[i][j] - Math.sqrt(dist), 2.0D);
                sum += Math.pow(d[i][j], 2.0D);
            }
        }
        return result / sum;
    }
}
