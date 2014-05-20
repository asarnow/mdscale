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
 * Time: 2:25 AM
 */
public class MDS {

    public static double[][] fullmds(double[][] d, int dim) {
        double[][] result = new double[dim][d.length];
        Data.randomize(result);
        double[] evals = new double[result.length];
        Data.squareEntries(d);
        Data.doubleCenter(d);
        Data.multiply(d, -0.5D);
        Data.eigen(d, result, evals);
        for (int i = 0; i < result.length; i++) {
            evals[i] = Math.sqrt(evals[i]);
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] *= evals[i];
            }
        }
        return result;
    }

    public static double[][] classicalScaling(double [][] distances, int dim) {
        double[][] d = Data.copyMatrix(distances);
        return fullmds(d, dim);
    }

    public static double[][] distanceScaling(double[][] d, int dim) {
        double[][] x = classicalScaling(d, dim);
        SMACOF smacof = new SMACOF(d, x);
        smacof.iterate(100, 3, 10*60000);
        return smacof.getPositions();
    }

    public static double[][] distanceScaling(double[][] d, int dim, int iter, int threshold, int timeout) {
        double[][] x = classicalScaling(d, dim);
        SMACOF smacof = new SMACOF(d, x);
        smacof.iterate(iter, threshold, timeout);
        return smacof.getPositions();
    }
}
