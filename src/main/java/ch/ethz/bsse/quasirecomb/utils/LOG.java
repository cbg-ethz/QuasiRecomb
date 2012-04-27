/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.utils;

/**
 *
 * @author toepfera
 */
public class LOG {

    public static double ZERO = Double.NaN;

    public static double add(double a, double b) {
        if (isNaN(a) && isNaN(b)) {
            return ZERO;
        } else if (isNaN(a)) {
            return b;
        } else if (isNaN(b)) {
            return a;
        } else {
            return a + b;
        }
    }

    public static double mult(double a, double b) {
        if (isNaN(a) || isNaN(b)) {
            return ZERO;
        } else {
            return a * b;
        }
    }

    public static double div(double a, double b) {
        if (isNaN(a) || isNaN(b)) {
            throw new IllegalAccessError("ZERO");
        } else if (isNaN(a) || isNaN(b)) {
            return ZERO;
        } else {
            return a / b;
        }
    }

    private static boolean isNaN(double a) {
        return Double.isNaN(a);
    }
}
