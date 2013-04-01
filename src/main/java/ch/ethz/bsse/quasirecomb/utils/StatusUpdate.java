/**
 * Copyright (c) 2011-2013 Armin Töpfer
 *
 * This file is part of QuasiRecomb.
 *
 * QuasiRecomb is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * QuasiRecomb is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * QuasiRecomb. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.quasirecomb.utils;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class StatusUpdate {

    private String oldOut = "";
    private double hammingCount = 0;
    private double PERCENTAGE = 0;
    private long start = System.currentTimeMillis();
    private final DateFormat df = new SimpleDateFormat("HH:mm:ss");
    private static final StatusUpdate INSTANCE = new StatusUpdate();

    public static StatusUpdate getINSTANCE() {
        return INSTANCE;
    }

    private StatusUpdate() {
        df.setTimeZone(TimeZone.getTimeZone("GMT"));
    }

    public void incPercentage() {
        PERCENTAGE += 100d / Globals.getINSTANCE().getREPEATS();
    }

    public void printBIC(int K, int bic) {
        System.out.print("\r                                                                                                                                                   ");
        if (Globals.getINSTANCE().isMODELSELECTION()) {
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        } else {
            System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + Math.round(PERCENTAGE * 1000) / 1000 + "%\t[BIC: " + (int) bic + "]                 ");
        }
    }

    public void printBIC(int K, int percentage, int bic) {
        System.out.print("\r                                                                                                                                                   ");
        if (Globals.getINSTANCE().isMODELSELECTION()) {
            System.out.print("\r" + time() + " Model selection [K " + K + "]:\t" + percentage + "%\t[BIC: " + (int) bic + "]                 ");
        } else {
            System.out.print("\r" + time() + " Model training  [K " + K + "]:\t" + percentage + "%\t[BIC: " + (int) bic + "]                 ");
        }
    }

    public void print(String s) {
        if (!oldOut.equals(s)) {
            this.oldOut = s;
            System.out.print("\r" + time() + " " + s);
        }
    }

    public void println(String s) {
        System.out.print("\n" + time() + " " + s);
    }

    public void printPercentage(int K, double read, double Kmin) {
        if (!Globals.getINSTANCE().isSILENT()) {
            if (!oldOut.equals(time())) {
                this.oldOut = time();
                System.out.print("\r                                                                                                                                                   ");
                System.out.print("\r" + time() + " Model " + (Globals.getINSTANCE().isMODELSELECTION() ? "selection" : "training ") + " [K " + (int) Kmin + "]:\t" + (Math.round(PERCENTAGE * 1000) / 1000) + "%");
            }
        }
    }

    public synchronized void incHamming(int inc) {
        hammingCount += inc * (100d / Globals.getINSTANCE().getHammingMax());
    }

    public synchronized void printHamming(int inc) {
        incHamming(inc);
        System.out.print("\r" + time() + " Computing:\t" + Math.round(hammingCount * 1000) / 1000 + "%");
    }

    public String time() {
        return df.format(new Date(System.currentTimeMillis() - start));
    }

    public double getPERCENTAGE() {
        return PERCENTAGE;
    }

    public void setPERCENTAGE(double PERCENTAGE) {
        this.PERCENTAGE = PERCENTAGE;
    }
}
