/**
 * Copyright (c) 2011-2012 Armin Töpfer
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

import ch.ethz.bsse.quasirecomb.informationholder.Read;
import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Plot {

    public static void plotCoverage(Read[] reads) {
        Multiset<Integer> coverage = HashMultiset.create();
        int rLength = reads.length;
        int j = 0;
        for (Read r : reads) {
            int amount = r.getCount();
            for (int i = r.getWatsonBegin(); i < r.getWatsonEnd(); i++) {
                coverage.add(i, amount);
            }
            for (int i = r.getCrickBegin(); i < r.getCrickEnd(); i++) {
                coverage.add(i, amount);
            }
            Globals.getINSTANCE().print("Plotting\t" + Math.round(1000 * ++j * 50d / rLength) / 1000 + "%");
        }
        XYSeries dataset = new XYSeries("Coverage");
        int alignmentLength = Globals.getINSTANCE().getALIGNMENT_END() - Globals.getINSTANCE().getALIGNMENT_BEGIN();
        for (int i = Globals.getINSTANCE().getALIGNMENT_BEGIN(); i < Globals.getINSTANCE().getALIGNMENT_END(); i++) {
            if (!coverage.contains(i)) {
                dataset.add((double) i, 1.0);
            } else {
                dataset.add((double) i, (double) coverage.count(i));
            }
            Globals.getINSTANCE().print("Plotting\t" + (50 + Math.round(1000 * (i - Globals.getINSTANCE().getALIGNMENT_BEGIN()) * 50d / alignmentLength) / 1000) + "%");
        }
        XYSeriesCollection collection = new XYSeriesCollection(dataset);
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "Coverage",
                "Position",
                "Coverage",
                collection,
                PlotOrientation.VERTICAL,
                false,
                false,
                false);

        final XYPlot plot = chart.getXYPlot();
//        final NumberAxis domainAxis = new NumberAxis("Position");
//        final NumberAxis rangeAxis = new LogarithmicAxis("Coverage (log)");
//        plot.setDomainAxis(domainAxis);
//        plot.setRangeAxis(rangeAxis);
        chart.setBackgroundPaint(Color.white);
        plot.setOutlinePaint(Color.black);
        plot.setBackgroundPaint(Color.white);
//        plot.setDomainGridlinePaint(Color.gray);
        plot.setRangeGridlinePaint(Color.decode("0xadadad"));

        try {
            ChartUtilities.saveChartAsPNG(new File(Globals.getINSTANCE().getSAVEPATH() + "coverage.png"),
                    chart, 1000, 500);
        } catch (IOException ex) {
            Logger.getLogger(Plot.class.getName()).log(Level.SEVERE, null, ex);
        }
        Globals.getINSTANCE().print("Plotting\t100%");
        System.out.println("");
    }

    public static void plot(double[] data, String name) {
        Multiset<Integer> coverage = HashMultiset.create();
        XYSeries dataset = new XYSeries(name);
        for (int i = 0; i < data.length; i++) {
            dataset.add((double) i, data[i]);
        }
        XYSeriesCollection collection = new XYSeriesCollection(dataset);
//        collection.addSeries(dataset);
        final JFreeChart chart = ChartFactory.createXYLineChart(
                name,
                "Position",
                "Coverage",
                collection,
                PlotOrientation.VERTICAL,
                false,
                false,
                false);

        final XYPlot plot = chart.getXYPlot();
//        final NumberAxis domainAxis = new NumberAxis("Position");
//        final NumberAxis rangeAxis = new LogarithmicAxis("Coverage (log)");
//        plot.setDomainAxis(domainAxis);
//        plot.setRangeAxis(rangeAxis);
        chart.setBackgroundPaint(Color.white);
        plot.setOutlinePaint(Color.black);
        plot.setBackgroundPaint(Color.white);
//        plot.setDomainGridlinePaint(Color.gray);
        plot.setRangeGridlinePaint(Color.decode("0xadadad"));

        try {
            ChartUtilities.saveChartAsPNG(new File(Globals.getINSTANCE().getSAVEPATH() + name + ".png"),
                    chart, 1000, 500);
        } catch (IOException ex) {
            Logger.getLogger(Plot.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
