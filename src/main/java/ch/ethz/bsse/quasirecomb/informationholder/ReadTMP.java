/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.informationholder;

/**
 *
 * @author toepfera
 */
public class ReadTMP {
    public String name;
    public double[] quality;
    public byte[] readBases;
    public int refStart;
    public boolean hasQuality;

    public ReadTMP(String name, double[] quality, byte[] readBases, int refStart, boolean hasQuality) {
        this.name = name;
        this.quality = quality;
        this.readBases = readBases;
        this.refStart = refStart;
        this.hasQuality = hasQuality;
    }

    
}
