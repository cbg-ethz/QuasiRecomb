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
package ch.ethz.bsse.quasirecomb.informationholder;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
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
