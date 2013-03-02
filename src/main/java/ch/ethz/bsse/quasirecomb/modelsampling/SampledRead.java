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
package ch.ethz.bsse.quasirecomb.modelsampling;

import java.util.List;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SampledRead {

    public byte[] watsonReadBases;
    public byte[] crickReadBases;
    public int watsonStart;
    public int watsonEnd;
    public int crickStart;
    public int crickEnd;

    public SampledRead(byte[] watsonReadBases, byte[] crickReadBases, int watsonStart, int watsonEnd, int crickStart, int crickEnd) {
        this.watsonReadBases = watsonReadBases;
        this.crickReadBases = crickReadBases;
        this.watsonStart = watsonStart;
        this.watsonEnd = watsonEnd;
        this.crickStart = crickStart;
        this.crickEnd = crickEnd;
    }

    public SampledRead(List<Byte> watsonReadBases, List<Byte> crickReadBases, int watsonStart, int watsonEnd, int crickStart, int crickEnd) {
        this.watsonReadBases = new byte[watsonReadBases.size()];
        for (int i = 0; i < watsonReadBases.size(); i++) {
            this.watsonReadBases[i] = watsonReadBases.get(i);
        }
        this.crickReadBases = new byte[crickReadBases.size()];
        for (int i = 0; i < crickReadBases.size(); i++) {
            this.crickReadBases[i] = crickReadBases.get(i);
        }
        this.watsonStart = watsonStart;
        this.watsonEnd = watsonEnd;
        this.crickStart = crickStart;
        this.crickEnd = crickEnd;
    }

    public SampledRead(List<Byte> watsonReadBases, int watsonStart, int watsonEnd) {
        this.watsonReadBases = new byte[watsonReadBases.size()];
        for (int i = 0; i < watsonReadBases.size(); i++) {
            this.watsonReadBases[i] = watsonReadBases.get(i);
        }
        this.watsonStart = watsonStart;
        this.watsonEnd = watsonEnd;
    }

    public SampledRead(byte[] watsonReadBases, int watsonStart, int watsonEnd) {
        this.watsonReadBases = watsonReadBases;
        this.watsonStart = watsonStart;
        this.watsonEnd = watsonEnd;
    }

    public byte[] getWatsonReadBases() {
        return watsonReadBases;
    }

    public byte[] getCrickReadBases() {
        return crickReadBases;
    }

    public int getWatsonStart() {
        return watsonStart;
    }

    public int getWatsonEnd() {
        return watsonEnd;
    }

    public int getCrickStart() {
        return crickStart;
    }

    public int getCrickEnd() {
        return crickEnd;
    }
}
