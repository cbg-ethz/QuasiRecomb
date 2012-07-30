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
package ch.ethz.bsse.quasirecomb.informationholder;

import java.util.Arrays;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Read {

    private byte[] sequence;
    private int begin;
    private int end;
    private int count;

    public Read(byte[] sequence, int begin, int end) {
        this.sequence = sequence;
        this.begin = begin;
        this.end = end;
    }
    public Read(byte[] sequence, int begin, int end, int count) {
        this.sequence = sequence;
        this.begin = begin;
        this.end = end;
        this.count = count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public int getBegin() {
        return begin;
    }

    public void incCount() {
        count++;
    }
    public int getCount() {
        return count;
    }

    public int getEnd() {
        return end;
    }

    public byte[] getSequence() {
        return sequence;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Arrays.hashCode(this.sequence);
        hash = 29 * hash + this.begin;
        hash = 29 * hash + this.end;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        return obj != null && obj.getClass() == this.getClass() && obj.hashCode() == this.hashCode();
    }
}
