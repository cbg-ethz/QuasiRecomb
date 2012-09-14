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

    private byte[] watsonSequence;
    private int watsonBegin;
    private int watsonEnd;
    private int count;
    private byte[] crickSequence;
    private int crickBegin;
    private int crickEnd = -1;
    private int end;

    public Read(byte[] sequence, int begin, int end) {
        this.watsonSequence = sequence;
        this.watsonBegin = begin;
        this.watsonEnd = end;
    }

    public Read(byte[] sequence, int begin, int end, int count) {
        this.watsonSequence = sequence;
        this.watsonBegin = begin;
        this.watsonEnd = end;
        this.count = count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public int getBegin() {
        return this.watsonBegin;
    }

    public void incCount() {
        count++;
    }

    public int getCount() {
        return count;
    }

    public int getEnd() {
        if (this.crickEnd == -1) {
            return end;
        } else {
            return this.crickEnd;
        }
    }

    public byte[] getSequence() {
        return this.watsonSequence;
    }
    
    public boolean isPaired() {
        return this.crickSequence != null;
    }
    
    public void setPairedEnd(byte[] sequence, int begin, int end) {
        this.crickSequence = sequence;
        this.crickBegin = begin;
        this.crickEnd = end;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Arrays.hashCode(this.watsonSequence);
        hash = 29 * hash + Arrays.hashCode(this.crickSequence);
        hash = 29 * hash + this.watsonBegin;
        hash = 29 * hash + this.watsonEnd;
        hash = 29 * hash + this.crickBegin;
        hash = 29 * hash + this.crickEnd;
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        return obj != null && obj.getClass() == this.getClass() && obj.hashCode() == this.hashCode();
    }
}
