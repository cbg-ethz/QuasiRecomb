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
package ch.ethz.bsse.quasirecomb.informatioholder;

import java.io.Serializable;
import java.util.List;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class OptimalResultsList implements Serializable{
    private List<OptimalResult> ors;

    public OptimalResultsList(List<OptimalResult> ors) {
        this.ors = ors;
    }

    public List<OptimalResult> getOrs() {
        return ors;
    }

    public void setOrs(List<OptimalResult> ors) {
        this.ors = ors;
    }
}
