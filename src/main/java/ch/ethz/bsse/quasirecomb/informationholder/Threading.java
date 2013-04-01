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

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

/**
 * Threading stores objects necessary for a parallel environment.
 *
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Threading {

    private static final Threading INSTANCE = new Threading();
    private static final BlockingQueue<Runnable> BLOCKING_QUEUE = new ArrayBlockingQueue<>(Runtime.getRuntime().availableProcessors() - 1);
    private static final RejectedExecutionHandler REH = new ThreadPoolExecutor.CallerRunsPolicy();
    private static ExecutorService EXECUTOR = refreshExecutor();

    private Threading() {
    }

    public static Threading getINSTANCE() {
        return INSTANCE;
    }
    
    private static ExecutorService refreshExecutor() {
//        return Executors.newSingleThreadExecutor();
        return new ThreadPoolExecutor(Runtime.getRuntime().availableProcessors() - 1, Runtime.getRuntime().availableProcessors() - 1, 5L, TimeUnit.MINUTES, BLOCKING_QUEUE, REH);
    }

    public static void renewExecutor() {
        EXECUTOR = refreshExecutor();
    }

    public ExecutorService getExecutor() {
        return EXECUTOR;
    }
}
