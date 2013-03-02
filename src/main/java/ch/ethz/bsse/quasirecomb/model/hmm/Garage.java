/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ch.ethz.bsse.quasirecomb.model.hmm;

import ch.ethz.bsse.quasirecomb.informationholder.Globals;
import ch.ethz.bsse.quasirecomb.informationholder.TempJHMMStorage;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author XLR
 */
public class Garage {

    protected Map<Integer, TempJHMMStorage> garage = new ConcurrentHashMap<>();
    protected final List<Integer> available = new ArrayList<>();

    protected void clearGarage(int L, int K, int n) {
        if (Globals.getINSTANCE().isSTORAGE()) {
            available.clear();
            garage.clear();
            for (int i = 0; i < Globals.getINSTANCE().getCpus(); i++) {
                available.add(i);
                garage.put(i, new TempJHMMStorage(L, K, n, i));
            }
        }
    }

    protected TempJHMMStorage mergeGarage() throws IllegalStateException {
        if (available.size() != Globals.getINSTANCE().getCpus()) {
            throw new IllegalStateException("Not all storages have been returned");
        }
        Iterator<TempJHMMStorage> iterator = this.garage.values().iterator();
        TempJHMMStorage store = iterator.next();
        while (iterator.hasNext()) {
            store.merge(iterator.next());
        }
        return store;
    }

    public TempJHMMStorage getStorage() {
        synchronized (this.available) {
            while (!available.iterator().hasNext()) {
                try {
                    notify();
                    TimeUnit.MILLISECONDS.sleep(10);
                    System.err.println("sleep");
                } catch (InterruptedException ex) {
                    Logger.getLogger(JHMM.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
//            System.out.println("GET " + s++);
            Integer i = available.iterator().next();
            available.remove(i);
            return garage.get(i);
        }
    }

    public void free(int id) {
        synchronized (this.available) {
            this.available.add(id);
        }
    }
}
