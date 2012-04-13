package ch.ethz.bsse.quasirecomb.align;

/**
 *
 * @author toepfera
 */
public class Fastq2 {
    private byte[] read;
    private byte[] quality;
    private String description;
    private String identifier;

    public Fastq2(byte[] read, byte[] quality, String description, String identifier) {
        this.read = read;
        this.quality = quality;
        this.description = description;
        this.identifier = identifier;
    }

    public Fastq2() {
    }

    public byte[] getRead() {
        return read;
    }

    public byte[] getQuality() {
        return quality;
    }

    public String getDescription() {
        return description;
    }

    public String getIdentifier() {
        return identifier;
    }
    
    public boolean isEmpty() {
        return read == null;
    }

    public void setRead(byte[] read) {
        this.read = read;
    }

    public void setQuality(byte[] quality) {
        this.quality = quality;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
}
