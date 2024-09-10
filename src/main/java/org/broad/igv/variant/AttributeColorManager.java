package org.broad.igv.variant;

import org.broad.igv.renderer.ColorScale;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class AttributeColorManager {

    private static final Map<Key, ColorScale> attributes = new TreeMap<>();

    public static PaletteColorTable DEFAULT_BOOLEAN_COLORS = new PaletteColorTable(Color.GRAY);
    static {
        DEFAULT_BOOLEAN_COLORS.put("true", Color.BLACK);
        DEFAULT_BOOLEAN_COLORS.put("false", Color.RED);
    }

    public enum Type{
        INFO, FORMAT
    }

    private record Key(Type type, String name) implements Comparable<Key> {

        final static Comparator<Key> COMPARATOR = Comparator.comparing(Key::type)
                .thenComparing(Key::name);

        @Override
        public int compareTo(Key o) {
            return COMPARATOR.compare(this, o);
        }
    };

    public static PaletteColorTable getBooleanColorTable(Type type, String id){
        return (PaletteColorTable) attributes.computeIfAbsent(new Key(type, id), k -> DEFAULT_BOOLEAN_COLORS);
    }

    public static PaletteColorTable getColorTable(Type type, String id){
        return (PaletteColorTable) attributes.computeIfAbsent(new Key(type, id), k -> new PaletteColorTable(ColorUtilities.getPalette("Set 1")));
    }

    static {

        PaletteColorTable svtypeColors = new PaletteColorTable(ColorUtilities.getPalette("Pastel 1"));
        //from VCF 4.5 spec section 1.4.5
        List.of("DEL",
                "INS",
                "DUP",
                "INV",
                "CNV",
                "CNV:TR",
                "DUP:TANDEM",
                "DEL:ME",
                "INS:ME")
                .forEach(svtypeColors::get);
        attributes.put(new Key(Type.INFO, "SVTYPE"), svtypeColors);

        PaletteColorTable defaultBooleanColors = new PaletteColorTable(Color.GRAY);
        defaultBooleanColors.put("true", Color.BLACK);
        defaultBooleanColors.put("false", Color.RED);


    }


}
