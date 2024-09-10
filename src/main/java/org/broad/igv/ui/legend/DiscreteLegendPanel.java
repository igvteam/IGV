package org.broad.igv.ui.legend;

import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.PaletteColorTable;

import java.awt.*;
import java.util.Map;

public class DiscreteLegendPanel extends LegendPanel {
    protected PaletteColorTable colorTable;

    public DiscreteLegendPanel(PaletteColorTable colorTable){
        this.colorTable = colorTable;
    }

    protected void reloadPreferences() {
    }

    protected void persistCurrentPreferences() {
    }

    @Override
    public ColorTable getColorScale() {
        return colorTable;
    }

    @Override
    protected void resetPreferencesToDefault() {
    }

    /**
     * Open the user preferences dialog
     */
    @Override
    public void edit() {}

    @Override
    public void paintLegend(Graphics2D g2D) {

        if (colorTable == null) {
            return;
        }

        g2D.setFont(FontManager.getFont(10));

        FontMetrics fm = g2D.getFontMetrics();
        int dh = fm.getHeight() / 2 + 3;

        int x = 0;
        int lineHeight = 12;
        int y = lineHeight;
        int colCount = 0;

        for (Map.Entry<String, Color> entry : colorTable.entrySet()) {

            String mutType = entry.getKey();
            String label = mutType.replace("_", " ");
            int labelWidth = (int) fm.getStringBounds(label, g2D).getWidth();

            g2D.setColor(entry.getValue());
            g2D.fillRect(x, y, 10, 10);
            g2D.setColor(Color.BLACK);
            g2D.drawRect(x, y, 10, 10);
            g2D.drawString(label, x + 20, y + dh);
            x += labelWidth + 40;
            colCount++;

            if (colCount % 5 == 0) {
                y += lineHeight + 5;
                x = 0;
            }
        }

    }
}
