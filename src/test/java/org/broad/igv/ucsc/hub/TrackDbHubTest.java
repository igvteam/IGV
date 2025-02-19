package org.broad.igv.ucsc.hub;

import org.junit.Test;

import java.util.Map;

import static org.junit.Assert.*;

public class TrackDbHubTest {

    //metadata differentiation="10 hour" treatment=X donor=A lab="List Meta Lab" data_set_id=ucscTest1 access=group assay=long-RNA-seq enriched_in=exon life_stage=postpartum species="Homo sapiens" ucsc_db=hg38

    @Test
    public void parseMetadata() {

        String metadata = "differentiation=\"10 hour\" treatment=X donor=A lab=\"List Meta Lab\" data_set_id=ucscTest1 access=group assay=long-RNA-seq enriched_in=exon life_stage=postpartum species=\"Homo sapiens\" ucsc_db=hg38";
        Map<String, String> metadataMap = TrackDbHub.parseMetadata(metadata);
        assertEquals(11, metadataMap.size());
        assertEquals("X", metadataMap.get("Treatment"));
        assertEquals("10 hour", metadataMap.get("Differentiation"));
        assertEquals("A", metadataMap.get("Donor"));
        assertEquals("group", metadataMap.get("Access"));
        assertEquals("hg38", metadataMap.get("Ucsc_db"));
        assertEquals("Homo sapiens", metadataMap.get("Species"));
        assertEquals("postpartum", metadataMap.get("Life_stage"));
        assertEquals("long-RNA-seq", metadataMap.get("Assay"));
        assertNotNull(metadataMap);
    }
}