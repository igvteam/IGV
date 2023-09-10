/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.feature.tribble.reader;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.FileInputStream;
import java.util.List;

import static org.junit.Assert.*;

/**
 * @author jrobinso
 * Date: 1/23/14
 * Time: 4:24 PM
 */
public class BamIndexTest {

    @Test
    public void testParse() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/refGene.hg19.bed.gz.tbi";

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        BamIndex tabixIndex = new BamIndex();
        tabixIndex.parse(bytes, null);
        assertNotNull(tabixIndex.sequenceIndexMap);


    }

    @Test
    public void testTBIIndex() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/refGene.hg19.bed.gz.tbi";
        int refID = 26;
        int beg = 55369194;
        int end = 55369443;

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        BamIndex tbiIndex = new BamIndex();
        tbiIndex.parse(bytes, null);

        List<BamIndex.Chunk> chunks = tbiIndex.chunksForRange(refID, beg, end);
        assertNotNull(tbiIndex.sequenceIndexMap);
        assertEquals(chunks.size(), 1);
        assertEquals(1640062, chunks.get(0).left.block);

    }

    @Test
    public void testCSIIndex() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/csi-test.vcf.gz.csi";
        int refID = 0;
        int beg = 1226000;
        int end = 1227000;

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        CSIIndex csiIndex = new CSIIndex();
        csiIndex.parse(bytes, null);

        List<CSIIndex.Chunk> chunks = csiIndex.chunksForRange(refID, beg, end);
        assertTrue(chunks.size() > 0);

    }

    @Test
    public void testCSIQuery() throws Exception {

        String path = TestUtils.DATA_DIR + "tabix/csi-test.vcf.gz.csi";
        int refID = 0;
        int beg = 1226000;
        int end = 1227000;

        FileInputStream fis = new FileInputStream(path);
        byte[] bytes = fis.readAllBytes();

        CSIIndex csiIndex = new CSIIndex();
        csiIndex.parse(bytes, null);

        List<CSIIndex.Chunk> chunks = csiIndex.chunksForRange(refID, beg, end);
        assertEquals(chunks.size(),  1);
        assertEquals(74500, chunks.get(0).left.block);
        assertEquals(61246, chunks.get(0).left.offset);
        assertEquals(94804, chunks.get(0).right.block);
    }


}
