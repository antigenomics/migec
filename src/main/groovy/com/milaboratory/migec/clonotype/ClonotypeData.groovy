/*
 * Copyright (c) 2014-2024, OOO «MiLaboratory»
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE LICENSOR HAS NO
 * OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE LICENSOR MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
 * PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * TRADEMARK OR OTHER RIGHTS.
 */

package com.milaboratory.migec.clonotype

import com.milaboratory.migec.Util
import com.milaboratory.migec.dalign.SegmentSearcher
import com.milaboratory.migec.segment.Allele

import static com.milaboratory.migec.Util.BLANK_FIELD

class ClonotypeData {
    Allele vAllele, jAllele
    Allele[] dAlleles
    String cdr3Seq
    int cdrFrom, cdrTo,
        vEndRel, jStartRel,
        dStartRel = -1, dEndRel = -1
    boolean rc
    int[] counts = new int[4] // counters: good events/total events/good reads/total reads with this V-J-CDR3
    String aaSeq

    ClonotypeData(Allele vAllele, Allele jAllele,
                  String cdr3Seq,
                  int cdrFrom, int cdrTo, int vEnd, int jStart,
                  boolean rc) {
        this.vAllele = vAllele
        this.jAllele = jAllele
        this.cdr3Seq = cdr3Seq
        this.cdrFrom = cdrFrom
        this.cdrTo = cdrTo
        this.vEndRel = vEnd - cdrFrom
        this.jStartRel = jStart - cdrFrom
        this.rc = rc
    }

    boolean isCanonical() {
        aaSeq = aaSeq ?: Util.translate(cdr3Seq)
        aaSeq.contains("?") || aaSeq.contains("*") || aaSeq.endsWith("F") || aaSeq.endsWith("W")
    }

    String getdAllele() {
        (dAlleles && dAlleles.length > 0) ? dAlleles.collect { it.alleleId }.join(",") : BLANK_FIELD
    }

    String getSignature() {
        [cdr3Seq, aaSeq = aaSeq ?: Util.translate(cdr3Seq),
         vAllele.alleleId, jAllele.alleleId, dAllele,
         vEndRel, dStartRel, dEndRel, jStartRel].join("\t") // for tabular output
    }

    String getShortSignature() {
        vAllele.alleleId + "|" + jAllele.alleleId + "|" + cdrFrom + "|" + cdr3Seq // for hash
    }

    void appendDInfo(SegmentSearcher dSegmentSearcher) {
        if (vAllele.hasD) {
            def result = dSegmentSearcher.scan(cdr3Seq, vEndRel, jStartRel)
            if (result) {
                dStartRel = result.segmentStart
                dEndRel = result.segmentEnd
                dAlleles = result.alleles as Allele[]
            }
        }
    }

    @Override
    boolean equals(Object obj) {
        this.shortSignature == (obj as ClonotypeData).shortSignature
    }

    @Override
    int hashCode() {
        this.shortSignature.hashCode()
    }
}
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
