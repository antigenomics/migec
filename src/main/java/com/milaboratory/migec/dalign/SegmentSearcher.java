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

package com.milaboratory.migec.dalign;

import com.milaboratory.migec.segment.Allele;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class SegmentSearcher {
    private final int minHitSize;
    private final List<HitTracker> hitTrackers;

    public SegmentSearcher(List<Allele> alleles) {
        this(alleles, 2);
    }

    public SegmentSearcher(List<Allele> alleles, int minHitSize) {
        this.hitTrackers = new LinkedList<>();
        for (Allele allele : alleles) {
            hitTrackers.add(new HitTracker(allele, minHitSize));
        }
        this.minHitSize = minHitSize;
    }

    public SearchResult scan(String cdr3seq, int vEndRel, int jStartRel) {
        if (jStartRel - vEndRel - 1 < minHitSize)
            return null;

        String subSeq = cdr3seq.substring(vEndRel + 1, jStartRel);
        List<Allele> hitAlleles = new ArrayList<>();

        // iterate from largest window to smallest one
        for (int i = subSeq.length(); i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < subSeq.length() - i; j++) {
                String kmer = subSeq.substring(j, j + i);

                boolean hasHit = false;
                for (HitTracker hitTracker : hitTrackers) {
                    if (hitTracker.hasHit(kmer)) {
                        hasHit = true;
                        hitAlleles.add(hitTracker.getAllele());
                    }
                }

                if (hasHit) {
                    int from = vEndRel + 1 + j;
                    return new SearchResult(from,
                            from + i - 1,
                            i,
                            hitAlleles
                    );
                }
            }
        }

        return null;
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
