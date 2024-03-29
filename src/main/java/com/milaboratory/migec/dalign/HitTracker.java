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

import com.milaboratory.migec.Util;
import com.milaboratory.migec.segment.Allele;

import java.util.HashSet;
import java.util.Set;

public class HitTracker {
    private final Allele allele;
    private final Set<String> kmers = new HashSet<>();

    public HitTracker(Allele allele, int minHitSize) {
        this.allele = allele;
        for (int i = minHitSize; i < allele.getSeq().length(); i++) {
            for (int j = 0; j < allele.getSeq().length() - i; j++) {
                String kmer = allele.getSeq().substring(j, j + i);
                kmers.add(kmer);
                kmers.add(Util.revCompl(kmer));
            }
        }
    }

    public boolean hasHit(String queryKmer) {
        return kmers.contains(queryKmer);
    }

    public Allele getAllele() {
        return allele;
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
