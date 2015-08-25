" " " this file was created in november 2014
as part of a de novo search for EPIC loci in 
the chaetognath species Pterosagitta draco
property of dr. Ferdinand Marlétaz
" " "

#!/usr/bin/env python
import sys
from collections import defaultdict
flanks=defaultdict(list)
for line in open(sys.argv[1]):
    bline=line.rstrip().split('\t')
    name=bline[0].rsplit('_',1)[0]
    side=bline[0].rsplit('_',1)[1]
    #start,end=sort(int(bline[8]),int(bline[9]))
    flanks[name].extend([int(bline[8]),int(bline[9])])

for intr,flk in flanks.items():
    if len(flk)==4:
        spos=sorted(flk)
        #print flk,'->',spos
        print intr,abs(spos[2]-spos[1])
    #else:
    #    print intr,','.join(map(str,flk))
    
