" " " this file was created in november 2014
as part of a de novo search for EPIC loci in 
the chaetognath species Pterosagitta draco
property of dr. Ferdinand Marlétaz
" " "

#!/usr/bin/env python
import sys
import re
from collections import defaultdict



def reverse(ali,taxa,clust):
    alen=len(ali[taxa[0]])
    #print alen
    positions=['']*alen
    for tax in taxa:
        seq=ali[tax]
        #print tax, len(seq)
        for i,res in enumerate(seq):
            positions[i]+=res
    #print taxa
    n_int=0
    onset=50
    maxgaps=20
    #We selection introns with at least 30 flanking positions out of 50
    #print ','.join(taxa)
    msk=[tax for tax in taxa if tax.startswith('Lgi')]
    id_msk=''.join(msk[0]) if len(msk)>0 else 'NA'
    for i,pos in enumerate(positions):
        if ''.join(set(pos))=='(':
            #print '('*len(taxa)
            items=dict((e,positions[i+1].count(e)) for e in set(positions[i+1]))
            sum_pres=sum([items[v] for v in ['0','1','2'] if v in items])
            sum_tot=sum(items.values())
            if sum_pres>sum_tot-5:
                cons_left=ali['cons'][i-onset:i]
                cons_right=ali['cons'][i+3:i+onset+3]
                cons_left_sc=cons_left.count('+')
                cons_right_sc=cons_right.count('+')                
                seq_left=ali[id_msk][i-onset:i].replace(']',')').split(')')[-1] if ')' in ali[id_msk][i-onset:i].replace(']',')') else ali[id_msk][i-onset:i]
                
                seq_right=ali[id_msk][i+3:i+onset+3].replace('[','(').split('(')[0]
                gap_left=cons_left.count('-')
                gap_right=cons_right.count('-')
                if len(seq_left.replace('-',''))>=onset-maxgaps and len(seq_right.replace('-',''))>=onset-maxgaps:
                    if gap_left<maxgaps and gap_right<maxgaps:
                        print "{0}\t{1}\t{2}/{3}\t{4}/{5}\t{6} / {7}\t{8} / {9}".format(clust,i,sum_pres,sum_tot,cons_left_sc,cons_right_sc,cons_left,cons_right,seq_left,seq_right)
                        out.write('>{0}_{1}_left\n{2}\n'.format(clust,i,seq_left))
                        out.write('>{0}_{1}_right\n{2}\n'.format(clust,i,seq_right))
                #print '\n'.join(positions[i-10:i+11])
            n_int+=1
            
        #print i,pos[0:50]
    #print n_int
    


gene=''
aliSet=defaultdict(str)
taxa=[]
elt=['cons','insert','sites','intron','sfilt']
out=open('flanks.fa','w')
for line in open(sys.argv[1]):
    if line.startswith('./METFAM'):
        clust=line.rstrip().split('/')[3]
        gene=''
        if len(aliSet)>0:
            #print "\n",clust
            #for tax in taxa:
                #print "{0}\t{1}".format(tax,len(aliSet[tax]))
            rev=reverse(aliSet,taxa,clust.split('.')[0])
        
        
        aliSet=defaultdict(str)
        taxa=[]
       
    alin=re.search(r'([^\s]+)(\s+)(.+)\n', line)
    if alin:
        name=line[0:40].split()[0]
        seq=line.rstrip()[40:]
        #print name,seq[0:40]
        aliSet[name]+=seq
        if not name in taxa and not name in elt:
            taxa.append(name)
        #print alin.groups()



#Bfl-G461809                             
#cons                                    
        
