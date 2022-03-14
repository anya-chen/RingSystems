import csv
'''
This script is used to get the refined coconut as NPs and also get the removed part as one file.
'''
source2remove = ["supernatural2","zincnp","ibs2019mar_nc","fooddb","gnps","conmednp","chembl_np",\
                 "lichendatabase","specsnp","bitterdb","exposome-explorer","Indofinechemical",\
                 "chemspidernp","swmd"]

coconut_sourceNP ='/home/ychen/projects/nps_ringsys/20210608_coconut_mongodb/coconut.sourceNP.csv'
refinedCoconut = '/home/ychen/projects/nps_ringsys/20210608_coconut_mongodb/refinedCoconut.sourceNP.csv'
removedCoconut = '/home/ychen/projects/nps_ringsys/20210608_coconut_mongodb/removedCoconut.sourceNP.csv'

count = 0
with open(refinedCoconut,'w') as of:
    with open(removedCoconut,'w') as of2:
        with open(coconut_sourceNP,'r') as f:
            moleculeCsv = csv.reader(f, delimiter=',')
            header = next(moleculeCsv)
            #print(header)
            of.write('counter\t')
            of.write( '\t'.join(map(str, header)))
            of.write('\n')
            of2.write('counter\t')
            of2.write( '\t'.join(map(str, header)))
            of2.write('\n')
            for row in moleculeCsv:
                count += 1
                if row[0] not in source2remove:
                    of.write(str(count)+'\t')
                    of.write('\t'.join(map(str,row)))
                    of.write('\n')
                else:
                    of2.write(str(count)+'\t')
                    of2.write('\t'.join(map(str,row)))
                    of2.write('\n')
        f.close
    of2.close
of.close() 
 
