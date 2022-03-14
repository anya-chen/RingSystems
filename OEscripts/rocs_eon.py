import os
import glob
import multiprocessing as mp
from datetime import datetime


def run_rocs(query,prefix):
    os.system('/data/shared/software/openeye2020/openeye/arch/redhat-RHEL8-x64/rocs/rocs -query '+query+' -dbase /home/ychen/projects/nps_ringsys/20210901_3dComparison/assign_charges/zinc_charges_assigned/zinc_rings_all.oeb.gz -prefix '+prefix+' -eon_input true -eon_input_size 500')

    
def run_eon(database,prefix):
    os.system('/data/shared/software/openeye2020/openeye/arch/redhat-RHEL8-x64/eon/eon -dbase '+database+' -charges existing -fixpka_query false -fixpka_dbase false -prefix '+prefix)


    
if __name__ == "__main__":
    start_time = datetime.now()
    
    queryFolder = '/data/local/ringsys/202111_3dComparison/assign_charges/queries_charges_assigned/'
    rocsFolder = '/data/local/ringsys/202111_3dComparison/rocs/'
    
    jobList = []
    for file in os.listdir(queryFolder):
        query = queryFolder+file
        prefix = rocsFolder + '_'.join(i for i in file.split('_')[0:2])  
        jobList.append((query,prefix))
    
    print('Length of joblist: {}'.format(len(jobList)))
    
    pool = mp.Pool(mp.cpu_count()-8)
    print('Start ROCS......')
    results = pool.starmap(run_rocs, jobList)
    pool.close()
    pool.join()
    
    
    eonFolder = '/data/local/ringsys/eon/'
    files = glob.glob(os.path.join(rocsFolder,"*eon_input*"))
    jobList = []
    for file in files:
        filename = os.path.basename(file)
        prefix = eonFolder + filename.split('.')[0]
        jobList.append((file,prefix))
    print('Length of joblist: {}'.format(len(jobList)))
    
    pool = mp.Pool(mp.cpu_count()-8)
    print('Start EON......')
    results = pool.starmap(run_eon, jobList)
    pool.close()
    pool.join()
    
    print('Finished in:')
    print(datetime.now()-start_time)
