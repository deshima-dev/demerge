from merge_to_dfits import MergeToDfits
import os
import glob
import multiprocessing as mp

#cosmos_dirs = ["/Users/deshima/hirado/ASTE2017/data/ASTE2017/LT119_FB3.3G_49ch/cosmos_20171112044900"]
cosmos_dirs = glob.glob("/Users/deshima/hirado/ASTE2017/data/ASTE2017/LT119_FB3.3G_49ch/cosmos_*")
ddb_fits = "/Users/deshima/workspace/takekoshi/ddb/toptica_20180612/DDB_20180619.fits.gz"
dfitsdict = "/Users/deshima/work/analysis/deshima/aste/dfits_dict.yaml"
cabinlog = "/Users/junyasuzuki/data/skydip/cabin.db" 

#outdir = "/Users/deshima/workspace/takekoshi/dfits/20180615/dfits_20180619"
outdir = "/Users/deshima/workspace/takekoshi/dfits/20180615/dfits_20180704"


def mp_wrapper(cosmos_dir):
    id = cosmos_dir[-14:]
    obstinst = os.path.join(cosmos_dir,"{}.obs".format(id))
    antennalog = os.path.join(cosmos_dir,"{}.ant".format(id))
    weatherlog = os.path.join(cosmos_dir,"{}.wea".format(id))
    rout_data = os.path.join(cosmos_dir,"reduced_{}.fits".format(id))
    print(obstinst,antennalog,weatherlog,rout_data)
    try:
        mtd = MergeToDfits(ddb_fits, dfitsdict, obstinst,
                           antennalog, rout_data, weatherlog, cabinlog)
        dfits_hdus = mtd.dfits
        dfits_hdus.writeto(os.path.join(outdir,"dfits_{}.fits.gz".format(id)))
        print("Suceeded. ID: {}".format(id))
    except:
        print("Failed. ID: {}".format(id))


#for cosmos_dir in cosmos_dirs:
ncores = mp.cpu_count()-2
p = mp.Pool(ncores)
p.map(mp_wrapper,cosmos_dirs)
p.close()

print("Processing finished")

