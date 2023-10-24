import glob
import pandas as pd
from astropy.io import fits
from tqdm import tqdm

def load_data():
    data_path = "/home/gbtdata/"
    all_sessions = glob.glob(f"{data_path}*")
    t_start = 60091

    save_df_dict = {
        "DMJD": [],
        "RA": [],
        "DEC": [],
        "VELDEF": [],
        "VFRAME": [],
        "RVSYS": [],
        "LO1FREQ": [],
        "RESTFREQ": [],
        "SIDEBAND": [],
        "IFFREQ": [],
        "LOMULT": [],
        "LOOFFSET": [],
        "DMJD_VEL": [],
        "VEL": [],
        "VDOT": [],
        "VDOTDOT": [],
    }

    for a in tqdm(all_sessions):        
        lo1a_files = glob.glob(f"{a}/LO1A/*.fits")

        for l in lo1a_files:
            go_file = glob.glob(f"{a}/GO/*.fits")[0]
            ldata = fits.open(l)
            gdata = fits.open(go_file)
            dmjd = float(ldata[3].data['DMJD'][0])

            if dmjd >= t_start:
                save_df_dict["DMJD"].append(dmjd)
                save_df_dict["RA"].append(ldata[3].data['RA'])
                save_df_dict["DEC"].append(ldata[3].data['DEC'])
                save_df_dict["VELDEF"].append(gdata[0].header['VELDEF'])
                save_df_dict["VFRAME"].append(ldata[3].data['VFRAME'][0])
                save_df_dict["RVSYS"].append(ldata[3].data['RVSYS'][0])
                save_df_dict["LO1FREQ"].append(ldata[3].data['LO1FREQ'][0])
                save_df_dict["RESTFREQ"].append(ldata[0].header['RESTFRQ'])
                save_df_dict["SIDEBAND"].append(ldata[0].header['SIDEBAND'])
                save_df_dict["IFFREQ"].append(ldata[0].header['IFFREQ'])
                save_df_dict["LOMULT"].append(ldata[0].header['SIDEBAND'])
                save_df_dict["LOOFFSET"].append(ldata[0].header['LOOFFSET'])
                save_df_dict["DMJD_VEL"].append(ldata[2].data['DMJD'][0])
                save_df_dict["VEL"].append(ldata[2].data['VELOCITY'][0])
                save_df_dict["VDOT"].append(ldata[2].data['VDOT'][0])
                save_df_dict["VDOTDOT"].append(ldata[2].data['VDOTDOT'][0])
            
            ldata.close()
            gdata.close()
            del dmjd
            del ldata
            del gdata

    df = pd.DataFrame(save_df_dict)
    df.to_csv('/home/sandboxes/vcatlett/repos/github/problem-data-fixes/2023-ephemeris/tests/local-data/true_values.csv', index=False)

load_data()