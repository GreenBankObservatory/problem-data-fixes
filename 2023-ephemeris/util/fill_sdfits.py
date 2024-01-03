import subprocess, os, datetime
import pandas as pd
from glob import glob
from tqdm import tqdm

def csv2flag(fpath):
    banks = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    session_code = fpath.split("/")[-1]
    data = pd.read_csv(f"{fpath}/row_shifts.csv")
    data = data.sort_values(by="SCAN")

    for b in banks:
        savepath = f"{fpath}/SDFITS/{session_code}.raw.vegas/{session_code}.raw.vegas.{b}.flag"

        indx_file = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/{session_code}/SDFITS/{session_code}.raw.vegas/{session_code}.raw.vegas.{b}.index"

        if os.path.exists(indx_file):
            indx_data = pd.read_csv(indx_file, skiprows=10, delim_whitespace = True)
            save_dict = {
                "scannum": list(indx_data["SCAN"].values),
                "intnum": list(indx_data["INT"].values),
                "plnum": list(indx_data["PLNUM"].values),
                "ifnum": list(indx_data["IFNUM"].values),
                "fdnum": list(indx_data["FDNUM"].values),
                "spurs": ["0" for si in range(len(list(indx_data["SCAN"].values)))],
            }

            if os.path.exists(savepath):
                os.remove(savepath)

            with open(savepath, "w+") as flag_file:
                t_now = datetime.datetime.now()
                if t_now.hour < 10:
                    t_now_hour = f"0{t_now.hour}"
                else:
                    t_now_hour = f"{t_now.hour}"
                if t_now.minute < 10:
                    t_now_minute = f"0{t_now.minute}"
                else:
                    t_now_minute = f"{t_now.minute}"
                if t_now.second < 10:
                    t_now_sec = f"0{t_now.second}"
                else:
                    t_now_sec = f"{t_now.second}"
                    
                flag_file.write("[header]\n")
                flag_file.write(f"created = Wed Dec 20 {t_now_hour}:{t_now_minute}:{t_now_sec} 2023\n")
                flag_file.write("version = 1.0\n")
                flag_file.write("created_by = sdfits\n")
                flag_file.write("modified_by = 2023_ephemeris_fix\n")
                flag_file.write("[flags]\n")
                flag_file.write("#RECNUM,SCAN,INTNUM,PLNUM,IFNUM,FDNUM,BCHAN,ECHAN,IDSTRING\n")
                for i in range(len(save_dict["scannum"])):
                    scannum = save_dict["scannum"][i]
                    intnum = save_dict["intnum"][i]
                    plnum = save_dict["plnum"][i]
                    ifnum = save_dict["ifnum"][i]
                    fdnum = save_dict["fdnum"][i]
                    data_matched = data[data["SCAN"]==scannum]
                    spurs = []
                    if len(data_matched["spurs"]) > 0:
                        spurs_matched = data_matched["spurs"].values[-1]
                        for s in spurs_matched.split(","):
                            if s not in spurs:
                                spurs.append(str(int(int(s)-1)))
                    spurs = (",").join(spurs)
                    line2add = f"*|{scannum}|{intnum}|{plnum}|{ifnum}|{fdnum}|{spurs}|{spurs}|VEGAS_SPUR\n"
                    flag_file.write(line2add)

bash_script = "/home/sandboxes/vcatlett/repos/github/ephem/2023-ephemeris/util/fill_sdfits.bash"
dirs = glob("/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/AGBT23*")

for d in tqdm(dirs):
    session_code = d.split("/")[-1]

    if not os.path.exists(f"{d}/GO"):
        source_path = f"/home/gbtdata/{session_code}/GO"
        out = subprocess.run([
        "cp",
        "-r", 
        str(source_path),
        f"{d}/"
        ], 
        capture_output=True)
        print(out.stdout)
    
    if not os.path.exists(f"{d}/Antenna"):
        source_path = f"/home/gbtdata/{session_code}/Antenna"
        out = subprocess.run([
        "cp",
        "-r", 
        str(source_path),
        f"{d}/"
        ], 
        capture_output=True)
        print(out.stdout)

    out = subprocess.run([
        "bash",
        bash_script, 
        session_code,
        ], 
        capture_output=True)
    print(out.stdout)
    csv2flag(d)