import pandas as pd
import numpy as np
import datetime, os

session_code = "AGBT23A_288"
nums = ["01"]
#session_code = "AGBT22B_065"
#nums = ["43", "44", "45", "46", "47", "48"]
banks = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

for n in nums:
    fpath = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/{session_code}_{n}"
    data = pd.read_csv(f"{fpath}/row_shifts.csv")
    data = data.sort_values(by="SCAN")

    for b in banks:
        savepath = f"{fpath}/{session_code}_{n}.raw.vegas.{b}.flag"

        indx_file = f"/stor/scratch/vcatlett/problem-data-temp/2023-ephemeris/modified/{session_code}_{n}/SDFITS/{session_code}_{n}.raw.vegas/{session_code}_{n}.raw.vegas.{b}.index"

        if os.path.exists(indx_file):
            indx_data = pd.read_csv(indx_file, skiprows=10, delim_whitespace = True)
            #print(indx_data)
            save_dict = {
                "scannum": list(indx_data["SCAN"].values),
                "intnum": list(indx_data["INT"].values),
                "plnum": list(indx_data["PLNUM"].values),
                "ifnum": list(indx_data["IFNUM"].values),
                "fdnum": list(indx_data["FDNUM"].values),
                "spurs": ["0" for si in range(len(list(indx_data["SCAN"].values)))],
            }
            #print(save_dict)
            #breakpoint()

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
                    print(line2add)
                    #breakpoint()
                    flag_file.write(line2add)