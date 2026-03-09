import re
import pandas as pd
import datetime
import sys

data_keys = ['Mass', 'X', 'Y', 'Z', 'Lxx', 'Lyy', 'Lzz', 'Lxy', 'Lyz', 'Lxz']

def read_MOI(inputfile, save_on=False):
    try:
        # inputfile = "MOI 061222.txt"
        with open(inputfile, "r") as f:
            text = f.read()
    except FileNotFoundError as e:
        print(f"{e}. Wrong input filename")
        sys.exit(1)

    data_keys = ['Mass', 'X', 'Y', 'Z', 'Lxx', 'Lyy', 'Lzz', 'Lxy', 'Lyz', 'Lxz']
    data = {key: [float(v) for v in re.findall(f'{key} = (-?\w+.\w+)', text)] for key in data_keys}
    df = pd.DataFrame.from_dict(data)
    if save_on:
        today = datetime.date.today().strftime("%Y-%m-%d")
        outputfile = f'MOI_data_{today}.csv'
        df.to_csv(outputfile, sep=';', decimal='.', header=True, index=True)

    return df


def read_MOI_data(filenames):
    df = pd.DataFrame(columns=data_keys)
    # df = df.iloc[0:0]
    for ind, f in enumerate(filenames):
        ddf = read_MOI(f)
        if ind == 0:
            df = ddf
            continue
        df = pd.concat([df, ddf], ignore_index=True)

    return df

def get_data(filenames):
    df = read_MOI_data(filenames)
    mass = df["Mass"].values
    rcgs = df[["X", "Y", "Z"]].values
    tensors = df[['Lxx', 'Lyy', 'Lzz', 'Lxy', 'Lyz', 'Lxz']].values
    return mass, rcgs, tensors