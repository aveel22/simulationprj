import pandas as pd
from Tools.process_MOI import read_MOI, data_keys, get_data


if __name__ == "__main__":
    f2 = "../Data/5х0ml.txt"
    dir = "../Data/"
    filenames = ["5x1200ml.txt", "5х0ml.txt"]
    file_list = [dir + file for file in filenames]
    df = pd.DataFrame(columns=data_keys)
    # df = df.iloc[0:0]
    for ind, f in enumerate(file_list):
        ddf = read_MOI(f)
        if ind == 0:
            df = ddf
            continue
        df = pd.concat([df, ddf], ignore_index=True)

    print(df.head())

    print("Extract info")
    m, r, t = get_data(file_list)
    print("Mass data:")
    print(m)
    print("Center of mass data:")
    print(r)
    print("Inertia data:")
    print(t)

