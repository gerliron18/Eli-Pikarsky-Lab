import pandas as pd

def main():
    df = pd.read_csv('C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/RAW_DATA/EGR1_targets/EGR1_targets.csv', sep=",", names=['Ensambel', 'Gene_name'])
    df_mat = pd.read_csv('C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/RAW_DATA/GSE162694_raw_counts.csv', sep=",")

    ensambel_dics = dict(zip(df.Ensambel, df.Gene_name))
    for x in df_mat['ID_REF']:
        if x in ensambel_dics.keys():
            df_mat.loc[df_mat['ID_REF'] == x, 'ID_REF'] = ensambel_dics[x]
        else:
            continue

    df_mat.to_csv('C:/Users/User/Documents/Computitional Biology/Year5/Pikarsky_Lab/NASH_option1/RAW_DATA/EGR1_targets/GSE162694_raw_counts_EGR1.csv', index=False)


if __name__ == '__main__':
    main()

