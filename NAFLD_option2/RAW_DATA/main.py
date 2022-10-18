import re

import pandas as pd

AFTER_SURGERY = ['GSM1178997', 'GSM1179000', 'GSM1179009', 'GSM1179011',
                 'GSM1179013', 'GSM1179014', 'GSM1179019', 'GSM1179020',
                 'GSM1179022', 'GSM1179023', 'GSM1179028', 'GSM1179029',
                 'GSM1179031', 'GSM1179032', 'GSM1179034', 'GSM1179039',
                 'GSM1179040', 'GSM1179041', 'GSM1179042']

def main():
    df = pd.read_csv('GPL11532-32230.txt', sep="\t", names=['ID', 'Ensambel'])
    df_mat = pd.read_csv('GSE48452_series_matrix.txt', sep="\t")

    # df['Ensambel'] = df['Ensambel'].str.replace('f.', 'ENST[^0-9]+', regex=True)
    # df.Ensambel.replace(value="ENST[0-9]+", regex=True)
    # df['Ensambel'].replace(to_replace='[nN]ENST[nN]', value='ENST[^0-9]+', regex=True)
    # df.Ensambel.str.extract('(.*(ENST[^0-9]+)).*', re.IGNORECASE, expand=False)

    # pat = re.compile('ENST[^0-9]+', re.IGNORECASE)
    #
    # df.assign(Ensambel=[re.split(pat, x, 1)[0] for x in df.Ensambel])

    df['Ensambel'] = df['Ensambel'].replace('^.*(ENSG[0-9]{11}).*$', r'\1', regex=True)

    # ensambel_dics = df.set_index('Ensambel').to_dict()
    # df_mat['ID_REF'] = df_mat['ID_REF'].apply(lambda x: ensambel_dics[x])

    ensambel_dics = dict(zip(df.ID, df.Ensambel))
    for x in df_mat['ID_REF']:
        if x in ensambel_dics.keys():
            df_mat.loc[df_mat['ID_REF'] == x, 'ID_REF'] = ensambel_dics[x]
        else:
            continue

    # df_mat['ID_REF'] = df_mat['ID_REF'].replace(ensambel_dics, inplace=True)

    # for k in ensambel_dics:
    #     try:
    #         df_mat['ID_REF'] = df_mat['ID_REF'].apply(ensambel_dics[k])
    #     except KeyError:
    #         continue

    df_mat.drop(columns=AFTER_SURGERY, axis=1, inplace=True)

    df_mat.to_csv('arranged_before_surgery_genes.csv', index=False)


if __name__ == '__main__':
    main()

