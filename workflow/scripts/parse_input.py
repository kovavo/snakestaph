'''
import inital sequence files into snakemake directories and name them accordingly.
concatenate and compress them if needed
input should be manifest file
output should be path in snakemake directory
'''

import glob
import os
import pandas as pd

import pandas as pd

def extract_platform_and_path(df):
    # Check if any of the platform columns are present in the DataFrame
    platform_columns = ['illumina_forward', 'Illumina_reverse', 'nanopore', 'iontorrent']
    if not any(col in df for col in platform_columns):
        # If no platform columns are present, raise an error
        raise ValueError("No platform information available")

    # Initialize an empty DataFrame to store the extracted platform and path data
    new_df = pd.DataFrame()

    # Iterate through the platform columns
    for platform in platform_columns:
        # Check if the current platform column exists in the DataFrame and if it has any non-null values
        if platform in df and df[platform].notna().any():
            # Filter the DataFrame to include only rows with non-null values in the current platform column
            non_none_data = df[df[platform].notna()]

            # Extract the 'reads' column and the current platform column
            platform_data = non_none_data[['reads', platform]].rename(columns={platform: 'path'}, inplace=False)

            # Insert a new column named 'platform' at index 1
            platform_data.insert(1, 'platform', platform)

            # Append the extracted platform and path data to the new DataFrame
            new_df = new_df.append(platform_data)

    # Return the extracted platform and path data
    return new_df

def extract_reads(df, reads_col='reads'):
    platform_col = df.melt(id_vars=reads_col).dropna()
    print(platform_col)


def recur_dictify(frame):
    '''
    No idea how it works but it does its job at casting input csv into hierarchical nested dictionary
    https://stackoverflow.com/questions/19798112/convert-pandas-dataframe-to-a-nested-dict
    '''
    if len(frame.columns) == 1:
        if frame.values.size == 1: return frame.values[0][0]
        return frame.values.squeeze()
    grouped = frame.groupby(frame.columns[0])
    d = {k: recur_dictify(g.iloc[:,1:]) for k,g in grouped}
    return d

def flatten(nested_list: list, type='list') -> list:
    '''
    Flatten nested list
    '''

    if not isinstance(nested_list, (list, tuple)):
        return [nested_list]
    else:
        flat_list = []
        for item in nested_list:
            items = flatten(item)
            flat_list = flat_list + items
        if type == 'list':
            return flat_list
        elif type == 'tuple':
            return tuple(flat_list)

def parse_csv(csv, colnames, mappings, colpath='path'):
    # load csv table with sequence files, unix style regexs are allowed - multiple files are globbed and inserted as list
    df = pd.read_csv(csv, dtype=str, names=colnames)
    df[colpath] = df[colpath].apply(lambda x: flatten([os.path.abspath(path) for path in glob.glob(x.strip())]))
    #rewrite column values or create new based on mappings
    for src_colname in mappings.keys():
        if src_colname in df.columns:
            for dis_colname, dictionary in mappings[src_colname].items():
                df[dis_colname] = df[src_colname].apply(lambda x: dictionary[x.strip()])

    return df.explode(colpath, ignore_index=True)

def collapse_df(df, groupby, exclude={}, to_list='path', explode=[], index="name"):

    #group datafame but escape excluded columns
    excluded = pd.DataFrame()
    explode = []
    for key in exclude:
        for value in set(df[key]):
            for new in df.columns:
                new_colname = f'{key}_{value}_{new}'
                excluded[new_colname] = df.apply(lambda x: x[new] if x[key] == value else None, axis=1)
                if new == exclude[key]:
                    explode.append(new_colname)

    df = pd.concat([df, excluded], axis=1)
    df = df.groupby(groupby, dropna=False).agg(lambda x: sorted(list(set(list(filter(None,x))))))

    #create index using groupby values
    df = df.reset_index()
    df.index = df[groupby].apply(lambda x: "_".join(x), axis=1)
    df['groupby'] = df[groupby].apply(lambda x: "_".join(x), axis=1)

    #explode
    for column in explode:
        df = df.explode(column)

    #now join all values into single string except path and replace all NaN to prevent it from messing with join()
    df = df.fillna('').apply(list)
    df[df.columns.drop(to_list)] = df.applymap(lambda x: "-".join([str(i) for i in flatten([x])]))[df.columns.drop(to_list)]

    return df

def expand_dataframe(pattern, df, **kwargs):
    '''
    Expand pattern filled from dataframe columns, if needed to pass wildcards, use double curly brackets
    :param pattern: string or dictionary pattern with dataframe column names in the curly braces
    :param df: pandas dataframe, consider filtering before passing
    :param kwargs: any other pattern not contained in dataframe, will override dataframe column name
    :return: list of strings or dictionary
    '''

    expanded = []
    for i, j in df.iterrows():
        locals().update(j)
        locals().update(kwargs)
        if isinstance(pattern, str):
            expanded.append(eval(f'f"{pattern}"'))
        elif isinstance(pattern, dict):
            expanded = {}
            for key, value in pattern.items():
                print(key,value)
                print(eval(f'f"{key}"'), eval(f'f"{value}"'))
                expanded[eval(f'f"{key}"')] = eval(f'f"{value}"')
    return expanded