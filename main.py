import os
import pandas as pd
import sqlite3

directories = ['./MOUSE BRAIN PROTEOME HRMS']

prefix = "mouse_brain_-_"
suffix = "_7_weeks"


def remove_starting_from(word, from_str):
    pos = word.find(from_str)
    result = ""
    for index in range(0, len(word) + 1):
        if index == pos:
            break
        result += word[index]

    return result


def excel_files_in_directory(directory):
    filenames = []
    for filename in os.listdir(directory):
        filename_lower = filename.lower()
        if filename_lower.startswith("supplementary") \
                or filename_lower.startswith(".") \
                or not filename_lower.endswith("xlsx"):
            continue

        filenames.append(directory + "/" + filename)

    return filenames


def column_names_as_single_string(dataframe):
    result = ""
    for column in dataframe.columns:
        result += column
        if column != dataframe.columns[-1]:
            result += ","

    return result


def table_exists(connection, table_name):
    result = connection.execute(f"SELECT count(*) FROM sqlite_master WHERE type='table' AND name='{table_name}'")
    return result.fetchone()[0] > 0


def clean_hrms_descriptions(dataframe):
    dataframe['description'] = dataframe['description'].apply(lambda x: remove_starting_from(x, "OS="))


def number_of_entries(connection, table_name):
    result = connection.execute(f"SELECT count(*) FROM {table_name}")
    return result.fetchone()[0]


def open_files_from_directory(directory):
    """
    Opens any excel files in the directory and returns DataFrames for each file
    :param directory: string
    :return: list of dataframes per excel file
    """
    filenames = excel_files_in_directory(directory)

    dataframes = []
    for file in filenames:
        brain_part = remove_starting_from(file.lower().removeprefix(directory.lower() + "/" + prefix), suffix)
        df = pd.read_excel(file)
        df['brain_part'] = brain_part
        df.columns = [
            name.replace('#', '').strip().replace('[', '').replace(']', '').replace(' ', '_').replace('.', '').lower()
            for name in df.keys()]
        dataframes.append(df)

    return dataframes


hrms = open_files_from_directory(directories[0])

# Clean description column
for brain_df in hrms:
    clean_hrms_descriptions(brain_df)

db_con = sqlite3.connect("biotechnology.db")

if not table_exists(connection=db_con, table_name="hrms"):
    for brain_df in hrms:
        brain_df.to_sql(name='hrms', if_exists='append', con=db_con, index=False)
        db_con.commit()

print(f"Table HRMS has {number_of_entries(connection=db_con, table_name='hrms')} entries")

unique_accessions = []
res = db_con.execute("SELECT DISTINCT accession FROM hrms")

for entry in res.fetchall():
    accession_id = entry[0]
    unique_accessions.append(accession_id)
    res2 = db_con.execute(f"SELECT description, brain_part, coverage FROM hrms WHERE accession='{accession_id}'")
    for protein in res2.fetchall():
        print(f"Found unique protein {accession_id} ({protein[0]}) in {protein[1]} with coverage={protein[2]}")

print(f"We have {len(unique_accessions)} unique accessions")

db_con.close()
