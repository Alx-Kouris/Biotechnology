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


def unique_proteins_of_table(connection, table_name, print_unique_info):
    unique_accessions = []
    res = connection.execute(f"SELECT DISTINCT accession FROM {table_name}")

    for entry in res.fetchall():
        accession_id = entry[0]
        unique_accessions.append(accession_id)
        if print_unique_info:
            res2 = connection.execute(
                f"SELECT description, brain_part, coverage FROM {table_name} WHERE accession='{accession_id}'")
            for protein in res2.fetchall():
                print(f"Found unique protein {accession_id} ({protein[0]}) in {protein[1]} with coverage={protein[2]}")

    print(f"We have {len(unique_accessions)} unique accessions")
    return unique_accessions


def annotate_dataframes(dataframes_to_annotate, gene_ontology_df, dataframes_accession_id, gene_ontology_accession_id):
    """
    Annotate dataframes with extra columns found in gene_ontology based on the column name for matching between dataframes
    :param dataframes_to_annotate: list of dataframes to annotate
    :param gene_ontology_df: gene ontology dataframe
    :param dataframes_accession_id: column name in dataframes that will try to match in gene_ontology dataframe
    :param gene_ontology_accession_id: column name in gene ontology that contains the accession id
    :return: list of annotated dataframes
    """
    for df in dataframes_to_annotate:
        for _, row in gene_ontology_df.iterrows():
            updated = df[dataframes_accession_id] == row[gene_ontology_accession_id]
            df.loc[updated, 'gene_names'] = row['Gene Names']
            df.loc[updated, 'biological_process'] = row['Gene Ontology (biological process)']
            df.loc[updated, 'gene_ontology'] = row['Gene Ontology (GO)']
            df.loc[updated, 'molecular_function'] = row['Gene Ontology (molecular function)']
            df.loc[updated, 'cellular_component'] = row['Gene Ontology (cellular component)']


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

gene_ontology = pd.read_excel("uniprot-gene-ontology.xlsx")

annotate_dataframes(dataframes_to_annotate=hrms,
                    gene_ontology_df=gene_ontology,
                    dataframes_accession_id="accession",
                    gene_ontology_accession_id="Entry")

if not table_exists(connection=db_con, table_name="hrms"):
    for brain_df in hrms:
        brain_df.to_sql(name='hrms', if_exists='append', con=db_con, index=False)
        db_con.commit()

print(f"Table HRMS has {number_of_entries(connection=db_con, table_name='hrms')} entries")

unique_proteins_of_table(db_con, "hrms", False)

db_con.close()
