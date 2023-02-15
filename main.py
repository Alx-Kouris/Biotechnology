import os
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt
import seaborn as sns

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


def extract_mouse_protein(sentence):
    for word in sentence.split(" "):
        if "_MOUSE" not in word:
            continue
        return string_between_brackets(word)


def string_between_brackets(word):
    return word[word.find("[") + 1:word.find("]")]


def add_identifier_column(dataframe, column):
    dataframe['protein_identifier'] = dataframe[column].apply(lambda x: extract_mouse_protein(x))


def clean_hrms_descriptions(dataframe):
    dataframe['description'] = dataframe['description'].apply(lambda x: remove_starting_from(x, "OS="))


def number_of_entries(connection, table_name):
    result = connection.execute(f"SELECT count(*) FROM {table_name}")
    return result.fetchone()[0]


def brain_parts_of_table(connection, table_name):
    brain_parts = []
    query = f"SELECT DISTINCT brain_part FROM {table_name}"
    res = connection.execute(query)
    for part in res.fetchall():
        if len(part[0]) == 0:
            continue
        if part[0] not in brain_parts:
            brain_parts.append(part[0])

    return brain_parts


def common_proteins_of_table(connection, table_name, part=""):
    common_proteins = []
    query = f"SELECT DISTINCT protein_identifier FROM {table_name} GROUP BY protein_identifier HAVING " + (
        f"brain_part='{part}' AND" if len(part) > 0 else "") + " COUNT(*) > 1"

    res = connection.execute(query)
    for entry in res.fetchall():
        common_proteins.append(entry[0])

    return common_proteins


def unique_proteins_of_table(connection, table_name, part="", print_unique_info=False):
    unique_accessions = []
    query = f"SELECT protein_identifier FROM {table_name} GROUP BY protein_identifier HAVING " + (
        f"brain_part='{part}' AND" if len(part) > 0 else "") + " COUNT(*) == 1"

    res = connection.execute(query)

    for entry in res.fetchall():
        accession_id = entry[0]
        unique_accessions.append(accession_id)
        if print_unique_info:
            res2 = connection.execute(
                f"SELECT description, brain_part, coverage FROM {table_name} WHERE accession='{accession_id}'")
            for protein in res2.fetchall():
                print(f"Found unique protein {accession_id} ({protein[0]}) in {protein[1]} with coverage={protein[2]}")

    return unique_accessions


def annotate_dataframes(dataframes_to_annotate, gene_ontology_df, dataframes_accession_id, gene_ontology_accession_id,
                        add_protein_identifier=False):
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
            updated = df[dataframes_accession_id] == str(row[gene_ontology_accession_id]).split(" ")[0]
            df.loc[updated, 'gene_names'] = row['Gene Names']
            df.loc[updated, 'biological_process'] = row['Gene Ontology (biological process)']
            df.loc[updated, 'gene_ontology'] = row['Gene Ontology (GO)']
            df.loc[updated, 'molecular_function'] = row['Gene Ontology (molecular function)']
            df.loc[updated, 'cellular_component'] = row['Gene Ontology (cellular component)']
            if add_protein_identifier:
                df.loc[updated, 'protein_identifier'] = row['Entry Name']


def cleanup_column_names(dataframe):
    dataframe.columns = [
        name.replace('#', '').strip().replace('[', '').replace(']', '').replace(' ', '_').replace('.', '').lower()
        for name in dataframe.keys()]


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
        cleanup_column_names(dataframe=df)
        dataframes.append(df)

    return dataframes


###################################################################
#
# LOAD GENE ONTOLOGY
#
###################################################################

gene_ontology = pd.read_excel("uniprot-gene-ontology.xlsx")

# Open SQL connection
db_con = sqlite3.connect("biotechnology.db")

###################################################################
#
# HRMS
#
###################################################################

hrms_table_name = "hrms"

if not table_exists(connection=db_con, table_name=hrms_table_name):
    hrms = open_files_from_directory('./MOUSE BRAIN PROTEOME HRMS')
    for brain_df in hrms:
        add_identifier_column(dataframe=brain_df, column='description')
        clean_hrms_descriptions(brain_df)

    # Annotate with gene ontology before storing in database
    annotate_dataframes(dataframes_to_annotate=hrms,
                        gene_ontology_df=gene_ontology,
                        dataframes_accession_id="protein_identifier",
                        gene_ontology_accession_id="Entry Name")

    for brain_df in hrms:
        brain_df.to_sql(name=hrms_table_name, if_exists='append', con=db_con, index=False)
        db_con.commit()

print(f"Table {hrms_table_name} has {number_of_entries(connection=db_con, table_name=hrms_table_name)} entries")

unique_proteins_of_table(db_con, hrms_table_name)

###################################################################
#
# TARASLIA
#
###################################################################

match_parts = {
    "OB": "olfactory_balb",
    "ΟΒ": "olfactory_balb",
    "HT": "hipothalamus",
    "MD": "medulla",
    "MB": "mid_brain",
    "OD": "olfactory_balb",
    "HC": "hipocampus",
    "CB": "cerebellum",
    "CC": "cortex"
}


def normalize_taraslia_dataframe(taraslia_df):
    taraslia_df.rename(columns={"accession_name": "protein_identifier",
                                "protein_name": "description",
                                "protein_mw": "mw_kda",
                                "pi-value": "calc_pi"}, inplace=True)
    taraslia_df['brain_part'] = taraslia_df['brain_part'].apply(lambda x: x.strip().upper().split(", "))
    taraslia_df = taraslia_df.explode('brain_part', ignore_index=True)
    taraslia_df['brain_part'] = taraslia_df['brain_part'].apply(lambda x: match_parts[x.strip()])
    return taraslia_df


taraslia_table_name = "taraslia"

if not table_exists(connection=db_con, table_name=taraslia_table_name):
    # Open only the first file
    taraslia = pd.read_excel("./MOUSE BRAIN 2DGE PROTEINS/TARASLIA et al TABLE 1.xls", skiprows=[0])

    # Cleanup
    cleanup_column_names(taraslia)
    taraslia = normalize_taraslia_dataframe(taraslia)

    # Annotate with gene ontology before storing in database
    annotate_dataframes(dataframes_to_annotate=[taraslia],
                        gene_ontology_df=gene_ontology,
                        dataframes_accession_id="protein_identifier",
                        gene_ontology_accession_id="Entry Name")

    taraslia.to_sql(name=taraslia_table_name, if_exists='append', con=db_con, index=False)

print(f"Table {taraslia_table_name} has {number_of_entries(connection=db_con, table_name=taraslia_table_name)} entries")

unique_proteins_of_table(connection=db_con, table_name=taraslia_table_name)

###################################################################
#
# JUNG
#
###################################################################

jung_match_parts = {
    "Olfactory": "olfactory_balb",
    "STR": "striatum",
    "CTX1": "cortex",
    "CTX2": "cortex",
    "TH1": "thalamus",
    "TH2": "thalamus",
    "PA": "cortex",
    "HY1": "hipothalamus",
    "HY2": "hipothalamus",
    "CXS": "cortex",
    "HP": "hipocampus",
    "MB": "mid_brain",
    "PO": "pons",
    "RHP": "hipocampus",
    "CB1": "cerebellum",
    "CB2": "cerebellum",
    "MY": "medulla"
}

jung_table_name = "jung"

if not table_exists(connection=db_con, table_name=jung_table_name):
    jung_df = pd.DataFrame()
    for key in jung_match_parts.keys():
        print(f"Reading sheet {key} from Jung")
        jung = pd.read_excel("./JUNG/jung.xlsx",
                             key, skiprows=[0, 1])
        jung = jung.drop_duplicates(subset='Gene ID')  # keep a gene only once from each sample
        jung = jung.drop(columns=['Sample', 'Recovered Sequence ( under score delimited)'])
        jung = jung.reset_index()

        jung['brain_part'] = jung_match_parts[key]
        jung_df = pd.concat([jung_df, jung])

    jung_df.reset_index()
    jung_df.drop(columns=['index'], inplace=True)
    jung_df = jung_df.drop_duplicates(subset=['Gene ID', 'brain_part'])
    print(jung_df['brain_part'].value_counts())
    cleanup_column_names(jung_df)
    jung_df.rename(columns={"full_name": "description"}, inplace=True)

    # Annotate with gene ontology before storing in database
    annotate_dataframes(dataframes_to_annotate=[jung_df],
                        gene_ontology_df=gene_ontology,
                        dataframes_accession_id="symbol",
                        gene_ontology_accession_id="Gene Names",
                        add_protein_identifier=True)

    jung_df.to_sql(name=jung_table_name, if_exists='append', con=db_con, index=False)

print(f"Table {jung_table_name} has {number_of_entries(connection=db_con, table_name=jung_table_name)} entries")

unique_proteins_of_table(connection=db_con, table_name=jung_table_name)

###################################################################
#
# SHARMA
#
###################################################################

sharma_table_name = "sharma"


def sharma_protein_identifiers(dataframe_to_annotate, gene_ontology_df):
    dataframe_to_annotate['protein_identifier'] = ""
    for index, row in dataframe_to_annotate.iterrows():
        for protein_id in row['majority_protein_ids'].split(";"):
            if "-" in protein_id:
                protein_id = remove_starting_from(protein_id, "-")

            gene_entry = gene_ontology_df.loc[gene_ontology_df['Entry'] == protein_id]
            if len(gene_entry) > 0:
                dataframe_to_annotate.loc[index, 'protein_identifier'] = gene_entry['Entry Name'].values[0]


def determine_brain_part(dataframe):
    dataframe['brain_part'] = ""
    for index, row in dataframe.iterrows():
        for column_name in dataframe.columns.values:
            if row[column_name] == "+":
                brain_part = "cortex" if column_name == "motor_cortex" or column_name == "prefrontal_cortex" else column_name
                dataframe.loc[index, 'brain_part'] += brain_part + ","


if not table_exists(connection=db_con, table_name=sharma_table_name):
    sharma = pd.read_excel("./SHARMA/sharma.xlsx", skiprows=[0])
    cleanup_column_names(sharma)
    sharma.rename(columns={"protein_names": "description",
                           "sequence_coverage_%": "coverage",
                           "mol_weight_kda": "mw_kda",
                           "olfactory_bulb": "olfactory_balb",
                           "hippocampus": "hipocampus"}, inplace=True)
    determine_brain_part(sharma)
    sharma['brain_part'] = sharma['brain_part'].apply(lambda x: x[:-1].split(",") if x.endswith(",") else x.split(","))
    sharma = sharma.explode('brain_part', ignore_index=True)
    sharma_protein_identifiers(sharma, gene_ontology)

    # Annotate with gene ontology before storing in database
    annotate_dataframes(dataframes_to_annotate=[sharma],
                        gene_ontology_df=gene_ontology,
                        dataframes_accession_id="protein_identifier",
                        gene_ontology_accession_id="Entry Name",
                        add_protein_identifier=True)

    sharma.drop(columns=['brainstem', 'cerebellum',
                         'corpus_callosum', 'motor_cortex', 'olfactory_balb', 'optic_nerve',
                         'prefrontal_cortex', 'striatum', 'thalamus', 'hipocampus'], inplace=True)
    sharma.to_sql(name=sharma_table_name, if_exists='append', con=db_con, index=False)

print(f"Table {sharma_table_name} has {number_of_entries(connection=db_con, table_name=sharma_table_name)} entries")

unique_proteins_of_table(connection=db_con, table_name=sharma_table_name)

###################################################################
#
# PLOTS
#
###################################################################
studies = [hrms_table_name, taraslia_table_name, jung_table_name, sharma_table_name]

all_reported_brain_parts = []
for study in studies:
    reported_brain_parts = brain_parts_of_table(db_con, study)
    for bp in reported_brain_parts:
        if len(bp[0]) == 0:
            continue
        if bp not in all_reported_brain_parts:
            all_reported_brain_parts.append(bp)


# Unique|Common proteins per brain_part per study
for study in studies:
    unique_proteins_per_part = []
    common_proteins_per_part = []

    for brain_part in all_reported_brain_parts:
        unique_proteins_per_part.append(len(unique_proteins_of_table(connection=db_con, table_name=study, part=brain_part)))
        common_proteins_per_part.append(len(common_proteins_of_table(connection=db_con, table_name=study, part=brain_part)))

    plt.bar(all_reported_brain_parts, unique_proteins_per_part, color='r')
    plt.bar(all_reported_brain_parts, common_proteins_per_part, bottom=unique_proteins_per_part, color='b')
    plt.xlabel("Brain Parts")
    plt.ylabel("# Proteins")
    plt.legend(["Unique", "Common"])
    plt.xticks(rotation='vertical')
    plt.title(study.upper())
    plt.show()

# Number of common proteins with previous studies
common_proteins_with_studies = []
for study in studies:
    if study == hrms_table_name:
        continue

    query = f"SELECT COUNT(protein_identifier) FROM {hrms_table_name} WHERE protein_identifier in (SELECT protein_identifier FROM {study})"
    res = db_con.execute(query)
    common_proteins_with_studies.append(res.fetchone()[0])

plt.bar(studies[1:], common_proteins_with_studies, color='maroon')
plt.xlabel("Previous studies")
plt.ylabel("# Common Proteins")
plt.title("Common Proteins with previous studies")
plt.show()

db_con.close()
