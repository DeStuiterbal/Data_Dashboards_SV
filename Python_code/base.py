#!/home/user/PycharmProjects/Data_Dashboards_SV/venv python3
import time

import panel as pn
import pandas as pd
import hvplot.pandas
import configparser
import numpy as np
import matplotlib
import os


# ------------------- PANEL -----------------------

def setup_panel(data, gene_options=0):
    """

    :param gene_options:
    :param data:
    :return:
    """
    #
    pn.extension(design="material", sizing_mode="stretch_width")

    gene_input = pn.widgets.TextInput(name="Gene Name:", placeholder="Type first letters here")

    gene_options = list(set(data.gene_loc["Symbol"]))

    gene_selector = pn.widgets.Select(name="Gene:", options=gene_options)

    return gene_input, gene_selector


def get_pro_region(data, gene):
    """

    :param data:
    :param gene:
    :return:
    """
    input_gene_info = data.gene_loc[data.gene_loc["Symbol"] == gene]

    promoter = [input_gene_info["Begin"].tolist()[0], input_gene_info["End"].tolist()[0]]

    if input_gene_info["Orientation"][0] == "plus":
        promoter[0] -= 1000
    else:
        promoter[1] += 1000

    return promoter


def visualise_panel(data, gene_input, gene_selector):
    row = pn.Row(gene_input)
    row.show()


# ------------------- CONFIG -------------------------

def read_config(config_loc):
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config


# ------------------- FILE PROCESSING -----------------------

def get_abs_paths(loc, ext):
    """

    :param loc:
    :param ext:
    :return:
    """
    return [os.path.abspath(os.path.join(loc, path)) for path in os.listdir(loc)
            if path.endswith(ext)]


def read_roi_files(abs_files):
    """

    :param abs_files:
    :return:
    """
    return {file[-29:-20]: (pd.read_csv(file, sep="\t", header=0, usecols=[0, 1, 2], names=["chr", "pos", "gene"]))
            for file in abs_files}


def read_all_files(abs_files):
    """

    :param abs_files:
    :return:
    """
    return {file[-29:-20]: (pd.read_csv(file, sep="\t")) for file in abs_files}


def read_gene_location(path):
    """

    :param path
    :return:
    """
    gene_loc_file = pd.read_csv(path, sep="\t",
                                usecols=["Begin", "End", "Chromosome", "Orientation", "Symbol", "Gene ID"])
    return gene_loc_file


def combine_df_dict_to_df(df_dict):
    """

    :param df_dict:
    :return:
    """
    # define empty dataframe to concat to
    full_roi_data = pd.DataFrame(columns=["chr", "pos", "gene", "barcode"])

    for barcode in df_dict:

        # create a copy to manipulate
        barcodex = df_dict[barcode].copy()

        # create a new column to get the patient barcode in te combined dataframe
        barcodex.insert(len(barcodex.columns), "barcode", barcode)

        # concat the data on top of the full_roi_data dataframe
        full_roi_data = pd.concat([full_roi_data, barcodex], join="inner", ignore_index=True, sort=False)

    return full_roi_data


class Data:
    # cache the data using panel:
    @pn.cache
    def __init__(self, config):
        # read the data and set it in self directory's
        # all_paths = get_abs_paths(config['LOCATIONS']['analysis'], "ALL.csv")
        roi_paths = get_abs_paths(config['LOCATIONS']['analysis'], "ROIs.csv")
        # all_dict = read_all_files(all_paths)
        roi_dict = read_roi_files(roi_paths)
        # all_data = combine_df_dict_to_df(all_dict)
        roi_data = combine_df_dict_to_df(roi_dict)
        self.gene_loc = read_gene_location(config['LOCATIONS']['ncbi'])
        self.roi_idf = roi_data.interactive()
        # self.all_idf = all_data.interactive()


def main():
    config = read_config("config.ini")
    info = Data(config)
    [v1, v2] = setup_panel(data=info)
    visualise_panel(info, v1, v2)


if __name__ == "__main__":
    main()
