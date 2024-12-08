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

class Panel:
    def __init__(self, data):
        """

        :param data:
        :return:
        """
        pn.extension(design="material", sizing_mode="stretch_width")

        self.gene_input = pn.widgets.TextInput(name="Gene Name:", placeholder="Type first letters here")

        gene_options = self.gene_input.param.watch(self.input_change, ['value'], onlychanged=False)


        self.gene_selector = pn.widgets.Select(name="Gene:", options=gene_options, value="")

        if len(gene_options) > 50:
            self.gene_selector.disabled = True
        else:
            self.gene_selector.disabled = False

        self.get_pro_region(data, )

        # self.all_pipeline =
        self.roi_pipeline = data.roi_idf[
            data.roi_idf["pos"] ==
            data.roi_idf["chr"] ==
        ]

    @staticmethod
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

        chrom = input_gene_info['Chromosome'].tolist()[0]

        return promoter, chrom

    @staticmethod
    def filter_on_pro(data, promoter):
        """

        :param data:
        :param promoter:
        :return:
        """
        return data.all_idf[(data.all_idf['pos'] >= promoter[0]) & (data.all_idf['pos'] <= promoter[1])]

    @staticmethod
    def filter_genes(data, value: str):
        """

        :param data:
        :param value:
        :return:
        """
        return [gene for gene in list(set(data.gene_loc["Symbol"])) if gene.startswith(value)]


    def input_change(self, value):
        try:
            pass
        except ValueError:
            pass


    def show_coords(self, idf, region):
        pass

    def visualise_panel(self):
        show_pipe = pn.bind(self.show_coords, self.roi_pipeline)
        row = pn.Row(self.gene_input, pn.Column(self.gene_selector))
        row.show()


# ------------------- CONFIG -------------------------

def read_config(config_loc):
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config





class Data:
    # cache the data using panel:
    @pn.cache
    def __init__(self, config):
        # read the data and set it in self directory's
        # all_paths = self.get_abs_paths(config['LOCATIONS']['analysis'], "ALL.csv")
        roi_paths = self.get_abs_paths(config['LOCATIONS']['analysis'], "ROIs.csv")
        # all_dict = self.read_all_files(all_paths)
        roi_dict = self.read_roi_files(roi_paths)
        # all_data = self.combine_df_dict_to_df(all_dict)
        roi_data = self.combine_df_dict_to_df(roi_dict)
        self.gene_loc = self.read_gene_location(config['LOCATIONS']['ncbi'])
        self.roi_idf = roi_data.interactive()
        # self.all_idf = self.all_data.interactive()

    # ------------------- FILE PROCESSING -----------------------

    @staticmethod
    def get_abs_paths(loc, ext):
        """

        :param loc:
        :param ext:
        :return:
        """
        return [os.path.abspath(os.path.join(loc, path)) for path in os.listdir(loc)
                if path.endswith(ext)]

    @staticmethod
    def read_roi_files(abs_files):
        """

        :param abs_files:
        :return:
        """
        return {file[-29:-20]: (pd.read_csv(file, sep="\t", header=0, usecols=[0, 1, 2], names=["chr", "pos", "gene"]))
                for file in abs_files}

    @staticmethod
    def read_all_files(abs_files):
        """

        :param abs_files:
        :return:
        """
        return {file[-29:-20]: (pd.read_csv(file, sep="\t")) for file in abs_files}

    @staticmethod
    def read_gene_location(path):
        """

        :param path
        :return:
        """
        gene_loc_file = pd.read_csv(path, sep="\t",
                                    usecols=["Begin", "End", "Chromosome", "Orientation", "Symbol", "Gene ID"])
        return gene_loc_file

    @staticmethod
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


def main():
    config = read_config("config.ini")
    info = Data(config)
    method = Panel(info)
    method.visualise_panel()


if __name__ == "__main__":
    main()
