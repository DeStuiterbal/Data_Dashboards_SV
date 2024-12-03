#!/home/user/PycharmProjects/Data_Dashboards_SV/venv python3

import panel as pn
import pandas as pd
import configparser
import numpy as np
import matplotlib
import os


def setup_panel(data):
    # panel
    pn.extension(design="material", sizing_mode="stretch_width")

    chr_select = pn.widgets.Select(name="Chr_N", )

    pos_slider_start = pn.widgets.IntSlider(name="Position_start", )
    pos_slider_end = pn.widgets.IntSlider(name="Position_end", )

    return chr_select, pos_slider_start, pos_slider_end


def panel_vis(data, chr_select, pos_slider_start, pos_slider_end):

    interactive_plot = pn.bind()

# No config is used at the moment


def read_config(config_loc):
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config


def get_abs_paths(loc, ext):
    return [os.path.abspath(os.path.join(loc, path)) for path in os.listdir(loc)
            if path.endswith(ext)]


def read_roi_files(abs_files):
    return {file[-29:-20]: (pd.read_csv(file, sep="\t", header=0, usecols=[0, 1, 2], names=["chr", "pos", "gene"]))
            for file in abs_files}


def read_all_files(abs_files):
    return {file[-29:-20]: (pd.read_csv(file, sep="\t")) for file in abs_files}


def combine_df_dict_to_df(df_dict):
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
    def __init__(self, path):
        # read the data and set it in self directory's
        # all_paths = get_abs_paths(path, "ALL.csv")
        roi_paths = get_abs_paths(path, "ROIs.csv")
        # all_dict = read_all_files(all_paths)
        roi_dict = read_roi_files(roi_paths)
        # self.all_data = combine_df_dict_to_df(all_dict)
        self.roi_data = combine_df_dict_to_df(roi_dict)


def main():
    # config = read_config("config.ini")
    info = Data("Data/analysis")
    # [chr_s, pos_ss, pos_se] = setup_panel(info)
    # panel_vis(info, chr_s, pos_ss, pos_se)
    print(info.roi_data)


if __name__ == "__main__":
    main()
