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

# def read_config(config_loc):
#     config = configparser.ConfigParser()
#     config.read(config_loc)
#     return config


def get_abs_paths(loc, ext):
    return [os.path.abspath(os.path.join(loc, path)) for path in os.listdir(loc)
            if path.endswith(ext)]


def read_roi_files(abs_files):
    return {file[-29:-20]: (pd.read_csv(file, sep="\t", header=0, usecols=[0, 1, 2], names=["chr", "pos", "gene"]))
            for file in abs_files}


def read_all_files(abs_files):
    return {file[-29:-20]: (pd.read_csv(file, sep="\t")) for file in abs_files}


class Data:
    # cache the data using panel:
    @pn.cache
    def __init__(self, path):
        # read the data and set it in self directory's
        # all_paths = get_abs_paths(path, "ALL.csv")
        roi_paths = get_abs_paths(path, "ROIs.csv")
        # self.all_dict = read_all_files(all_paths)
        self.roi_dict = read_roi_files(roi_paths)


def main():
    # config = read_config("config.ini")
    info = Data("Data/analysis")
    [chr_s, pos_ss, pos_se] = setup_panel(info)
    panel_vis(info, chr_s, pos_ss, pos_se)


if __name__ == "__main__":
    main()
