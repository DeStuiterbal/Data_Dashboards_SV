#!/home/user/PycharmProjects/Data_Dashboards_SV/venv python3
import time
import configparser
import param
import pandas as pd
import panel as pn
from dataclass import Data
from viewing import View


# ------------------- CONFIG -------------------------


def read_config(config_loc):
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config


def main():
    config = read_config("config.ini")
    info = Data(config)
    View(data=info.roi_data, genes=info.gene_loc).show()


if __name__ == "__main__":
    main()
