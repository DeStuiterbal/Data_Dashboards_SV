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
    """
    Parses the config file
    :param config_loc: relative path to the config.ini file
    :return: parsed config file
    """
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config


def main():
    config = read_config("config.ini")
    info = Data(config)
    View(data=info.all_data, genes=info.gene_loc).show()


if __name__ == "__main__":
    main()
