#!/home/user/PycharmProjects/Data_Dashboards_SV/venv python3

import panel
import pandas
import configparser


def setup_panel():
    # panel
    panel.extension(design="material", sizing_mode="stretch_width")


class Data:
    # cache the data using panel:
    @panel.cache
    def __init__(self, config):
        # read the data and set it in self directory's
        self.all_frame = {
            pandas.read_csv(config["FILE"])
        }
        self.roi_frame = {
            pandas.read_csv(config["FILE"])
        }


def read_config(config_loc):
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config


def main():
    config = read_config("config.ini")
    setup_panel()
    Data(config)


if __name__ == "__main__":
    main()
