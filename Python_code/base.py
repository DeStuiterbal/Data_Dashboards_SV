#!/home/user/PycharmProjects/Data_Dashboards_SV/venv python3

import panel
import pandas
import configparser


def read_config(config_loc):
    config = configparser.ConfigParser()
    config.read(config_loc)
    return config


def main():
    config = read_config('config.ini')


if __name__ == "__main__":
    main()
