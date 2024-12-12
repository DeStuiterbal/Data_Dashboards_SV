import pandas as pd
import os
import panel as pn

class Data:
    # cache the data using panel:
    @pn.cache()
    def __init__(self, config):
        # read the data and set it in self directory's
        # all_paths = Data.get_abs_paths(config['LOCATIONS']['analysis'], "ALL.csv")
        roi_paths = Data.get_abs_paths(config['LOCATIONS']['analysis'], "ROIs.csv")
        # all_dict = Data.read_all_files(all_paths)
        roi_dict = Data.read_roi_files(roi_paths)
        # self.all_data = Data.combine_df_dict_to_df(all_dict)
        self.roi_data = Data.combine_df_dict_to_df(roi_dict)
        self.gene_loc = Data.read_gene_location(config['LOCATIONS']['ncbi'])

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

    def __str__(self):
        return """
        
        """


if __name__ == "__main__":
    print(Data)
