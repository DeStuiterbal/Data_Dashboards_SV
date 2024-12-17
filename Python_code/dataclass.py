import pandas as pd
import os
import panel as pn


class Data:
    # cache the data using panel:
    @pn.cache()
    def __init__(self, config):
        """

        :param config: parsed configuration file with relative paths to the data.
        """
        # read the data and set it in self directory's
        all_paths = Data.get_abs_paths(config['LOCATIONS']['analysis'], "ALL.csv")
        all_dict = Data.read_all_files(all_paths)
        # self.all_data = Data.combine_df_dict_to_df(all_dict)
        self.all_data = all_dict["barcode13"]
        self.gene_loc = Data.read_gene_location(config['LOCATIONS']['ncbi'])

    @staticmethod
    def get_abs_paths(loc, ext):
        """

        :param loc: relative path to file directory
        :param ext: extension of desired files (ROI or ALL) + .***
        :return: absolute paths of files in the directory with the desired extension
        """
        return [os.path.abspath(os.path.join(loc, path)) for path in os.listdir(loc)
                if path.endswith(ext)]

    @staticmethod
    def read_roi_files(abs_files):
        """

        :param abs_files: absolute file paths of all the ROI files to be read
        :return: dictionary of pandas dataframes of the files with the barcodes as keys
        """
        # file indexing is to obtain te barcode of the filename, it is consistent with negative numbers
        return {file[-29:-20]: (pd.read_csv(file, sep="\t", header=0, usecols=[0, 1, 2], names=["chr", "pos", "gene"]))
                for file in abs_files}

    @staticmethod
    def read_all_files(abs_files):
        """

        :param abs_files: absolute file paths of all the ALL files to be read
        :return: dictionary of pandas dataframes of the files with the barcodes as keys
        """
        # file indexing is to obtain te barcode of the filename, it is consistent with negative numbers
        return {file[-28:-19]: (pd.read_csv(file, sep="\t")) for file in abs_files}

    @staticmethod
    def read_gene_location(path):
        """

        :param path: path to the ncbi file containing all the gene locations
        :return: pandas dataframe of gene locations from the path, containing the columns used in the code
        """
        gene_loc_file = pd.read_csv(path, sep="\t",
                                    usecols=["Begin", "End", "Chromosome", "Orientation", "Symbol", "Gene ID"])
        return gene_loc_file

    @staticmethod
    def combine_df_dict_to_df(df_dict):
        """

        :param df_dict: dictionary with pandas dataframes and as key the barcode ids
        :return: a pandas dataframe with an additional column for the barcode
        """
        # define empty dataframe to concat to
        full_data = None

        for barcode in df_dict:
            # create a copy to manipulate
            barcodex = df_dict[barcode].copy()

            # create a new column to get the patient barcode in te combined dataframe
            barcodex.insert(len(barcodex.columns), "barcode", barcode)

            # concat the data on top of the full_roi_data dataframe only if full_data already has data
            full_data = (barcodex.copy() if full_data is None
                         else pd.concat([full_data, barcodex], join="inner", ignore_index=True, sort=False))
        return full_data

    def __str__(self):
        return """
        data reading class
        """


if __name__ == "__main__":
    print(Data)
