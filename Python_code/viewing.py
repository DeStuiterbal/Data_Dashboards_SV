import panel as pn
from panel.viewable import Viewer
import param
import pandas as pd
import matplotlib


class View(Viewer):
    pn.extension(design="material", sizing_mode="stretch_width")

    data = param.DataFrame()
    genes = param.DataFrame()

    gene_input = param.String(default="")

    # define input for all the chromosomes
    chromosome = param.ListSelector(objects=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                             "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"])

    # set default to empty string because values are determined by other input values
    gene_name = param.ObjectSelector(objects=[""], precedence=-1)

    def __init__(self, **params):
        super().__init__(**params)
        self.param.objects = self.data
        self.param.objects = self.genes

        self.gene_selectables = list(set(self.genes["Symbol"]))

        # promoter, chrom = self.get_pro_region(self.gene_input.value)

        # self.all_pipeline = self.data.all_idf[
        #     (self.data.all_idf["pos"] in range(promoter[0], promoter[1])) &
        #     (self.data.all_idf["chr"] == self.Chromosome.value)
        # ]
        # self.roi_pipeline = self.data.roi_idf[
        #     (self.data.roi_idf["pos"] in range(promoter[0], promoter[1])) &
        #     (self.data.roi_idf["chr"] == self.Chromosome.value)
        # ]

    @param.depends("genes")
    def get_pro_region(self, gene):
        """

        :param gene:
        :return:
        """
        input_gene_info = self.genes[self.genes["Symbol"] == gene]

        promoter = [input_gene_info["Begin"].tolist()[0], input_gene_info["End"].tolist()[0]]

        if input_gene_info["Orientation"][0] == "plus":
            promoter[0] -= 1000
        else:
            promoter[1] += 1000

        chrom = input_gene_info['Chromosome'].tolist()[0]

        return promoter, chrom

    @param.depends("data")
    def filter_on_pro(self, promoter):
        """

        :param promoter:
        :return:
        """
        return self.data.all_idf[(self.data['pos'] >= promoter[0]) & (self.data['pos'] <= promoter[1])]

    @param.depends("genes", "gene_input")
    def filter_by_gene_input(self):
        """
        filter on list of genes based on typed input, keep strings that match beginning of input
        :return: list of strings that either match or begin with the input string
        """
        filtered = [gene for gene in self.gene_selectables if gene.startswith(self.gene_input)]
        if len(filtered) != 0:
            self.gene_selectables = filtered

    @param.depends("genes", "chromosome")
    def filter_by_chrom(self):
        if self.chromosome is None or self.chromosome == []:
            self.gene_selectables = list(set(self.genes["Symbol"]))
        else:
            self.gene_selectables = list(set(self.genes[self.genes["Chromosome"].isin(self.chromosome)]["Symbol"]))

    @param.depends("gene_name", "genes", "chromosome", "gene_input")
    def set_options(self):
        self.filter_by_chrom()
        self.filter_by_gene_input()
        if len(self.gene_selectables) > 1000:
            self.param.gene_name.objects = [""]
        else:
            self.param.gene_name.objects = self.gene_selectables

    @param.depends("data", "genes", "chromosome", "gene_name", "gene_input")
    def show_coords(self):
        df = self.data
        return df

    @param.depends("gene_input", "chromosome", "gene_name")
    def __panel__(self):
        self.set_options()
        return pn.Column(pn.Row(pn.widgets.TextInput.from_param(self.param.gene_input),
                                (pn.widgets.MultiChoice.from_param(self.param.chromosome))),
                         pn.widgets.Select.from_param(self.param.gene_name),
                         pn.widgets.Tabulator(self.show_coords, page_size=10, pagination="remote"))

    def __str__(self):
        return """
        Viewer class
        """


if __name__ == "__main__":
    print(View)
