import matplotlib.pyplot as plt
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
    gene_name = param.ObjectSelector(objects=[None], allow_None=True)

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

    @param.depends("gene_name", "genes")
    def get_pro_region(self):
        """

        :return: promoter and chromosome location of the given gene
        """
        input_gene_info = self.genes[self.genes["Symbol"] == self.gene_name]

        # get the beginning and end location of the gene
        promoter = [input_gene_info["Begin"].tolist()[0], input_gene_info["End"].tolist()[0]]

        # add or subtract from the promoter depending on the orientation of the gene
        if input_gene_info["Orientation"].tolist()[0] == "plus":
            promoter[0] -= 1000
        else:
            promoter[1] += 1000

        # get the chromosome of the gene
        chrom = input_gene_info['Chromosome'].tolist()[0]

        return promoter, chrom

    @param.depends("data")
    def filter_on_pro(self, promoter, chrom):
        """

        :param chrom:
        :param promoter:
        :return:
        """
        return self.data[(self.data['start'] >= promoter[0]) & (self.data['end'] <= promoter[1]) &
                         (self.data['chr'] == "chr" + chrom)]

    @param.depends("genes", "gene_input")
    def filter_by_gene_input(self):
        """
        filter on list of genes based on typed input, keep strings that match beginning of input
        """
        filtered = [gene for gene in self.gene_selectables if gene.startswith(self.gene_input)]
        if len(filtered) != 0:
            self.gene_selectables = filtered

    @param.depends("genes", "chromosome")
    def filter_by_chrom(self):
        """
        filter selectable genes bases on selected chromosome
        """
        if self.chromosome is None or self.chromosome == []:
            self.gene_selectables = list(set(self.genes["Symbol"]))
        else:
            self.gene_selectables = list(set(self.genes[self.genes["Chromosome"].isin(self.chromosome)]["Symbol"]))

    @param.depends("gene_name", "genes", "chromosome", "gene_input")
    def set_options(self):
        self.filter_by_chrom()
        self.filter_by_gene_input()
        if len(self.gene_selectables) > 1000:
            self.param.gene_name.objects = [None]
        else:
            self.param.gene_name.objects = [None] + self.gene_selectables

    @param.depends("data", "gene_name")
    def update_data(self):
        if self.gene_name is None:
            return None
        else:
            [promoter, chrom] = self.get_pro_region()
            df = self.filter_on_pro(promoter, chrom)
            return df

    @param.depends("data", "genes", "chromosome", "gene_name", "gene_input")
    def show_coords(self):
        df = self.update_data()
        if df is None:
            return df
        else:

            return

    @param.depends("gene_input", "chromosome", "gene_name")
    def __panel__(self):
        self.set_options()
        self.update_data()
        return pn.Column(pn.Row((pn.widgets.MultiChoice.from_param(self.param.chromosome))),
                         pn.Row(pn.widgets.TextInput.from_param(self.param.gene_input),
                                pn.widgets.Select.from_param(self.param.gene_name)))

    def __str__(self):
        return """
        Viewer class
        """


if __name__ == "__main__":
    print(View)
