import matplotlib.pyplot as plt
import panel as pn
from panel.viewable import Viewer, Viewable
import param
import pandas as pd


class View(Viewer):
    pn.extension(design="material", sizing_mode="stretch_width")

    data = param.DataFrame()
    genes = param.DataFrame()

    gene_input = param.String(default="")

    # set default to empty string because values are determined by other input values
    gene_name = param.ObjectSelector(objects=[None], allow_None=True)

    # define input for all the chromosomes
    chromosome = param.ListSelector(objects=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                             "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"])

    # define range on chromosome
    chrom_range_down = pn.widgets.IntInput(value=0, start=0, name="start chromosome region")
    chrom_range_up = pn.widgets.IntInput(value=0, start=0, name="end chromosome region")

    def __init__(self, **params):
        super().__init__(**params)
        self.param.objects = self.data
        self.param.objects = self.genes

        self.gene_selectables = list(set(self.genes["Symbol"]))

        # set chrom bounds
        chrom_bound_up = max(self.data["end"] + 5)
        self.chrom_range_down.end = chrom_bound_up
        self.chrom_range_up.end = chrom_bound_up

        # self.all_pipeline = self.data.idf[
        #     (self.data['start'] >= self.chrom_range_down.value) &
        #     (self.data['end'] <= self.chrom_range_up.value) &
        #     (self.data['chr'].isin(["chr" + chrom for chrom in self.chromosome]))
        # ]

    @param.depends("gene_name", "genes", "chromosome", "chrom_range_down", "chrom_range_up")
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
        self.chromosome = input_gene_info['Chromosome'].tolist()

        self.chrom_range_down.value = promoter[0]
        self.chrom_range_up.value = promoter[1]

    @param.depends("data", "chromosome", "chrom_range_down", "chrom_range_up")
    def filter_on_pro(self):
        """

        :return:
        """
        return self.data[(self.data['start'] >= self.chrom_range_down.value) &
                         (self.data['end'] <= self.chrom_range_up.value) &
                         (self.data['chr'].isin(["chr" + chrom for chrom in self.chromosome]))]

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
            self.gene_selectables.sort()
            self.param.gene_name.objects = [None] + self.gene_selectables

    @param.depends("data", "gene_name")
    def update_data(self):
        if self.gene_name is not None:
            self.get_pro_region()
            return self.filter_on_pro()
        elif self.chromosome & (self.chrom_range_down.value - self.chrom_range_up.value) > 1000:
            return self.filter_on_pro()
        return None

    @staticmethod
    def create_barcode_counts(df, all_barcode):
        counts = []
        for barcode in all_barcode:
            counts.append(df["barcode"].to_list().count(barcode))

        fig, ax = plt.subplots(figsize=(4, 3))
        ax.bar(all_barcode, counts)
        ax.set(title="Barcode Counts")

        plt.setp(ax.get_xticklabels(), rotation=15, horizontalalignment='center', fontsize='x-small')

        return fig

    @param.depends("data", "chrom_range_down", "chrom_range_up")
    def show_coords(self):
        df = self.update_data()
        if df is not None:
            x_range = range(self.chrom_range_down.value, self.chrom_range_up.value)
            all_barcode = list(set(df["barcode"]))
            all_barcode.sort()

            barcode_bar = self.create_barcode_counts(df, all_barcode)

            return barcode_bar

    @param.depends("gene_input", "chromosome", "gene_name", "chrom_range_down", "chrom_range_up")
    def __panel__(self) -> Viewable:
        self.set_options()
        bar = self.show_coords()
        return pn.Column(pn.Row((pn.widgets.MultiChoice.from_param(self.param.chromosome))),
                         pn.Row(self.chrom_range_down, self.chrom_range_up),
                         pn.Row(pn.widgets.TextInput.from_param(self.param.gene_input),
                                pn.widgets.Select.from_param(self.param.gene_name)),
                         pn.Row(pn.pane.Matplotlib(bar)))

    def __str__(self):
        return """
        Viewer class
        """


if __name__ == "__main__":
    print(View)
