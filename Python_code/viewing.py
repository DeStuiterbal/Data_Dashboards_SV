import matplotlib
import matplotlib.pyplot as plt
import panel as pn
from panel.viewable import Viewer, Viewable
import plotly.express as px
import param
import pandas as pd
import seaborn as sns


class View(Viewer):
    pn.extension(design="material", sizing_mode="stretch_width")
    matplotlib.use("agg")
    pn.extension("plotly")

    data = param.DataFrame()
    genes = param.DataFrame()

    gene_input = param.String(default="", label="Input gene name")

    # set default to empty string because values are determined by other input values
    gene_name = param.ObjectSelector(objects=[None], allow_None=True, label="Select gene")

    # define input for all the chromosomes
    chromosome = param.ListSelector(objects=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                             "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"],
                                    label="Selected chromosome(s)", default=[])
    clear_button = param.Action(label="Clear all")

    def __init__(self, **params):
        super().__init__(**params)
        self.param.objects = self.data
        self.param.objects = self.genes

        self.gene_selectables = list(set(self.genes["Symbol"]))

        # chrom range is to give the ranges for each chromosome
        self.chrom_range = {}
        self.clear_button = self.clear_input

    @param.depends("clear_button")
    def clear_input(self, *args):
        """
        set all input to default values
        """
        self.gene_input = self.param.gene_input.default
        self.gene_name = self.param.gene_name.default
        self.chromosome = self.param.chromosome.default
        self.chrom_range = {}

    @param.depends("gene_name", "genes")
    def get_pro_region(self):
        """
        obtain the promoter based on the ncbi data
        :return: promoter and chromosome location of the given gene
        """
        self.chrom_range = {}
        input_gene_info = self.genes[self.genes["Symbol"] == self.gene_name]
        # get the chromosomes of the genes
        self.chromosome = input_gene_info['Chromosome'].tolist()

        for i, chrom in enumerate(self.chromosome):
            # get the beginning and end location of the gene
            promoter = [input_gene_info["Begin"].tolist()[i], input_gene_info["End"].tolist()[i]]

            # add or subtract from the promoter depending on the orientation of the gene
            if input_gene_info["Orientation"].tolist()[i] == "plus":
                promoter[0] -= 1000
            else:
                promoter[1] += 1000

            self.chrom_range.update({chrom: promoter})

    @param.depends("data", "chromosome")
    def filter_on_pro(self):
        """
        filter on data obtained by get_pro_region()
        :return: filtered data
        """
        data = [self.data[(self.data['start'] >= self.chrom_range[chrom][0]) &
                          (self.data['end'] <= self.chrom_range[chrom][1]) &
                          (self.data['chr'] == ("chr" + chrom))] for chrom in self.chromosome]
        return data

    @param.depends("genes", "gene_input", "gene_name")
    def filter_by_gene_input(self):
        """
        filter on list of genes based on typed input, keep strings that match beginning of input
        """
        filtered = [gene for gene in self.gene_selectables if gene.startswith(self.gene_input)]
        if len(filtered) != 0:
            self.gene_selectables = filtered

    @param.depends("genes", "chromosome", "gene_name")
    def filter_by_chrom(self):
        """
        filter selectable genes bases on selected chromosome
        """
        if self.chromosome is None or self.chromosome == []:
            self.gene_selectables = list(set(self.genes["Symbol"]))
        else:
            selectable = []
            for chromosome in self.chromosome:
                if selectable:
                    temp = list(set(self.genes[self.genes["Chromosome"] == chromosome]["Symbol"]))
                    selectable = [gene for gene in selectable if gene in temp]
                else:
                    selectable = list(set(self.genes[self.genes["Chromosome"] == chromosome]["Symbol"]))
            self.gene_selectables = selectable

    @param.depends("gene_name", "genes", "chromosome", "gene_input")
    def set_options(self):
        """
        filter selectable genes based on selected chromosomes
        """
        self.filter_by_chrom()
        self.filter_by_gene_input()
        if len(self.gene_selectables) > 200:
            self.param.gene_name.objects = [None]
        else:
            self.gene_selectables.sort()
            self.param.gene_name.objects = [None] + self.gene_selectables

    @param.depends("data", "gene_name", "chromosome")
    def update_data(self):
        """
        update the data
        :return: dataframes generated
        """
        if self.gene_name is not None:
            self.get_pro_region()
            out = self.filter_on_pro()
            return out
        return None

    @staticmethod
    def create_barcode_counts(df, all_barcode):
        """
        create barplot with counts for each barcode
        :param df: filtered pandas dataframe
        :param all_barcode: all barcodes in the dataframe
        :return: matplotlib barplot figure
        """
        counts = []
        for barcode in all_barcode:
            counts.append(df["barcode"].to_list().count(barcode))

        fig, ax = plt.subplots(figsize=(3, 3))
        ax.bar(all_barcode, counts)
        ax.set(title="Barcode Counts")

        plt.setp(ax.get_xticklabels(), rotation=15, horizontalalignment='center', fontsize='x-small')
        plt.close(fig)
        return fig

    def create_dataframe_points(self, df, all_barcode, chrom):
        """
        create heatmap for all the points in range
        :param df: filtered pandas dataframe
        :param all_barcode: list of all barcodes in the df
        :param chrom: chromosome on which the range sits
        :return: matplotlib heatmap figure
        """
        heat_df = pd.DataFrame(index=all_barcode, columns=range(self.chrom_range[chrom][0], self.chrom_range[chrom][1]))
        for barcode in all_barcode:
            heat_df.loc[barcode] = [1.0 if col in df[df["barcode"] == barcode]["start"].tolist() else 0.0 for col in
                                    heat_df.columns.tolist()]
        heat_df = heat_df.astype(float)
        print(heat_df)
        fig = px.imshow(heat_df)
        return fig

    @staticmethod
    def create_density(df):
        """
        create density plot
        :param df: filtered pandas dataframe
        :return: matplotlib density plot figure
        """
        fig, density_ax = plt.subplots(figsize=(10, 4))
        sns.kdeplot(data=df, x="start", hue="barcode")
        density_ax.set(title="Methylation density", xlabel="Position (bp)")
        plt.close(fig)
        return fig

    @param.depends("data", "chromosome")
    def show_coords(self):
        """
        manage all visualisations
        :return: list of lists with each list containing all of a kind of visualisations.
        """
        dfs = self.update_data()
        barcode_bar = []
        density_fig = []
        heatmap_fig = []
        if dfs is not None:
            for i, chromosome in enumerate(self.chromosome):
                df = dfs[i]
                # skip iteration if df is empty
                if df.empty:
                    continue
                # get all barcodes from dataframe
                all_barcode = list(set(df["barcode"]))
                all_barcode.sort()

                # append plots
                density_fig.append(self.create_density(df))
                barcode_bar.append(self.create_barcode_counts(df, all_barcode))
                heatmap_fig.append(self.create_dataframe_points(df, all_barcode, chromosome))

        figures = [barcode_bar, density_fig, heatmap_fig]

        return figures

    @param.depends("gene_input", "chromosome", "gene_name")
    def __panel__(self) -> Viewable:
        """
        visualize all objects
        :return: visualisation
        """
        self.set_options()
        figs = self.show_coords()
        visualise = pn.Column(pn.Row((pn.widgets.MultiChoice.from_param(self.param.chromosome))),
                              pn.Row(pn.widgets.TextInput.from_param(self.param.gene_input),
                                     pn.widgets.Select.from_param(self.param.gene_name)),
                              pn.Row(pn.widgets.Button.from_param(self.param.clear_button)))
        if self.chromosome is not None:
            for i, pic in enumerate(figs[1]):
                visualise.append(pn.Row(pn.pane.Str("Chromosome " + self.chromosome[i], styles={"font-size": "20pt"},
                                                    align="center")))
                visualise.append(pn.Row(pn.pane.Matplotlib(figs[0][i])))
                visualise.append(pn.Row(pn.pane.Matplotlib(figs[1][i])))
                visualise.append(pn.Row(pn.pane.Plotly(figs[2][i])))
                visualise.append(pn.Row(pn.layout.Divider()))
        return visualise

    def __str__(self):
        return """
        Viewer class
        """


if __name__ == "__main__":
    print(View)
