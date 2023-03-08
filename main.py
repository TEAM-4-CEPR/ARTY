from __future__ import print_function


__author__ = "Ilango Guy"
__version__ = "0.0.3"
__licence__ = "CEPR"
__doc__ = "Window for data analysis"
__draw__ = "Cezard Adeline"


"""
Change from 0.0.2
@ A. Cezard made some cool draw to the GUI ! 
GO now working
Change from 0.0.1
@ This script no longer use PydesEq2 (cause I can't fucking export the DE matrix , and tbh tf is a pickle format????) , now use DESEQ2 from R and pass all args to R. 
Now all can be done from the main menu (Y)
"""


"""
Library
"""
import sys
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
import pandas as pd
from goatools.anno.genetogo_reader import Gene2GoReader
import mygene
from genes_ncbi_10090_proteincoding import GENEID2NT as GeneID2nt_mus
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from scipy import stats
import pandas as pd
import argparse
import re
import os
import numpy as np
import pickle as pkl
from pandas.api.types import is_numeric_dtype
from pathlib import Path
from tkinter import *
from tkinter.filedialog import *
import plotly.express as px
import plotly.graph_objects as go
from sklearn.impute import SimpleImputer
from diffexpr.py_deseq import py_DESeq2
from scipy import stats
from sklearn.decomposition import PCA
"""
@ TODO : TEST tools and debug, curate code .... , I should have done some POO ... 
"""
def PCA3d():
	"""
	__doc__ = Plot a 3D PCA from chosen tables (count and design)  
	"""
	df = pd.read_table(path_to_data.get() , sep ="\t")
	id_list =  list(df["id"])
	df=df.set_index("id")
	df=df.T
	
	design = pd.read_table(target_data.get(), sep = ",")
	print(df)
	pval_dict = {}
	for i in df:
		shapiro_test = stats.shapiro(df[i])	
		pval_dict[i] = shapiro_test.pvalue

	print(pval_dict)
	val_list= []
	for j in df.columns:
		val_list.append(j)
	X = df[val_list]

	pca = PCA(n_components=3)
	components = pca.fit_transform(X)
	total_var = pca.explained_variance_ratio_.sum() * 100
	print(components)
	fig = px.scatter_3d(
    components, x=0, y=1, z=2, color=design["sample"].to_list(),
    title=f'Total Explained Variance: {total_var:.2f}%',
    labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
)
	fig.write_html("pca3d.html")


def PCA2d():
	"""
	__doc__ = Plot a 2D pCA from Chosen tables (count and design)
	"""
	df = pd.read_table(path_to_data.get() , sep ="\t")
	id_list =  list(df["id"])
	df=df.set_index("id")
	df=df.T
	
	design = pd.read_table(target_data.get(), sep = ",")
	print(df)
	pval_dict = {}
	for i in df:
		shapiro_test = stats.shapiro(df[i])	
		pval_dict[i] = shapiro_test.pvalue

	print(pval_dict)
	val_list= []
	for j in df.columns:
		val_list.append(j)
	X = df[val_list]

	pca = PCA(n_components=2)
	components = pca.fit_transform(X)
	print(components)
	fig = px.scatter(components, x=0, y=1, color = design["sample"].to_list())
	fig.write_html("pca.html")


def  Differential_Analysis():
	"""
	__doc__ = use diffexpr (Python implementation of DesEQ2) to run a DEA on
 		count matrix. It also create the design matrix
	"""
	df = pd.read_table(path_to_data.get())
	print(df.head())
	sample_df = pd.DataFrame({'samplename': df.columns}) \
        .query('samplename != "id"')\
        .assign(sample = lambda d: d.samplename.str.extract('([AB])_', expand=False)) \
        .assign(replicate = lambda d: d.samplename.str.extract('_([123])', expand=False)) 
	sample_df.index = sample_df.samplename
	sample_df.to_csv("design.csv")
	print(df.head())
	dds = py_DESeq2(count_matrix = df,
               design_matrix = sample_df,
               design_formula = '~ replicate + sample',
               gene_column = 'id') 
    
	dds.run_deseq() 
	dds.get_deseq_result(contrast = ['sample','B','A'])
	res = dds.deseq_result 
	print(res.head())
	dds.normalized_count()
	dds.comparison 
	lfc_res = dds.lfcShrink(coef=4, method='apeglm')
	print(lfc_res)
	lfc_res["id"] = df["id"]
	print(lfc_res)
	lfc_res.to_csv("DE.csv")
	


def heatmap():
	"""
	__doc__ = Create a Heatmap (what else ? ) 
	"""
	try:
		os.mkdir("HEATOMAP")
	except OSError as error:
		print(error)
	df = pd.read_csv(path_to_data.get() , sep ="\t")
	id_list =  list(df["id"])
	print(df)
	feature_name = []
	df.pop(df.columns[0])
	print(df)
	val_list = []
	for j in df.columns:
		val_list.append(j)
	print(len(val_list) , len(id_list))
	print(val_list, id_list)
	fig = px.imshow(df, labels=dict(x="ID", y="Gene", color='Expression'), y = id_list,\
    x=val_list, width=700, height=2800)
	fig.update_xaxes(side='top')
	fig.write_html("./HEATOMAP/heatmap.html")



def volcanoPlot():
	"""
	__doc__ = Create a volcano plot , automatically create groups : 					['Significatively_down_regulated','Significatively_up_regulated', 
'Non_significatively_up_regulated', 'non_significatively_down_regulated']
	"""
	df = pd.read_csv(path_to_data.get(), sep=",")
	conditions = [
    (df['log2FoldChange'] < 0) & (df["padj"] <= 0.05),
    (df['log2FoldChange'] > 0) & (df['padj'] <= 0.05),
    (df['log2FoldChange'] > 0) & (df['padj'] > 0.05),
    (df['log2FoldChange'] < 0) & (df['padj'] > 0.05)]
	values = ['Significatively_down_regulated', 'Significatively_up_regulated', 'Non_significatively_up_regulated', 'non_significatively_down_regulated']
	df['info'] = np.select(conditions, values)
	df.loc[df["info"] == '0', "info"] = "Other"
	groups = df['info'].unique()
	print(df['info'])
	print(df["id"])
	data = []
	for group in groups:
		df_group = df[df['info'] == group]
		trace = go.Scatter(x=df_group['log2FoldChange'], 
                        y=-np.log10(df_group['padj']),
                        mode='markers',
                        name=group)
		data.append(trace)

	layout = go.Layout(title='Grouping')
	fig = go.Figure(data=data, layout=layout)
	fig.show()
	

	
def goterm():
	"""
	__doc__ = use GOatool to retrieve all ontology from up and downregulated gene/prot(/metabolite <- SHIT)
	"""
	obo_fname = download_go_basic_obo()
	file_gene2go = download_ncbi_associations()
	obodag = GODag("go-basic.obo")
	mg = mygene.MyGeneInfo()
	objanno = Gene2GoReader(file_gene2go, taxids=[10090])
	ns2assoc = objanno.get_ns2assc()
	for nspc, id2gos in ns2assoc.items():
		print("{NS} {N:,} annotated mouse genes".format(NS=nspc, N=len(id2gos)))
	goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_mus.keys(),  # List of mouse protein-coding genes
        ns2assoc,  # geneid/GO associations
        obodag,  # Ontologies
        propagate_counts=False,
        alpha=0.05,  # default significance cut-off
        methods=['fdr_bh'])  # defult multipletest correction method
    
	df = pd.read_csv(path_to_data.get() , sep = ",")
	df_down  = df[df["log2FoldChange"] < 0 ] 
	df_up = df[df["log2FoldChange"] > 0 ]
	print(df_down.shape)
	print(df_up.shape)
	xli_up = df_up['id'].to_list()
	xli_down = df_down['id'].to_list()
	#print(df['id'].to_list())
	out_up = mg.querymany(xli_up, scopes='symbol', fields='entrezgene', species='mouse')
	out_down = mg.querymany(xli_down, scopes='symbol', fields='entrezgene', species='mouse')
	entrez_id_up = []
	entrez_id_down = []

	for i in out_up:
		#print(i)
		try:
			entrez_id_up.append(int(i["_id"]))
		except (KeyError , ValueError) as error:
			pass
	for i in out_down:
		#print(i)
		try:
			entrez_id_down.append(int(i["_id"]))
		except (KeyError , ValueError) as error:
			pass
	
	geneids_study_up = entrez_id_up
	
	goea_results_all = goeaobj.run_study(geneids_study_up)
	goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
	goeaobj.wr_xlsx("GO_up.xlsx", goea_results_sig)
	geneids_study_down = entrez_id_down
	
	goea_results_all = goeaobj.run_study(geneids_study_down)
	goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
	goeaobj.wr_xlsx("GO_down.xlsx", goea_results_sig)
    

def get_data_frame():
	"""
	browser from tk
	"""
	filename = askopenfilename(title="Find text file", filetypes=[('txt files', '.txt'), ('all files', '.*')])
	path_to_data.set(filename)
def get_target_frame():
	"""
	browser from tk
	"""
	filename = askopenfilename(title="Find text file", filetypes=[('txt files', '.txt'), ('all files', '.*')])
	target_data.set(filename)

fenetre = Tk()
fenetre.title('Nom Classe')
fenetre.geometry("960x540")
fenetre.configure(bg='black')
bgimg = PhotoImage(file = "interface.png")
limg = Label(fenetre , i = bgimg)
limg.place(x=0 , y =0)
window_label = Label(fenetre, text="Proteomics/Transcriptomics analysis")
window_label.pack()

# Frame
Frame_super  = Frame(fenetre ,bg ="orange", borderwidth = 2 , relief = GROOVE)
Frame_super.pack(side=TOP, anchor=NW)
Frame_data = Frame(Frame_super,bg ="orange")
Frame_data.pack()

Frame_target = Frame(Frame_super,bg ="orange")
Frame_target.pack()


Frame_button = Frame(fenetre, borderwidth = 2 , relief = GROOVE)
Frame_button.place(in_=fenetre, anchor="c", relx=.5, rely=.5)

#In Frame_data
path_label = Label(Frame_data, bg ="orange", text="Select dataframe")
path_label.pack()
path_to_data = StringVar()
path_to_data.set("Path_to_the_counts")
path_t = Entry(Frame_data, textvariable=path_to_data, width=30)
path_t.pack()
sample_data_image = PhotoImage(file = "sample_data.png")
path_button = Button(Frame_data ,  image =sample_data_image , command = get_data_frame)
path_button.pack()

#In Frame_design
target_label = Label(Frame_target,  bg ="orange",text="Select design matrix")
target_label.pack()
target_data = StringVar()
target_data.set("write col names of features")
target_t = Entry(Frame_target, textvariable=target_data, width=30)
target_t.pack()
design_matrix_image = PhotoImage(file = "design.png")
target_button = Button(Frame_target , image = design_matrix_image,  command = get_target_frame)
target_button.pack()

"""
value = StringVar()
value.set("Number of condition")
entree = Entry(fenetre, textvariable=value, width=30)
entree.pack()
"""
pca_2d_image = PhotoImage(file = "pca2d.png")
pca_3d_image = PhotoImage(file = "pca3d.png")
go_image = PhotoImage(file = "Go.png")
heatmap_image = PhotoImage(file = "heatmap.png")
de_image = PhotoImage(file = "De.png")
gsea_image = PhotoImage(file = "gsea.png")
volcano_image = PhotoImage(file = "volcano.png")




pca3d_button = Button(Frame_button, image = pca_3d_image, command = PCA3d)
pca2d_button = Button(Frame_button , image = pca_2d_image, command = PCA2d)
GO = Button(Frame_button , image = go_image, command = goterm)
de = Button(Frame_button , image = de_image , command = Differential_Analysis)
heatomap = Button(Frame_button, image =heatmap_image , command = heatmap)
volcoco = Button(Frame_button, image =volcano_image, command = volcanoPlot)
exit=Button(Frame_button, text="Fermer", command=fenetre.quit)


GO.pack()
de.pack()
pca2d_button.pack()
pca3d_button.pack()
heatomap.pack()
volcoco.pack()
exit.pack()
fenetre.mainloop()
