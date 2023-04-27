import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random

np.random.seed(12)
random.seed(12)

nb_patients = np.load("files/nb_patients.npy")

#==================================================================================================================
def compute_confusion_matrix(file, first_run, second_run):
    patno_cluster_a = first_run
    patno_cluster_b = second_run
    a = 5
    b = 5
    print(a,b)
    confusion_matrix = np.matrix(np.zeros((a,b)))
    for cla in range(a):
        for clb in range(b):
            confusion_matrix[cla,clb] = len(set(patno_cluster_a[patno_cluster_a['group'] == cla+1].index) & set(patno_cluster_b[patno_cluster_b['group'] == clb+1].index))
    plt.figure()
    sns.heatmap(confusion_matrix.T, square = True, annot=True, xticklabels = ["cluster "+str(ii) for ii in range(1,a+1)], yticklabels = ["cluster "+str(ii) for ii in range(1,b+1)], fmt = '.0f' )
    plt.xlabel("4779 somamers")
    plt.ylabel("919 somamers")
    plt.title('Cluster intersection, nb_patients = '+ str(nb_patients))
    plt.margins(x=0.1, y=0.1)
    plt.savefig(file, bbox_inches='tight', dpi = 300)
    return confusion_matrix
#==================================================================================================================

#1st: PD-PROD all proteins
first_run = pd.read_csv("/shared-data/research/genomics/projects/ppmi_analysis/proteomics_analysis/proteomic_VAE/VAEs/PD_PROD_patients/optimized_clustering/02_patno_type_cluster_runid=191.csv").set_index('PATNO')
#2nd: PD-PROD most variable proteins
run_id = 154
second_run = pd.read_csv("files/res/02_patno_type_cluster_runid="+str(run_id)+".csv", index_col=0)
second_run = second_run.replace(3,30)
second_run = second_run.replace(4,40)
second_run = second_run.replace(40,3)
second_run = second_run.replace(30,4)


file = 'figs/most_var_prot_vs_all_prot_cluster_intersection.png'
confusion_matrix = compute_confusion_matrix(file, first_run, second_run)

