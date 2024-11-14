# Rshiny_CZIArterialCellAtlas

Developing the RShiny webpage for the utility of bench scientist to perform the analysis. 

The webpage has been divided into 3 tabs.

## 1. UMAP

The UMAP tab allow the user to make UMAP, featureplots, dotplots and violin plots for different site of the artery as well as for the different genes with publication downloadable images. See the figure below:

<img width="1732" alt="Screenshot 2024-11-14 at 08 06 15" src="https://github.com/user-attachments/assets/3aa7f468-32bd-4c36-8b96-3af8d54fadca">

## 2. DEG

The DEG tab will allow the user to perform the differential for different tab items:

#### 1. Site
In this user can perform the differential across different site. The user can combine different site with combination like ARCH, and CARTOID can combined and compared ILIAC and DESC. 

#### 2. Cluster
In this user can perform the differential across different seurat cluster. The user can combine different celltypes and compare.

#### 3. Cell Type
In this user can perform the differential across different celltypes. The user can combine different celltypes and compare.

#### 4. Site Cluster
In this user can perform the differential between the site in combination with the cluster like ARCH_3 vs DESC_3. The user can combine different site clusters and compare.

#### 5. Site Celltype
In this user can perform the differential between the site in combination with the celltype like ARCH_Endo vs DESC_Macro. The user can combine different site clusters and compare.

The differential gene expression table is easily downloadable. Please see the figure below:

<img width="1764" alt="Screenshot 2024-11-14 at 08 20 33" src="https://github.com/user-attachments/assets/ebebfd74-de40-4379-8fb4-7215f4f136e0">


## 3. Module Score

Module Score tab include the combined score for the genes that you need to calculate the aggregate score. Upload the gene list with the header. This page will calculate the module score and plot the VlnPlot across all the sites and featureplot. All the image are downloadable. See the figure below:

<img width="1761" alt="Screenshot 2024-11-14 at 08 10 10" src="https://github.com/user-attachments/assets/b1b23f04-1968-4435-a7d5-de8661bdd5a8">
