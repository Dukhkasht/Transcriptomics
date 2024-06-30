# Identification and characterization of differentially expressed genes in Type 2 Diabetes and male infertility using in silico approach
Type 2 diabetes (T2D) and male infertility are both prevalent conditions with 
rising trends. Despite their recognized connection, the literature lacks in-depth pathway analysis specific to T2D and male infertility. The goal of this 
study was to explore the underlying biological mechanisms and potential 
intervention targets for diabetes-related male infertility.
T2D datasets (GSE38642, GSE25724, and GSE20966) and male infertility 
datasets (GSE6872 and GSE45887) were retrieved from the GEO database. 
After batch correction and normalization, we performed GSEA for both 
conditions separately. WGCNA identified key T2D and male infertility 
modules. Intersection of these key modules revealed 15 common genes: 
SKIC3, CSRNP3, OPA1, C11orf58, EIF4G2, MATR3, TUSC3, PCLO, BTG3, 
SNAP25, CLCN4, SALL2, SLC17A6, NKX2-2, and ISL1. Functional 
enrichment analysis indicated that RNA degradation, synaptic vesicle cycle, 
and oxidative phosphorylation pathways may regulate or affect both 
diseases. EIF4G2, ISL1, OPA1, PCLO, and SNAP25 emerged as potential 
therapeutic targets based on their involvement in enriched pathways and 
functional annotation.
This comprehensive analysis provides insights into the shared molecular mechanisms of T2D and male infertility, highlighting potential pathways and therapeutic targets for further research

## Methodology 
![image](https://github.com/Dukhkasht/Transcriptomics/assets/126152802/75e69260-65b5-44a0-9fd3-83f0e275acf1)

Scripts related to these analyses are present in [here](/Dissertation/Scripts)

The data for these analyses are retrieved from GEO website and the datasets (two datasets for male infertility and 3 for type 2 diabetes)

## Result

### Batch effect correction
- **For diabetes dataset**

<div style="display: flex;">
  <img src="https://github.com/Dukhkasht/Transcriptomics/assets/126152802/35715bc4-ecd5-463b-b5a3-220089fa5bb7" alt="Image 1" style="width: 300px; margin-right: 10px;">
  <img src="https://github.com/Dukhkasht/Transcriptomics/assets/126152802/bdbecf63-860b-4bfb-ae00-093bd9918a6e" alt="Image 2" style="width: 300px; margin-left: 10px;">
</div>

- **For male infertility dataset**

  <div style="display: flex;">
  <img src="https://github.com/Dukhkasht/Transcriptomics/assets/126152802/e234e5de-29cd-4ba4-8583-dab6f1cc88c4" alt="Image 1" style="width: 300px; margin-right: 10px;">
  <img src="https://github.com/Dukhkasht/Transcriptomics/assets/126152802/89704103-b649-4c1e-8824-31a5ea672676" alt="Image 2" style="width: 300px; margin-left: 10px;">
</div>

### WGCNA result
- **Type 2 diabetes**
![modulerelationt2d2](https://github.com/Dukhkasht/Transcriptomics/assets/126152802/cb109569-2504-4113-9d53-3162ebe28319)

- **Male Infertility**

![modulerelationmf](https://github.com/Dukhkasht/Transcriptomics/assets/126152802/e4f611ce-2d23-4a21-bf5b-806be0cd0393)

### Shared genes in both diseases
![venn](https://github.com/Dukhkasht/Transcriptomics/assets/126152802/959c5a16-4423-434d-b6c2-1b210784581d)


