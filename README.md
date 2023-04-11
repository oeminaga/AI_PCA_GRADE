# Artificial Intelligence Unravels Interpretable Malignancy Grades of Prostate Cancer on Histology Images

## CONTENT
* You need a source-code editor (e.g., Visual Studio Code from https://code.visualstudio.com)
* Python ≥3.8 (Either built-in or from https://www.python.org)
* Jupyter notebook (How to install: https://jupyter.org/install)

* Packages installation:
</br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/install_packages.sh</br>

* Extraction of TMA core images
   * Download QuPath </br> https://qupath.github.io
   * The instruction given by Andrew Janowczyk</br>http://www.andrewjanowczyk.com/de-array-a-tissue-microarray-tma-using-qupath-and-python/</br>

* Neural Architecture Search:
</br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/NAS_PlexusNet.ipynb</br>

* Model training and Inference on TMA images
</br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/main_notebook_for_model_development.ipynb</br>

* Inference on whole-slide images (the ImageScope's xml annotation file including the demarcation of cancer lesions is required)
   * Aperio ImageScope (Image Viewer and Annotion Tool): </br>https://www.leicabiosystems.com/us/digital-pathology/manage/aperio-imagescope/</br>
   * SCN files: </br>https://github.com/oeminaga/AI_PCA_GRADE/tree/main/PROCESS_SCN_WSI
   * SVS files: </br>https://github.com/oeminaga/AI_PCA_GRADE/tree/main/PROCESS_SVS_WSI</br>
   * TIFF and other formats: ⚙️🚧

* Statistical analyses -abstract version for clarity- </br>
  * Survival analyses and comparsion of cox models:
    </br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/Analysis/evaluation.R</br>
      - Please use R Studio to open this file. To install, please follow the instruction given in the link [https://posit.co/download/rstudio-desktop/]
  * Distribution of 64 representation features for our novel model:</br>
    https://github.com/oeminaga/AI_PCA_GRADE/blob/main/Analysis/FeatureDistributionAnalysis.ipynb</br>

Should you have issues, please open a thread in the issue section.

Please cite the following paper when you use these scripts
[PAPER]
