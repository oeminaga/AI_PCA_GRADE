# Artificial Intelligence Unravels Interpretable Malignancy Grades of Prostate Cancer on Histology Images

THE SCRIPTS 📝 ARE:
</br>

* Packages installation:
</br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/install_packages.sh</br>

* Extraction of TMA spot images
   * Download QuPath </br> https://qupath.github.io
   * Please follow the instruction given by Andrew Janowczyk</br>http://www.andrewjanowczyk.com/de-array-a-tissue-microarray-tma-using-qupath-and-python/</br>

* Neural Architecture Search:
</br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/NAS_PlexusNet.ipynb</br>

* Model training and Inference on TMA images
</br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/main_notebook_for_model_development.ipynb</br>

* Inference on whole-slide images (the ImageScope's xml annotation file for cancer lesions is required)
   * Aperio ImageScope (Image Viewer and Annotion Tool): </br>https://www.leicabiosystems.com/us/digital-pathology/manage/aperio-imagescope/</br>
   * SCN files: </br>https://github.com/oeminaga/AI_PCA_GRADE/tree/main/PROCESS_SCN_WSI
   * SVS files: </br>https://github.com/oeminaga/AI_PCA_GRADE/tree/main/PROCESS_SVS_WSI</br>
   * TIFF and other formats: ⚙️🚧

* Statistical analyses -abstract version for clarity- </br>
  * Survival analyses and comparsion of cox models:
    </br>https://github.com/oeminaga/AI_PCA_GRADE/blob/main/Analysis/evaluation.R</br>
  * Distribution of 64 representation features for our novel model:</br>
    https://github.com/oeminaga/AI_PCA_GRADE/blob/main/Analysis/FeatureDistributionAnalysis.ipynb</br>

Please cite the following paper when you use these script
[PAPER]
