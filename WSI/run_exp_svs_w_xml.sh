echo "INFO: Generate a list for svs files with the corresponding xml annotation files"
python get_file_list.py
echo "INFO: Read the files and Run the Prediction"
python main_10x_50_50.py
echo "INFO: Merge the csv output files into a single files"
python merge.py
echo "INFO: DONE -> EXIST"
