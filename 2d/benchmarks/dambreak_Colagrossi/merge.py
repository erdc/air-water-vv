

import os
import glob
import pandas as pd
os.chdir("\Users\vis.msk\Desktop\new flume_corrected\New set of runs\Twophaseflow\case2\d2\5")
extension = 'csv'
all_filenames = [i for i in glob.glob('*.'.format(extension))]
#combine all files in the list
combined = pd.concat([pd.read_csv('all_filenames')])
#export to csv
combined.csv("combined.csv", index=False, encoding='utf-8-sig')


