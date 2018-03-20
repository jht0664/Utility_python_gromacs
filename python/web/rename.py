from os import listdir
from os.path import isfile, join
import os

onlyfiles = [f for f in listdir('/home/htjung/card/') if isfile(join('/home/htjung/card/', f))]
onlyfiles.sort()
i_num = 0
for isheet in range(int(len(onlyfiles)/8)):
	temp_num = i_num
	os.rename(onlyfiles[i_num],str(temp_num)+'.jpg')
	os.rename(onlyfiles[i_num+2],str(temp_num+1)+'.jpg')
	os.rename(onlyfiles[i_num+4],str(temp_num+2)+'.jpg')
	os.rename(onlyfiles[i_num+6],str(temp_num+3)+'.jpg')
	os.rename(onlyfiles[i_num+5],str(temp_num+4)+'.jpg')
	os.rename(onlyfiles[i_num+7],str(temp_num+5)+'.jpg')
	os.rename(onlyfiles[i_num+1],str(temp_num+6)+'.jpg')
	os.rename(onlyfiles[i_num+3],str(temp_num+7)+'.jpg')
	i_num += 8

