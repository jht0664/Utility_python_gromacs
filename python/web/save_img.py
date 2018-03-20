import shutil
import requests

basic_url = 'http://www.sejonghakdang.org/storage/subTextbook/card'
end_url = '.jpg'
exit_loop = False
for num0 in range(33):
	print("current {}".format(num0))
	exit_loop = False
	if num0 < 10:
		num00 = '0' + str(num0)
	else:
		num00 = str(num0)
	for num1 in range(1,100):
		if exit_loop:
			break
		for num2 in range(1,3):
			response = requests.get(basic_url+num00+'_'+str(num1)+'_'+str(num2)+end_url, stream=True)
			if '404' in str(response):
				del response
				exit_loop = True
				break
			if num1 < 10:
				num11 = '0' + str(num1)
			else:
				num11 = str(num1)
			with open('card'+num00+'_'+num11+'_'+str(num2)+end_url, 'wb') as out_file:
				shutil.copyfileobj(response.raw, out_file)
			del response

