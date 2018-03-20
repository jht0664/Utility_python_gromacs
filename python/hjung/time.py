## timer
import time

# initialize starting time
# output: process tiem and wall time
# Example: start_process, start_wall = init_time()

def init():
	print("time.init: ")
	return time.process_time(), time.perf_counter() # process time, wall time

# print difference of times
# input: process, prof are from the function, init().
# Example: end_print(process, prof)
def end_print(process, prof):
	print("time.printf: ")
	elapsed_proc = time.process_time() - process
	elapsed_prof = time.perf_counter() - prof
	print(" process time = %s sec." %(int(elapsed_proc)))
	print(" performance time = %s sec." %(int(elapsed_prof)))
	sleep = int(elapsed_prof-elapsed_proc)
	print(" sleep time = %s sec." %sleep)
	return

# initialize processing bar
# output: process tiem and wall time
# Example: mod_frame = process_init()
def process_init():
	print("time.process_init: ")
	return 10 # initial mode number = 10

# print process bar
# input: itime, current time
#        ftime, final time
#        mod_time, frequency of printing.
# output: mod_time, new or old printing frequency 
# Example: mod_frame = process_print(itime+1, ftime, mod_frame)
def process_print(itime, ftime, mod_time):
	if itime%mod_time == 0:
		print("... {0} th frame reading ({1:.0%}) ...".format(itime,itime/ftime))
		if (itime/mod_time)%10 == 0:
			mod_time = mod_time*10
	return mod_time