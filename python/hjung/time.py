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
