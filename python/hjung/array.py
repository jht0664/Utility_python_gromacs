# merge all arraylists to one arraylist in order to remove time info.
# input: x_t is array along time (t1, t2, ...)
#		[[x1(t0), x2(t0), x3(t0), ...], [x1(t1), x2(t1), x3(t1)], ...]
# output: [x1(t0), x2(t0), x3(t0), ..., x1(t1), x2(t1), x3(t1), ...] 
# Example: number_1d = merge_t_to_1d(number_t_1d)
def merge_t_to_1d(x_t):
	# internally merge arraylist in a list
	import itertools
	# merge 
	return [i for i in itertools.chain.from_iterable(x_t)]

