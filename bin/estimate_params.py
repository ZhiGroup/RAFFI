import math
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
import operator as op
import gzip
import numpy as np

temp_x_vals = []
w_rho_vals = []


def get_total_length(mapping_file):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    f = open(mapping_file)
    min_val = 1000000
    max_val = 0
    for line in f:
        vals  = line.split()
        if (float(vals[1]) < min_val):
            min_val = float(vals[1])
        if (float(vals[1]) > max_val):                
            max_val = float(vals[1])

    f.close()
    return max_val - min_val


def compute_mafs(vcf_input,max_window_size):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    f = gzip.open(vcf_input, "rt")
    entries_started = False
    num_haps = 0
    num_sites = 0
    expected_vals = []

    site_counter = 0
    mafs = []
    for line in f:
        if ('#CHROM' in line):
            entries_started = True
        elif (entries_started):
            _values = line.split()
            _pos = _values[1]
            alt = _values[4]
            if (len(alt.split(',')) > 1):
                continue
            i = 2
            while(i < len(_values) and _values[i] != 'GT'):
                i += 1
            i += 1
            tags = _values[7]
            i = 9
            num_ones = 0
            num_zeros = 0
            while (i < len(_values)):
                vals = _values[i].split('|')
                if (vals[0] == '1'):
                    num_ones = num_ones + 1
                elif (vals[0] == '0'):
                    num_zeros = num_zeros + 1
                if (vals[1] == '1'):
                    num_ones = num_ones + 1
                elif (vals[1] == '0'):
                    num_zeros = num_zeros + 1
                i = i + 1
            v =  min(num_zeros,num_ones)/float(num_zeros+num_ones)
            num_haps = num_zeros + num_ones
            mafs.append(v)
    f.close()
    num_sites = len(mafs)
    x_vals = [0]
    y_vals = [0]
    for o in range(1,max_window_size):
        window_size = o
        expected_vals = []
        for i in range(0,len(mafs),window_size):
            expected_maf = 0
            _sum = 0.0
            for j in range(0,window_size):
                if (i + j  >= len(mafs)):
                    continue
                _sum += mafs[i+j]
            for j in range(0,window_size):
                if (i + j >= len(mafs)):
                    continue
                if (_sum == 0):
                    expected_maf = expected_maf
                else:
                    expected_maf += (mafs[i+j] * mafs[i+j]/_sum)
            expected_vals.append(expected_maf)

        
        _a = np.array(expected_vals)

        pa = np.percentile(_a,1)
        rho_v = pa*pa + (1-pa)*(1-pa)
        x_vals.append(o)
        y_vals.append(rho_v)

    temp_x_vals = list(x_vals)
    w_rho_vals = list(y_vals)

    return num_haps,num_sites

def w_rho_func(_w):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    if (_w >= len(w_rho_vals)):
        return w_rho_vals[len(w_rho_vals)-1]
    else:
        return w_rho_vals[int(_w)]


def ncr(n,r):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    f = math.factorial
    return f(n) / f(r) / f(n-r)


def fp(e,N,L,rho,r,c,w):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    sum = 0.0
    #rho = w_rho_func(w)
    for i in range(0,c):
        sum += (ncr(r,i)* ((rho**(L/w))**i) * ( (1-rho**(L/w))**(r-i) ))
    return 1 - sum


def tp(er,N,L,rho,r,c,w):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    sum = 0.0
    for i in range(0,c):
        sum += (ncr(r,i) * (math.exp(-(er*L)/w)**i) * ((1-math.exp(-(L*er)/w))**(r-i)) )
    return 1 - sum


def compute_w(error_rate,N,L,rho,r,c,max_w=300):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    fun = lambda w: 0.5* N*(N-1)*fp(error_rate,N,L,w_rho_func(w),r,c,w) / tp(error_rate,N,L,w_rho_func(int(w)),r,c,w)
    lambda_val = 0.5* N*(N-1)
    w_min = -1
    _started = False
    for w in range(1,max_w):
        _tp = tp(error_rate,N,L,w_rho_func(int(w)),r,c,w)
        _fp = fp(error_rate,N,L,w_rho_func(w),r,c,w)
        if (round(_tp,2) -  0.5* N*(N-1)* round(_fp,2)) == 1:
            if (_started):
                continue
            else:
                w_min = w
                _started = True
        elif _started:
            return ((w_min +w-1) / 2)
    return -1


def main(vcf_path, map_path, error_rate, target_length, num_runs, num_successes):
    global math, lowess, op, gzip, np, get_total_length, compute_mafs, w_rho_func, fp, tp, main, compute_w, ncr, temp_x_vals, w_rho_vals

    vcf_input = vcf_path
    mapping_input = map_path
    min_length_cM = float(target_length)

    min_window_size = 2
    max_window_size = 300
    num_haps, total_sites = compute_mafs(vcf_input,max_window_size)
    tl = get_total_length(mapping_input)
    min_length_SNPs = ( min_length_cM * total_sites ) / tl
    
    rho_initial = 0.9

    w = compute_w(error_rate,num_haps,min_length_SNPs,rho_initial,num_runs,num_successes,max_window_size)
    if w == -1:
        w_rho_vals = lowess(w_rho_vals, temp_x_vals,frac=0.1,return_sorted=False)
        w = compute_w(error_rate,num_haps,min_length_SNPs,rho_initial,num_runs,num_successes,max_window_size)

    return max(w, min_window_size)


if __name__ == "__main__":
	import sys 
	vcf_path = sys.argv[1]
	map_path = sys.argv[2]
	error_rate = 0.001
	target_length = 5
	num_runs = 3
	num_successes = 1
	if len(sys.argv) > 3:
		error_rate = float(sys.argv[3])
	if (len(sys.argv) > 4):
		target_length = float(sys.argv[4])

	if (len(sys.argv)> 5):
		num_runs = int(sys.argv[5])
 
	if (len(sys.argv)> 6):
		num_successes = int(sys.argv[6])
 
	window_size = main(vcf_path, map_path, error_rate, target_length, num_runs, num_successes)

	print (window_size)
