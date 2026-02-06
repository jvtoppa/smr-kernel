import math
z = 813
n = 1995
sigma = 256
total_bits = z * (max(math.ceil(math.log2(n)),1) + max(math.ceil(math.log2(n / z)), 1) + max(math.ceil(math.log2(sigma)), 1) + 2)

print(total_bits)