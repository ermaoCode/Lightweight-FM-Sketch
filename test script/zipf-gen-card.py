import sys
import numpy as np
import pandas as pd
import time

# print ("脚本名：", sys.argv[0])

if len(sys.argv) != 2:
	print ("1 parameter required!")
	sys.exit()

size = int(sys.argv[1])
print ("Total number: ", size)

a = 2
arr = np.random.zipf(a, size)
df = pd.DataFrame(arr)
df.to_csv("zipf.csv")
strarr = df[0].apply(lambda x: str(x) + '\n')
f = open('zipf.txt', 'w')
f.writelines(strarr)
f.close()

# print ("Generate completed!")

# start = time.clock()
se = df.iloc[:, 0]
nums = []

def count_cord(x):
	if x not in nums:
		nums.append(x)
	return x

se.apply(count_cord)
# end = time.clock()
# print ("Time used:")
# print (end - start)
print ("Cardinality:")
print (len(nums))