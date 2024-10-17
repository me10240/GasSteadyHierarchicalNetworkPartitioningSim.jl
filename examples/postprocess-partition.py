import json

import os
dirname = "./data/Texas7k_Gas/"
partition_file = dirname + "partition-test-script.json"
with open(partition_file, "r") as read_file:
    partition_data = json.load(read_file)

num_partitions = partition_data["num_partitions"]
num_small_blocks = 0
biggest_block = 0
block_size_array = []
for i in range(1, num_partitions+1):
    i_str = str(i)
    var = len(partition_data[i_str])
    if var == 2:
        num_small_blocks += 1
    biggest_block = max(biggest_block, var)
    if var not in block_size_array:
        block_size_array.append(var)

print(num_partitions, "\n", biggest_block, "\n", num_small_blocks, "\n", block_size_array)
    


