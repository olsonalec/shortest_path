import numpy as np
import random
import time

#######################################

original = list(np.zeros(10000))
for i in range(len(original)):
    original[i] = random.randint(0, 123456)

copy = list(np.zeros(10000))
for i in range(len(copy)):
    copy[i] = random.randint(0, 123456)

start = time.time()
smallest = float('inf')
total = 0
for i in range(len(original)):
    if (original[i] not in copy):
        total += original[i]
    # if (original[i] < smallest) and (original[i] not in copy):
        # if original[i] not in copy:
        # smallest = original[i]
end = time.time()

print(smallest)
print(end - start)
print()

#######################################

class Number:
    def __init__(self, value):
        self.value = value

my_list = list(np.zeros(20000))
for i in range(len(my_list)):
    my_list[i] = Number(random.randint(0, 500000000))

start = time.time()
smallest = float('inf')
for i in range(len(my_list)):
    if my_list[i].value < smallest:
        smallest = my_list[i].value
end = time.time()

print(smallest)
print(end - start)